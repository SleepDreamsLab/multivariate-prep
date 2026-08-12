function [EEG, KeepTime] = run_filter(EEG, opts)
% RUN_FILTER  Resample, remove DC offset, and suppress line noise on an EEG struct.
%
% USAGE:
%   EEG = run.run_filter(EEG)
%   [EEG, KeepTime] = run.run_filter(EEG, removeDC=true, zapline=true)
%
% INPUTS:
%   EEG   — EEGLAB EEG struct (already imported by the caller)
%
% OPTIONAL NAME-VALUE:
%   noteegchannels  channel indices to drop before processing  (default 257:264)
%   targetsrate     resample to this rate in Hz; 0 = skip      (default 125)
%   removeDC        apply DC-removal filter                    (default true)
%   zapline        apply Zapline-plus line-noise removal      (default true)
%   KeepTime        struct of prior timings to merge into the output (default [])
%   cleanline       apply CleanLine after Zapline               (default true)
%   JsonFile        full path for a JSON sidecar with processing timings;
%                   '' = skip                                  (default '')
%
% Methods secton:
%
% Continuous EEG (256 channels, 250 Hz) was high-pass filtered to remove DC offset and then 
% cleaned of power-line artifacts using Zapline-plus (Klug & Kloosterman, 2022), an extension 
% of the Zapline algorithm (de Cheveigné, 2020). Zapline separates the data into a 
% line-frequency-dominated and a residual component using a comb filter, isolates the artifactual 
% subspace via DSS applied to the line frequency and its harmonics, and removes it by spatial 
% projection, thereby avoiding the spectral distortion introduced by notch filtering. The target 
% frequency was fixed at 50 Hz; at a 250 Hz sampling rate the fundamental and its second harmonic 
% (100 Hz) were removed jointly. To accommodate non-stationarity of the line artifact across the 
% [8-h] recordings, data were processed in fixed 300-s chunks, with the individual noise peak 
% re-estimated within each chunk (search window 50 ± 3 Hz) and the number of removed components 
% determined adaptively per chunk by iterative outlier detection on the DSS component scores 
% (minimum 1 component). The detection threshold was adapted iteratively until the cleaned 
% spectrum showed neither residual noise above nor over-correction below the surrounding-noise 
% criterion. Residual artifacts at 50 and 100 Hz were subsequently attenuated using 
% a custom parallelised reimplementation CleanLine (Mullen, 2012; EEGLAB), leveraging 
% multi-taper sinusoidal regression, .
% Within 4-s windows advanced in 2-s steps, the Thomson F-test for a deterministic sinusoid 
% was evaluated at 50 and 100 Hz using 7 Slepian tapers (time–bandwidth product 4, 
% resolution bandwidth 2 Hz) and an FFT zero-padded to 4096 points (0.061 Hz bin spacing). 
% Sinusoids reaching p < 0.01 were fitted by least squares and subtracted; target frequencies 
% not reaching significance in a given window and channel were left unmodified. Fitted 
% components were crossfaded sigmoidally across the 2-s overlap between adjacent windows.
% Unlike the default CleanLine implementation, which subtracts an estimated sinusoid at each 
% target frequency regardless of significance, subtraction here was conditional on the F-test.

arguments
    EEG                  struct
    opts.noteegchannels  (1,:) double  = 257:264
    opts.targetsrate     (1,1) double  = 0
    opts.removeDC        (1,1) logical = true
    opts.zapline         (1,1) logical = true
    opts.KeepTime        struct        = struct()
    opts.cleanline       (1,1) logical = true
    opts.JsonFile        (1,1) string  = ""
    opts.noisefreqs                    = 50
    opts.adaptiveNremove (1,1) logical = true
    opts.fixedNremove    (1,1) double  = 1
    opts.chunkLength     (1,1) double  = 300
    opts.plotResults     (1,1) logical = true
end

KeepTime = opts.KeepTime;

%%% Drop non-EEG channels
EEG = pop_select(EEG, 'nochannel', intersect(1:EEG.nbchan, opts.noteegchannels));

%%% Optional downsampling
if opts.targetsrate > 0 && EEG.srate ~= opts.targetsrate
    D = tic; fprintf('Resampling %d → %d Hz ...\n', EEG.srate, opts.targetsrate)
    EEG = pop_resample(EEG, opts.targetsrate);
    KeepTime.Downsample = toc(D);
end

%%% Build filters
fprintf('Building filters (srate = %d Hz) ...\n', EEG.srate)
EEG_DCFilter_NumDen = filterbank(EEG.srate, 'DC_RCSquareFilt');

%%% DC removal
if opts.removeDC
    D = tic; fprintf('\nDC removal ...\n')
    EEG.data = filtfilt(EEG_DCFilter_NumDen(1,:), EEG_DCFilter_NumDen(2,:), double(EEG.data'))';
    KeepTime.DCRemoval = toc(D);
end

%%% Zapline
if opts.zapline
    D = tic; fprintf('\nZapline plus ...\n')
    [EEG.data, zaplineConfig, analyticsResults] = clean_data_with_zapline_plus( ...
        double(EEG.data), EEG.srate, ...
        'noisefreqs',      opts.noisefreqs, ...
        'adaptiveNremove', opts.adaptiveNremove, ...
        'fixedNremove',    opts.fixedNremove, ...
        'chunkLength',     opts.chunkLength, ...
        'plotResults',     opts.plotResults);
    EEG.etc.zapline.config    = zaplineConfig;
    EEG.etc.zapline.analytics = analyticsResults;    
    KeepTime.Zapline = toc(D);
    fprintf('ZapLine-plus: %.2f min\n', KeepTime.Zapline / 60);
end

%%% Cleanline
if opts.cleanline 

%     % base CleanLine directly on what ZapLine actually found and treated
%     linefreqs = zaplineConfig.noisefreqs;
%     fprintf('ZapLine-plus detected and cleaned: %s Hz\n', mat2str(round(linefreqs,2)));
% 
%     % only pass freqs that had a real noise ratio to CleanLine
%     keepIdx = analyticsResults.ratioNoiseRaw > 1.5;   % adjust threshold based on what you see
%     linefreqs = zaplineConfig.noisefreqs(keepIdx);
%     fprintf('Only keep those with noise ratio >1.5: %s Hz\n', mat2str(round(linefreqs,2)));

    D = tic; fprintf('\nClean line ...\n')
    % EEG.data = double(EEG.data);
    EEG = cleanline_fast(EEG, 'linefreqs', [opts.noisefreqs opts.noisefreqs*2], ...
        'winsize', 4, 'winstep', 2, 'sigtest', true, 'pad', 2);

    % EEG = pop_cleanline(EEG, ...
    %     'chanlist', 1:EEG.nbchan, ...
    %     'linefreqs', 50, ...
    %     'winsize', 4, ...
    %     'winstep', 2, ...
    %     'computepower', false);
    KeepTime.Cleanline = toc(D);

end

%%% Write JSON sidecar
if opts.JsonFile ~= ""
    sidecarjson(KeepTime, char(opts.JsonFile));
end
end
