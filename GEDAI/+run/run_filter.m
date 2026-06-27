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

arguments
    EEG                  struct
    opts.noteegchannels  (1,:) double  = 257:264
    opts.targetsrate     (1,1) double  = 0
    opts.removeDC        (1,1) logical = true
    opts.zapline         (1,1) logical = true
    opts.KeepTime        struct        = struct()
    opts.cleanline       (1,1) logical = true
    opts.JsonFile        (1,1) string  = ""
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
        'noisefreqs', 'line', ...
        'maxfreq', floor(EEG.srate/2) - 5, ...   % let zapline search up to Nyquist itself, not just default 99 Hz        
        'plotResults', 0);
    KeepTime.Zapline = toc(D);
    fprintf('ZapLine-plus: %.2f min\n', KeepTime.Zapline / 60);
end

%%% Cleanline
if opts.cleanline 

    % base CleanLine directly on what ZapLine actually found and treated
    linefreqs = zaplineConfig.noisefreqs;
    fprintf('ZapLine-plus detected and cleaned: %s Hz\n', mat2str(round(linefreqs,2)));

    % only pass freqs that had a real noise ratio to CleanLine
    keepIdx = analyticsResults.ratioNoiseRaw > 1.5;   % adjust threshold based on what you see
    linefreqs = zaplineConfig.noisefreqs(keepIdx);
    fprintf('Only keep those with noise ratio >1.5: %s Hz\n', mat2str(round(linefreqs,2)));

    D = tic; fprintf('\nClean line ...\n')
    EEG.data = double(EEG.data);
    EEG = pop_cleanline(EEG, ...
        'bandwidth', 2, ...
        'chanlist', 1:EEG.nbchan, ...
        'computepower', 1, ...
        'linefreqs', linefreqs, ...
        'normSpectrum', 0, ...
        'p', 0.01, ...
        'pad', 2, ...
        'plotfigures', 0, ...
        'scanforlines', 1, ...
        'sigtype', 'Channels', ...
        'tau', 100, ...
        'verb', 1, ...
        'winsize', 4, ...
        'winstep', 2);
    KeepTime.Cleanline = toc(D);

end

%%% Write JSON sidecar
if opts.JsonFile ~= ""
    sidecarjson(KeepTime, char(opts.JsonFile));
end
end
