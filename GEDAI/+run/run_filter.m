function [EEG, KeepTime] = run_filter(eegFile, opts)
% RUN_FILTER  Import EEG from file, resample, remove DC offset, and suppress line noise.
%
% USAGE:
%   EEG = run.run_filter(eegFile)
%   [EEG, KeepTime] = run.run_filter(eegFile, removeDC=true, removeLN=true)
%
% INPUTS:
%   eegFile   — full path to an EEG file readable by eeg_import (.vhdr etc.)
%
% OPTIONAL NAME-VALUE:
%   noteegchannels  channel indices to drop before processing  (default 257:264)
%   targetsrate     resample to this rate in Hz; 0 = skip      (default 125)
%   removeDC        apply DC-removal filter                    (default true)
%   removeLN        apply Zapline-plus line-noise removal      (default true)
%   JsonFile        full path for a JSON sidecar with processing timings;
%                   '' = skip                                  (default '')

arguments
    eegFile              (1,1) string
    opts.noteegchannels  (1,:) double  = 257:264
    opts.targetsrate     (1,1) double  = 125
    opts.removeDC        (1,1) logical = true
    opts.removeLN        (1,1) logical = true
    opts.JsonFile        (1,1) string  = ""
end

KeepTime = struct();

%%% Import EEG
D = tic; fprintf('\nEEG import ...\n')
EEG = eeg_import(char(eegFile));
KeepTime.EEGimport = toc(D);

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
if opts.removeLN
    D = tic; fprintf('\nZapline plus ...\n')
    [EEG.data, ~, ~] = clean_data_with_zapline_plus( ...
        EEG.data, EEG.srate, ...
        'noisefreqs', 'line', ...
        'plotResults', 0);
    KeepTime.Zapline = toc(D);
    fprintf('ZapLine-plus: %.2f min\n', KeepTime.Zapline / 60);
end

%%% Write JSON sidecar
if opts.JsonFile ~= ""
    sidecarjson(KeepTime, char(opts.JsonFile));
end
end
