function [EEG, KeepTime] = run_prep(EEG, opts)
% RUN_PREP  Drop non-EEG channels, resample, excise zero patches, remove DC offset.
%
%   The steps here are shared, verbatim, by bidsfun_detect_badchans and
%   bidsfun_hp_zap_cleanline: whatever bad-channel detection sees has to be exactly
%   what filtering later sees, or the two stages are answering different questions.
%   Bad-channel detection and line-noise removal are NOT here - each lives directly in
%   the bidsfun that owns it.
%
% USAGE:
%   [EEG, KeepTime] = run.run_prep(EEG)
%   [EEG, KeepTime] = run.run_prep(EEG, targetsrate=125, KeepTime=KeepTime)
%
% OPTIONAL NAME-VALUE:
%   noteegchannels    channel indices to drop before processing      (default 257:300)
%   targetsrate       resample to this rate in Hz; 0 = skip          (default 0)
%   removeDC          apply DC-removal filter                        (default true)
%   zeropatchseconds  cut out all-zero patches (amplifier crash padding) longer than
%                     this many seconds; 0 = skip                    (default 5)
%   KeepTime          struct of prior timings to merge into the output (default struct())
%   ramsaver          apply the DC-removal filter one channel at a time instead of to
%                     the whole data matrix at once, trading speed for a lower peak
%                     memory footprint                                (default false)
%
% See also: run.excise_zero_patches, run.restore_zero_patches

arguments
    EEG                     struct
    opts.noteegchannels     (1,:) double  = 257:300
    opts.targetsrate        (1,1) double  = 0
    opts.removeDC           (1,1) logical = true
    opts.zeropatchseconds   (1,1) double  = 5
    opts.KeepTime           struct        = struct()
    opts.ramsaver           (1,1) logical = false
end

KeepTime = opts.KeepTime;

%%% Optional downsampling
if opts.targetsrate > 0 && EEG.srate ~= opts.targetsrate
    D = tic; fprintf('Resampling %d → %d Hz ...\n', EEG.srate, opts.targetsrate)
    EEG = pop_resample(EEG, opts.targetsrate);
    KeepTime.Downsample = toc(D);
end

%%% Cut out all-zero patches (amplifier crash padding) so no filter ever sees them.
%%% Must happen before the DC filter: filtfilt would smear them into large transients.
if opts.zeropatchseconds > 0
    EEG = run.excise_zero_patches(EEG, 'minseconds', opts.zeropatchseconds);
end

%%% DC removal
if opts.removeDC
    fprintf('Building filters (srate = %d Hz) ...\n', EEG.srate)
    EEG_DCFilter_NumDen = filterbank(EEG.srate, 'DC_RCSquareFilt');

    D = tic; fprintf('\nDC removal ...\n')
    if opts.ramsaver
        % Overwrite each channel's row in place - EEG.data already exists at full size,
        % so this never holds a second full-size copy alongside it.
        for iCh = 1:EEG.nbchan
            EEG.data(iCh,:) = filtfilt(EEG_DCFilter_NumDen(1,:), EEG_DCFilter_NumDen(2,:), double(EEG.data(iCh,:)));
        end
    else
        EEG.data = filtfilt(EEG_DCFilter_NumDen(1,:), EEG_DCFilter_NumDen(2,:), double(EEG.data'))';
    end
    KeepTime.DCRemoval = toc(D);
end

end
