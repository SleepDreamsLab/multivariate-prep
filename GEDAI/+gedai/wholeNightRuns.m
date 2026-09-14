function runs = wholeNightRuns(mode, modeBB)
% WHOLENIGHTRUNS  GEDAI run configuration for RunMode 'WholeNight'.
%
%   runs = gedai.wholeNightRuns()
%   runs = gedai.wholeNightRuns(mode, modeBB)
%
% gedai.defaultRuns sets a different cleaning strength per sleep stage (auto for
% N3/N2/REM, auto+ for Wake), which is exactly what 'WholeNight' cannot use: one pass
% covers every stage, and bidsfun_gedai's resolveStageMode refuses to silently pick one
% of two disagreeing modes. This returns the same configuration with a single strength.
%
% The strength that is safe to use for the whole night is the SLEEP one. Wake is then
% handled by the sliding threshold (bidsfun_gedai's 'thresholdwindow'), not by a
% stronger mode - which is the point: cleaning follows the data rather than the
% hypnogram, so the scoring never enters the cleaning, and stage contrasts stop being
% confounded by having been processed differently.
%
% Measured on sub-drop0001 ses-t1 (9.2 h, 243 ch) against the stage-specific output:
%   slow waves 100.0 % of reference (99.6 % at the 10th percentile), K-complexes
%   99.9 %, spindles 100.1 %, N3 delta 714 vs 717 uV^2, N2 slightly cleaner.
%   The cost is wake: ~9 % more residual variance and ~31 % more residual wake delta.
%
% Do NOT reach for a stronger mode to close that wake gap. With a global threshold
% 'auto+' is what destroys slow waves (N3 delta 546 vs 717 uV^2 over the full night,
% slow waves down to 68 % of reference at the 10th percentile), and even with the
% sliding threshold it costs K-complexes (81 % of reference) and boundary-adjacent slow
% waves. If wake really must be cleaner, raise modeBB to 'auto' and check the
% K-complexes - that setting cleans wake better than the stage-specific reference but
% drops K-complexes to 96 % (absolute median 93 vs 115 uV).
%
% IF YOU ARE ANALYSING SLEEP-TO-WAKE TRANSITIONS, set the aggregation to 'min':
%
%   bidsfun_gedai(..., 'bandopts', struct('thresh_window_aggregate','min'))
%
% The default 3-window moving average lets the wake side's threshold reach back over the
% sleep before an awakening. Measured over 16 genuine sleep->wake transitions, retained
% delta 90 s before the first wake epoch was 0.296 of input with 'mean' against 0.563
% for the stage-specific reference - i.e. the sliding threshold is WORSE than stage
% cleaning there. With 'min' it is 0.539, back at parity, and sigma is the best of any
% variant across the whole transition. This is not a tuning preference; without it the
% method is actively unsuitable for transition analyses.
%
%   mode    strength for the wavelet bands. Default 'auto'.
%   modeBB  strength for the broadband pass.  Default 'auto-'.
%           The broadband pass is where whole-night cleaning does its damage: its cut
%           sits ~12x deeper than the stage-specific one on N3 epochs, against ~2.3x
%           for the slow wavelet bands and ~1x for the fast ones. It runs first, on the
%           full-band signal, where N3's variance IS the slow waves.
%
% See also gedai.defaultRuns, bidsfun_gedai.

arguments
    mode   (1,:) char {mustBeMember(mode,   {'auto-','auto','auto+'})} = 'auto'
    modeBB (1,:) char {mustBeMember(modeBB, {'auto-','auto','auto+'})} = 'auto-'
end

    N3 = -3; N2 = -2; REM = 0; Wake = 1;
    stages = [N3, N2, REM, Wake];

    runs = {struct( ...
        'GEDAIMode_dict',             dictionary(stages, repmat(string(mode),   1, 4)), ...
        'GEDAIModeBB_dict',           dictionary(stages, repmat(string(modeBB), 1, 4)), ...
        'GEDAILowCutOffFreq',         0.1, ...
        'GEDAIEpochSize',             12, ...
        'GEDAIBroadbandEpochSize',    10, ...
        'boost1',                     1, ...
        'boost2',                     1, ...
        'broadbandOnly',              false, ...
        'percentileThreshold',        98, ...
        'WeightKC',                   0, ...
        'BBMinThreshold',             -2, ...
        'computeSENSAI',              false, ...
        'GEDAIEnovaChannelThreshold', Inf, ...
        'ICAtype',                    'none')};
end
