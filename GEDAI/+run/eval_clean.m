function eval_clean(EEGraw, EEGclean, StageScoring, opts)
% GEDAI.EVAL_CLEAN  Evaluate GEDAI cleaning quality: spectral + topographic + epoch comparison.
%
%   gedai.eval_clean(EEGraw, EEGclean, StageScoring)
%   gedai.eval_clean(EEGraw, EEGclean, StageScoring, Name, Value)
%
%   Inputs
%   ------
%   EEGraw     : Raw continuous EEG fed into GEDAI (EEGraw output of gedai.stage_specific).
%   EEGclean     : GEDAI-cleaned EEG (EEGclean output of gedai.stage_specific).
%   StageScoring : Per-epoch sleep-stage digits. Pass Scoring(epochIdx).
%
%   Optional name-value pairs
%   -------------------------
%   EpochLength  : Epoch duration in seconds (default 30).
%   WelchWindow  : Welch window length in seconds (default 4).
%   WelchOverlap : Welch overlap fraction (default 0.5).
%   FreqLim      : [fMin fMax] for power plots (default [0 srate/2]).
%   FreqScale    : 'linear', 'log', or 'both' (default).
%   SavePath     : Full path + base filename; suffix + '.png'/'.mat' appended per figure/cache file.
%   EpochsToPlot : Epoch indices for overlay plots; default = first epoch of each stage.
%   cachepower   : If true, cache run.run_pwelch output via smartcache; if false, run uncached (default false).
%   cachefooof   : If true, cache run.run_fooof output via smartcache; if false, run uncached (default false).
%
%   Plot toggles (all default true)
%   --------------------------------
%   PlotCharacteristics   : GEDAI cleaning characteristics summary.
%   PlotPsdPerStage       : Clean vs. raw PSD per sleep stage (single channel).
%   PlotPsdPerStageChans  : Clean vs. raw PSD per sleep stage, all channels.
%   PlotPsdOverview       : Two-tile PSD coloured by stage.
%   PlotTopoBandPower     : Topographic band power, raw vs. clean.
%   PlotTopoBandStage     : Topographic band power by stage.
%   PlotEpochOverlay      : Per-epoch signal overlay.
%   PlotTimefreq          : Time-frequency PSD with FOOOF exponent overlay.
%   PlotExponentByStage   : FOOOF aperiodic exponent by sleep stage.
%   PlotSlopesTimecourse  : FOOOF aperiodic exponent time course.
%
%   FOOOF is only computed if at least one of PlotTimefreq, PlotExponentByStage,
%   PlotSlopesTimecourse is true.

arguments
    EEGraw
    EEGclean
    StageScoring  {mustBeVector}
    opts.EpochLength  = 30
    opts.WelchWindow  = 4
    opts.WelchOverlap = 0.5
    opts.FreqLim      = [0 EEGraw.srate/2]
    opts.FreqScale    = 'both'
    opts.SavePath     = ''
    opts.EpochsToPlot = []
    opts.FooofMeanSmooth    = 2       % Hz — mean smoothing window
    opts.FooofMedianSmooth  = 3       % Hz — median smoothing window
    opts.FooofPeakWidthLims = [0.5 20]
    opts.FooofAperiodicMode = 'fixed' % 'fixed' or 'knee'
    opts.TopoBandLims = []              % nBands x 2 — [lo hi] per band; NaN row = auto
    opts.refresh = false;
    opts.net = 'EGI256';
    opts.cachepower (1,1) logical = false
    opts.cachefooof (1,1) logical = false
    opts.PlotCharacteristics  (1,1) logical = true
    opts.PlotPsdPerStage      (1,1) logical = true
    opts.PlotPsdPerStageChans (1,1) logical = true
    opts.PlotPsdOverview      (1,1) logical = true
    opts.PlotTopoBandPower    (1,1) logical = true
    opts.PlotTopoBandStage    (1,1) logical = true
    opts.PlotEpochOverlay     (1,1) logical = true
    opts.PlotTimefreq         (1,1) logical = true
    opts.PlotExponentByStage  (1,1) logical = true
    opts.PlotSlopesTimecourse (1,1) logical = true
end

%%% --- Interpolate ---
%%% Both datasets go back to the full montage before anything below runs: bad channels
%%% were dropped from the cleaned data upstream, and every plot function here compares
%%% the two through a single channel index or a single chanlocs, so they have to share a
%%% montage. Only call pop_interp when something is actually missing - it returns early
%%% otherwise, but the call still drags eeg_checkset behind it. On a full night this is
%%% the single most memory-hungry step of the evaluation, so patches/eeg_interp.m
%%% replaces EEGLAB's spherical interpolation with a chunked, class-preserving version
%%% (same output, ~4x less peak RAM) - keep 'patches' ahead of EEGLAB on the path.
if isfield(EEGraw, 'urchanlocs') && ~isempty(EEGraw.urchanlocs)
    urchanlocs = EEGraw.urchanlocs;
    if numel(EEGraw.chanlocs) < numel(urchanlocs)
        fprintf('gedai.eval_clean: interpolating %d channel(s) in raw ...\n', ...
            numel(urchanlocs) - numel(EEGraw.chanlocs));
        EEGraw = pop_interp(EEGraw, urchanlocs, 'spherical');
    end
    
    urchanlocs = EEGclean.urchanlocs;
    if numel(EEGclean.chanlocs) < numel(urchanlocs)
        fprintf('gedai.eval_clean: interpolating %d channel(s) in clean ...\n', ...
            numel(urchanlocs) - numel(EEGclean.chanlocs));
        EEGclean = pop_interp(EEGclean, urchanlocs, 'spherical');
    end
end

%%% --- Rename 10-20 channels ---
EEGclean    = chans1020(EEGclean, false, 'add_eog', 0, 'net', opts.net, 'chanprefix', 'E');
EEGraw      = chans1020(EEGraw, false, 'add_eog', 0, 'net', opts.net, 'chanprefix', 'E');

%%% --- Compute Welch power spectra (cached only if opts.cachepower) ---
%%% Power matrices are chans x freqs x epochs.
fprintf('gedai.eval_clean: computing Welch power (clean) ...\n');
if opts.cachepower
    [PwrClean, FreqsClean] = smartcache( ...
        @() run.pwelch_fast(EEGclean, opts.EpochLength, ...
            opts.WelchWindow, opts.WelchOverlap), ...
                fullfile([opts.SavePath '_' 'PSDclean' '.mat']), ...
                opts.refresh , {'Power', 'Freqs'});
else
    [PwrClean, FreqsClean] = run.pwelch_fast(EEGclean, opts.EpochLength, ...
        opts.WelchWindow, opts.WelchOverlap);
end

fprintf('gedai.eval_clean: computing Welch power (raw) ...\n');
if opts.cachepower
    [PwrRaw, FreqsRaw] = smartcache( ...
        @() run.pwelch_fast(EEGraw, opts.EpochLength, ...
            opts.WelchWindow, opts.WelchOverlap), ...
                fullfile([opts.SavePath '_' 'PSDraw' '.mat']), ...
                opts.refresh , {'Power', 'Freqs'});
else
    [PwrRaw, FreqsRaw] = run.pwelch_fast(EEGraw, opts.EpochLength, ...
        opts.WelchWindow, opts.WelchOverlap);
end

%%% --- Channel index used for PSD line plots ---
goodLabels = {EEGclean.chanlocs.labels};
FzIdx = find(strcmpi(goodLabels, 'Fz'), 1);
if isempty(FzIdx)
    warning('gedai.eval_clean: Fz not found; using first channel for PSD plots.');
    FzIdx = 1;
end

%%% --- Guard epoch-count mismatch from pwelch_fast rounding ---
stageScoring = StageScoring(:)';
if size(PwrClean, 3) ~= numel(stageScoring)
    warning('Sleep scoring vector length does not match # epochs');
end

%%% --- Figures ---
%%% The order of the blocks below is the order the PNGs land on disk, and
%%% gedai.lastEvalFigure mirrors it so callers can tell a finished recording from one
%%% that merely started. Adding, moving or renaming a plot here means updating that
%%% table too, or resume checks will key on the wrong file.

%%% --- Figures (no FOOOF) ---
if opts.PlotCharacteristics
    % etc field capitalization varies by GEDAI version; try both, swallow the one that doesn't exist
    try
        gedai.evalplots.gedai_characteristics(EEGclean.etc.GEDAI, 'Srate', EEGclean.srate, 'SavePath', opts.SavePath);
    end
    try
        gedai.evalplots.gedai_characteristics(EEGclean.etc.gedai, 'Srate', EEGclean.srate, 'SavePath', opts.SavePath);
    end
end

if opts.PlotPsdPerStage
    evalplots.psd_per_stage(PwrClean, PwrRaw, FreqsClean, FreqsRaw, stageScoring, FzIdx, ...
        'FreqLim', opts.FreqLim, 'FreqScale', opts.FreqScale, 'SavePath', opts.SavePath);
end

if opts.PlotPsdPerStageChans
    evalplots.psd_per_stage_chans(PwrClean, PwrRaw, FreqsClean, FreqsRaw, stageScoring, ...
        'FreqLim', opts.FreqLim, 'FreqScale', opts.FreqScale, 'SavePath', opts.SavePath);
end

if opts.PlotPsdOverview
    evalplots.psd_overview(PwrClean, PwrRaw, FreqsClean, FreqsRaw, stageScoring, FzIdx, ...
        'FreqLim', opts.FreqLim, 'FreqScale', opts.FreqScale, 'SavePath', opts.SavePath);
end

if opts.PlotTopoBandPower
    evalplots.topo_band_power(PwrClean, PwrRaw, FreqsClean, FreqsRaw, stageScoring, ...
        EEGraw.chanlocs, 'CLims', opts.TopoBandLims, 'SavePath', opts.SavePath);
end

if opts.PlotTopoBandStage
    evalplots.topo_band_stage(PwrClean, PwrRaw, FreqsClean, FreqsRaw, stageScoring, ...
        EEGclean.chanlocs, 'SavePath', opts.SavePath);
end

if opts.PlotEpochOverlay
    %%% --- Append EOG + EMG channels for epoch overlay ---
    EEGclean    = chans1020(EEGclean, 0, 'add_eog', 1, 'net', opts.net);
    EEGraw      = chans1020(EEGraw, 0, 'add_eog', 1, 'net', opts.net);

    evalplots.epoch_overlay(EEGraw, EEGclean, stageScoring(1:30/opts.EpochLength:end), ...
        'EpochLength', 30, 'EpochsToPlot', opts.EpochsToPlot, ...
        'SavePath', opts.SavePath);
end

%%% --- FOOOF: only run if a FOOOF-dependent plot is requested ---
needFooof = opts.PlotTimefreq || opts.PlotExponentByStage || opts.PlotSlopesTimecourse;
if needFooof
    %%% --- Pre-smooth power spectra for FOOOF (done once, reused across frequency ranges) ---
    %%% oscip/FOOOF expect chans x epochs x freqs, so permute out of the
    %%% chans x freqs x epochs layout used everywhere else here.
    pwrRawFz   = permute(PwrRaw(FzIdx, :, :),   [1 3 2]);
    pwrCleanFz = permute(PwrClean(FzIdx, :, :), [1 3 2]);

    powerRawSmooth   = oscip.smooth_spectrum_median(pwrRawFz,   FreqsRaw,   opts.FooofMedianSmooth);
    powerRawSmooth   = oscip.smooth_spectrum(powerRawSmooth,   FreqsRaw,   opts.FooofMeanSmooth);
    powerCleanSmooth = oscip.smooth_spectrum_median(pwrCleanFz, FreqsClean, opts.FooofMedianSmooth);
    powerCleanSmooth = oscip.smooth_spectrum(powerCleanSmooth, FreqsClean, opts.FooofMeanSmooth);

    %%% --- FOOOF loop over frequency ranges ---
    fooofRanges = {[2 30], [30 45], [2 45]};
    slopesRaw   = cell(1, numel(fooofRanges));
    slopesClean = cell(1, numel(fooofRanges));
    frLabels    = cell(1, numel(fooofRanges));

    for iRange = 1:numel(fooofRanges)
        fr      = fooofRanges{iRange};
        frLabel = sprintf('%d–%d Hz', fr(1), fr(2));
        frTag   = sprintf('fooof%d-%d', fr(1), fr(2));

        fprintf('eval_clean: running FOOOF (raw, %s) ...\n', frLabel);
        fprintf('eval_clean: running FOOOF (clean, %s) ...\n', frLabel);
        if opts.cachefooof
            FooofRaw = smartcache( ...
                @() run.run_fooof(powerRawSmooth, FreqsRaw, [], ...
                    'FrequencyRange', fr, ...
                    'MeanSmoothSpan', opts.FooofMeanSmooth, ...
                    'MedianSmoothSpan', opts.FooofMedianSmooth, ...
                    'PeakWidthLimits', opts.FooofPeakWidthLims, ...
                    'AperiodicMode', opts.FooofAperiodicMode), ...
                        fullfile([opts.SavePath '_FOOOFraw_' frTag '.mat']), ...
                        false, {'FOOOF'});

            FooofClean = smartcache( ...
                @() run.run_fooof(powerCleanSmooth, FreqsClean, [], ...
                    'FrequencyRange', fr, ...
                    'MeanSmoothSpan', opts.FooofMeanSmooth, ...
                    'MedianSmoothSpan', opts.FooofMedianSmooth, ...
                    'PeakWidthLimits', opts.FooofPeakWidthLims, ...
                    'AperiodicMode', opts.FooofAperiodicMode), ...
                        fullfile([opts.SavePath '_FooofClean_' frTag '.mat']), ...
                        false, {'FOOOF'});
        else
            FooofRaw = run.run_fooof(powerRawSmooth, FreqsRaw, [], ...
                'FrequencyRange', fr, ...
                'MeanSmoothSpan', opts.FooofMeanSmooth, ...
                'MedianSmoothSpan', opts.FooofMedianSmooth, ...
                'PeakWidthLimits', opts.FooofPeakWidthLims, ...
                'AperiodicMode', opts.FooofAperiodicMode);

            FooofClean = run.run_fooof(powerCleanSmooth, FreqsClean, [], ...
                'FrequencyRange', fr, ...
                'MeanSmoothSpan', opts.FooofMeanSmooth, ...
                'MedianSmoothSpan', opts.FooofMedianSmooth, ...
                'PeakWidthLimits', opts.FooofPeakWidthLims, ...
                'AperiodicMode', opts.FooofAperiodicMode);
        end

        slopesRaw{iRange}   = FooofRaw.Exponents;
        slopesClean{iRange} = FooofClean.Exponents;
        frLabels{iRange}    = frLabel;
    end

    %%% --- Figures (FOOOF) ---
    if opts.PlotTimefreq
        saveSuffix = [opts.SavePath '_' frTag]; % last range in the loop (broadest)
        evalplots.timefreq(PwrClean, PwrRaw, FreqsClean, FreqsRaw, stageScoring, FzIdx, ...
            'ChanLabel', goodLabels{FzIdx}, ...
            'EpochLength', opts.EpochLength, 'FreqLim', opts.FreqLim, ...
            'ExponentsClean', FooofClean.Exponents, ...
            'ExponentsRaw',   FooofRaw.Exponents, ...
            'FooofLabel', frLabel, ...
            'SavePath', saveSuffix);
    end

    if opts.PlotExponentByStage
        evalplots.exponent_by_stage(slopesClean, slopesRaw, stageScoring, ...
            'FooofLabel', frLabels, 'SavePath', opts.SavePath);
    end

    if opts.PlotSlopesTimecourse
        evalplots.slopes_timecourse(slopesRaw, slopesClean, fooofRanges, stageScoring, ...
            'EpochLength', opts.EpochLength, 'SavePath', opts.SavePath);
    end
end
end
