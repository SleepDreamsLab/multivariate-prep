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
%   SavePath     : Full path + base filename; suffix + '.png' appended per figure.
%   EpochsToPlot : Epoch indices for overlay plots; default = first epoch of each stage.
%
%   Produces four figures via gedai.evalplots sub-functions:
%     psd_per_stage  — clean vs. raw PSD per sleep stage
%     psd_overview   — two-tile PSD coloured by stage
%     topo_band_power — topographic band power, raw vs. clean
%     epoch_overlay  — per-epoch signal overlay

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
    opts.EOGMontage   = [54 248 1 230]  % [L1 L2 R1 R2] EGI-256 bipolar EOG; [] = skip
    opts.TopoBandLims = []              % nBands x 2 — [lo hi] per band; NaN row = auto
    opts.refresh = false;
    opts.net = 'EGI256';
end

%%% --- Pathing ---
basepath = fileparts(opts.SavePath);

%%% --- Interpolate ---
tic;
EEGraw      = pop_interp(EEGraw, EEGraw.urchanlocs, 'spherical');
EEGclean    = pop_interp(EEGclean, EEGraw.urchanlocs, 'spherical');
toc

%%% --- Rename 10-20 channels ---
EEGclean    = chans1020(EEGclean, false, 'add_eog', 0, 'net', opts.net);
EEGraw      = chans1020(EEGraw, false, 'add_eog', 0, 'net', opts.net);

%%% --- Compute Welch power spectra ---
fprintf('gedai.eval_clean: computing Welch power (clean) ...\n');
[PwrClean, FreqsClean] = smartcache( ...
    @() run.run_pwelch(EEGclean, opts.EpochLength, ...
        opts.WelchWindow, opts.WelchOverlap), ...
            fullfile([opts.SavePath '_' 'PSDclean' '.mat']), ...
            opts.refresh , {'Power', 'Freqs'});

fprintf('gedai.eval_clean: computing Welch power (raw) ...\n');
[PwrRaw, FreqsRaw] = smartcache( ...
    @() run.run_pwelch(EEGraw, opts.EpochLength, ...
        opts.WelchWindow, opts.WelchOverlap), ...
            fullfile([opts.SavePath '_' 'PSDraw' '.mat']), ...
            opts.refresh , {'Power', 'Freqs'});

% [PwrClean, Freqs] = run.run_pwelch(EEGclean, opts.EpochLength, ...
%     opts.WelchWindow, opts.WelchOverlap);
% 
% [PwrRaw, ~] = run.run_pwelch(EEGraw, opts.EpochLength, ...
%     opts.WelchWindow, opts.WelchOverlap);

%%% --- Channel index used for PSD line plots ---
goodLabels = {EEGclean.chanlocs.labels};
FzIdx = find(strcmpi(goodLabels, 'Fz'), 1);
if isempty(FzIdx)
    warning('gedai.eval_clean: Fz not found; using first channel for PSD plots.');
    FzIdx = 1;
end

%%% --- Guard epoch-count mismatch from run_pwelch rounding ---
stageScoring = StageScoring(:)';
if size(PwrClean, 2) ~= numel(stageScoring)
    warning('Sleep scoring vector length does not match # epochs');
end

% %%% --- Fix chanlocs for topoplot ---
% EEGclean = qol.bids_fixchanlocs(EEGclean);

%%% --- Figures (no FOOOF) ---
try
gedai.evalplots.gedai_characteristics(EEGclean.etc.GEDAI, 'Srate', EEGclean.srate, 'SavePath', opts.SavePath);
end
try
gedai.evalplots.gedai_characteristics(EEGclean.etc.gedai, 'Srate', EEGclean.srate, 'SavePath', opts.SavePath);
end


evalplots.psd_per_stage(PwrClean, PwrRaw, FreqsClean, FreqsRaw, stageScoring, FzIdx, ...
    'FreqLim', opts.FreqLim, 'FreqScale', opts.FreqScale, 'SavePath', opts.SavePath);

evalplots.psd_per_stage_chans(PwrClean, PwrRaw, FreqsClean, FreqsRaw, stageScoring, ...
    'FreqLim', opts.FreqLim, 'FreqScale', opts.FreqScale, 'SavePath', opts.SavePath);

evalplots.psd_overview(PwrClean, PwrRaw, FreqsClean, FreqsRaw, stageScoring, FzIdx, ...
    'FreqLim', opts.FreqLim, 'FreqScale', opts.FreqScale, 'SavePath', opts.SavePath);

evalplots.topo_band_power(PwrClean, PwrRaw, FreqsClean, FreqsRaw, stageScoring, ...
    EEGraw.chanlocs, 'CLims', opts.TopoBandLims, 'SavePath', opts.SavePath);

evalplots.topo_band_stage(PwrClean, PwrRaw, FreqsClean, FreqsRaw, stageScoring, ...
    EEGclean.chanlocs, 'SavePath', opts.SavePath);

%%% --- Append EOG + EMG channels for epoch overlay ---
EEGclean    = chans1020(EEGclean, 0, 'add_eog', 1, 'net', opts.net);
EEGraw      = chans1020(EEGraw, 0, 'add_eog', 1, 'net', opts.net);

%%% --- Epoch overlay ---
evalplots.epoch_overlay(EEGraw, EEGclean, stageScoring(1:30/opts.EpochLength:end), ...
    'EpochLength', 30, 'EpochsToPlot', opts.EpochsToPlot, ...
    'SavePath', opts.SavePath);

%%% --- Pre-smooth power spectra for FOOOF (done once, reused across frequency ranges) ---
powerRawSmooth   = oscip.smooth_spectrum_median(PwrRaw(FzIdx, :, :),   FreqsRaw,   opts.FooofMedianSmooth);
powerRawSmooth   = oscip.smooth_spectrum(powerRawSmooth,   FreqsRaw,   opts.FooofMeanSmooth);
powerCleanSmooth = oscip.smooth_spectrum_median(PwrClean(FzIdx, :, :), FreqsClean, opts.FooofMedianSmooth);
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

    %%% --- Raw FOOOF ---
    fprintf('eval_clean: running FOOOF (raw, %s) ...\n', frLabel);
    [FooofRaw] = smartcache( ...
        @() run.run_fooof(powerRawSmooth, FreqsRaw, [], ...
            'FrequencyRange', fr, ...
            'MeanSmoothSpan', opts.FooofMeanSmooth, ...
            'MedianSmoothSpan', opts.FooofMedianSmooth, ...
            'PeakWidthLimits', opts.FooofPeakWidthLims, ...
            'AperiodicMode', opts.FooofAperiodicMode), ...
                fullfile([opts.SavePath '_FOOOFraw_' frTag '.mat']), ...
                false, {'FOOOF'});

    %%% --- Clean FOOOF ---
    fprintf('eval_clean: running FOOOF (clean, %s) ...\n', frLabel);
    [FooofClean] = smartcache( ...
        @() run.run_fooof(powerCleanSmooth, FreqsClean, [], ...
            'FrequencyRange', fr, ...
            'MeanSmoothSpan', opts.FooofMeanSmooth, ...
            'MedianSmoothSpan', opts.FooofMedianSmooth, ...
            'PeakWidthLimits', opts.FooofPeakWidthLims, ...
            'AperiodicMode', opts.FooofAperiodicMode), ...
                fullfile([opts.SavePath '_FooofClean_' frTag '.mat']), ...
                false, {'FOOOF'});       

% [PwrRaw, Freqs] = smartcache( ...
%     @() run.run_pwelch(EEGraw, opts.EpochLength, ...
%         opts.WelchWindow, opts.WelchOverlap), ...
%             fullfile([opts.SavePath '_' 'PSDraw' '.mat']), ...
%             false, {'Power', 'Freqs'});    

%     FooofClean = run.run_fooof(pwr, Freqs, [], ...
%         'FrequencyRange', fr, ...
%         'MeanSmoothSpan', opts.FooofMeanSmooth, ...
%         'MedianSmoothSpan', opts.FooofMedianSmooth, ...
%         'PeakWidthLimits', opts.FooofPeakWidthLims, ...
%         'AperiodicMode', opts.FooofAperiodicMode);

    %%% --- Save slopes ---
    slopesRaw{iRange}   = FooofRaw.Exponents;
    slopesClean{iRange} = FooofClean.Exponents;
    frLabels{iRange}    = frLabel;
end

%%% --- Figures (FOOOF) ---
saveSuffix = [opts.SavePath '_' frTag];
evalplots.timefreq(PwrClean, PwrRaw, FreqsClean, FreqsRaw, stageScoring, FzIdx, ...
    'ChanLabel', goodLabels{FzIdx}, ...
    'EpochLength', opts.EpochLength, 'FreqLim', opts.FreqLim, ...
    'ExponentsClean', FooofClean.Exponents, ...
    'ExponentsRaw',   FooofRaw.Exponents, ...
    'FooofLabel', frLabel, ...
    'SavePath', saveSuffix);

evalplots.exponent_by_stage(slopesClean, slopesRaw, stageScoring, ...
    'FooofLabel', frLabels, 'SavePath', opts.SavePath);

evalplots.slopes_timecourse(slopesRaw, slopesClean, fooofRanges, stageScoring, ...
    'EpochLength', opts.EpochLength, 'SavePath', opts.SavePath);
end
