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
    opts.egi256_to_1020 = struct( ...
    'Fp1',  37,  'Fpz',  26,  'Fp2',  18, 'F7',   47,  'F3',   36,  'Fz',   21,  'F4', 224, 'F8',   2, ...
    'T7',   69,  'C3',   59, 'C4', 183, 'T8', 202, 'P7',   96,  'P3',   87,  'Pz',  101,  'P4', 153, 'P8', 170, ...
    'O1',  116,  'O2',  150 ...
    );
    opts.refresh = false;
end

%%% --- Pathing ---
basepath = fileparts(opts.SavePath);

%%% --- Select non-EOG channels --- (for the old datasets with EOG)
cleanLabels = {EEGclean.chanlocs.labels};
eogMask     = cellfun(@(s) ~isempty(regexpi(s, 'EOG')), cleanLabels);
goodChans   = find(~eogMask);
if isempty(goodChans)
    warning('gedai.eval_clean: no non-EOG channels found; using all channels.');
    goodChans = 1:numel(cleanLabels);
end

EEGclean = pop_select(EEGclean, 'channel', goodChans);
EEGraw = pop_select(EEGraw, 'channel', goodChans);


%%% --- Rename 10-20 channels ---
egi256indices = struct2array(opts.egi256_to_1020);
labels1020 = fieldnames(opts.egi256_to_1020);
for iCh = 1:numel(egi256indices)
    EEGclean.chanlocs(egi256indices(iCh)).labels = labels1020{iCh};
    EEGraw.chanlocs(egi256indices(iCh)).labels = labels1020{iCh};
end


%%% --- Compute Welch power spectra ---
fprintf('gedai.eval_clean: computing Welch power (clean) ...\n');
[PwrClean, Freqs] = smartcache( ...
    @() run.run_pwelch(EEGclean, opts.EpochLength, ...
        opts.WelchWindow, opts.WelchOverlap), ...
            fullfile([opts.SavePath '_' 'PSDclean' '.mat']), ...
            opts.refresh , {'Power', 'Freqs'});

fprintf('gedai.eval_clean: computing Welch power (raw) ...\n');
[PwrRaw, Freqs] = smartcache( ...
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


gedai.evalplots.psd_per_stage(PwrClean, PwrRaw, Freqs, stageScoring, FzIdx, ...
    'FreqLim', opts.FreqLim, 'FreqScale', opts.FreqScale, 'SavePath', opts.SavePath);

gedai.evalplots.psd_overview(PwrClean, PwrRaw, Freqs, stageScoring, FzIdx, ...
    'FreqLim', opts.FreqLim, 'FreqScale', opts.FreqScale, 'SavePath', opts.SavePath);

gedai.evalplots.topo_band_power(PwrClean, PwrRaw, Freqs, stageScoring, ...
    EEGclean.chanlocs, 'CLims', opts.TopoBandLims, 'SavePath', opts.SavePath);

gedai.evalplots.topo_band_stage(PwrClean, PwrRaw, Freqs, stageScoring, ...
    EEGclean.chanlocs, 'SavePath', opts.SavePath);

%%% --- Append EOG + EMG channels for epoch overlay ---
m = opts.EOGMontage;
if numel(m) == 4 && max(m) <= EEGraw.nbchan
    EEGraw.data(end+1,:) = EEGraw.data(m(1),:) - EEGraw.data(m(3),:);
    EEGraw.chanlocs(end+1).labels = 'EOG1';
    EEGraw.data(end+1,:) = EEGraw.data(m(4),:) - EEGraw.data(m(2),:);
    EEGraw.chanlocs(end+1).labels = 'EOG2';
    EEGclean.data(end+1,:) = EEGclean.data(m(1),:) - EEGclean.data(m(3),:);
    EEGclean.chanlocs(end+1).labels = 'EOG1';
    EEGclean.data(end+1,:) = EEGclean.data(m(4),:) - EEGclean.data(m(2),:);
    EEGclean.chanlocs(end+1).labels = 'EOG2';
end
EEGraw.nbchan = size(EEGraw.data, 1);
EEGclean.nbchan = size(EEGclean.data, 1);

%%% --- Epoch overlay ---
gedai.evalplots.epoch_overlay(EEGraw, EEGclean, stageScoring(1:30/opts.EpochLength:end), ...
    'EpochLength', 30, 'EpochsToPlot', opts.EpochsToPlot, ...
    'SavePath', opts.SavePath);

%%% --- Pre-smooth power spectra for FOOOF (done once, reused across frequency ranges) ---
powerRawSmooth   = oscip.smooth_spectrum_median(PwrRaw(FzIdx, :, :),   Freqs, opts.FooofMedianSmooth);
powerRawSmooth   = oscip.smooth_spectrum(powerRawSmooth,   Freqs, opts.FooofMeanSmooth);
powerCleanSmooth = oscip.smooth_spectrum_median(PwrClean(FzIdx, :, :), Freqs, opts.FooofMedianSmooth);
powerCleanSmooth = oscip.smooth_spectrum(powerCleanSmooth, Freqs, opts.FooofMeanSmooth);

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
    fprintf('gedai.eval_clean: running FOOOF (raw, %s) ...\n', frLabel);
    [FooofRaw] = smartcache( ...
        @() run.run_fooof(powerRawSmooth, Freqs, [], ...
            'FrequencyRange', fr, ...
            'MeanSmoothSpan', opts.FooofMeanSmooth, ...
            'MedianSmoothSpan', opts.FooofMedianSmooth, ...
            'PeakWidthLimits', opts.FooofPeakWidthLims, ...
            'AperiodicMode', opts.FooofAperiodicMode), ...
                fullfile([opts.SavePath '_FOOOFraw_' frTag '.mat']), ...
                false, {'FOOOF'});

    %%% --- Clean FOOOF ---
    fprintf('gedai.eval_clean: running FOOOF (clean, %s) ...\n', frLabel);
    [FooofClean] = smartcache( ...
        @() run.run_fooof(powerCleanSmooth, Freqs, [], ...
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
gedai.evalplots.timefreq(PwrClean, PwrRaw, Freqs, stageScoring, FzIdx, ...
    'ChanLabel', goodLabels{FzIdx}, ...
    'EpochLength', opts.EpochLength, 'FreqLim', opts.FreqLim, ...
    'ExponentsClean', FooofClean.Exponents, ...
    'ExponentsRaw',   FooofRaw.Exponents, ...
    'FooofLabel', frLabel, ...
    'SavePath', saveSuffix);

gedai.evalplots.exponent_by_stage(slopesClean, slopesRaw, stageScoring, ...
    'FooofLabel', frLabels, 'SavePath', opts.SavePath);

gedai.evalplots.slopes_timecourse(slopesRaw, slopesClean, fooofRanges, stageScoring, ...
    'EpochLength', opts.EpochLength, 'SavePath', opts.SavePath);
end
