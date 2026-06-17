function KCpower(EEG, opts)
% GEDAI.EVALPLOTS.KCPOWER  Topographic band power during K-complex windows.
%
%   gedai.evalplots.KCpower(EEG, KCperiods)
%   gedai.evalplots.KCpower(EEG, KCperiods, Name, Value)
%
%   Inputs
%   ------
%   EEG       : EEGLAB struct (continuous).
%   KCperiods : Nx2 matrix of [start stop] times in seconds.
%
%   Optional name-value pairs
%   -------------------------
%   CLims    : nBands x 2 matrix of [lo hi] per band; NaN/missing rows → auto.
%   SavePath : base path for saving; suffix + '.png' appended.

arguments
    EEG
    opts.EpochLength  = 4
    opts.WelchWindow  = 4
    opts.WelchOverlap = 2
    opts.CLims        = []   % nBands×2 [lo hi]; NaN/missing rows → auto
    opts.FocusChans   = {'Fz', 'Fz', 'Oz', 'Pz', 'Fz'}  % 1×5 cell; channel whose mean band power sets the upper clim
end

bands = {
    'SWA',   [0.5  4];
    'Theta', [4    8];
    'Alpha', [8   12];
    'Sigma', [12  16];
    'Beta',  [16  30];
};

nBands = size(bands, 1);
nFrequencies = floor(opts.WelchWindow*EEG.srate/2) + 1;

% Pad EEG
epochPnts   = EEG.srate*opts.EpochLength;
nEpochs     = ceil(EEG.pnts / epochPnts);
EEG.data(:, EEG.pnts+1:nEpochs*epochPnts) = 0;

% Epoch data
EEG.event(strcmpi({EEG.event.type}, 'boundary')) = [];
EEG = eeg_checkset(EEG);
EEG = eeg_regepochs(EEG, 'recurrence', opts.EpochLength, 'limits', [0 opts.EpochLength]);

% Compute power
PowerKC = nan(EEG.nbchan, nFrequencies, EEG.trials);
for epo = 1:EEG.trials
    [power, freqs] = pwelch(EEG.data(:,:,epo)', hanning(opts.WelchWindow*EEG.srate), ...
        opts.WelchOverlap*EEG.srate, opts.WelchWindow*EEG.srate, EEG.srate);
    PowerKC(:, :, epo) = power';
end

% Color limits: user-supplied where provided; otherwise lower=data-driven,
% upper=mean band power of the reference channel (opts.FocusChans{iBand}).
clims = nan(nBands, 2);
for iBand = 1:nBands
    if size(opts.CLims, 1) >= iBand && all(isfinite(opts.CLims(iBand, :)))
        clims(iBand, :) = opts.CLims(iBand, :);
    else
        freqInds = freqs >= bands{iBand,2}(1) & freqs <= bands{iBand,2}(2);
        vals     = mean(PowerKC(:, freqInds, :), [2 3]);
        finVals  = vals(isfinite(vals));
        if isempty(finVals), continue; end

        loLim = floor(min(finVals) * 11) / 10;

        refChan = opts.FocusChans{iBand};
        chanIdx = find(strcmpi({EEG.chanlocs.labels}, refChan), 1);
        if ~isempty(chanIdx)
            refVal = mean(PowerKC(chanIdx, freqInds, :), 'all');
            hiLim  = ceil(refVal * 11) / 10;
        else
            hiLim  = ceil(max(finVals) * 11) / 10;
        end

        clims(iBand, :) = [loLim, hiLim];
    end
end

%% Open figure
fig = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [2 2 11 25], ...
    'Name', 'KC power topographies');
tiledlayout(nBands, 2, 'TileSpacing', 'tight', 'Padding', 'tight');


for iBand = 1:nBands
    bandName  = bands{iBand, 1};
    bandMin   = bands{iBand, 2}(1);
    bandMax   = bands{iBand, 2}(2);

    % Power
    freqInds    = freqs >= bandMin & freqs <= bandMax;
    bandmean    = mean(PowerKC(:,freqInds,:), [2, 3]);
    bandmedian  = median(PowerKC(:,freqInds,:), [2, 3]);
    bandmax     = max(PowerKC(:,freqInds,:), [], [2, 3]);

    bandTitle = sprintf('%s (%.1f-%d)', bandName, bandMin, bandMax);
    cl        = clims(iBand, :);
    hasCL     = all(isfinite(cl));
    maplimArg = {};
    if hasCL, maplimArg = {'maplimits', cl}; end

    % Topoplot
    nexttile()
    topoplot(bandmean, EEG.chanlocs, maplimArg{:}, 'electrodes', 'on', 'numcontour', 0);
    cb = colorbar();
    if hasCL, cb.Limits = cl; cb.Ticks = linspace(cl(1), cl(2), 5); end
    title(sprintf('%s\n(Mean)', bandTitle), 'FontSize', 10, 'FontWeight', 'bold');

    nexttile()
    topoplot(bandmedian, EEG.chanlocs, maplimArg{:}, 'electrodes', 'on', 'numcontour', 0);
    cb = colorbar();
    if hasCL, cb.Limits = cl; cb.Ticks = linspace(cl(1), cl(2), 5); end
    title(sprintf('%s\n(Median)', bandTitle), 'FontSize', 10, 'FontWeight', 'bold');

%     nexttile()
%     topoplot(bandmax, EEG.chanlocs, maplimArg{:}, 'electrodes', 'on', 'numcontour', 0);
%     cb = colorbar();
%     if hasCL, cb.Limits = cl; cb.Ticks = linspace(cl(1), cl(2), 5); end
%     title(sprintf('%s\n(Max)', bandTitle), 'FontSize', 10, 'FontWeight', 'bold');

end

colormap(custom_cmap())
set(gcf, 'color', 'w')
end
