function topo_band_power(PwrClean, PwrRaw, FreqsClean, FreqsRaw, StageScoring, Chanlocs, opts)
% GEDAI.EVAL.TOPO_BAND_POWER  Topographic band power: raw vs. clean, per sleep stage.
%
%   gedai.eval.topo_band_power(PwrClean, PwrRaw, Freqs, StageScoring, Chanlocs)
%   gedai.eval.topo_band_power(..., Name, Value)
%
%   Inputs
%   ------
%   PwrClean     : chans x epochs x freqs.
%   PwrRaw       : chans x epochs x freqs.
%   Freqs        : frequency vector (Hz).
%   StageScoring : per-epoch sleep-stage codes.
%   Chanlocs     : EEGLAB chanlocs struct (already position-corrected).
%
%   Optional name-value pairs
%   -------------------------
%   Bands     : nBands x 3 cell — {label, [fMin fMax], stageCode}.
%               Default covers SWA/Sigma/Alpha/Theta for N2, N3, Wake, REM.
%               Color limits are derived from the clean data (min/max across
%               channels, floor/ceil to 1 decimal) and applied to both raw and clean.
%   SavePath  : base path for saving; suffix + '.png' appended.
%   CLims     : nBands x 2 matrix of [lo hi] per band. NaN (or missing) rows
%               fall back to data-driven limits (default: all auto).

arguments
    PwrClean
    PwrRaw
    FreqsClean   {mustBeVector}
    FreqsRaw     {mustBeVector}
    StageScoring {mustBeVector}
    Chanlocs
    opts.Bands    = {}
    opts.CLims    = []
    opts.SavePath = ''
end

if isempty(opts.Bands)
    bands = {
        'SWA 0.5-4 Hz',   [0.5  4],  -3;
        'SWA 0.5-4 Hz',   [0.5  4],  -2;
        'Sigma 12-16 Hz', [12  16],  -2;
        'Alpha 8-12 Hz',  [8   12],   1;
        'Theta 4-8 Hz',   [4    8],   0;
        };
else
    bands = opts.Bands;
end

nBands = size(bands, 1);
nChans = numel(Chanlocs);

rawBands   = nan(nChans, nBands);
cleanBands = nan(nChans, nBands);

for iBand = 1:nBands
    freqRange = bands{iBand, 2};
    stageCode = bands{iBand, 3};

    if stageCode == 0
        stageMask = StageScoring == 0 | StageScoring == 5;
    else
        stageMask = StageScoring == stageCode;
    end
    fBandMaskRaw   = FreqsRaw   >= freqRange(1) & FreqsRaw   <= freqRange(2);
    fBandMaskClean = FreqsClean >= freqRange(1) & FreqsClean <= freqRange(2);

    if ~any(stageMask) || (~any(fBandMaskRaw) && ~any(fBandMaskClean)), continue; end

    if any(fBandMaskRaw)
        rawBands(:, iBand)   = squeeze(mean(mean( ...
            log10(PwrRaw(:,   stageMask, fBandMaskRaw)   + eps), 3), 2));
    end
    if any(fBandMaskClean)
        cleanBands(:, iBand) = squeeze(mean(mean( ...
            log10(PwrClean(:, stageMask, fBandMaskClean) + eps), 3), 2));
    end
end

% Color limits: user-supplied where provided, data-driven (clean min/max) otherwise
clims = nan(nBands, 2);
for iBand = 1:nBands
    if size(opts.CLims, 1) >= iBand && all(isfinite(opts.CLims(iBand, :)))
        clims(iBand, :) = opts.CLims(iBand, :);
    else
        vals = cleanBands(:, iBand);
        vals = vals(isfinite(vals));
        if ~isempty(vals)
            clims(iBand, :) = [floor(min(vals)*10)/10, ceil(max(vals)*10)/10];
        end
    end
end

fig = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [2 2 15 6*nBands], ...
    'Name', 'GEDAI evaluation — band power topographies');
tiledlayout(nBands, 2, 'TileSpacing', 'tight', 'Padding', 'tight');

axRaw   = gobjects(1, nBands);
axClean = gobjects(1, nBands);
cbRaw   = gobjects(1, nBands);
cbClean = gobjects(1, nBands);

% topoplot tries to set axes Position, which raises a harmless warning inside TiledChartLayout
warnState = warning('off', 'all');

for iBand = 1:nBands
    bandName  = bands{iBand, 1};
    stageCode = bands{iBand, 3};
    cl        = clims(iBand, :);

    if any(isnan(rawBands(:, iBand))), continue; end
    if any(~isfinite(cl)), continue; end

    bandTitle = sprintf('%s  (%s)', bandName, stage_name(stageCode));

    % Left column: raw
    axRaw(iBand) = nexttile((iBand-1)*2 + 1);
    topoplot(rawBands(:, iBand), Chanlocs, 'maplimits', cl, 'electrodes', 'on', 'numcontour', 0, 'conv', 'on');
    set(axRaw(iBand), 'CLim', cl);
    cbRaw(iBand) = colorbar;
    cbRaw(iBand).FontSize = 8;
    cbRaw(iBand).Limits = cl;
    cbRaw(iBand).Ticks  = linspace(cl(1), cl(2), 5);
    if iBand == 1
        title(axRaw(iBand), sprintf('Raw\n%s', bandTitle), 'FontSize', 10, 'FontWeight', 'bold');
    else
        title(axRaw(iBand), bandTitle, 'FontSize', 10, 'FontWeight', 'bold');
    end

    % Right column: clean
    axClean(iBand) = nexttile((iBand-1)*2 + 2);
    topoplot(cleanBands(:, iBand), Chanlocs, 'maplimits', cl, 'electrodes', 'on', 'numcontour', 0, 'conv', 'on');
    set(axClean(iBand), 'CLim', cl);
    cbClean(iBand) = colorbar;
    cbClean(iBand).FontSize = 8;
    cbClean(iBand).Limits = cl;
    cbClean(iBand).Ticks  = linspace(cl(1), cl(2), 5);
    if iBand == nBands, cbClean(iBand).Label.String = 'log_{10}(\muV^2/Hz)'; end
    if iBand == 1
        title(axClean(iBand), sprintf('Clean\n%s', bandTitle), 'FontSize', 10, 'FontWeight', 'bold');
    else
        title(axClean(iBand), bandTitle, 'FontSize', 10, 'FontWeight', 'bold');
    end
end

warning(warnState);

% Re-apply per-axes colormap and CLim after topoplot sets the figure-level colormap
for iBand = 1:nBands
    cl = clims(iBand, :);
    if any(~isfinite(cl)), continue; end
    if isgraphics(axRaw(iBand))
        colormap(axRaw(iBand), custom_cmap());
        set(axRaw(iBand), 'CLim', cl);
        if isgraphics(cbRaw(iBand))
            cbRaw(iBand).Limits = cl;
            cbRaw(iBand).Ticks  = linspace(cl(1), cl(2), 5);
        end
    end
    if isgraphics(axClean(iBand))
        colormap(axClean(iBand), custom_cmap());
        set(axClean(iBand), 'CLim', cl);
        if isgraphics(cbClean(iBand))
            cbClean(iBand).Limits = cl;
            cbClean(iBand).Ticks  = linspace(cl(1), cl(2), 5);
        end
    end
end
set(gcf, 'Color', 'w');

save_fig(fig, opts.SavePath, 'topo_band_power');
end
