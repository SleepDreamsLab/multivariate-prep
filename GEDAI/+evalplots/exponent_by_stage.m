function exponent_by_stage(ExponentsClean, ExponentsRaw, StageScoring, opts)
% GEDAI.EVALPLOTS.EXPONENT_BY_STAGE  Aperiodic exponent per sleep stage: raw vs. clean.
%
%   gedai.evalplots.exponent_by_stage(ExponentsClean, ExponentsRaw, StageScoring)
%   gedai.evalplots.exponent_by_stage(..., Name, Value)
%
%   For each sleep stage: a split half-violin (raw on left, clean on right)
%   with individual epoch values connected by lines showing the per-epoch change.
%   Accepts either a single exponent vector or a cell of vectors (one per
%   frequency range); in the latter case a tiled row is drawn per range.
%
%   Inputs
%   ------
%   ExponentsClean : 1 × nEpochs vector, or 1 × nRanges cell of such vectors.
%   ExponentsRaw   : same format as ExponentsClean.
%   StageScoring   : per-epoch sleep-stage codes.
%
%   Optional name-value pairs
%   -------------------------
%   FooofLabel : string or cell of strings — one label per frequency range.
%   SavePath   : base path for saving; suffix + '.png' appended.

arguments
    ExponentsClean
    ExponentsRaw
    StageScoring {mustBeVector}
    opts.FooofLabel = {}
    opts.SavePath   = ''
end

% Normalise to cells so the rest of the code is uniform
if ~iscell(ExponentsClean), ExponentsClean = {ExponentsClean}; end
if ~iscell(ExponentsRaw),   ExponentsRaw   = {ExponentsRaw};   end
if ischar(opts.FooofLabel) || isstring(opts.FooofLabel)
    opts.FooofLabel = {char(opts.FooofLabel)};
end

nRanges = numel(ExponentsClean);

figName = 'GEDAI evaluation — exponent by sleep stage';
fig = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [2 2 22 13*nRanges], 'Name', figName);
tiledlayout(nRanges, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

for iRange = 1:nRanges
    ax = nexttile;
    frLabel = '';
    if numel(opts.FooofLabel) >= iRange
        frLabel = opts.FooofLabel{iRange};
    end
    draw_tile(ax, ExponentsClean{iRange}, ExponentsRaw{iRange}, StageScoring, frLabel);
end

save_fig(fig, opts.SavePath, 'exponent_by_stage');
end

% -------------------------------------------------------------------------
function draw_tile(ax, ExponentsClean, ExponentsRaw, StageScoring, frLabel)

stageOrder  = [  1,        0,        -1,       -2,       -3    ];
stageLabels = {'Wake',   'REM',    'N1',     'N2',     'N3'   };
stageColors = [
    0.96  0.60  0.15;   % Wake  — sunset amber
    0.55  0.50  0.90;   % REM   — twilight lavender
    0.35  0.60  0.90;   % N1    — early-night blue
    0.15  0.40  0.75;   % N2    — deep-night blue
    0.05  0.10  0.50;   % N3    — midnight navy
];

expClean = ExponentsClean(:)';
expRaw   = ExponentsRaw(:)';
scoring  = StageScoring(:)';
n        = min([numel(expClean), numel(expRaw), numel(scoring)]);
expClean = expClean(1:n);
expRaw   = expRaw(1:n);
scoring  = scoring(1:n);

presentIdx  = arrayfun(@(s) any(scoring == s | (s == 0 & any(scoring == 5))), stageOrder);
stageOrder  = stageOrder(presentIdx);
stageLabels = stageLabels(presentIdx);
stageColors = stageColors(presentIdx, :);
nStages     = numel(stageOrder);

maxViolinW  = 0.28;
dotOff      = 0.15;
dotSz       = 18;
colDecrease = [0.85 0.15 0.15];

axes(ax); hold on;

for iStage = 1:nStages
    d    = stageOrder(iStage);
    col  = stageColors(iStage, :);
    xPos = iStage;

    if d == 0
        mask = scoring == 0 | scoring == 5;
    else
        mask = scoring == d;
    end
    rawVals   = expRaw(mask);
    cleanVals = expClean(mask);

    valid     = ~isnan(rawVals) & ~isnan(cleanVals);
    rawVals   = rawVals(valid);
    cleanVals = cleanVals(valid);
    if numel(rawVals) < 3, continue; end

    jitter = (rand(1, numel(rawVals)) - 0.5) * 0.018;
    xRaw   = xPos - dotOff + jitter;
    xClean = xPos + dotOff + jitter;

    % Paired epoch lines
    for iEp = 1:numel(rawVals)
        if cleanVals(iEp) >= rawVals(iEp)
            lineCol = [col,         0.28];
        else
            lineCol = [colDecrease, 0.28];
        end
        line([xRaw(iEp), xClean(iEp)], [rawVals(iEp), cleanVals(iEp)], ...
            'Color', lineCol, 'LineWidth', 0.7, 'HandleVisibility', 'off');
    end

    % Half violin — raw (extends left)
    [densR, yiR] = ksdensity(rawVals);
    densR_norm   = densR / max(densR) * maxViolinW;
    patch([(xPos - dotOff - densR_norm), repmat(xPos - dotOff, 1, numel(yiR))], ...
          [yiR, fliplr(yiR)], col, 'FaceAlpha', 0.30, 'EdgeColor', 'none', ...
          'HandleVisibility', 'off');
    medR     = median(rawVals);
    [~, iMR] = min(abs(yiR - medR));
    line([xPos - dotOff - densR_norm(iMR), xPos - dotOff], [medR, medR], ...
        'Color', col * 0.75, 'LineWidth', 2.0, 'HandleVisibility', 'off');

    % Half violin — clean (extends right)
    [densC, yiC] = ksdensity(cleanVals);
    densC_norm   = densC / max(densC) * maxViolinW;
    patch([(xPos + dotOff + densC_norm), repmat(xPos + dotOff, 1, numel(yiC))], ...
          [yiC, fliplr(yiC)], col, 'FaceAlpha', 0.68, 'EdgeColor', 'none', ...
          'HandleVisibility', 'off');
    medC     = median(cleanVals);
    [~, iMC] = min(abs(yiC - medC));
    line([xPos + dotOff, xPos + dotOff + densC_norm(iMC)], [medC, medC], ...
        'Color', col * 0.75, 'LineWidth', 2.0, 'HandleVisibility', 'off');

    % Scatter dots
    scatter(xRaw,   rawVals,   dotSz, col, 'filled', ...
        'MarkerFaceAlpha', 0.40, 'MarkerEdgeColor', 'none', 'HandleVisibility', 'off');
    scatter(xClean, cleanVals, dotSz, col, 'filled', ...
        'MarkerFaceAlpha', 0.82, 'MarkerEdgeColor', 'none', 'HandleVisibility', 'off');
end

set(ax, 'XTick', 1:nStages, 'XTickLabel', stageLabels, ...
    'FontSize', 11, 'Box', 'off', 'TickDir', 'out');
ylabel(ax, 'Aperiodic exponent', 'FontSize', 11);
if ~isempty(frLabel)
    title(ax, frLabel, 'FontSize', 11, 'FontWeight', 'normal');
end
xlim(ax, [0.5, nStages + 0.5]);
grid(ax, 'on');
set(ax, 'GridAlpha', 0.12, 'GridLineStyle', ':', 'XGrid', 'off');

if strcmp(frLabel, '30–45 Hz')
    ylim(ax, [-8 10]);
else
    ylim(ax, [-2 4]);
end

% Legend (first tile only — drawn via invisible patches for alpha encoding)
patch(nan, nan, [0.55 0.55 0.55], 'FaceAlpha', 0.30, 'EdgeColor', 'none', ...
    'DisplayName', 'before GEDAI');
patch(nan, nan, [0.55 0.55 0.55], 'FaceAlpha', 0.68, 'EdgeColor', 'none', ...
    'DisplayName', 'after GEDAI');
legend(ax, 'Location', 'best', 'FontSize', 9, 'Box', 'off');
end
