function fig = amica_lt_by_stage(Lht, V, scoringDigits, opts)
% ICA.PLOT.AMICA_LT_BY_STAGE  Per-model AMICA log-likelihood/posterior-odds vs.
% hypnogram/sleep stage.
%
%   ica.plot.amica_lt_by_stage(Lht, V, scoringDigits)
%   ica.plot.amica_lt_by_stage(Lht, V, scoringDigits, Label=root, SavePath=fullfile(outDir, id))
%
%   Lht           – mod.Lht from loadmodout15: num_models x nsamples,
%                   per-model posterior log-likelihood per sample — unlike
%                   mod.Lt, which is already marginalized across models and
%                   hides per-model fit differences.
%   V             – mod.v from loadmodout15: num_models x nsamples, log10
%                   posterior odds of each model vs. the combined mixture,
%                   i.e. v(h,t) = log10 P(model h | sample t), so v <= 0.
%                   Converted to per-sample probability (10.^V) and averaged
%                   per epoch in probability space for display as 0-100%
%                   (100/num_models = evenly split) — averaging V itself
%                   (log-odds) before exponentiating would give the geometric
%                   mean instead, which is noisier and doesn't sum to 100%
%                   across models.
%   scoringDigits – per-epoch sleep-stage codes from scoreloader (N3=-3 .. W=1).
%   Srate         – sample rate of Lht/V, in Hz (default 250/4, i.e. the AMICA
%                   decimated rate used for this pipeline).
%   EpochLength   – scoring epoch length in seconds (default 30).
%   ModelNames    – cell of per-model names (default {'Model 1', ...}).
%   YLim          – y-limits for the log-likelihood panels (default [-800 -200]).
%   Label         – text (e.g. the AMICA output folder) annotated on the figure.
%   SavePath      – base path for saving; '_lt_by_stage.png' appended.
%                   '' (default) = don't save.

arguments
    Lht
    V
    scoringDigits (:,1) double
    opts.Srate       (1,1) double = 250/4
    opts.EpochLength (1,1) double = 30
    opts.ModelNames  cell         = {}
    opts.YLim        (1,2) double = [-800 -200]
    opts.Label       char         = ''
    opts.SavePath    char         = ''
end

if isvector(Lht), Lht = Lht(:)'; end
if isvector(V),   V   = V(:)';   end
nModels = size(Lht, 1);
if isempty(opts.ModelNames)
    opts.ModelNames = arrayfun(@(h) sprintf('Model %d', h), 1:nModels, 'uni', 0);
end

modelColors = [
    0.20 0.36 0.60;   % Model 1 — slate blue
    0.85 0.33 0.20;   % Model 2 — burnt orange
    0.30 0.65 0.35;   % Model 3 — green
    0.55 0.35 0.75];  % Model 4 — purple
if nModels > size(modelColors, 1)
    modelColors = lines(nModels);
else
    modelColors = modelColors(1:nModels, :);
end

epochLen = opts.EpochLength * opts.Srate;
nEpochs  = floor(size(Lht, 2) / epochLen);

epochMean = @(X) cell2mat(arrayfun(@(h) ...
    mean(reshape(X(h, 1:nEpochs*epochLen), epochLen, nEpochs), 1), (1:size(X,1))', 'uni', 0));
Lht_epoch = epochMean(Lht);
Pct_epoch = 100 * epochMean(10.^V);   % average the per-sample posterior probability, not exp(mean(log10-odds))

stageVals  = [-3 -2 -1 0 1];
stageNames = {'N3','N2','N1','REM','W'};
stageLabel = scoringDigits(1:nEpochs);
validIdx   = ismember(stageLabel, stageVals);
stageCat   = categorical(stageLabel(validIdx), stageVals, stageNames, 'Ordinal', true);

% With a single model, the posterior-probability panel is trivially 100% for
% every epoch (nothing to compare against), so skip it and use one less row.
nRows = 4;
if nModels == 1
    nRows = 3;
end

fig = figure('Color', 'w');
set(fig, 'Position', [300 150 900 225*nRows])

% --- hypnogram ---
ax1 = subplot(nRows,1,1);
plot(scoringDigits, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1); hold on;
plot(scoringDigits, 'o', 'MarkerFaceColor', [0.2 0.4 0.8], ...
    'MarkerEdgeColor', 'none', 'MarkerSize', 4);
ylabel('Sleep stage');
title('Hypnogram');
box off;
set(gca, 'YTick', -3:1, 'YTickLabel', {'N3','N2','N1','REM','W'});
xlim([1 nEpochs]);

% --- per-model log-likelihood per epoch ---
ax2 = subplot(nRows,1,2);
hold on;
wakeIdx = find(scoringDigits(1:nEpochs) == 1);
yMin = min(Lht_epoch(:)); yMax = max(Lht_epoch(:));
for k = 1:numel(wakeIdx)
    xc = wakeIdx(k);
    patch([xc-0.5 xc+0.5 xc+0.5 xc-0.5], [yMin yMin yMax yMax], ...
        [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.6, 'HandleVisibility', 'off');
end
for h = 1:nModels
    plot(Lht_epoch(h,:), '-', 'Color', modelColors(h,:), 'LineWidth', 1.2, ...
        'DisplayName', opts.ModelNames{h});
end
xlabel('Epoch (30 s)');
ylabel('Mean log-likelihood');
title('Higher (less negative) = better fit to that model''s learned sources');
box off;
xlim([1 nEpochs]);
ylim(opts.YLim)
legend('show', 'Box', 'off', 'Location', 'southeast', 'Orientation', 'horizontal', 'FontSize', 8);

axLinked = [ax1 ax2];

if nModels > 1
    % --- per-model posterior probability per epoch ---
    ax3 = subplot(nRows,1,3);
    hold on;
    pctMin = 0;                   % theoretical min: P(model h | sample) -> 0
    pctMax = 100;                 % theoretical max: P(model h | sample) = 1 (one model certain)
    for k = 1:numel(wakeIdx)
        xc = wakeIdx(k);
        patch([xc-0.5 xc+0.5 xc+0.5 xc-0.5], [pctMin pctMin pctMax pctMax], ...
            [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.6, 'HandleVisibility', 'off');
    end
    yline(100/nModels, ':', 'Color', [.45 .45 .45], 'HandleVisibility', 'off');   % evenly-split reference
    for h = 1:nModels
        plot(Pct_epoch(h,:), '-', 'Color', modelColors(h,:), 'LineWidth', 1.2, ...
            'DisplayName', opts.ModelNames{h});
    end
    xlabel('Epoch (30 s)');
    ylabel('Posterior probability (%)');
    title('Higher = epoch more likely generated by that model');
    box off;
    xlim([1 nEpochs]);
    ylim([pctMin pctMax]);
    legend('show', 'Box', 'off', 'Location', 'southeast', 'Orientation', 'horizontal', 'FontSize', 8);

    axLinked   = [axLinked ax3];
    boxplotRow = 4;
else
    boxplotRow = 3;
end

linkaxes(axLinked, 'x');

% --- per-model log-likelihood by stage ---
validEpochs = find(validIdx);
stageIdx    = double(stageCat);   % 1..numel(stageNames), matches stageCat's ordinal coding

ax4 = subplot(nRows,1,boxplotRow);
hold(ax4, 'on'); box(ax4, 'off');
if nModels > 1
    % Manual x-offsets: boxchart's GroupByColor doesn't expose spacing
    % directly, so cluster each stage's per-model boxes tightly together
    % and leave a clear gap to the next stage's cluster.
    clusterWidth = 0.6;                              % fraction of each stage's 1-unit slot used by its box cluster
    boxWidth     = clusterWidth / nModels * 0.85;     % boxes packed close within a category
    offsets      = linspace(-clusterWidth/2 + boxWidth/2, clusterWidth/2 - boxWidth/2, nModels);
    for h = 1:nModels
        boxchart(ax4, stageIdx + offsets(h), Lht_epoch(h, validEpochs)', ...
            'BoxWidth', boxWidth, 'BoxFaceColor', modelColors(h,:), ...
            'MarkerColor', modelColors(h,:), 'MarkerStyle', '.', ...
            'DisplayName', opts.ModelNames{h});
    end
    legend('show', 'Box', 'off', 'Location', 'best', 'FontSize', 8);
else
    % Single model: colour boxes by sleep stage instead of by model.
    stageColors = [
        0.20 0.20 0.55;   % N3
        0.30 0.45 0.75;   % N2
        0.45 0.65 0.85;   % N1
        0.80 0.40 0.40;   % REM
        0.55 0.55 0.55];  % W
    for i = 1:numel(stageNames)
        idx = stageCat == stageNames{i};
        boxchart(ax4, stageIdx(idx), Lht_epoch(1, validEpochs(idx))', ...
            'BoxFaceColor', stageColors(i,:), 'MarkerColor', stageColors(i,:), 'MarkerStyle', '.');
    end
end
set(ax4, 'XTick', 1:numel(stageNames), 'XTickLabel', stageNames);
xlim(ax4, [0.5 numel(stageNames) + 0.5]);
xlabel('Sleep stage');
ylabel('Mean log-likelihood');
if nModels > 1
    title('AMICA Lht by stage (per model)');
else
    title('AMICA Lht by stage');
end
ylim(opts.YLim)

if ~isempty(opts.Label)
    annotation(fig, 'textbox', [0.01 0.005 0.3 0.02], 'String', opts.Label, ...
        'FontSize', 6, 'EdgeColor', 'none', 'Interpreter', 'none');
end

save_fig(fig, opts.SavePath, 'lt_by_stage');
end
