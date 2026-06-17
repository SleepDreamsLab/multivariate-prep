function gedai_characteristics(GedaiStages, opts)
% GEDAI.EVALPLOTS.GEDAI_CHARACTERISTICS  Plot GEDAI cleaning characteristics per sleep stage.
%
%   gedai.evalplots.gedai_characteristics(EEGgedai.etc.gedai)
%   gedai.evalplots.gedai_characteristics(..., Name, Value)
%
%   Input
%   -----
%   GedaiStages : struct array from EEGgedai.etc.gedai — one element per stage.
%                 Required fields: sleepStage, SENSAI_score, SENSAI_score_per_band,
%                 artifact_threshold_per_band, mean_ENOVA, ENOVA_per_epoch, time_spent_min.
%
%   Optional name-value pairs
%   -------------------------
%   Srate    : EEG sampling rate (Hz); used to label frequency bands via
%              [Fs/2^(f+1), Fs/2^f]. Falls back to band indices if omitted.
%   SavePath : base path for saving; suffix + '.png' appended.

arguments
    GedaiStages
    opts.Srate   (1,1) double = 0
    opts.SavePath             = ''
end

nStages = numel(GedaiStages);
nBands  = numel(GedaiStages(1).SENSAI_score_per_band);

%%% --- Band-axis tick labels ---
if opts.Srate > 0
    fmt      = @(x) strtrim(sprintf('%g', round(x, 2)));
    bandLbls = cell(1, nBands);
    bandLbls{1} = 'BB';
    for iBand = 2:nBands
        f = iBand - 1;
        bandLbls{iBand} = sprintf('%s–%s Hz', fmt(opts.Srate/2^(f+1)), fmt(opts.Srate/2^f));
    end
else
    bandLbls = arrayfun(@(f) num2str(f), 1:nBands, 'UniformOutput', false);
end

%%% --- Collect per-stage data ---
sensai      = arrayfun(@(s) s.SENSAI_score,   GedaiStages);
meanEnova   = arrayfun(@(s) s.mean_ENOVA,     GedaiStages);
timeSpent   = arrayfun(@(s) s.time_spent_min, GedaiStages);

sensaiBands = zeros(nStages, nBands);
threshBands = zeros(nStages, nBands);
for i = 1:nStages
    sensaiBands(i, :) = GedaiStages(i).SENSAI_score_per_band;
    threshBands(i, :) = GedaiStages(i).artifact_threshold_per_band;
end

enovaEpochs = {GedaiStages.ENOVA_per_epoch};

cols = zeros(nStages, 3);
lbls = cell(1, nStages);
for i = 1:nStages
    sc = GedaiStages(i).sleepStage(:)';
    if numel(sc) == 1
        cols(i, :) = stage_color(sc);
        lbls{i}    = stage_name(sc);
    else
        c = zeros(numel(sc), 3);
        n = cell(1, numel(sc));
        for k = 1:numel(sc)
            c(k, :) = stage_color(sc(k));
            n{k}    = stage_name(sc(k));
        end
        cols(i, :) = mean(c, 1);
        lbls{i}    = strjoin(n, '+');
    end
end

%%% --- Figure ---
fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 28 22], ...
    'Name', 'GEDAI evaluation — cleaning characteristics');
tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

%%% --- SENSAI score ---
ax = nexttile(1);
hold on;
for i = 1:nStages
    bar(i, sensai(i), 'FaceColor', cols(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.85);
    text(i, sensai(i) + 3, sprintf('%.1f', sensai(i)), ...
        'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', 'k');
end
set(ax, 'XTick', 1:nStages, 'XTickLabel', lbls, 'FontSize', 10, ...
    'Box', 'off', 'TickDir', 'out');
ylabel('Score', 'FontSize', 10);
title('SENSAI score', 'FontSize', 11, 'FontWeight', 'bold');
ylim([0 100]); grid on; set(ax, 'GridAlpha', 0.15, 'GridLineStyle', ':', 'XGrid', 'off');

%%% --- Time spent ---
ax = nexttile(2);
hold on;
for i = 1:nStages
    bar(i, timeSpent(i), 'FaceColor', cols(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.85);
    text(i, timeSpent(i), sprintf('%.0f', timeSpent(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 8, 'Color', 'k');
end
totalTime = sum(timeSpent);
bar(nStages+1, totalTime, 'FaceColor', [0.4 0.4 0.4], 'EdgeColor', 'none', 'FaceAlpha', 0.85);
text(nStages+1, totalTime, sprintf('%.0f', totalTime), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 8, 'Color', 'k');
set(ax, 'XTick', 1:nStages+1, 'XTickLabel', [lbls, {'Total'}], 'FontSize', 10, ...
    'Box', 'off', 'TickDir', 'out');
ylabel('Minutes', 'FontSize', 10);
title('Time spent', 'FontSize', 11, 'FontWeight', 'bold');
grid on; set(ax, 'GridAlpha', 0.15, 'GridLineStyle', ':', 'XGrid', 'off');

%%% --- SENSAI score per band ---
ax = nexttile(3);
hold on;
for i = 1:nStages
    plot(1:nBands, sensaiBands(i,:), '-o', 'Color', cols(i,:), ...
        'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerFaceColor', cols(i,:), ...
        'DisplayName', lbls{i});
end
set(ax, 'XTick', 1:nBands, 'XTickLabel', bandLbls, 'XTickLabelRotation', 45, ...
    'FontSize', 10, 'Box', 'off', 'TickDir', 'out');
ylabel('Score', 'FontSize', 10);
title('SENSAI score per band', 'FontSize', 11, 'FontWeight', 'bold');
ylim([0 100]);
legend('Location', 'best', 'FontSize', 8, 'Box', 'off');
grid on; set(ax, 'GridAlpha', 0.15, 'GridLineStyle', ':');

%%% --- Artifact threshold per band ---
ax = nexttile(4);
hold on;
for i = 1:nStages
    plot(1:nBands, threshBands(i,:), '-o', 'Color', cols(i,:), ...
        'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerFaceColor', cols(i,:), ...
        'HandleVisibility', 'off');
end
set(ax, 'XTick', 1:nBands, 'XTickLabel', bandLbls, 'XTickLabelRotation', 45, ...
    'FontSize', 10, 'Box', 'off', 'TickDir', 'out');
ylabel('Threshold', 'FontSize', 10);
title('Artifact threshold per band', 'FontSize', 11, 'FontWeight', 'bold');
ylim([-6 8]); grid on; set(ax, 'GridAlpha', 0.15, 'GridLineStyle', ':');

%%% --- Mean ENOVA ---
ax = nexttile(5);
hold on;
for i = 1:nStages
    bar(i, meanEnova(i), 'FaceColor', cols(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.85);
    text(i, meanEnova(i) + 0.03, sprintf('%.2f', meanEnova(i)), ...
        'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', 'k');
end
set(ax, 'XTick', 1:nStages, 'XTickLabel', lbls, 'FontSize', 10, ...
    'Box', 'off', 'TickDir', 'out');
ylabel('ENOVA', 'FontSize', 10);
title('Mean ENOVA', 'FontSize', 11, 'FontWeight', 'bold');
ylim([0 1]); grid on; set(ax, 'GridAlpha', 0.15, 'GridLineStyle', ':', 'XGrid', 'off');

%%% --- ENOVA per epoch (violin) ---
ax = nexttile(6);
hold on;
maxViolinW = 0.35;
for i = 1:nStages
    vals = enovaEpochs{i}(:);
    vals = vals(isfinite(vals));
    if numel(vals) < 3, continue; end

    [dens, yi] = ksdensity(vals, 'BoundaryCorrection', 'reflection', ...
        'Support', [0 Inf]);
    dens_norm  = dens / max(dens) * maxViolinW;

    % Full symmetric violin
    patch([i + dens_norm, i - fliplr(dens_norm)], [yi, fliplr(yi)], ...
        cols(i,:), 'FaceAlpha', 0.55, 'EdgeColor', 'none', 'HandleVisibility', 'off');

    % Median line
    med = median(vals);
    [~, iM] = min(abs(yi - med));
    line([i - dens_norm(iM), i + dens_norm(iM)], [med med], ...
        'Color', cols(i,:) * 0.7, 'LineWidth', 2, 'HandleVisibility', 'off');
end
set(ax, 'XTick', 1:nStages, 'XTickLabel', lbls, 'FontSize', 10, ...
    'Box', 'off', 'TickDir', 'out');
ylabel('ENOVA', 'FontSize', 10);
title('ENOVA per epoch', 'FontSize', 11, 'FontWeight', 'bold');
ylim([0 1]); grid on; set(ax, 'GridAlpha', 0.15, 'GridLineStyle', ':', 'XGrid', 'off');

save_fig(fig, opts.SavePath, 'gedai_characteristics');
end
