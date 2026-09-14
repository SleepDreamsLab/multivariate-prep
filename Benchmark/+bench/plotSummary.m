function files = plotSummary(metrics, details, chanlocs, figDir, baseName)
% PLOTSUMMARY  Summary figures for one recording and one GEDAI run of the benchmark.
%
%   files = bench.plotSummary(metrics, details, chanlocs, figDir, baseName)
%
%   metrics    rows of the benchmark metrics table for one run and one scale
%   details    struct with one field per condition ('signals', 'eog', ...) holding the
%              detail output of bench.scoreSignals / bench.scoreArtifacts
%   chanlocs   EEG.chanlocs of the data channels (artifact topographies)
%   figDir     output folder
%   baseName   file name prefix, e.g. <fileID>_desc-<desc>_run-1
%
%   Writes, when the data are there:
%     <baseName>_summary.png   SER and gain per class, detector recall, ARR and in-band
%                              residual excess per artifact, each grouped by stage
%     <baseName>_events.png    median-SER example event and event-locked means per class
%     <baseName>_spectra.png   PSD ratio contaminated/clean before and after GEDAI
%     <baseName>_topo.png      per-channel ARR and in-band residual per artifact
%   Figures are created invisible and closed after saving.

arguments
    metrics table
    details struct
    chanlocs struct
    figDir char
    baseName char
end

if ~isfolder(figDir), mkdir(figDir); end
files = {};
artifactTypes = intersect({'eog', 'emg', 'ecg'}, fieldnames(details), 'stable');

%%% ---- Summary bars ----
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 1400 850]);
layout = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(layout, strrep(baseName, '_', '\_'));
groupedBars(nexttile(layout), metrics, 'SER_dB', 'SER (dB)', 'Injected graphoelements: signal-to-error');
groupedBars(nexttile(layout), metrics, 'gain', 'gain', 'Fraction of injected field preserved');
groupedBars(nexttile(layout), metrics, 'recall_retention', 'recall out / in', 'Detector recall retention');
groupedBars(nexttile(layout), metrics, 'ARR_dB', 'ARR (dB)', 'Injected artifacts: artifact-to-residue');
groupedBars(nexttile(layout), metrics, 'residual_gain', 'residual gain', 'Fraction of artifact left');
groupedBars(nexttile(layout), metrics, 'psd_residual_excess', 'residual excess', 'In-band excess power left');
files{end+1} = saveFigure(fig, fullfile(figDir, [baseName '_summary.png']));

%%% ---- Example events ----
if isfield(details, 'signals') && ~isempty(fieldnames(details.signals.example))
    classNames = fieldnames(details.signals.example);
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 1400 300 * numel(classNames)]);
    layout = tiledlayout(fig, numel(classNames), 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    for k = 1:numel(classNames)
        ex = details.signals.example.(classNames{k});
        ax = nexttile(layout);
        plot(ax, ex.t, ex.X, 'Color', [0.75 0.75 0.75]); hold(ax, 'on');
        plot(ax, ex.t, ex.Xhat, 'Color', [0.35 0.35 0.35]);
        plot(ax, ex.t, ex.S, 'k', 'LineWidth', 1.2);
        plot(ax, ex.t, ex.Shat, 'r', 'LineWidth', 1.2);
        title(ax, sprintf('%s at %s, median-SER event (%.1f dB)', strrep(classNames{k}, '_', ' '), ...
            ex.channel, ex.SER_dB));
        ylabel(ax, '\muV');
        if k == 1, legend(ax, {'input', 'GEDAI output', 'injected', 'recovered'}, 'Location', 'best'); end

        erp = details.signals.erp.(classNames{k});
        ax = nexttile(layout);
        plot(ax, erp.t, erp.S, 'k', 'LineWidth', 1.2); hold(ax, 'on');
        plot(ax, erp.t, erp.Shat, 'r', 'LineWidth', 1.2);
        title(ax, sprintf('event-locked mean at peak channel (n = %d)', erp.n));
        if k == numel(classNames), xlabel(ax, 's from event onset'); end
    end
    files{end+1} = saveFigure(fig, fullfile(figDir, [baseName '_events.png']));
end

%%% ---- PSD ratio spectra and topographies ----
if ~isempty(artifactTypes)
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 450 * numel(artifactTypes) 380]);
    layout = tiledlayout(fig, 1, numel(artifactTypes), 'TileSpacing', 'compact', 'Padding', 'compact');
    for k = 1:numel(artifactTypes)
        spectra = details.(artifactTypes{k}).spectra;
        ax = nexttile(layout); hold(ax, 'on');
        colours = lines(numel(spectra));
        for s = 1:numel(spectra)
            plot(ax, spectra(s).f, spectra(s).ratioIn_dB, '--', 'Color', colours(s, :));
            plot(ax, spectra(s).f, spectra(s).ratioOut_dB, '-', 'Color', colours(s, :), 'LineWidth', 1.2, ...
                'DisplayName', spectra(s).stage);
        end
        yline(ax, 0, 'k:');
        xlim(ax, [0 min(100, spectra(1).f(end))]);
        xlabel(ax, 'Hz'); ylabel(ax, 'dB over clean');
        title(ax, sprintf('%s: dashed before, solid after GEDAI', upper(artifactTypes{k})));
        legend(ax, findobj(ax, 'LineStyle', '-', '-not', 'DisplayName', ''), 'Location', 'best');
    end
    files{end+1} = saveFigure(fig, fullfile(figDir, [baseName '_spectra.png']));

    if exist('topoplot', 'file') == 2 && isfield(chanlocs, 'X')
        fig = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 350 * numel(artifactTypes) 600]);
        for k = 1:numel(artifactTypes)
            topo = details.(artifactTypes{k}).topo;
            subplot(2, numel(artifactTypes), k);
            topoplot(topo.psd_contamination_dB, chanlocs, 'electrodes', 'off');
            colorbar; title(sprintf('%s added (dB, in-band)', upper(artifactTypes{k})));
            subplot(2, numel(artifactTypes), numel(artifactTypes) + k);
            topoplot(topo.psd_residual_dB, chanlocs, 'electrodes', 'off');
            colorbar; title('left after GEDAI (dB)');
        end
        files{end+1} = saveFigure(fig, fullfile(figDir, [baseName '_topo.png']));
    end
end
end

% -------------------------------------------------------------------------
function groupedBars(ax, metrics, metricName, yLabel, titleText)
rows = metrics(strcmp(metrics.metric, metricName), :);
if isempty(rows)
    axis(ax, 'off');
    title(ax, [titleText ' (not run)']);
    return
end
classes = unique(rows.class, 'stable');
stages  = unique(rows.stage, 'stable');
values  = nan(numel(classes), numel(stages));
for i = 1:height(rows)
    values(strcmp(classes, rows.class{i}), strcmp(stages, rows.stage{i})) = rows.value(i);
end
nClasses = numel(classes);
if nClasses == 1
    values(2, :) = NaN;          % a single row would be drawn as ungrouped bars
end
bar(ax, values);
xticks(ax, 1:nClasses);
xticklabels(ax, classes);
xlim(ax, [0.5, nClasses + 0.5]);
legend(ax, stages, 'Location', 'best');
ylabel(ax, yLabel);
title(ax, titleText);
set(ax, 'TickLabelInterpreter', 'none');
grid(ax, 'on');
end

% -------------------------------------------------------------------------
function file = saveFigure(fig, file)
exportgraphics(fig, file, 'Resolution', 130);
close(fig);
fprintf('[FIG] %s\n', file);
end
