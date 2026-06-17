function epoch_overlay(EEGstage, EEGclean, StageScoring, opts)
% GEDAI.EVAL.EPOCH_OVERLAY  Per-epoch raw vs. clean signal overlay.
%
%   gedai.eval.epoch_overlay(EEGstage, EEGclean, StageScoring)
%   gedai.eval.epoch_overlay(..., Name, Value)
%
%   Inputs
%   ------
%   EEGstage     : raw continuous EEG (full channel set, including EOG).
%   EEGclean     : GEDAI-cleaned EEG (same channel set).
%   StageScoring : per-epoch sleep-stage codes.
%
%   Optional name-value pairs
%   -------------------------
%   EpochLength  : epoch duration in seconds (default 30).
%   EpochsToPlot : epoch indices to plot. Default = first epoch of each stage.
%   ChanOffset   : vertical spacing between channels in µV (default 150).
%   SavePath     : base path for saving; suffix + '.png' appended.

arguments
    EEGstage
    EEGclean
    StageScoring {mustBeVector}
    opts.EpochLength  (1,1) double = 30
    opts.EpochsToPlot double       = []
    opts.ChanOffset   (1,1) double = 150
    opts.SavePath                  = ''
end

epochSamples = round(opts.EpochLength * EEGclean.srate);
stageScoring = StageScoring(:)';
uniqueStages = unique(stageScoring);

if isempty(opts.EpochsToPlot)
    epochsToPlot = arrayfun(@(d) find(stageScoring == d, 1, 'first'), uniqueStages);
else
    epochsToPlot = opts.EpochsToPlot(:)';
end

allLabels    = {EEGclean.chanlocs.labels};
eegTargets   = {'FP1', 'FP2', 'F3', 'F4', 'C3', 'C4', 'P3', 'P4', 'O1', 'O2'};
displayChans = [];
for iT = 1:numel(eegTargets)
    idx = find(strcmpi(allLabels, eegTargets{iT}), 1);
    if ~isempty(idx), displayChans(end+1) = idx; end  %#ok<AGROW>
end
eogChans     = find(cellfun(@(s) ~isempty(regexpi(s, 'EOG')), allLabels));
displayChans = [displayChans, eogChans];

if isempty(displayChans)
    warning('gedai.eval: no target or EOG channels found; using all channels.');
    displayChans = 1:numel(allLabels);
end

colorRaw   = [0.85 0.40 0.50];
colorClean = [0.08 0.22 0.55];
lw_raw     = 1.0;
lw_clean   = 1.2;

nPlotChans = numel(displayChans);
plotNames  = allLabels(displayChans);
offsets    = (nPlotChans-1:-1:0) * opts.ChanOffset;

for iEp = epochsToPlot
    if iEp < 1 || iEp > numel(stageScoring), continue; end

    epSamps = (iEp-1)*epochSamples + (1:epochSamples);
    if epSamps(end) > size(EEGclean.data, 2), continue; end

    t_s      = (0:epochSamples-1) / EEGclean.srate;
    rawEp    = EEGstage.data(displayChans, epSamps);
    cleanEp  = EEGclean.data(displayChans, epSamps);
    stageLbl = stage_name(stageScoring(iEp));

    fig = figure('Color', 'w', 'Units', 'centimeters', ...
        'Position', [2 2 55 24], ...
        'Name', sprintf('Epoch %d (%s) — signal overlay', iEp, stageLbl));
    tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    nexttile; hold on;

    for iCh = 1:nPlotChans
        plot(t_s, rawEp(iCh,:)   + offsets(iCh), 'Color', colorRaw,   ...
            'LineWidth', lw_raw,   'HandleVisibility', 'off');
        plot(t_s, cleanEp(iCh,:) + offsets(iCh), 'Color', colorClean, ...
            'LineWidth', lw_clean, 'HandleVisibility', 'off');
    end
    yline(offsets, ':', 'HandleVisibility', 'off');

    % Amplitude scale bar — centred at the bottom channel's baseline
    scaleAmp = 100;                    % µV
    sx       = t_s(end) * 0.985;
    sy0      = offsets(end) - scaleAmp/2;
    sy1      = offsets(end) + scaleAmp/2;
    capLen   = t_s(end) * 0.004;
    line([sx sx],                   [sy0 sy1], 'Color', 'k', 'LineWidth', 2,   'HandleVisibility', 'off');
    line([sx-capLen sx+capLen], [sy0 sy0],     'Color', 'k', 'LineWidth', 2,   'HandleVisibility', 'off');
    line([sx-capLen sx+capLen], [sy1 sy1],     'Color', 'k', 'LineWidth', 2,   'HandleVisibility', 'off');
    text(sx, (sy0+sy1)/2, sprintf(' %d \\muV', scaleAmp), ...
        'FontSize', 9, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');

    plot(nan, nan, 'Color', colorRaw,   'LineWidth', lw_raw,   'DisplayName', 'raw');
    plot(nan, nan, 'Color', colorClean, 'LineWidth', lw_clean, 'DisplayName', 'clean');

    legend('Location', 'northeast', 'FontSize', 9, 'Box', 'off');
    set(gca, 'YTick', fliplr(offsets), 'YTickLabel', fliplr(plotNames), ...
        'FontSize', 9, 'Box', 'off', 'TickDir', 'out');
    xlim([0 t_s(end)]);
    xtickformat(gca, '%g s');
    title(sprintf('Epoch %d — %s', iEp, stageLbl), 'FontSize', 11, 'FontWeight', 'bold');
    ylim([opts.ChanOffset * -1, nPlotChans * opts.ChanOffset]);

    save_fig(fig, opts.SavePath, sprintf('epoch_%03d_%s', iEp, stageLbl));
end
end
