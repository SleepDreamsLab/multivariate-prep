function slopes_timecourse(SlopesRaw, SlopesClean, FreqRanges, StageScoring, opts)
% GEDAI.EVALPLOTS.SLOPES_TIMECOURSE  Aperiodic exponent time course: raw vs. clean.
%
%   gedai.evalplots.slopes_timecourse(SlopesRaw, SlopesClean, FreqRanges, StageScoring)
%   gedai.evalplots.slopes_timecourse(..., Name, Value)
%
%   Inputs
%   ------
%   SlopesRaw   : 1 x nRanges cell — each cell is 1 x nEpochs exponent vector (raw).
%   SlopesClean : 1 x nRanges cell — each cell is 1 x nEpochs exponent vector (clean).
%   FreqRanges  : 1 x nRanges cell — each cell is [fMin fMax] (Hz).
%   StageScoring: per-epoch sleep-stage codes (used for hypnogram row).
%
%   Optional name-value pairs
%   -------------------------
%   EpochLength : epoch duration in seconds (default 30).
%   SavePath    : base path for saving; suffix + '.png' appended.

arguments
    SlopesRaw    cell
    SlopesClean  cell
    FreqRanges   cell
    StageScoring {mustBeVector}
    opts.EpochLength (1,1) double = 30
    opts.SavePath                 = ''
end

nRanges    = numel(FreqRanges);
colorRaw   = [0.65 0.65 0.65];
colorClean = [0.08 0.22 0.55];

nEpochs    = max(cellfun(@numel, SlopesRaw));
T_min      = (0:nEpochs-1) * opts.EpochLength / 60;
xLim       = [0  T_min(end) + opts.EpochLength/60];
scoringTF  = StageScoring(:)';
if numel(scoringTF) > nEpochs, scoringTF = scoringTF(1:nEpochs); end

fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 32 4 + 6*nRanges], ...
    'Name', 'GEDAI evaluation — aperiodic slopes');
tiledlayout(nRanges + 1, 1, 'TileSpacing', 'tight', 'Padding', 'compact');

%%% --- Hypnogram row ---
ax_hyp = nexttile;
hold on;

codes = [-3, -2, -1,  0,  1,  5];
yPos  = [  1,  2,  3,  5,  4,  5];
yVals = nan(size(scoringTF));
for iEp = 1:numel(scoringTF)
    idx = find(codes == scoringTF(iEp), 1);
    if ~isempty(idx), yVals(iEp) = yPos(idx); end
end

remMask = scoringTF == 0 | scoringTF == 5;
for iEp = 1:numel(scoringTF)
    if remMask(iEp)
        t0 = T_min(iEp);
        t1 = t0 + opts.EpochLength / 60;
        patch([t0 t1 t1 t0], [4.5 4.5 5.5 5.5], [1 0.85 0.85], ...
            'EdgeColor', 'none', 'FaceAlpha', 0.5);
    end
end

stairs(T_min, yVals, 'Color', [0.15 0.15 0.15], 'LineWidth', 1);
xlim(xLim); ylim([0.5 5.5]);
set(ax_hyp, 'YTick', 1:5, 'YTickLabel', {'N3','N2','N1','W','R'}, ...
    'FontSize', 8, 'Box', 'off', 'TickDir', 'out', 'XTick', []);

%%% --- Slope rows ---
for iRange = 1:nRanges
    fr      = FreqRanges{iRange};
    frLabel = sprintf('%d–%d Hz', fr(1), fr(2));
    raw     = SlopesRaw{iRange}(:)';
    clean   = SlopesClean{iRange}(:)';

    ax = nexttile;
    hold on;
    plot(T_min(1:numel(raw)),   raw,   '-', 'Color', colorRaw,   'LineWidth', 1.2, 'DisplayName', 'before GEDAI');
    plot(T_min(1:numel(clean)), clean, '-', 'Color', colorClean, 'LineWidth', 1.5, 'DisplayName', 'after GEDAI');
    xlim(xLim);
    ylabel('1/f exponent', 'FontSize', 9);
    title(frLabel, 'FontSize', 9, 'FontWeight', 'normal');
    set(ax, 'FontSize', 9, 'Box', 'off', 'TickDir', 'out');
    grid on; set(ax, 'GridAlpha', 0.15, 'GridLineStyle', ':');
    if iRange == 1
        legend('Location', 'northeast', 'FontSize', 8, 'Box', 'off');
    end
    if iRange == nRanges
        xlabel('Time (min)', 'FontSize', 9);
    else
        set(ax, 'XTick', []);
    end
end

linkaxes(findobj(fig, 'Type', 'axes'), 'x');

save_fig(fig, opts.SavePath, 'slopes_timecourse');
end
