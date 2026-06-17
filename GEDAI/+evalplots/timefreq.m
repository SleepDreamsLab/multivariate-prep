function timefreq(PwrClean, PwrRaw, Freqs, StageScoring, ChanIdx, opts)
% GEDAI.EVALPLOTS.TIMEFREQ  Time-frequency spectrogram: raw vs. clean.
%
%   gedai.evalplots.timefreq(PwrClean, PwrRaw, Freqs, StageScoring, ChanIdx)
%   gedai.evalplots.timefreq(..., Name, Value)
%
%   Inputs
%   ------
%   PwrClean     : chans x epochs x freqs — GEDAI-cleaned power.
%   PwrRaw       : chans x epochs x freqs — raw power.
%   Freqs        : frequency vector (Hz).
%   StageScoring : per-epoch sleep-stage codes (used for hypnogram strip).
%   ChanIdx      : channel index within PwrClean/PwrRaw (default 1).
%
%   Optional name-value pairs
%   -------------------------
%   ChanLabel       : channel name shown in subplot titles (default 'C3').
%   EpochLength     : epoch duration in seconds (default 30).
%   FreqLim         : [fMin fMax] displayed frequency range (default [0 40]).
%   CLim            : [cLow cHigh] colorbar limits in log10 µV²/Hz (default [-3 3]).
%   ExponentsClean  : 1 × nEpochs aperiodic exponent from FOOOF on clean data.
%   ExponentsRaw    : 1 × nEpochs aperiodic exponent from FOOOF on raw data.
%   SavePath        : base path for saving; suffix + '.png' appended.

arguments
    PwrClean
    PwrRaw
    Freqs        {mustBeVector}
    StageScoring {mustBeVector}
    ChanIdx      (1,1) double = 1
    opts.ChanLabel               = 'C3'
    opts.EpochLength (1,1) double = 30
    opts.FreqLim     (1,2) double = [0 40]
    opts.CLim        (1,2) double = [-3  3]
    opts.ExponentsClean double    = []
    opts.ExponentsRaw   double    = []
    opts.FooofLabel               = ''
    opts.SavePath                = ''
end

nEpochs = size(PwrRaw, 2);
fMask   = Freqs >= opts.FreqLim(1) & Freqs <= opts.FreqLim(2);
T_min   = (0:nEpochs-1) * opts.EpochLength / 60;   % epoch start-times in minutes
xLim    = [0  T_min(end) + opts.EpochLength/60];

P_raw_log   = log10(squeeze(PwrRaw(ChanIdx,   :, fMask))' + eps);   % freqs × epochs
P_clean_log = log10(squeeze(PwrClean(ChanIdx, :, fMask))' + eps);

scoringTF = StageScoring(:)';
if numel(scoringTF) > nEpochs, scoringTF = scoringTF(1:nEpochs); end

expClean = opts.ExponentsClean(:)';
expRaw   = opts.ExponentsRaw(:)';
if numel(expClean) > nEpochs, expClean = expClean(1:nEpochs); end
if numel(expRaw)   > nEpochs, expRaw   = expRaw(1:nEpochs);   end

% Layout: 12 rows × 2 cols
%   Row  1       (1×2)  : hypnogram strip
%   Rows 2-3     (2×2)  : aperiodic exponent time series
%   Rows 4-12    (9×1)  : spectrograms (raw | clean)
figName = 'GEDAI evaluation — time-frequency';
if ~isempty(opts.FooofLabel)
    figName = [figName '  (' opts.FooofLabel ')'];
end
fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 32 16], 'Name', figName);
tiledlayout(14, 2, 'TileSpacing', 'tight', 'Padding', 'compact');

%%% --- Rows 1-3: hypnogram ---
ax_hyp = nexttile(1, [3 2]);

% Map stage codes to y positions: R=5 (top) … N3=1 (bottom)
codes = [-3, -2, -1,  0,  1,  5];
yPos  = [  1,  2,  3,  5,  4,  5];
yVals = nan(size(scoringTF));
for iEp = 1:numel(scoringTF)
    idx = find(codes == scoringTF(iEp), 1);
    if ~isempty(idx), yVals(iEp) = yPos(idx); end
end

% Red shaded background for REM epochs
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

%%% --- Rows 4-7: aperiodic exponent ---
ax_exp = nexttile(7, [4 2]);
hold on;
colorRaw   = [0.65 0.65 0.65];
colorClean = [0.08 0.22 0.55];
if ~isempty(expRaw)
    plot(T_min(1:numel(expRaw)),   expRaw,   '-', 'Color', colorRaw,   'LineWidth', 1.2, 'DisplayName', 'before GEDAI');
end
if ~isempty(expClean)
    plot(T_min(1:numel(expClean)), expClean, '-', 'Color', colorClean, 'LineWidth', 1.5, 'DisplayName', 'after GEDAI');
end
xlim(xLim);
ylabel('1/f exp.', 'FontSize', 9);
if ~isempty(opts.FooofLabel)
    title(ax_exp, opts.FooofLabel, 'FontSize', 9, 'FontWeight', 'normal');
end
set(ax_exp, 'FontSize', 9, 'Box', 'off', 'TickDir', 'out');
xlabel(ax_exp, 'Time (min)', 'FontSize', 9);
legend('Location', 'northeast', 'FontSize', 8, 'Box', 'off');
grid on; set(ax_exp, 'GridAlpha', 0.15, 'GridLineStyle', ':');

%%% --- Rows 8-14: spectrograms ---
nContours = 40;

ax_raw = nexttile(15, [7 1]); hold on;
contourf(T_min, Freqs(fMask), P_raw_log, nContours, 'LineColor', 'none');
clim(opts.CLim);
colormap(ax_raw, custom_cmap());
cb = colorbar; cb.Label.String = 'log_{10}(\muV^2/Hz)'; cb.FontSize = 9;
xlabel('Time (min)', 'FontSize', 10);
ylabel('Frequency (Hz)', 'FontSize', 10);
title(sprintf('before GEDAI  (%s)', opts.ChanLabel), 'FontSize', 11, 'FontWeight', 'bold');
set(ax_raw, 'FontSize', 10, 'Box', 'off', 'TickDir', 'out', 'YDir', 'normal');

ax_clean = nexttile(16, [7 1]); hold on;
contourf(T_min, Freqs(fMask), P_clean_log, nContours, 'LineColor', 'none');
clim(opts.CLim);
colormap(ax_clean, custom_cmap());
cb = colorbar; cb.Label.String = 'log_{10}(\muV^2/Hz)'; cb.FontSize = 9;
xlabel('Time (min)', 'FontSize', 10);
ylabel('Frequency (Hz)', 'FontSize', 10);
title(sprintf('after GEDAI  (%s)', opts.ChanLabel), 'FontSize', 11, 'FontWeight', 'bold');
set(ax_clean, 'FontSize', 10, 'Box', 'off', 'TickDir', 'out', 'YDir', 'normal');

linkaxes([ax_hyp, ax_exp, ax_raw, ax_clean], 'x');

save_fig(fig, opts.SavePath, 'timefreq');
end
