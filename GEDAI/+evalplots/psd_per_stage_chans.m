function psd_per_stage_chans(PwrClean, PwrRaw, FreqsClean, FreqsRaw, StageScoring, opts)
% GEDAI.EVALPLOTS.PSD_PER_STAGE_CHANS  Per-stage PSD: per-channel epoch-means (transparent)
%   overlaid with grand mean across all channels (thick line), raw vs. clean.
%
%   Inputs
%   ------
%   PwrClean     : chans x freqs x epochs — GEDAI-cleaned power.
%   PwrRaw       : chans x freqs x epochs — raw power.
%   FreqsClean   : frequency vector for PwrClean (Hz).
%   FreqsRaw     : frequency vector for PwrRaw (Hz).
%   StageScoring : per-epoch sleep-stage codes.
%
%   Optional name-value pairs
%   -------------------------
%   FreqLim   : [fMin fMax] (default [0 40]).
%   FreqScale : 'linear', 'log', or 'both' (default).
%   SavePath  : base path for saving; suffix + '.png' appended.

arguments
    PwrClean
    PwrRaw
    FreqsClean   {mustBeVector}
    FreqsRaw     {mustBeVector}
    StageScoring {mustBeVector}
    opts.FreqLim   (1,2) double = [0 40]
    opts.FreqScale             = 'both'
    opts.SavePath              = ''
end

uniqueStages = unique(StageScoring);
nStages      = numel(uniqueStages);
fMaskClean   = FreqsClean >= opts.FreqLim(1) & FreqsClean <= opts.FreqLim(2);
fMaskRaw     = FreqsRaw   >= opts.FreqLim(1) & FreqsRaw   <= opts.FreqLim(2);

[scaleList, logXLim, pow2ticks] = freq_scale_setup(opts.FreqLim, opts.FreqScale);
nScales = numel(scaleList);

fig = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [2 2 50 8*nScales], ...
    'Name', 'GEDAI evaluation — channel PSD per sleep stage');
tiledlayout(nScales, nStages, 'TileSpacing', 'compact', 'Padding', 'compact');

for iScale = 1:nScales
    sc = scaleList{iScale};
    for iStage = 1:nStages
        d    = uniqueStages(iStage);
        mask = StageScoring == d;
        nexttile; hold on;
        if ~any(mask), continue; end

        % Per-channel epoch-averages: [nChans x nFreqSel]
        chanMeansRaw   = squeeze(mean(log10(PwrRaw(:,   fMaskRaw,   mask) + eps), 3));
        chanMeansClean = squeeze(mean(log10(PwrClean(:, fMaskClean, mask) + eps), 3));
        if isvector(chanMeansRaw),   chanMeansRaw   = chanMeansRaw(:)';   end
        if isvector(chanMeansClean), chanMeansClean = chanMeansClean(:)'; end

        % Grand mean across channels: [1 x nFreqSel]
        grandMeanRaw   = mean(chanMeansRaw,   1);
        grandMeanClean = mean(chanMeansClean, 1);

        % Per-channel lines (transparent)
        plot(FreqsRaw(fMaskRaw),     chanMeansRaw',   '-', 'Color', [0.75 0.75 0.75 0.25], ...
            'LineWidth', 0.5, 'HandleVisibility', 'off');
        plot(FreqsClean(fMaskClean), chanMeansClean', '-', 'Color', [stage_color(d), 0.25], ...
            'LineWidth', 0.5, 'HandleVisibility', 'off');

        % Grand mean (thick)
        plot(FreqsRaw(fMaskRaw),     grandMeanRaw,   '-', 'Color', [0.45 0.45 0.45], ...
            'LineWidth', 1.6, 'DisplayName', 'before');
        plot(FreqsClean(fMaskClean), grandMeanClean, '-', 'Color', 'k', ...
            'LineWidth', 2.0, 'DisplayName', 'after');

        if iStage == 1
            ylabel('Power (log_{10} \muV^2/Hz)', 'FontSize', 10);
        end
        if iScale == 1
            title(sprintf('%s  (n=%d)', stage_name(d), sum(mask)), ...
                'FontSize', 11, 'FontWeight', 'bold');
        end
        xlabel('Frequency (Hz)', 'FontSize', 10);
        legend('Location', 'northeast', 'FontSize', 9, 'Box', 'off');
        set(gca, 'FontSize', 10, 'Box', 'off', 'TickDir', 'out');
        grid on; set(gca, 'GridAlpha', 0.15, 'GridLineStyle', ':');
        xlim(opts.FreqLim); ylim([-8 4]);
        apply_freq_scale(gca, sc, logXLim, pow2ticks);
    end
end

save_fig(fig, opts.SavePath, 'psd_per_stage_chans');
end
