function psd_per_stage(PwrClean, PwrRaw, FreqsClean, FreqsRaw, StageScoring, ChanIdx, opts)
% GEDAI.EVAL.PSD_PER_STAGE  Clean vs. raw mean PSD plotted per sleep stage.
%
%   gedai.eval.psd_per_stage(PwrClean, PwrRaw, FreqsClean, FreqsRaw, StageScoring, ChanIdx)
%   gedai.eval.psd_per_stage(..., Name, Value)
%
%   Inputs
%   ------
%   PwrClean     : chans x epochs x freqs — GEDAI-cleaned power.
%   PwrRaw       : chans x epochs x freqs — raw power.
%   FreqsClean   : frequency vector for PwrClean (Hz).
%   FreqsRaw     : frequency vector for PwrRaw (Hz).
%   StageScoring : per-epoch sleep-stage codes.
%   ChanIdx      : channel index (within PwrClean/PwrRaw) used for the PSD line.
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
    ChanIdx      (1,1) double = 1
    opts.FreqLim   (1,2) double = [0 40]
    opts.FreqScale           = 'both'
    opts.SavePath            = ''
end

uniqueStages = unique(StageScoring);
nStages      = numel(uniqueStages);
fMaskClean   = FreqsClean >= opts.FreqLim(1) & FreqsClean <= opts.FreqLim(2);
fMaskRaw     = FreqsRaw   >= opts.FreqLim(1) & FreqsRaw   <= opts.FreqLim(2);

[scaleList, logXLim, pow2ticks] = freq_scale_setup(opts.FreqLim, opts.FreqScale);
nScales = numel(scaleList);

fig = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [2 2 50 8*nScales], ...
    'Name', 'GEDAI evaluation — power per sleep stage');
tiledlayout(nScales, nStages, 'TileSpacing', 'compact', 'Padding', 'compact');

for iScale = 1:nScales
    sc = scaleList{iScale};
    for iStage = 1:nStages
        d    = uniqueStages(iStage);
        mask = StageScoring == d;
        nexttile; hold on;
        if ~any(mask), continue; end

        meanClean = squeeze(mean(log10(PwrClean(ChanIdx, mask, fMaskClean)), 2));
        meanRaw   = squeeze(mean(log10(PwrRaw(ChanIdx,   mask, fMaskRaw)),   2));

        % Individual epoch lines — plotted first so means sit on top
        epochsRaw   = squeeze(log10(PwrRaw(ChanIdx,   mask, fMaskRaw)));
        epochsClean = squeeze(log10(PwrClean(ChanIdx, mask, fMaskClean)));
        if isvector(epochsRaw),   epochsRaw   = epochsRaw(:)';   end
        if isvector(epochsClean), epochsClean = epochsClean(:)'; end
        plot(FreqsRaw(fMaskRaw),     epochsRaw',   '-', 'Color', [0.65 0.65 0.65 0.3], ...
            'LineWidth', 0.4, 'HandleVisibility', 'off');
        plot(FreqsClean(fMaskClean), epochsClean', '-', 'Color', [stage_color(d), 0.3], ...
            'LineWidth', 0.4, 'HandleVisibility', 'off');

        plot(FreqsRaw(fMaskRaw),     meanRaw,   '-', 'Color', [0.65 0.65 0.65], ...
            'LineWidth', 1.6, 'DisplayName', 'before GEDAI');
        plot(FreqsClean(fMaskClean), meanClean, '-',  'Color', 'k', ...
            'LineWidth', 2.0, 'DisplayName', 'after GEDAI');

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

save_fig(fig, opts.SavePath, 'psd_per_stage');
end
