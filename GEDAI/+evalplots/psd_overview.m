function psd_overview(PwrClean, PwrRaw, Freqs, StageScoring, ChanIdx, opts)
% GEDAI.EVAL.PSD_OVERVIEW  Two-tile PSD: raw (left) vs. clean (right), stages as coloured lines.
%
%   gedai.eval.psd_overview(PwrClean, PwrRaw, Freqs, StageScoring, ChanIdx)
%   gedai.eval.psd_overview(..., Name, Value)
%
%   Inputs
%   ------
%   PwrClean     : chans x epochs x freqs — GEDAI-cleaned power.
%   PwrRaw       : chans x epochs x freqs — raw power.
%   Freqs        : frequency vector (Hz).
%   StageScoring : per-epoch sleep-stage codes.
%   ChanIdx      : channel index used for the PSD lines.
%
%   Optional name-value pairs
%   -------------------------
%   FreqLim   : [fMin fMax] (default [0 40]).
%   FreqScale : 'linear', 'log', or 'both' (default).
%   SavePath  : base path for saving; suffix + '.png' appended.

arguments
    PwrClean
    PwrRaw
    Freqs        {mustBeVector}
    StageScoring {mustBeVector}
    ChanIdx      (1,1) double = 1
    opts.FreqLim   (1,2) double = [0 40]
    opts.FreqScale           = 'both'
    opts.SavePath            = ''
end

uniqueStages = unique(StageScoring);
nStages      = numel(uniqueStages);
fMask        = Freqs >= opts.FreqLim(1) & Freqs <= opts.FreqLim(2);

[scaleList, logXLim, pow2ticks] = freq_scale_setup(opts.FreqLim, opts.FreqScale);
nScales = numel(scaleList);

fig = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [2 2 24 8*nScales], ...
    'Name', 'GEDAI evaluation — PSD by sleep stage');
tiledlayout(nScales, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

tileLabels = {'before GEDAI', 'after GEDAI'};
PwrBoth    = {PwrRaw, PwrClean};

for iScale = 1:nScales
    sc = scaleList{iScale};
    for iTile = 1:2
        Pwr = PwrBoth{iTile};
        nexttile; hold on;
        for iStage = 1:nStages
            d    = uniqueStages(iStage);
            mask = StageScoring == d;
            if ~any(mask), continue; end
            meanPwr = squeeze(mean(log10(Pwr(ChanIdx, mask, fMask)), 2));
            plot(Freqs(fMask), meanPwr, '-', 'Color', stage_color(d), ...
                'LineWidth', 1.2, 'DisplayName', stage_name(d));
        end
        if iTile == 1
            ylabel('Power (log_{10} \muV^2/Hz)', 'FontSize', 10);
        end
        if iScale == 1
            title(tileLabels{iTile}, 'FontSize', 11, 'FontWeight', 'bold');
        end
        xlabel('Frequency (Hz)', 'FontSize', 10);
        legend('Location', 'southwest', 'FontSize', 9, 'Box', 'off');
        set(gca, 'FontSize', 10, 'Box', 'off', 'TickDir', 'out');
        grid on; set(gca, 'GridAlpha', 0.15, 'GridLineStyle', ':');
        xlim(opts.FreqLim); ylim([-8 3]);
        apply_freq_scale(gca, sc, logXLim, pow2ticks);
    end
end

save_fig(fig, opts.SavePath, 'psd_overview');
end
