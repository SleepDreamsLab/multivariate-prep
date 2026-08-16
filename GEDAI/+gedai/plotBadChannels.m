function plotBadChannels(corr, znoise, chanlocs, savefile, flatprop, params)
% PLOTBADCHANNELS  Topoplots of the bad channel criteria.
%   plotBadChannels(corr, znoise, chanlocs, savefile)
%   plotBadChannels(corr, znoise, chanlocs, savefile, flatprop)
%   plotBadChannels(corr, znoise, chanlocs, savefile, flatprop, params)
%
%   Every electrode is drawn as a small rose dot; the ones a criterion marks are
%   overdrawn in red. Pass flatprop (from gedai.detectFlatChannels) to add the flat-line
%   panel, and params (EEG.etc.badchans.params) so the markers use the thresholds
%   actually in force - a marker that disagrees with the mask is worse than no marker.
%   Without params the historical defaults are used, so old callers still work.
    if nargin < 5, flatprop = []; end
    if nargin < 6, params   = struct(); end

    %%% All three panels mark what was actually removed, so every red dot on this figure
    %%% corresponds to a channel missing from the data. That means the line-noise panel
    %%% keys on noiseThreshold, not on the reporting threshold: with noiseThreshold = Inf
    %%% the criterion removes nothing and the panel is deliberately unmarked. The
    %%% channels.tsv still reports which channels are line-noisy, via
    %%% noiseReportThreshold - that is the place to look for them.
    corrTh   = fieldOr(params, 'corrThreshold',     0.7);
    brokenTh = fieldOr(params, 'maxBrokenTime',     0.5);
    noiseTh  = fieldOr(params, 'noiseThreshold',    4);
    flatTh   = fieldOr(params, 'flatMaxBrokenTime', 0.5);

    lowcorrprop = sum(corr < corrTh, 2) ./ size(corr, 2);
    nTiles      = 2 + ~isempty(flatprop);

    figure();
    tiledlayout(1, nTiles, 'TileSpacing', 'compact', 'Padding', 'compact');

    %%% Panels follow the order the criteria are decided in: flat first (a dead channel
    %%% is invisible to the other two), then correlation, then line noise.
    if ~isempty(flatprop)
        nexttile();
        drawTopoPanel(flatprop, chanlocs, flatprop > flatTh, [0 1]);
        title({'Prop. of recording', sprintf('flat  (marked > %g)', flatTh)});
    end

    nexttile();
    drawTopoPanel(lowcorrprop, chanlocs, lowcorrprop > brokenTh, [0 .8]);
    title({'Prop. of recording', sprintf('with corr < %g  (marked > %g)', corrTh, brokenTh)});

    nexttile();
    drawTopoPanel(znoise, chanlocs, znoise > noiseTh, [0 10]);
    if isfinite(noiseTh)
        title({'Line noise', sprintf('(removed at z > %g)', noiseTh)});
    else
        title({'Line noise', '(reported only, removes nothing)'});
    end

    colormap('gray');
    set(gcf, 'Color', 'w', 'Units', 'centimeters', 'Position', [2 2 10*nTiles 10]);
    gedai.printFigure(gcf, savefile);
    close
end

% -------------------------------------------------------------------------
function v = fieldOr(s, f, d)
% Field value, or the default when the caller passed no params (or an older struct).
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
