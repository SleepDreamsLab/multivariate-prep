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

    corrTh   = fieldOr(params, 'corrThreshold',        0.7);
    brokenTh = fieldOr(params, 'maxBrokenTime',        0.5);
    noiseTh  = fieldOr(params, 'noiseReportThreshold', 4);
    flatTh   = fieldOr(params, 'flatMaxBrokenTime',    0.5);

    lowcorrprop = sum(corr < corrTh, 2) ./ size(corr, 2);
    nTiles      = 2 + ~isempty(flatprop);

    figure();
    tiledlayout(1, nTiles, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile();
    drawTopo(lowcorrprop, chanlocs, lowcorrprop > brokenTh, [0 .8]);
    title({'Prop. of recording', sprintf('with corr < %g  (marked > %g)', corrTh, brokenTh)});

    nexttile();
    drawTopo(znoise, chanlocs, znoise > noiseTh, [0 10]);
    title({'Line noise', sprintf('(marked z > %g, not removed)', noiseTh)});

    if ~isempty(flatprop)
        nexttile();
        drawTopo(flatprop, chanlocs, flatprop > flatTh, [0 1]);
        title({'Prop. of recording', sprintf('flat  (marked > %g)', flatTh)});
    end

    colormap('gray');
    set(gcf, 'Color', 'w', 'Units', 'centimeters', 'Position', [2 2 10*nTiles 10]);
    print(gcf, savefile, '-dpng', '-r100');
    close
end

% -------------------------------------------------------------------------
function drawTopo(vals, chanlocs, hilite, cl)
% One panel: every electrode as a small rose dot, the ones over threshold in red.
% emarker2 is only added when something is above threshold - topoplot does not accept an
% empty highlight list.
    args = {'numcontour', 0, 'electrodes', 'on', 'emarker', {'.', [1 0.45 0.6], 3, 1}};
    idx  = find(hilite);
    if ~isempty(idx)
        args = [args, {'emarker2', {idx, '.', 'r', 10}}];
    end
    topoplot(vals, chanlocs, args{:});
    colorbar(); caxis(cl);
end

% -------------------------------------------------------------------------
function v = fieldOr(s, f, d)
% Field value, or the default when the caller passed no params (or an older struct).
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
