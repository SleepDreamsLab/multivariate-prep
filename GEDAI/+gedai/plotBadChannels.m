function plotBadChannels(corr, znoise, chanlocs, savefile, flatprop)
% PLOTBADCHANNELS  Topoplots of the bad channel criteria.
%   plotBadChannels(corr, znoise, chanlocs, savefile)
%   plotBadChannels(corr, znoise, chanlocs, savefile, flatprop)
%
%   Red dots mark the channels each criterion flags. Pass flatprop (from
%   gedai.detectFlatChannels) to add the flat-line panel; omit it for two panels.
    if nargin < 5, flatprop = []; end

    lowcorrprop = sum(corr < .7, 2) ./ size(corr, 2);
    nTiles = 2 + ~isempty(flatprop);

    figure();
    tiledlayout(1, nTiles, 'TileSpacing', 'compact', 'Padding', 'compact');
    nexttile();
    topoplot(lowcorrprop, chanlocs, 'numcontour', 0, 'emarker2', {find(lowcorrprop > .5), '.', 'r', 10});
    colorbar(); caxis([0 .8]); title({'Prop. of recording', 'with low correlation'});
    nexttile();
    topoplot(znoise, chanlocs, 'numcontour', 0, 'emarker2', {find(znoise > 4), '.', 'r', 10});
    colorbar(); caxis([0 4]); title('Line noise');
    if ~isempty(flatprop)
        nexttile();
        topoplot(flatprop, chanlocs, 'numcontour', 0, 'emarker2', {find(flatprop > .5), '.', 'r', 10});
        colorbar(); caxis([0 1]); title({'Prop. of recording', 'flat'});
    end
    colormap('gray');
    set(gcf, 'Color', 'w', 'Units', 'centimeters', 'Position', [2 2 10*nTiles 10]);
    print(gcf, savefile, '-dpng', '-r100');
    close
end
