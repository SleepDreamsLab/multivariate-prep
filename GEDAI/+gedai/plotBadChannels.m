function plotBadChannels(corr, znoise, chanlocs, savefile)
    lowcorrprop = sum(corr < .7, 2) ./ size(corr, 2);
    figure();
    tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    nexttile();
    topoplot(lowcorrprop, chanlocs, 'numcontour', 0, 'emarker2', {find(lowcorrprop > .5), '.', 'r', 10});
    colorbar(); caxis([0 .8]); title({'Prop. of recording', 'with low correlation'});
    nexttile();
    topoplot(znoise, chanlocs, 'numcontour', 0, 'emarker2', {find(znoise > 4), '.', 'r', 10});
    colorbar(); caxis([0 4]); title('Line noise');
    colormap('gray');
    set(gcf, 'Color', 'w', 'Units', 'centimeters', 'Position', [2 2 20 10]);
    print(gcf, savefile, '-dpng', '-r150');
    close
end
