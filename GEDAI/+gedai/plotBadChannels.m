function plotBadChannels(corr, znoise, chanlocs, savefile, flatprop, params, savechans)
% PLOTBADCHANNELS  Topoplots of the bad channel criteria.
%   plotBadChannels(corr, znoise, chanlocs, savefile)
%   plotBadChannels(corr, znoise, chanlocs, savefile, flatprop)
%   plotBadChannels(corr, znoise, chanlocs, savefile, flatprop, params)
%   plotBadChannels(corr, znoise, chanlocs, savefile, flatprop, params, savechans)
%
%   Every electrode is drawn as a small rose dot; the ones a criterion marks are
%   overdrawn in red. Pass flatprop (from gedai.detectFlatChannels) to add the flat-line
%   panel, and params (EEG.etc.badchans.params) so the markers use the thresholds
%   actually in force - a marker that disagrees with the mask is worse than no marker.
%   Without params the historical defaults are used, so old callers still work.
%
%   savechans   channels the noise criterion flagged but that were spared (see the
%               module help of bidsfun_detect_badchans). Overdrawn in green, on the
%               noise panel only, instead of red - they no longer count as bad, and the
%               override should be visible rather than silent.
    if nargin < 5, flatprop  = []; end
    if nargin < 6, params    = struct(); end
    if nargin < 7, savechans = []; end

    %%% All three panels mark what was actually removed, so every red dot on this figure
    %%% corresponds to a channel missing from the data. The line-noise panel keys on
    %%% noiseThreshold, the same threshold clean_channels used for detection - if a caller
    %%% ever sets it to Inf the criterion removes nothing and the panel is left unmarked,
    %%% but znoise itself still lands in the .mat and channels.tsv either way, so a
    %%% pathological channel remains visible even when unmarked.
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
        drawTopoPanel(flatprop, chanlocs, flatprop > flatTh, [0 1], 'Prop. of recording flat');
        title({'Prop. of recording', sprintf('flat  (marked > %g)', flatTh)});
    end

    nexttile();
    drawTopoPanel(lowcorrprop, chanlocs, lowcorrprop > brokenTh, [0 .8], 'Prop. of recording low-corr');
    title({'Prop. of recording', sprintf('with corr < %g  (marked > %g)', corrTh, brokenTh)});

    nexttile();
    %%% Channels in savechans still exceed noiseTh - that is why they were flagged in
    %%% the first place - but they no longer count as bad (see the module help of
    %%% bidsfun_detect_badchans), so they are pulled out of the red set and shown in
    %%% green instead.
    noiseFail = znoise > noiseTh;
    savedMask = false(size(noiseFail));
    if ~isempty(savechans)
        savedMask(savechans) = true;
        noiseFail(savechans) = false;
    end
    drawTopoPanel(znoise, chanlocs, noiseFail, [0 5], 'Noise z-score (robust)', savedMask);
    if isfinite(noiseTh)
        title({'Line noise', sprintf('(removed at z > %g, green = recovered)', noiseTh)});
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
