function drawTopoPanel(vals, chanlocs, hilite, cl, cbLabel, hilite2)
% DRAWTOPOPANEL  One bad-channel topoplot panel: every electrode as a small rose dot,
%   the ones over threshold overdrawn in red, and an optional second set overdrawn in
%   green.
%
%   vals      one value per entry of chanlocs. NaN drops that electrode from both the
%             interpolation and the dots, which is how a channel removed by an earlier
%             round is shown as simply absent.
%   hilite    logical over the same channels, or indices; marked in red.
%   cl        colour limits.
%   cbLabel   colorbar label; '' (default) leaves the colorbar unlabelled.
%   hilite2   logical over the same channels, or indices; marked in green instead of
%             red - e.g. channels a later criterion recovers (default: none).
%
%   emarker2 is only added when something is above threshold - topoplot does not accept
%   an empty highlight list. Its indices are matched against the channels that survived
%   the NaN filter, so full-montage indices are safe to pass.
%
%   hilite2 is drawn by hand rather than through topoplot, which only supports one
%   highlighted subset per call (emarker2 is already spent on hilite). The position
%   formula - pol2cart on (theta, radius), no rotation - is topoplot's own for its
%   default nosedir ('+X'; see topoplot.m); it is reused here, not re-derived, so the
%   green dots land exactly on the electrodes topoplot itself would draw.
%
%   Th/Rd come from readlocs(chanlocs), not a direct [chanlocs.theta] extraction: any
%   channel with an empty theta/radius (no location) is silently DROPPED by a
%   comma-list concatenation, shifting every later index - the green dots ended up off
%   by exactly that amount. readlocs() instead NaN-pads such channels, keeping every
%   index aligned with chanlocs, which is what topoplot's own internal call to it does.
%
%   Private to +gedai so plotBadChannels and plotLineNoiseZ cannot drift apart.

if nargin < 5, cbLabel = ''; end
if nargin < 6, hilite2 = []; end

args = {'numcontour', 0, 'electrodes', 'on', 'emarker', {'.', [1 0.45 0.6], 3, 1}};
idx  = find(hilite);
if ~isempty(idx)
    args = [args, {'emarker2', {idx, '.', 'r', 10}}];
end
topoplot(vals, chanlocs, args{:});
cb = colorbar(); caxis(cl);   %#ok<CAXIS> - clim() is R2022a+, caxis works everywhere
if ~isempty(cbLabel)
    cb.Label.String = cbLabel;
end

idx2 = find(hilite2);
if ~isempty(idx2)
    [~, ~, Th, Rd] = readlocs(chanlocs);
    Th = pi/180 * Th;
    [gx, gy] = pol2cart(Th, Rd);
    hold on
    plot(gy(idx2), gx(idx2), '.', 'Color', [0 0.6 0], 'MarkerSize', 10);
    hold off
end
end
