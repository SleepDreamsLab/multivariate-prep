function drawTopoPanel(vals, chanlocs, hilite, cl)
% DRAWTOPOPANEL  One bad-channel topoplot panel: every electrode as a small rose dot,
%   the ones over threshold overdrawn in red.
%
%   vals      one value per entry of chanlocs. NaN drops that electrode from both the
%             interpolation and the dots, which is how a channel removed by an earlier
%             round is shown as simply absent.
%   hilite    logical over the same channels, or indices; marked in red.
%   cl        colour limits.
%
%   emarker2 is only added when something is above threshold - topoplot does not accept
%   an empty highlight list. Its indices are matched against the channels that survived
%   the NaN filter, so full-montage indices are safe to pass.
%
%   Private to +gedai so plotBadChannels and plotLineNoiseZ cannot drift apart.

args = {'numcontour', 0, 'electrodes', 'on', 'emarker', {'.', [1 0.45 0.6], 3, 1}};
idx  = find(hilite);
if ~isempty(idx)
    args = [args, {'emarker2', {idx, '.', 'r', 10}}];
end
topoplot(vals, chanlocs, args{:});
colorbar(); caxis(cl);   %#ok<CAXIS> - clim() is R2022a+, caxis works everywhere
end
