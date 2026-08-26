function fig = plotICTopos(winv, chanlocs, comps, titles, figName)
% PLOTICTOPOS  Grid of IC scalp topographies, one tile per component.
%
%   fig = viewer.plotICTopos(winv, chanlocs, comps, titles, figName)
%
%   winv      nChan x nComp mixing matrix (EEG.icawinv)
%   chanlocs  chanlocs entries matching winv's ROWS (i.e. EEG.chanlocs(icachansind))
%   comps     component indices to draw
%   titles    per-component title strings, indexed by component number
%   figName   window title
%
% Channels without usable coordinates are dropped before plotting rather
% than handed to topoplot: the viewer's decomposition can legitimately
% include derived channels (the bipolar EOG pair chans1020 appends, which
% viewer.extendICAToDerived folds into the mixing matrix) that have no
% scalp position at all, and topoplot cannot place them.
if nargin < 5 || isempty(figName), figName = 'IC topographies'; end
if isempty(comps), fig = gobjects(0); return; end
if exist('topoplot', 'file') ~= 2
    error('plotICTopos:noTopoplot', ...
        'EEGLAB''s topoplot is not on the path -- cannot draw IC topographies.');
end

hasLoc = false(1, numel(chanlocs));
for i = 1:numel(chanlocs)
    v = [chanlocs(i).X, chanlocs(i).Y, chanlocs(i).Z];
    hasLoc(i) = numel(v) == 3 && all(isfinite(v)) && any(v ~= 0);
end
if ~any(hasLoc)
    error('plotICTopos:noChanLocs', ...
        'None of the %d decomposition channels carry scalp coordinates.', numel(chanlocs));
end
if ~all(hasLoc)
    warning('plotICTopos:someChansUnplaced', ...
        'Omitting %d channel(s) with no scalp coordinates from the topographies: %s.', ...
        sum(~hasLoc), strjoin({chanlocs(~hasLoc).labels}, ', '));
end

nComp = numel(comps);
nCol  = ceil(sqrt(nComp));
nRow  = ceil(nComp / nCol);

fig = figure('Name', figName, 'NumberTitle', 'off', 'Color', 'w', ...
    'Units', 'normalized', 'Position', [0.20 0.20 0.60 0.65]);

for i = 1:nComp
    c  = comps(i);
    ax = subplot(nRow, nCol, i, 'Parent', fig);
    map = double(winv(hasLoc, c));
    topoplot(map, chanlocs(hasLoc), 'electrodes', 'off');
    % Symmetric colour limits: an IC map's sign is arbitrary, so anything
    % not centred on zero would imply a polarity the decomposition does
    % not actually carry.
    lim = max(abs(map));
    if lim > 0, clim(ax, [-lim lim]); end
    if numel(titles) >= c && ~isempty(titles{c})
        ttl = titles{c};
    else
        ttl = sprintf('IC %d', c);
    end
    title(ax, ttl, 'Interpreter', 'none', 'FontSize', 9);
end

colormap(fig, divergingMap());
set(fig, 'Color', 'w');   % topoplot repaints the figure background
end


function m = divergingMap()
% Blue-white-red, so zero reads as neutral. Built inline rather than pulled
% from a colormap toolbox that may not be installed alongside EEGLAB.
t = linspace(-1, 1, 256)';
m = [min(1, 1 + t), 1 - abs(t), min(1, 1 - t)] .^ 0.85;
end
