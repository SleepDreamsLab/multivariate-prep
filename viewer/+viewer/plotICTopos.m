function fig = plotICTopos(winv, chanlocs, comps, titles, figName, opts)
% PLOTICTOPOS  Grid of IC scalp topographies, one tile per component, tinted
% by ICLabel category and each carrying a good/bad toggle.
%
%   fig = viewer.plotICTopos(winv, chanlocs, comps, titles, figName, opts)
%
%   winv      nChan x nComp mixing matrix (EEG.icawinv)
%   chanlocs  chanlocs entries matching winv's ROWS (i.e. EEG.chanlocs(icachansind))
%   comps     component indices to draw
%   titles    per-component title strings, indexed by component number
%   figName   window title
%   opts      .classOf     per-component ICLabel class index (0 = unknown)
%             .classNames  cellstr of class names
%             .isGood      per-component logical override for the initial state
%             .onChange    handle called with the full good/bad vector after
%                          every toggle, so the caller -- not this window --
%                          owns the verdict and it survives closing the window
%
% Each tile is tinted and outlined in its category's colour, so a wall of
% topographies can be read by category at a glance instead of tile by tile
% through the titles. The tint sits on the enclosing panel rather than on the
% axes, because topoplot repaints axes backgrounds itself.
%
% The good/bad state starts from ICLabel -- Brain is good, everything else
% bad -- and each tile's toggle flips it (green = good, red = bad). The
% button's caption IS the good/bad label, so it is not repeated in the title.
% With no ICLabel classification available every component starts GOOD:
% nothing is known to be artefactual, and defaulting to bad would amount to
% recommending you discard a decomposition nobody has judged yet.
%
% The current state lives on the figure and can be read back with
%   getappdata(fig, 'icIsGood')   % logical, indexed by component number
%   getappdata(fig, 'icComps')    % the components drawn
if nargin < 5 || isempty(figName), figName = 'IC topographies'; end
if nargin < 6, opts = struct(); end
if ~isfield(opts, 'classOf'),    opts.classOf    = []; end
if ~isfield(opts, 'classNames'), opts.classNames = {}; end
if ~isfield(opts, 'isGood'),     opts.isGood     = []; end
if ~isfield(opts, 'onChange'),   opts.onChange   = []; end
if isempty(comps), fig = gobjects(0); return; end
if exist('topoplot', 'file') ~= 2
    error('plotICTopos:noTopoplot', ...
        'EEGLAB''s topoplot is not on the path -- cannot draw IC topographies.');
end

comps = comps(:).';

% Channels without usable coordinates are dropped rather than handed to
% topoplot: the decomposition can legitimately include derived channels (the
% bipolar EOG pair chans1020 appends, which viewer.extendICAToDerived folds
% into the mixing matrix) that have no scalp position at all.
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

%%% --- Initial good/bad state, indexed by component number ---
% Sized to span the caller's whole vector, not just the components drawn, so
% onChange hands back something the caller can assign straight over its own
% state without the components it did not draw silently reverting.
nAll   = max([max(comps), numel(opts.isGood), numel(opts.classOf)]);
isGood = true(1, nAll);
if ~isempty(opts.isGood)
    n = numel(opts.isGood);
    isGood(1:n) = logical(opts.isGood(:).');
elseif ~isempty(opts.classNames) && ~isempty(opts.classOf)
    brainCls = find(strcmpi(opts.classNames, 'Brain'), 1);
    if ~isempty(brainCls)
        cls = zeros(1, nAll);
        n   = min(numel(opts.classOf), numel(cls));
        cls(1:n) = opts.classOf(1:n);
        isGood  = cls == brainCls;
    end
end

goodCol = [0.42 0.78 0.45];
badCol  = [0.90 0.36 0.36];

nComp = numel(comps);
nCol  = ceil(sqrt(nComp));
nRow  = ceil(nComp / nCol);

fig = figure('Name', figName, 'NumberTitle', 'off', 'Color', 'w', ...
    'Units', 'normalized', 'Position', [0.15 0.15 0.68 0.72]);

padX = 0.006; padY = 0.008;
tileW = (1 - padX*(nCol+1)) / nCol;
tileH = (1 - padY*(nRow+1)) / nRow;
btns  = gobjects(1, nComp);

for i = 1:nComp
    c   = comps(i);
    row = ceil(i / nCol);
    col = i - (row-1)*nCol;
    x   = padX + (col-1)*(tileW + padX);
    y   = 1 - row*(tileH + padY);

    [tint, edge] = categoryColour(classNameOf(c, opts));
    pan = uipanel(fig, 'Units', 'normalized', 'Position', [x y tileW tileH], ...
        'BackgroundColor', tint, 'BorderType', 'line', ...
        'HighlightColor', edge, 'BorderWidth', 2);

    ax  = axes('Parent', pan, 'Units', 'normalized', 'Position', [0.06 0.30 0.88 0.54]);
    map = double(winv(hasLoc, c));
    topoplot(map, chanlocs(hasLoc), 'electrodes', 'off');
    % Symmetric colour limits: an IC map's sign is arbitrary, so anything not
    % centred on zero would imply a polarity the decomposition does not carry.
    lim = max(abs(map));
    if lim > 0, clim(ax, [-lim lim]); end
    set(ax, 'Color', tint);

    if numel(titles) >= c && ~isempty(titles{c})
        ttl = titles{c};
    else
        ttl = sprintf('IC %d', c);
    end
    title(ax, ttl, 'Interpreter', 'none', 'FontSize', 9, 'FontWeight', 'bold');

    btns(i) = uicontrol(pan, 'Style', 'togglebutton', 'Units', 'normalized', ...
        'Position', [0.20 0.05 0.60 0.17], 'FontSize', 8, 'FontWeight', 'bold', ...
        'Value', isGood(c), 'Callback', @(o,e) onToggle(i, c));
    paintButton(i, c);
end

colormap(fig, divergingMap());
set(fig, 'Color', 'w');   % topoplot repaints the figure background
publish();


    function onToggle(i, c)
        isGood(c) = ~isGood(c);
        set(btns(i), 'Value', isGood(c));
        paintButton(i, c);
        publish();
    end

    function paintButton(i, c)
        if isGood(c)
            set(btns(i), 'String', 'good', 'BackgroundColor', goodCol, ...
                'ForegroundColor', [0 0 0]);
        else
            set(btns(i), 'String', 'bad', 'BackgroundColor', badCol, ...
                'ForegroundColor', [1 1 1]);
        end
    end

    function publish()
        setappdata(fig, 'icIsGood', isGood);
        setappdata(fig, 'icComps',  comps);
        if ~isempty(opts.onChange)
            opts.onChange(isGood);
        end
    end
end


function name = classNameOf(c, opts)
% ICLabel class name for one component, '' when unclassified.
name = '';
if isempty(opts.classNames) || isempty(opts.classOf) || c > numel(opts.classOf)
    return
end
k = opts.classOf(c);
if k >= 1 && k <= numel(opts.classNames)
    name = opts.classNames{k};
end
end


function [tint, edge] = categoryColour(name)
% One hue per ICLabel category: a saturated edge for the tile outline and a
% heavily whitened version of it for the fill, so the tint identifies the
% category without competing with the topography's own blue/red scale.
switch lower(strtrim(name))
    case 'brain',         edge = [ 60 160  75];
    case 'muscle',        edge = [214  90  40];
    case 'eye',           edge = [ 55 110 200];
    case 'heart',         edge = [200  60 140];
    case 'line noise',    edge = [185 155  30];
    case 'channel noise', edge = [130  95  70];
    case 'other',         edge = [140 140 140];
    otherwise,            edge = [175 175 175];   % unclassified
end
edge = edge / 255;
tint = 1 - (1 - edge) * 0.20;
end


function m = divergingMap()
% Blue-white-red, so zero reads as neutral. Built inline rather than pulled
% from a colormap toolbox that may not be installed alongside EEGLAB.
t = linspace(-1, 1, 256)';
m = [min(1, 1 + t), 1 - abs(t), min(1, 1 - t)] .^ 0.85;
end
