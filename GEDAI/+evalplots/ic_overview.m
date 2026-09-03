function fig = ic_overview(EEG, opts)
% EVALPLOTS.IC_OVERVIEW  One page of independent components: topography, PSD, spectrogram.
%
%   evalplots.ic_overview(EEG)
%   fig = evalplots.ic_overview(EEG, Name, Value, ...)
%
%   Draws the first 40 components of an attached decomposition on a single page,
%   eight across and five down. Each component gets a 3x3 block of its own:
%
%       top left, 2x2     scalp topography of its mixing weights (icawinv column)
%       top right, 2x1    its power spectrum, log10(power) against linear frequency
%       bottom, 1x3       its time-frequency map, over the whole recording
%
%   and a title reading "IC #7 - Muscle", coloured by the ICLabel class named in
%   it. The colours are ica.plot.amica_iclabel_bars' palette, so a class keeps
%   the same colour wherever it is drawn. A decomposition with no ICLabel
%   results still plots; its titles are the bare "IC #7" in neutral grey.
%
%   What the panels are for. The topography says where a component sits and the
%   spectrum says what it is made of - between them that is the ICLabel decision
%   restated in a form a human can argue with. The spectrogram is the part
%   neither of the other two can give: whether the component is a property of
%   the whole night or of ten minutes of it. A muscle component that is only
%   ever active while the subject is awake and a "brain" component that turns on
%   halfway through and never stops look identical in a topography, and quite
%   different here.
%
%   Time-frequency scaling. Component activations are in arbitrary units - ICA
%   fixes the scale of neither the weights nor the activations - so the numbers
%   are never meaningful on their own; only shapes and relative changes are. The
%   spectrum line is the per-window power averaged across time and then
%   log10-transformed; the spectrogram is log10-transformed pixel by pixel, with
%   no time-averaging and no normalisation against a per-component reference
%   such as each row's own median - so a component that is simply louder than
%   another still looks louder here, rather than being rescaled to hide it.
%   Frequency stays on a linear axis throughout; only power is logged.
%
%   Required
%   --------
%   EEG        EEGLAB struct carrying icaweights, icasphere, icawinv,
%              icachansind and chanlocs. ICLabel results are read from
%              EEG.etc.ic_classification.ICLabel when present.
%
%   Name-Value
%   ----------
%   Components  Component indices to draw. Default: the first GridSize(1) *
%               GridSize(2) of them, in order.
%   GridSize    [rows cols] of components per page. Default: [5 8].
%   Window      Welch/spectrogram window length in seconds. Default: 4.
%   FreqLim     [lo hi] frequency range of both the PSD and the spectrogram.
%               Default: [0.5 40].
%   CLim        [lo hi] spectrogram colour limits, in log10(power). [] (default)
%               auto-scales from the data: the 1st and 99th percentile of every
%               finite log10 power value on the page, rounded outward to one
%               decimal - fixed limits would need a different pair for every
%               decomposition, since the values are in arbitrary ICA units.
%   MaxCols     Upper bound on spectrogram columns per component. Default: 600.
%   Title       Page title. Default: 'Independent components'.
%   SavePath    Base path for saving; '_ic_overview.png' appended. '' = no save.
%   Suffix      Filename suffix. Default: 'ic_overview'.
%
%   Returns the figure handle, so a caller running a batch can close it before
%   the next recording starts allocating.

arguments
    EEG struct
    opts.Components (1,:) double = []
    opts.GridSize   (1,2) double {mustBePositive} = [5 8]
    opts.Window     (1,1) double {mustBePositive} = 4
    opts.FreqLim    (1,2) double = [0.5 40]
    opts.CLim       (1,:) double {mustBeNumeric} = []
    opts.MaxCols    (1,1) double {mustBePositive} = 600
    opts.Title            char   = 'Independent components'
    opts.SavePath         char   = ''
    opts.Suffix           char   = 'ic_overview'
end

if ~isempty(opts.CLim) && numel(opts.CLim) ~= 2
    error('evalplots:ic_overview:badCLim', 'CLim must be [] or a [lo hi] pair.');
end

%%% --- What there is to plot ---
if ~isfield(EEG, 'icaweights') || isempty(EEG.icaweights) ...
        || ~isfield(EEG, 'icawinv') || isempty(EEG.icawinv)
    error('evalplots:ic_overview:noICA', ...
        'This EEG carries no decomposition - nothing to plot.');
end
if ~isfield(EEG, 'icachansind') || isempty(EEG.icachansind)
    EEG.icachansind = 1:size(EEG.icasphere, 2);
end

nComp   = size(EEG.icaweights, 1);
nRows   = round(opts.GridSize(1));
nCols   = round(opts.GridSize(2));
comps   = opts.Components;
if isempty(comps), comps = 1:min(nComp, nRows * nCols); end
comps   = comps(comps >= 1 & comps <= nComp);
comps   = comps(1:min(numel(comps), nRows * nCols));
if isempty(comps)
    error('evalplots:ic_overview:noComponents', 'No component indices left to plot.');
end
%%% Only the rows that have something in them. A decomposition of eleven
%%% components on a five-row grid would otherwise be two rows of tiles above
%%% three rows of blank page.
nRows = max(1, ceil(numel(comps) / nCols));

%%% --- ICLabel class per component ---
%%% Top class only. The probability is deliberately not in the title: forty
%%% titles is already a lot of text at this size, and a class the classifier was
%%% unsure about is not something a number in 6 pt type is going to convey.
classIdx  = nan(1, nComp);
className = repmat({''}, 1, nComp);
if isfield(EEG, 'etc') && isfield(EEG.etc, 'ic_classification') ...
        && isfield(EEG.etc.ic_classification, 'ICLabel')
    L = EEG.etc.ic_classification.ICLabel;
    if isfield(L, 'classifications') && ~isempty(L.classifications)
        nHave = min(nComp, size(L.classifications, 1));
        [~, top] = max(L.classifications(1:nHave, :), [], 2);
        classIdx(1:nHave) = top(:)';
        for k = 1:nHave
            if top(k) <= numel(L.classes), className{k} = L.classes{top(k)}; end
        end
    end
end

%%% --- Component activations ---
%%% Computed here rather than leaning on EEG.icaact, which loadica deliberately
%%% leaves empty: materialising every component of a full night costs several GB
%%% and only forty of them are wanted. Done in groups over a chunked matrix
%%% multiply, so the peak is one group of activations plus one chunk of channel
%%% data instead of the whole decomposition applied at once.
fprintf('evalplots.ic_overview: deriving %d component activation(s) ...\n', numel(comps));
unmix = single(EEG.icaweights * EEG.icasphere);      % nComp x nChanICA

%%% Epoched data laid flat. This stage's recordings are continuous, but a .set
%%% saved with trials would otherwise be indexed with two subscripts into a
%%% three-dimensional array - which silently reads epochs as extra time and
%%% would then report the wrong duration on the x axis.
data = EEG.data;
if ~ismatrix(data), data = reshape(data, size(data, 1), []); end
nPts = size(data, 2);

Pwr   = cell(1, numel(comps));
Freqs = [];
Tsec  = [];

sampChunk = 2^18;
compGroup = 8;
for g0 = 1:compGroup:numel(comps)
    g1 = min(g0 + compGroup - 1, numel(comps));
    W  = unmix(comps(g0:g1), :);

    act = zeros(size(W, 1), nPts, 'single');
    for i0 = 1:sampChunk:nPts
        ii = i0:min(i0 + sampChunk - 1, nPts);
        act(:, ii) = W * single(data(EEG.icachansind, ii));
    end

    for k = 1:size(act, 1)
        [Pwr{g0 + k - 1}, Freqs, Tsec] = ic_spectrogram(act(k, :), EEG.srate, ...
            'Window', opts.Window, 'FreqLim', opts.FreqLim, 'MaxCols', opts.MaxCols);
    end
    clear act
end
clear data

%%% --- Spectrogram colour limits ---
%%% Data-driven when not given explicitly: the values are log10(power) in
%%% arbitrary ICA units, so a fixed pair would be right for one decomposition
%%% and wrong for the next. Pooled across every component on the page rather
%%% than fit per-tile, so the forty tiles stay on one comparable scale.
cLimUse = opts.CLim;
if isempty(cLimUse)
    allLog = cellfun(@(p) log10(p(:) + eps), Pwr, 'UniformOutput', false);
    allLog = vertcat(allLog{:});
    allLog = allLog(isfinite(allLog));
    if isempty(allLog)
        cLimUse = [-1 1];
    else
        cLimUse = prctile(allLog, [1 99]);
        cLimUse = [floor(cLimUse(1)*10)/10, ceil(cLimUse(2)*10)/10];
        if cLimUse(2) <= cLimUse(1), cLimUse(2) = cLimUse(1) + 1; end
    end
end

%%% Minutes for a nap, hours for a night - a night's x axis in minutes is four
%%% digits wide and unreadable in a tile this size.
tHours = Tsec / 3600;
if isempty(tHours) || tHours(end) <= 2
    tPlot = Tsec / 60;  tLabel = 'Time (min)';
else
    tPlot = tHours;     tLabel = 'Time (h)';
end
if isempty(tPlot), tXLim = [0 1]; else, tXLim = [tPlot(1) max(tPlot(end), tPlot(1) + eps)]; end

%%% --- The page ---
warnState = warning('off', 'all');
cleanupWarn = onCleanup(@() warning(warnState));

%%% Taller per component row than the old 2x2 block needed, now that each one
%%% is a 3x3 with its own spectrum column on top of the spectrogram.
fig = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [1 1 min(nCols * 6.5, 60) min(nRows * 8.5, 48)], ...
    'Name', opts.Title);
outer = tiledlayout(fig, nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
%%% A strip of the figure kept clear on the right for the shared colourbar. A
%%% colorbar with Layout.Tile = 'east' would dock to the *inner* layout its axes
%%% belongs to - one component's block - not to the page.
outer.OuterPosition = [0.004 0 0.936 1];
title(outer, opts.Title, 'FontSize', 12, 'FontWeight', 'bold');

topoChanlocs = EEG.chanlocs(EEG.icachansind);
fMask        = Freqs >= opts.FreqLim(1) & Freqs <= opts.FreqLim(2);
fTicks       = 0:10:opts.FreqLim(2);
fTicks       = fTicks(fTicks >= opts.FreqLim(1) & fTicks <= opts.FreqLim(2));
%%% The spectrum panel is a 2-row-tall, 1-col-wide tile - too narrow for
%%% fTicks' full set without MATLAB auto-rotating the labels to fit. Three
%%% ticks is what that width has room for horizontally.
fTicksPsd    = unique(round(linspace(opts.FreqLim(1), opts.FreqLim(2), 3)));

for iC = 1:numel(comps)
    c       = comps(iC);
    %%% Bottom of its own column, not of the page: with a part-filled last row
    %%% the columns past its end end one row higher, and those are the tiles
    %%% whose time axis a reader actually has in front of them.
    isLast  = iC + nCols > numel(comps);
    isFirst = mod(iC - 1, nCols) == 0;

    %%% 3x3: topo takes the top-left 2x2, the spectrum the top-right 2x1 column,
    %%% and the spectrogram the full bottom row. 'compact' rather than 'none':
    %%% at 'none' the topography's ear cartoons run into the spectrum's axis and
    %%% the spectrum's tick labels are covered by the spectrogram below it. The
    %%% space this costs is space the panels could not have used anyway.
    inner = tiledlayout(outer, 3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    inner.Layout.Tile = iC;

    if isempty(className{c}), lbl = sprintf('IC #%d', c);
    else,                     lbl = sprintf('IC #%d - %s', c, className{c});
    end
    title(inner, lbl, 'FontSize', 11, 'FontWeight', 'bold', ...
        'Color', iclabel_color(classIdx(c)));

    %%% --- Top left, 2x2: topography ---
    %%% The pale blue disc behind each head is a patch topoplot draws itself (its
    %%% BACKCOLOR skirt), not an axes colour, so it stays - and it is what every
    %%% other topography in this repo looks like.
    axTopo = nexttile(inner, 1, [2 2]);
    quietTopoplot(EEG.icawinv(:, c), topoChanlocs);
    colormap(axTopo, custom_cmap());

    %%% --- Top right, 2x1: power spectrum, log10(power) against linear Hz ---
    axPsd = nexttile(inner, 3, [2 1]);
    P     = Pwr{iC};
    if isempty(P)
        axis(axPsd, 'off');
    else
        %%% Averaged across time first, then logged - not the other way round,
        %%% which would be the log-mean rather than the mean-then-log this line
        %%% is meant to read as (and would be dragged around by the loudest
        %%% instant rather than showing the typical level).
        logPsd = log10(mean(P, 2, 'omitnan') + eps);
        plot(axPsd, Freqs(fMask), logPsd(fMask), '-', ...
            'Color', [0.10 0.20 0.45], 'LineWidth', 0.9);
        xlim(axPsd, opts.FreqLim);
        set(axPsd, 'FontSize', 6, 'Box', 'off', 'TickDir', 'out', 'XTick', fTicksPsd);
        ylabel(axPsd, 'log_{10}(power)', 'FontSize', 6);
    end

    %%% --- Bottom, 1x3: time-frequency, spanning the full width ---
    axTf = nexttile(inner, 7, [1 3]);
    if isempty(P) || size(P, 2) < 2
        axis(axTf, 'off');
        text(0.5, 0.5, 'no data', 'Parent', axTf, 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'FontSize', 7, 'Color', [0.6 0.6 0.6]);
    else
        lastTf = axTf;
        %%% Plain log10 of the power itself, pixel by pixel - no per-row median
        %%% division, so this stays on the same scale as the spectrum line above
        %%% and a louder component still looks louder.
        rel = log10(P + eps);
        imagesc(axTf, tPlot, Freqs, rel);
        set(axTf, 'YDir', 'normal');
        clim(axTf, cLimUse);
        colormap(axTf, custom_cmap());
        xlim(axTf, tXLim);
        ylim(axTf, opts.FreqLim);
        set(axTf, 'FontSize', 6, 'Box', 'off', 'TickDir', 'out', 'YTick', fTicks);
        %%% Axis labels only around the edge of the page: forty copies of
        %%% "Time (h)" is forty times the ink for the same information.
        if isFirst, ylabel(axTf, 'Hz', 'FontSize', 6); end
        if isLast,  xlabel(axTf, tLabel, 'FontSize', 6); else, set(axTf, 'XTickLabel', []); end
    end
end

%%% One colourbar for the page - every map is on the same log10(power) scale,
%%% so forty of them would say the same thing forty times. Hung on an invisible
%%% axes of its own rather than on one of the spectrograms: an axes inside a
%%% tiledlayout has its position managed by the layout, and this one has to sit
%%% in the strip reserved above.
if exist('lastTf', 'var') && isgraphics(lastTf)
    cbPos = [0.955 0.30 0.012 0.40];
    axCb  = axes(fig, 'Units', 'normalized', 'Position', cbPos, 'Visible', 'off');
    colormap(axCb, custom_cmap());
    clim(axCb, cLimUse);
    cb = colorbar(axCb, 'Position', cbPos);
    cb.Label.String = 'log_{10}(power)';
    cb.FontSize = 8;
end

%%% topoplot paints the figure EEGLAB's pale blue ([.93 .96 1]) on its way past,
%%% and it is the page background that ends up in the PNG. Put white back.
set(fig, 'Color', 'w');

save_fig(fig, opts.SavePath, opts.Suffix);
end

% -------------------------------------------------------------------------
function quietTopoplot(vals, locs) %#ok<INUSD>
% topoplot into the current axes, without the log it normally leaves behind.
%
% topoplot ends its contour block with a bare "warning on", undoing whatever the
% caller had suppressed, and from there warns five times per call that it cannot
% set the position of an axes a tiledlayout owns. It is right and it does not
% matter - the layout places the axes and the map is drawn correctly either way -
% but at forty components that is two hundred lines of batch log saying nothing,
% and there is no warning state a caller can set that survives topoplot turning
% them all back on. evalc is what is left: it captures the printing and nothing
% else, so errors still throw. vals and locs are read from inside the string.

evalc(['topoplot(vals, locs, ''electrodes'', ''off'', ' ...
       '''numcontour'', 0, ''conv'', ''on'');']);
end
