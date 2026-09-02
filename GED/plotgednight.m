function plotgednight(GED, EEG, scoring, opts)
% PLOTGEDNIGHT  Whole-night overview of a ged() result: hypnogram + activations.
%
%   plotgednight(GED, EEG, scoring)
%   plotgednight(GED, EEG, scoring, Name, Value, ...)
%
%   The time-resolved counterpart to plotged: where that one answers "is this
%   component real and what does it look like", this one answers "when during the
%   night is it active, and in which sleep stage". Layout, top to bottom:
%
%     1. the eigenspectrum, full width
%     2. the hypnogram, full width
%     3. one activation trace per selected component, stacked, full width, all
%        sharing the hypnogram's time axis (linked, so zooming one zooms all)
%     4. a block with the power spectra of all selected components overlaid in
%        one axis on the left, and the component maps beside it, up to four per
%        row; the spectrum spans the full height of that block
%
%   Do I need the EEG?
%   ------------------
%   Not always - GED.comp already holds the component time series. Pass the EEG
%   when either of these applies, which between them covers most real uses:
%
%     - you want components beyond the ncomps that ged() computed (10 by default);
%     - the GED was run on stage-selected data (as bidsfun_subcomp does), so
%       GED.comp has the dropped epochs cut out of it and no longer lines up
%       with a full-night hypnogram.
%
%   Given the EEG, the filters are re-projected onto it (the same mean-centring
%   and channel scaling ged() used), which puts the traces back on the
%   recording's own continuous time base - the one the hypnogram lives on.
%   Without it, GED.comp is used as-is and the time axis is only meaningful if
%   the GED ran on the whole recording.
%
%   Required
%   --------
%   GED       The struct returned by ged().
%   EEG       EEGLAB struct (or channels x time array) to project the filters
%             onto. [] to use GED.comp instead - see above.
%   scoring   One sleep-stage digit per epoch, as scoreloader returns them:
%             -3 N3, -2 N2, -1 N1, 0 Wake, 1 REM. [] to omit the hypnogram.
%
%   Name-Value
%   ----------
%   ncomps       How many components to show. Default: 5.
%   comps        Explicit component indices, overriding ncomps. E.g. [1 4].
%   epochlength  Scoring epoch length in seconds. Default: 30.
%   keptepochs   Epoch indices that actually entered the GED (bidsfun_subcomp
%                returns these). Shaded in the hypnogram, so it is visible which
%                part of the night the filters were built on. If the EEG passed
%                in is itself the stage-selected recording, the scoring is
%                subset with these before plotting, so the two still align.
%   acttype      How to draw each activation:
%                  'signal'    (default) the component time series itself, min/max
%                              decimated so a whole night stays honest at screen
%                              resolution. Best for minutes, not hours.
%                  'envelope'  its smoothed absolute amplitude as a line, which
%                              stays readable once the raw trace turns into a
%                              solid band
%                  'heat'      the same envelope as a colour strip - the easiest
%                              to read across 8 hours, with bursts of spindle or
%                              slow-wave activity showing up as bright bands
%   smoothsec    Envelope smoothing for 'envelope' and 'heat', in seconds. Default: 5.
%   freqlim      x-limits for the spectra, in Hz. Default: [0 min(45, srate/2)].
%   evalscale    y-scale of the eigenspectrum: 'auto' (default) switches to log
%                when one component dwarfs the rest, which is otherwise the case
%                where none of the others can be read; 'linear' or 'log' to force.
%   cmap         Colormap for the component maps and the heat strips: any
%                colormap function name or an n-by-3 matrix. Default: 'turbo'.
%   timeunit     'auto' (default), 'h', 'min' or 's' for the shared time axis.
%   linewidth    Line width for traces and spectra. Default: 1.6.
%   xwindow      Seconds of data to show at once. [] (default) fits the whole
%                recording into the axes; give it a duration - 30 for a screen of
%                sleep scoring - and the traces show that much at a time, with a
%                scrollbar underneath. The hypnogram scrolls with them, since it
%                shares their time axis.
%   maxpoints    Points drawn per activation trace. Default: 20000.
%   title        Text for the figure title, e.g. the fileID. Default: ''.
%
%   Example
%   -------
%     [failures, GEDs] = bidsfun_subcomp(BIDS, 'gedargs', {'peakfreq', 12});
%     EEG     = fast_eeg_import(theSetFile);       % the full night
%     scoring = scoreloader(theScoringFile);
%     plotgednight(GEDs(1).ged, EEG, scoring, ...
%         'comps', [1 2 3], 'keptepochs', GEDs(1).keptepochs, 'acttype', 'heat')
%
%   See also PLOTGED, GED, EXTRACTSLEEPEPOCHS.

arguments
    GED     struct
    EEG                  = []
    scoring       double = []
    opts.ncomps      (1,1) double = 5
    opts.comps             double = []
    opts.epochlength (1,1) double {mustBePositive} = 30
    opts.keptepochs        double = []
    opts.acttype     (1,:) char {mustBeMember(opts.acttype, ...
        {'signal', 'envelope', 'heat'})} = 'signal'
    opts.smoothsec   (1,1) double {mustBePositive} = 5
    opts.freqlim           double = []
    opts.evalscale   (1,:) char {mustBeMember(opts.evalscale, {'auto', 'linear', 'log'})} = 'auto'
    opts.cmap                     = 'turbo'
    opts.timeunit    (1,:) char {mustBeMember(opts.timeunit, {'auto', 'h', 'min', 's'})} = 'auto'
    opts.linewidth   (1,1) double {mustBePositive} = 1.6
    opts.xwindow           double = []
    opts.maxpoints   (1,1) double {mustBePositive} = 20000
    opts.title       (1,:) char = ''
end

srate = GED.info.srate;

comps = opts.comps;
if isempty(comps)
    comps = 1:min(opts.ncomps, numel(GED.evals));
end
comps = comps(:)';
if any(comps > numel(GED.evals))
    error('plotgednight:badComponent', ...
        'Asked for component %d but the GED only has %d.', max(comps), numel(GED.evals));
end
ncomps = numel(comps);

%% ------------------------------------------------- the component time series
%%% Either re-project the filters onto the EEG (continuous, full-night time base)
%%% or fall back to what ged() already computed.
if ~isempty(EEG)
    if isstruct(EEG)
        X = double(EEG.data);
    else
        X = double(EEG);
    end
    if size(X, 1) ~= size(GED.filters, 1)
        error('plotgednight:channelMismatch', ...
            'The EEG has %d channels but the filters were built on %d.', ...
            size(X, 1), size(GED.filters, 1));
    end
    X    = reshape(X, size(X, 1), []);              % epoched input -> continuous
    X    = scalechannels(X, GED.info.channorm);
    act  = GED.filters(:, comps)' * (X - mean(X, 2));
else
    if any(comps > size(GED.comp, 1))
        error('plotgednight:noSuchComponent', ...
            ['GED.comp only holds %d component(s); pass the EEG so component %d ' ...
             'can be projected, or re-run ged() with a larger ncomps.'], ...
            size(GED.comp, 1), max(comps));
    end
    act = GED.comp(comps, :);
    act = reshape(act, ncomps, []);
end
npnts = size(act, 2);
t     = (0:npnts - 1) / srate;

%% ------------------------------------------------------------ the hypnogram
%%% The scoring and the traces have to describe the same stretch of recording.
%%% When they do not, the usual cause is a stage-selected EEG paired with the
%%% full night's scoring - which keptepochs can repair.
nEpochsData = floor(npnts / (opts.epochlength * srate));
if ~isempty(scoring)
    scoring = scoring(:)';
    if abs(numel(scoring) - nEpochsData) > 2
        if ~isempty(opts.keptepochs) && numel(opts.keptepochs) == nEpochsData
            scoring = scoring(opts.keptepochs);
            fprintf(['plotgednight: scoring subset to the %d kept epoch(s) to match the ' ...
                     'stage-selected recording.\n'], numel(scoring));
            opts.keptepochs = [];   % already applied; nothing left to shade
        else
            warning('plotgednight:scoringMismatch', ...
                ['the scoring covers %d epoch(s) but the data covers %d. The hypnogram ' ...
                 'and the traces will not line up - pass the full-night EEG, or the ' ...
                 'keptepochs that produced this recording.'], numel(scoring), nEpochsData);
        end
    end
end

%% ------------------------------------------------------------- time scaling
switch opts.timeunit
    case 'auto'
        if t(end) > 7200,   tscale = 3600; tlabel = 'Time (h)';
        elseif t(end) > 120, tscale = 60;  tlabel = 'Time (min)';
        else,                tscale = 1;   tlabel = 'Time (s)';
        end
    case 'h',   tscale = 3600; tlabel = 'Time (h)';
    case 'min', tscale = 60;   tlabel = 'Time (min)';
    case 's',   tscale = 1;    tlabel = 'Time (s)';
end
tplot = t / tscale;
xlims = [0 max(tplot(end), eps)];

freqlim = opts.freqlim;
if isempty(freqlim), freqlim = [0 min(45, srate / 2)]; end

%% ------------------------------------------------------------------ layout
%%% Rows: eigenspectrum, hypnogram, one per component, then the map/spectrum
%%% block. The block is at least two rows tall so the spectrum keeps a sane
%%% aspect ratio; with four or fewer maps they simply grow to fill it.
hasHypno  = ~isempty(scoring);
topoCols  = min(4, ncomps);
topoRows  = ceil(ncomps / topoCols);
blockRows = max(2, topoRows);
topoSpan  = floor(blockRows / topoRows);
specCols  = 2;
ncols     = specCols + topoCols;
nrows     = 1 + double(hasHypno) + ncomps + blockRows;

fig = figure('Color', 'w', 'Name', 'GED night overview', 'NumberTitle', 'off', ...
    'Position', [60 60 1500 min(1300, 240 + 105 * nrows)]);
tl = tiledlayout(nrows, ncols, 'TileSpacing', 'compact', 'Padding', 'compact');

palette  = gedpalette(ncomps);
cmapAxes = gobjects(0);   % re-coloured at the end; see the note down there

%% ------------------------------------------------------------ eigenspectrum
axE   = nexttile(1, [1 ncols]);
nspec = min(numel(GED.evals), max(20, max(comps)));
ev    = GED.evals(1:nspec);
plot(axE, 1:nspec, ev, 's-', 'Color', [.35 .35 .35], ...
    'MarkerFaceColor', [.35 .35 .35], 'MarkerSize', 4, 'LineWidth', 1.1); hold(axE, 'on')
if isfield(GED, 'perm') && ~isempty(GED.perm.crit95) && ~isnan(GED.perm.crit95)
    yline(axE, GED.perm.crit95, 'r--', 'p < .05');
end
for i = 1:ncomps
    plot(axE, comps(i), GED.evals(comps(i)), 'o', 'MarkerSize', 9, ...
        'MarkerFaceColor', palette(i, :), 'MarkerEdgeColor', 'k', 'LineWidth', 0.8);
end
xlabel(axE, 'Component'); title(axE, 'Eigenspectrum'); box(axE, 'off');
xlim(axE, [0.5 nspec + 0.5]);

%%% A dominant first component flattens every other one against the axis, which
%%% is exactly when you most want to see where the rest sit relative to the
%%% permutation threshold. Log scale when the spread demands it, and say so.
uselog = strcmp(opts.evalscale, 'log') || (strcmp(opts.evalscale, 'auto') && ...
    all(ev > 0) && max(ev) / median(ev) > 20);
if uselog
    set(axE, 'YScale', 'log');
    ylabel(axE, '\lambda (log)');
else
    ylabel(axE, '\lambda');
end

%% ---------------------------------------------------------------- hypnogram
axTime = gobjects(0);
row    = 2;
if hasHypno
    axH = nexttile((row - 1) * ncols + 1, [1 ncols]);
    drawhypnogram(axH, scoring, opts.epochlength / tscale, opts.keptepochs, opts.linewidth);
    xlim(xlims); set(axH, 'XTickLabel', []);
    axTime(end + 1) = axH;
    row = row + 1;
end

%% -------------------------------------------------------- activation traces
for i = 1:ncomps
    ax = nexttile((row - 1) * ncols + 1, [1 ncols]);
    if drawactivation(ax, tplot, act(i, :), srate, opts.acttype, palette(i, :), ...
            opts.maxpoints, opts.smoothsec)
        cmapAxes(end + 1) = ax; %#ok<AGROW>
    end
    ylabel(ax, sprintf('#%d', comps(i)), 'Color', palette(i, :), 'FontWeight', 'bold');
    xlim(ax, xlims);
    if i < ncomps
        set(ax, 'XTickLabel', []);
    else
        xlabel(ax, tlabel);
    end
    axTime(end + 1) = ax; %#ok<AGROW>
    row = row + 1;
end

%%% One time axis for the night: zoom or pan any of these and the rest follow,
%%% which is the whole point of stacking them. The hypnogram joins the x link but
%%% not the y one - its y axis is the sleep stages, not an amplitude. Hence
%%% linkprop and two separate groups: linkaxes keeps only one link per axis and
%%% the second call would silently undo the first.
actAxes = axTime(end - ncomps + 1:end);
links   = {};
if numel(axTime) > 1
    links{end + 1} = linkprop(axTime, 'XLim');
end
if ncomps > 1
    links{end + 1} = linkprop(actAxes, 'YLim');
    %%% A common amplitude range, set once - see the note in plotged.
    yl = get(actAxes, 'YLim');
    yl = vertcat(yl{:});
    ylim(actAxes(1), [min(yl(:, 1)) max(yl(:, 2))]);
end
setappdata(fig, 'gedAxisLinks', links);
addtimescroll(fig, tl, axTime, opts.xwindow / tscale, xlims);

%% ------------------------------------------------- spectra and component maps
R0 = row;
axS = nexttile((R0 - 1) * ncols + 1, [blockRows specCols]);
hold(axS, 'on');
labels = cell(1, ncomps);
for i = 1:ncomps
    [pxx, hz] = compspectrum(act(i, :), srate);
    plot(axS, hz, 10 * log10(pxx), 'Color', palette(i, :), 'LineWidth', opts.linewidth);
    labels{i} = sprintf('#%d (\\lambda = %.2f)', comps(i), GED.evals(comps(i)));
end
xlim(axS, freqlim);
xlabel(axS, 'Frequency (Hz)'); ylabel(axS, 'Power (dB)');
title(axS, 'Component spectra'); box(axS, 'off');
legend(axS, labels, 'Box', 'off', 'Location', 'northeast');

hastopo = hastopocoords(GED.info.chanlocs);
for i = 1:ncomps
    r   = R0 + (ceil(i / topoCols) - 1) * topoSpan;
    c   = specCols + mod(i - 1, topoCols) + 1;
    ax  = nexttile((r - 1) * ncols + c, [topoSpan 1]);
    if drawmap(ax, GED.maps(:, comps(i)), GED.info.chanlocs, hastopo, palette(i, :))
        cmapAxes(end + 1) = ax; %#ok<AGROW>
    end
    title(ax, sprintf('#%d', comps(i)), 'Color', palette(i, :));
end

if ~isempty(opts.title)
    sgtitle(opts.title, 'FontWeight', 'bold', 'Interpreter', 'none');
end

%%% Colours are asserted here, at the very end, and not where each axis is drawn:
%%% EEGLAB's topoplot pulls its own colormap and background out of icadefs and
%%% applies them to whatever figure it lands in, so anything set before the last
%%% topoplot call gets overwritten by it.
cmap = resolvecmap(opts.cmap);
for ax = cmapAxes
    colormap(ax, cmap);
end
set(fig, 'Color', 'w');
end

%% ========================================================================= %%
%  Local functions
%% ========================================================================= %%

function istopo = drawmap(ax, map, chanlocs, hastopo, barcolor)
% One component map, as a topography where the montage allows it and as a plain
% weight bar where it does not. topoplot is wrapped because a montage it dislikes
% should cost one panel, not the whole figure. Returns whether a topography was
% drawn, so the caller knows which axes still want the component colormap.

istopo = false;
if hastopo
    try
        topoplot(map, chanlocs, 'electrodes', 'off', 'numcontour', 0);
        istopo = true;
        return
    catch ME
        warning('plotgednight:topoplotFailed', ...
            'topoplot failed (%s); falling back to a weight bar.', ME.message);
        cla(ax);
    end
end
bar(ax, map, 'FaceColor', barcolor, 'EdgeColor', 'none');
axis(ax, 'tight'); box(ax, 'off');
end

% -------------------------------------------------------------------------
function drawhypnogram(ax, scoring, epochdur, keptepochs, lw)
% Standard hypnogram: Wake on top, then REM, N1, N2, N3 going down. REM is drawn
% over the trace in its own colour, the way sleep labs read it.

digits = [ 0  1 -1 -2 -3];          % W  REM  N1  N2  N3
ypos   = [ 5  4  3  2  1];
names  = {'N3', 'N2', 'N1', 'REM', 'W'};

y = nan(size(scoring));
for k = 1:numel(digits)
    y(scoring == digits(k)) = ypos(k);
end

x = (0:numel(scoring)) * epochdur;              % epoch edges
hold(ax, 'on');

%%% Which epochs fed the GED. Drawn first, as a background band, and run-length
%%% encoded so a night of alternating stages does not become 900 patches.
if ~isempty(keptepochs)
    mask = false(1, numel(scoring));
    mask(keptepochs(keptepochs >= 1 & keptepochs <= numel(scoring))) = true;
    for run = findruns(mask)
        patch(ax, x([run(1) run(2) + 1 run(2) + 1 run(1)]), [0.4 0.4 5.6 5.6], ...
            [0.85 0.90 0.97], 'EdgeColor', 'none');
    end
end

stairs(ax, x, [y y(end)], 'Color', [0.20 0.20 0.20], 'LineWidth', lw);

%%% REM on top of the trace, in red.
for run = findruns(scoring == 1)
    idx = run(1):run(2);
    stairs(ax, x([idx idx(end) + 1]), [y(idx) y(idx(end))], ...
        'Color', [0.84 0.19 0.15], 'LineWidth', lw + 0.6);
end

set(ax, 'YTick', ypos(end:-1:1), 'YTickLabel', names, 'YLim', [0.4 5.6]);
ylabel(ax, 'Stage'); box(ax, 'off');
end

% -------------------------------------------------------------------------
function runs = findruns(mask)
% Start/stop index pairs of each true run in a logical vector, as a 2 x n array
% so it can be walked with a for loop (each iteration yields one column).

mask = mask(:)';
d    = diff([false mask false]);
starts = find(d == 1);
stops  = find(d == -1) - 1;
runs   = [starts; stops];
end

% -------------------------------------------------------------------------
function X = scalechannels(X, mode)
% Mirror of ged()'s channel scaling, so a re-projection through the filters
% reproduces the component time series ged() would have computed.

switch mode
    case 'pooled'
        X = X - mean(X, 2);
        X = X / std(X(:));
    case 'zscore'
        X = (X - mean(X, 2)) ./ std(X, 0, 2);
end
end
