function plotged(GED, opts)
% PLOTGED  Plot the diagnostic figure for a ged() result (Cohen 2022, Fig. 4).
%
%   plotged(GED)
%   plotged(GED, Name, Value, ...)
%
%   Layout, top to bottom:
%
%     1. the eigenspectrum (with the permutation threshold, if computed) beside
%        the power spectra of all selected components, overlaid in one axis
%     2. the component maps, double height
%     3. one activation trace per component, stacked, sharing one time axis
%
%   This is exactly what ged(..., 'plot', true) draws internally; call it
%   directly on a GED struct you already have - e.g. one stored by
%   bidsfun_subcomp, which never plots on its own because it runs unattended over
%   a whole BIDS dataset. GED.info carries the srate and chanlocs it needs, so no
%   other arguments are required.
%
%   The activations come straight out of GED.comp, so no EEG is needed - but they
%   therefore span exactly the data ged() was given, and nothing more. When the
%   GED ran on a stage-selected recording that is a concatenation of the kept
%   epochs, and its time axis is not the recording's own. Use plotgednight to put
%   the traces back on a full-night time base beside a hypnogram.
%
%   Required
%   --------
%   GED    The struct returned by ged().
%
%   Name-Value
%   ----------
%   ncomps     How many components to show. Default: 5.
%   comps      Explicit list of component indices to show instead of 1:ncomps,
%              e.g. [1 4] to compare the top component against a later one you
%              picked by eye (3.7). Overrides ncomps.
%   freqlim    x-axis limits for the spectra, in Hz. Default: [0 min(45, srate/2)].
%   cmap       Colormap for the component maps and heat strips: any colormap
%              function name or an n-by-3 matrix. Default: 'hot'.
%   acttype    How to draw each activation:
%                'signal'    (default) the component time series itself, min/max
%                            decimated so nothing is aliased away at screen
%                            resolution
%                'envelope'  its smoothed absolute amplitude as a line - the
%                            readable choice once the raw trace turns into a
%                            solid band
%                'heat'      the same envelope as a colour strip
%              The same three modes as plotgednight.
%   showact    Draw the activation traces at all. Default: true.
%   xwindow    Seconds of data to show at once. [] (default) fits the whole
%              recording into the axes; give it a duration - 30 for a screen of
%              sleep scoring - and the traces show that much at a time, with a
%              scrollbar underneath to page through the rest. The activation
%              axes are linked in x and y, so scrolling and zooming move them
%              together and their amplitudes stay comparable.
%   maxpoints  Points drawn per activation trace. Default: 20000.
%   smoothsec  Envelope smoothing for 'envelope' and 'heat', in seconds.
%              Default: 5.
%
%   Example
%   -------
%     [failures, GEDs] = bidsfun_subcomp(BIDS, 'gedargs', {'peakfreq', 12});
%     plotged(GEDs(1).ged)
%     plotged(GEDs(1).ged, 'comps', [1 3], 'freqlim', [0 25])
%     plotged(GEDs(1).ged, 'acttype', 'envelope')
%     plotged(GEDs(1).ged, 'xwindow', 30)          % 30 s at a time, scrollable
%
%   See also PLOTGEDNIGHT, GED.

arguments
    GED    struct
    opts.ncomps (1,1) double = 5
    opts.comps        double = []
    opts.freqlim      double = []
    opts.cmap                = 'hot'
    opts.acttype   (1,:) char {mustBeMember(opts.acttype, ...
        {'signal', 'envelope', 'heat', 'line'})} = 'signal'
    opts.showact   (1,1) logical = true
    opts.xwindow         double = []
    opts.maxpoints (1,1) double {mustBePositive} = 20000
    opts.smoothsec (1,1) double {mustBePositive} = 5
end

comps = opts.comps;
if isempty(comps)
    comps = 1:min(opts.ncomps, numel(GED.evals));
end
comps  = comps(:)';
ncomps = numel(comps);

srate   = GED.info.srate;
freqlim = opts.freqlim;
if isempty(freqlim), freqlim = [0 min(45, srate / 2)]; end

%%% Activations are optional, and only available for the components ged() computed
%%% a time series for - GED.comp holds the first ncomps of them.
showact = opts.showact && ~isempty(GED.comp);
if showact && any(comps > size(GED.comp, 1))
    warning('plotged:missingActivations', ...
        ['GED.comp only holds %d component(s), so the activation traces are omitted. ' ...
         'Re-run ged() with a larger ncomps, or use plotgednight with the EEG to ' ...
         'project the missing ones.'], size(GED.comp, 1));
    showact = false;
end

%%% Two columns per component: that makes the grid divisible both by the number
%%% of maps in row 2 and by the two panels sharing row 1, whatever ncomps is.
ncols    = 2 * ncomps;
topoSpan = 2;                                   % maps get double height
nrows    = 1 + topoSpan + ncomps * double(showact);

fig = figure('Color', 'w', 'Name', 'GED diagnostics', 'NumberTitle', 'off', ...
    'Position', [80 80 1250 min(1250, 220 + 105 * nrows)]);
tl = tiledlayout(nrows, ncols, 'TileSpacing', 'compact', 'Padding', 'compact');

palette  = gedpalette(ncomps);
cmapAxes = gobjects(0);

%% ----------------------------------------- eigenspectrum and power spectra
%%% Eigenspectrum. The elbow says how many directions actually separate S from R;
%%% the dashed line is the permutation threshold, when one was computed.
nexttile(1, [1 ncomps]);
nspec = min(numel(GED.evals), max(20, max(comps)));
plot(1:nspec, GED.evals(1:nspec), 's-', 'Color', [.35 .35 .35], ...
    'MarkerFaceColor', [.35 .35 .35], 'MarkerSize', 4); hold on
if isfield(GED, 'perm') && ~isempty(GED.perm.crit95) && ~isnan(GED.perm.crit95)
    yline(GED.perm.crit95, 'r--', 'p < .05');
end
for i = 1:ncomps
    plot(comps(i), GED.evals(comps(i)), 'o', 'MarkerSize', 9, ...
        'MarkerFaceColor', palette(i, :), 'MarkerEdgeColor', 'k', 'LineWidth', 0.8);
end
xlabel('Component'); ylabel('\lambda (S:R ratio)'); title('Eigenspectrum'); box off
xlim([0.5 nspec + 0.5]);

%%% All spectra in one axis: overlaid they can be compared directly, which is the
%%% whole question when deciding which component carries the band of interest.
axS = nexttile(1 + ncomps, [1 ncomps]);
hold(axS, 'on');
labels = cell(1, ncomps);
for i = 1:ncomps
    c = min(comps(i), size(GED.comp, 1));
    [pxx, hz] = compspectrum(GED.comp(c, :), srate);
    plot(axS, hz, 10 * log10(pxx), 'Color', palette(i, :), 'LineWidth', 1.6);
    labels{i} = sprintf('#%d (\\lambda = %.2f)', comps(i), GED.evals(comps(i)));
end
xlim(axS, freqlim);
xlabel(axS, 'Frequency (Hz)'); ylabel(axS, 'Power (dB)');
title(axS, 'Component spectra'); box(axS, 'off');
legend(axS, labels, 'Box', 'off', 'Location', 'northeast');

%% ------------------------------------------------------------ component maps
hastopo  = hastopocoords(GED.info.chanlocs);
cmap     = resolvecmap(opts.cmap);
topoAxes = gobjects(0);
for i = 1:ncomps
    c  = comps(i);
    ax = nexttile(ncols + 2 * i - 1, [topoSpan 2]);
    if hastopo
        try
            topoplot(GED.maps(:, c), GED.info.chanlocs, 'electrodes', 'off', 'numcontour', 0);
            topoAxes(end + 1) = ax; %#ok<AGROW>
        catch ME
            warning('plotged:topoplotFailed', ...
                'topoplot failed (%s); falling back to a weight bar.', ME.message);
            cla(ax); bar(ax, GED.maps(:, c), 'k'); axis(ax, 'tight'); box(ax, 'off');
        end
    else
        bar(ax, GED.maps(:, c), 'k'); axis(ax, 'tight'); box(ax, 'off');
    end
    title(sprintf('#%d, \\lambda = %.2f', c, GED.evals(c)), 'Color', palette(i, :));
end

%% -------------------------------------------------------- activation traces
%%% Stacked and sharing one time axis, straight from GED.comp - so they span
%%% exactly the data ged() was handed; see the note in the help above.
if showact
    act = reshape(GED.comp(comps, :), ncomps, []);
    t   = (0:size(act, 2) - 1) / srate;
    if t(end) > 7200,    tscale = 3600; tlabel = 'Time (h)';
    elseif t(end) > 120, tscale = 60;   tlabel = 'Time (min)';
    else,                tscale = 1;    tlabel = 'Time (s)';
    end
    tplot  = t / tscale;
    xlims  = [0 max(tplot(end), eps)];
    axTime = gobjects(1, ncomps);
    for i = 1:ncomps
        row = 1 + topoSpan + i;
        ax  = nexttile((row - 1) * ncols + 1, [1 ncols]);
        if drawactivation(ax, tplot, act(i, :), srate, opts.acttype, palette(i, :), ...
                opts.maxpoints, opts.smoothsec)
            cmapAxes(end + 1) = ax; %#ok<AGROW>
        end
        ylabel(ax, sprintf('#%d', comps(i)), 'Color', palette(i, :), 'FontWeight', 'bold');
        xlim(ax, xlims);
        if i < ncomps, set(ax, 'XTickLabel', []); else, xlabel(ax, tlabel); end
        axTime(i) = ax;
    end

    %%% Link x and y across the traces: scrolling or zooming moves them together,
    %%% and a shared amplitude axis is what makes their heights comparable.
    %%% linkprop rather than linkaxes, because plotgednight has to link two
    %%% different sets of axes on two different properties and linkaxes would
    %%% overwrite one with the other; the handle has to be kept alive, hence the
    %%% appdata.
    if ncomps > 1
        setappdata(fig, 'gedAxisLinks', linkprop(axTime, {'XLim', 'YLim'}));
        %%% A common amplitude range, set once: each axis was scaled to its own
        %%% component, so linking alone would hand every trace whichever range
        %%% happened to be set last.
        yl = get(axTime, 'YLim');
        yl = vertcat(yl{:});
        ylim(axTime(1), [min(yl(:, 1)) max(yl(:, 2))]);
    end
    addtimescroll(fig, tl, axTime, opts.xwindow / tscale, xlims);
end

%%% Set after every topoplot call, not inside the loop: EEGLAB's topoplot applies
%%% its own colormap and background from icadefs to the figure it lands in, and
%%% would otherwise overwrite whatever the earlier tiles were given.
for ax = [topoAxes cmapAxes]
    colormap(ax, cmap);
end
set(fig, 'Color', 'w');
end
