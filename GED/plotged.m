function plotged(GED, opts)
% PLOTGED  Plot the diagnostic figure for a ged() result (Cohen 2022, Fig. 4).
%
%   plotged(GED)
%   plotged(GED, Name, Value, ...)
%
%   The eigenspectrum (with the permutation threshold, if computed), the
%   activation trace of each component, the component maps, and the power
%   spectrum of each component. This is exactly what ged(..., 'plot', true)
%   draws internally; call this directly on a GED struct you already have - e.g.
%   one stored by bidsfun_subcomp, which never plots on its own because it runs
%   unattended over a whole BIDS dataset. GED.info carries the srate and chanlocs
%   it needs, so no other arguments are required.
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
%   acttype    'line' (default) draws each activation as a signal, min/max
%              decimated for display; 'heat' draws its smoothed absolute
%              amplitude as a colour strip, which reads better over long
%              recordings. The same two modes as plotgednight.
%   showact    Draw the activation traces at all. Default: true.
%   maxpoints  Points drawn per activation trace. Default: 20000.
%   smoothsec  Envelope smoothing for acttype 'heat', in seconds. Default: 5.
%
%   Example
%   -------
%     [failures, GEDs] = bidsfun_subcomp(BIDS, 'gedargs', {'peakfreq', 12});
%     plotged(GEDs(1).ged)
%     plotged(GEDs(1).ged, 'comps', [1 3], 'freqlim', [0 25])

arguments
    GED    struct
    opts.ncomps (1,1) double = 5
    opts.comps        double = []
    opts.freqlim      double = []
    opts.cmap                = 'hot'
    opts.acttype   (1,:) char {mustBeMember(opts.acttype, {'line', 'heat'})} = 'line'
    opts.showact   (1,1) logical = true
    opts.maxpoints (1,1) double {mustBePositive} = 20000
    opts.smoothsec (1,1) double {mustBePositive} = 5
end

comps = opts.comps;
if isempty(comps)
    comps = 1:min(opts.ncomps, numel(GED.evals));
end
ncomps = numel(comps);

srate = GED.info.srate;
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

nrows = 1 + ncomps * double(showact) + 2;

figure('Color', 'w', 'Name', 'GED diagnostics', 'NumberTitle', 'off', ...
    'Position', [80 80 1200 min(1200, 300 + 95 * nrows)]);
tiledlayout(nrows, ncomps, 'TileSpacing', 'compact', 'Padding', 'compact');

palette  = gedpalette(ncomps);
cmapAxes = gobjects(0);

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

%%% Activation traces, stacked and sharing one time axis. Straight from GED.comp,
%%% so they span exactly the data ged() was handed - see the note in the help.
row = 2;
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
        ax = nexttile((row - 1) * ncomps + 1, [1 ncomps]);
        if drawactivation(ax, tplot, act(i, :), srate, opts.acttype, palette(i, :), ...
                opts.maxpoints, opts.smoothsec)
            cmapAxes(end + 1) = ax; %#ok<AGROW>
        end
        ylabel(ax, sprintf('#%d', comps(i)), 'Color', palette(i, :), 'FontWeight', 'bold');
        xlim(ax, xlims);
        if i < ncomps, set(ax, 'XTickLabel', []); else, xlabel(ax, tlabel); end
        axTime(i) = ax;
        row = row + 1;
    end
    if ncomps > 1, linkaxes(axTime, 'x'); end
end

hastopo = hastopocoords(GED.info.chanlocs);
cmap    = resolvecmap(opts.cmap);
topoAxes = gobjects(0);
for i = 1:ncomps
    c = comps(i);
    ax = nexttile((row - 1) * ncomps + i);
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

row = row + 1;
for i = 1:ncomps
    c = comps(i);
    nexttile((row - 1) * ncomps + i);
    [pxx, hz] = compspectrum(GED.comp(min(c, size(GED.comp, 1)), :), srate);
    plot(hz, 10 * log10(pxx), 'Color', palette(i, :), 'LineWidth', 1.6);
    xlim(freqlim); xlabel('Hz');
    if i == 1, ylabel('Power (dB)'); end
    box off
end

%%% Set after every topoplot call, not inside the loop: EEGLAB's topoplot applies
%%% its own colormap and background from icadefs to the figure it lands in, and
%%% would otherwise overwrite whatever the earlier tiles were given.
for ax = [topoAxes cmapAxes]
    colormap(ax, cmap);
end
set(gcf, 'Color', 'w');
end
