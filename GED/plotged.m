function plotged(GED, opts)
% PLOTGED  Plot the diagnostic figure for a ged() result (Cohen 2022, Fig. 4).
%
%   plotged(GED)
%   plotged(GED, Name, Value, ...)
%
%   Three rows: the eigenspectrum (with the permutation threshold, if computed),
%   the component maps, and the power spectrum of each component time series.
%   This is exactly what ged(..., 'plot', true) draws internally; call this
%   directly on a GED struct you already have - e.g. one stored by
%   bidsfun_subcomp, which never plots on its own because it runs unattended
%   over a whole BIDS dataset. GED.info carries the srate and chanlocs it needs,
%   so no other arguments are required.
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
%   cmap       Colormap for the component maps: any colormap function name or an
%              n-by-3 matrix. Default: 'hot'.
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
end

comps = opts.comps;
if isempty(comps)
    comps = 1:min(opts.ncomps, numel(GED.evals));
end
ncomps = numel(comps);

srate = GED.info.srate;
freqlim = opts.freqlim;
if isempty(freqlim), freqlim = [0 min(45, srate / 2)]; end

figure('Color', 'w', 'Name', 'GED diagnostics', 'NumberTitle', 'off');
tiledlayout(3, ncomps, 'TileSpacing', 'compact', 'Padding', 'compact');

%%% Eigenspectrum. The elbow says how many directions actually separate S from R;
%%% the dashed line is the permutation threshold, when one was computed.
nexttile([1 ncomps]);
nspec = min(numel(GED.evals), 20);
plot(1:nspec, GED.evals(1:nspec), 'ks-', 'MarkerFaceColor', 'k'); hold on
if isfield(GED, 'perm') && ~isempty(GED.perm.crit95) && ~isnan(GED.perm.crit95)
    yline(GED.perm.crit95, 'r--', 'p < .05');
end
for c = comps
    plot(c, GED.evals(c), 'ro', 'MarkerSize', 10, 'LineWidth', 1.5);
end
xlabel('Component'); ylabel('\lambda (S:R ratio)'); title('Eigenspectrum'); box off

hastopo = hastopocoords(GED.info.chanlocs);
cmap    = resolvecmap(opts.cmap);
topoAxes = gobjects(0);
for i = 1:ncomps
    c = comps(i);
    ax = nexttile;
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
    title(sprintf('#%d, \\lambda = %.2f', c, GED.evals(c)));
end

for i = 1:ncomps
    c = comps(i);
    nexttile;
    [pxx, hz] = compspectrum(GED.comp(c, :), srate);
    plot(hz, 10 * log10(pxx), 'k'); xlim(freqlim);
    xlabel('Hz');
    if i == 1, ylabel('Power (dB)'); end
    box off
end

%%% Set after every topoplot call, not inside the loop: EEGLAB's topoplot applies
%%% its own colormap and background from icadefs to the figure it lands in, and
%%% would otherwise overwrite whatever the earlier tiles were given.
for ax = topoAxes
    colormap(ax, cmap);
end
set(gcf, 'Color', 'w');
end
