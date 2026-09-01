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

hastopo = exist('topoplot', 'file') == 2 && ~isempty(GED.info.chanlocs);
for i = 1:ncomps
    c = comps(i);
    nexttile;
    if hastopo
        topoplot(GED.maps(:, c), GED.info.chanlocs, 'electrodes', 'off', 'numcontour', 0);
    else
        bar(GED.maps(:, c), 'k'); axis tight; box off
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
end

% -------------------------------------------------------------------------
function [pxx, hz] = compspectrum(x, srate)
% Welch spectrum of a component time series, with a plain-MATLAB fallback for
% installations without the Signal Processing Toolbox.

x   = x(:);
win = min(numel(x), round(4 * srate));
if exist('pwelch', 'file') == 2
    [pxx, hz] = pwelch(x, hann(win), [], [], srate);
    return
end
nseg  = max(1, floor(numel(x) / win));
taper = 0.5 - 0.5 * cos(2 * pi * (0:win - 1)' / win);
pxx   = zeros(floor(win / 2) + 1, 1);
for s = 1:nseg
    seg = x((s - 1) * win + (1:win)) .* taper;
    amp = abs(fft(seg)).^2;
    pxx = pxx + amp(1:floor(win / 2) + 1);
end
pxx = pxx / (nseg * srate * sum(taper.^2));
hz  = linspace(0, srate / 2, numel(pxx))';
end
