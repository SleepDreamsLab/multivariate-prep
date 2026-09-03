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
%              function name or an n-by-3 matrix. Default: 'turbo'.
%   acttype    How to draw each activation:
%                'signal'    (default) the component time series itself, min/max
%                            decimated so nothing is aliased away at screen
%                            resolution
%                'envelope'  its smoothed absolute amplitude as a line - the
%                            readable choice once the raw trace turns into a
%                            solid band
%                'heat'      the same envelope as a colour strip
%                'wavelet'   time-frequency power from complex Morlet wavelets,
%                            as a spectrogram - the one that shows which band a
%                            component is in at each moment, rather than only
%                            how strong it is. See the block below.
%                'multitaper' the same spectrogram, estimated instead by
%                            Thomson's multitaper method. Sharper in frequency
%                            and much steadier, blunter in time. See below.
%                'wavelet_signal', 'multitaper_signal'
%                            either spectrogram with the component's own trace
%                            drawn flush underneath it, at a fifth of the height.
%                            Worth the room whenever the question is whether a
%                            band is really oscillating or whether an artefact
%                            is being read as power - the two are hard to tell
%                            apart in a spectrogram alone and obvious in the
%                            trace. The figure grows to make room, so the
%                            spectrograms are not shortened to fit.
%              The first three are the same modes as plotgednight.
%   showact    Draw the activation traces at all. Default: true.
%   xwindow    Seconds of data to show at once. [] (default) fits the whole
%              recording into the axes; give it a duration - 30 for a screen of
%              sleep scoring - and the traces show that much at a time, with a
%              scrollbar underneath to page through the rest. The activation
%              axes are linked in x, so scrolling and zooming move them together;
%              for the trace modes they share a y axis as well, so their
%              amplitudes stay comparable.
%   maxpoints  Points drawn per activation trace. Default: 20000.
%   smoothsec  Envelope smoothing for 'envelope' and 'heat', in seconds.
%              Default: 5.
%
%   Which spectrogram?
%   ------------------
%   'wavelet' holds the number of cycles roughly fixed, so its resolution scales
%   with frequency: sharp in time and blunt in frequency at the top of the range,
%   the reverse at the bottom. 'multitaper' uses one window length and one
%   bandwidth for every row, so 1 Hz and 30 Hz are measured alike, and averages
%   several orthogonal tapers, which makes each estimate far steadier. Reach for
%   wavelets to time an event, multitapers to identify a frequency - a spindle's
%   exact peak, say, or whether two components sit in the same band.
%
%   Name-Value, acttype 'wavelet' and 'multitaper'
%   ----------------------------------------------
%   tffreqlim  [lo hi] of the frequency bank, in Hz. Default: [1 40], clipped to
%              the Nyquist frequency. Evenly spaced, and drawn on a linear
%              frequency axis to match.
%   tfnfreq    Number of frequencies in the bank. Default: 40, which at the
%              default limits puts a row every 1 Hz.
%   tfnorm     What the colours mean:
%                'db'      (default) dB relative to that frequency's median
%                          power over the recording. Without a per-frequency
%                          reference, 1/f leaves the bottom rows a solid bright
%                          bar and everything above them flat.
%                'percent' percent change from the same median
%                'raw'     10*log10 of the power itself, unreferenced
%              For 'db' and 'percent' the colour scale is symmetric about zero,
%              so the midpoint of the colormap marks a typical moment - which is
%              what a diverging map like bwr needs to read correctly.
%
%   Name-Value, acttype 'wavelet' only
%   ----------------------------------
%   tfcycles   [lo hi] wavelet cycles, spread linearly from the bottom of the
%              bank to the top. Default: [3 10] - fewer cycles buy time
%              precision at the low end, more buy frequency precision at the
%              high end.
%   tfsmooth   Smooth each row in time by this many of its own wavelet widths
%              before drawing it. Default: 1. A single power estimate scatters
%              over some 25 dB no matter what the signal is doing, and the
%              colours would otherwise show that scatter and nothing else; one
%              wavelet width discards no detail the wavelet could have resolved.
%              0 to see the transform unsmoothed. Multitaper needs no equivalent:
%              averaging K tapers is already what steadies it.
%
%   Name-Value, acttype 'multitaper' only
%   -------------------------------------
%   tfwindow   Taper window length, in seconds. Default: 4. Sets the time
%              resolution outright, and with tftw the frequency resolution.
%   tftw       Time-bandwidth product NW. Default: 3, giving a resolution of
%              2*tftw/tfwindow = 1.5 Hz at the defaults. Raising it buys more
%              usable tapers, and so less variance, over a wider band.
%   tftapers   Number of tapers. 0 (default) takes 2*tftw-1, past which the
%              sequences are no longer well concentrated.
%
%   The spectrogram is computed for the visible window only, and again whenever
%   the x limits change, so xwindow scrolls at full resolution rather than
%   paging through one coarse precomputed image. The reference and the colour
%   limits are fixed up front from snippets spread over the whole recording, so
%   two windows can be compared by eye.
%
%   Example
%   -------
%     [failures, GEDs] = bidsfun_subcomp(BIDS, 'gedargs', {'peakfreq', 12});
%     plotged(GEDs(1).ged)
%     plotged(GEDs(1).ged, 'comps', [1 3], 'freqlim', [0 25])
%     plotged(GEDs(1).ged, 'acttype', 'envelope')
%     plotged(GEDs(1).ged, 'xwindow', 30)          % 30 s at a time, scrollable
%     plotged(GEDs(1).ged, 'acttype', 'wavelet')   % whole night, time-frequency
%     plotged(GEDs(1).ged, 'acttype', 'wavelet', 'xwindow', 300, ...
%         'tffreqlim', [0.5 20])                   % 5 min at a time, scrollable
%     plotged(GEDs(1).ged, 'acttype', 'multitaper')
%     plotged(GEDs(1).ged, 'acttype', 'multitaper', 'tfwindow', 6, 'tftw', 4, ...
%         'cmap', slanCM('bwr'))                   % sharper in frequency
%     plotged(GEDs(1).ged, 'acttype', 'wavelet_signal', 'xwindow', 30)
%                                                  % spectrogram over its trace
%
%   See also PLOTGEDNIGHT, GED, TFMORLET, TFMULTI.

arguments
    GED    struct
    opts.ncomps (1,1) double = 5
    opts.comps        double = []
    opts.freqlim      double = []
    opts.cmap                = 'turbo'
    opts.acttype   (1,:) char {mustBeMember(opts.acttype, ...
        {'signal', 'envelope', 'heat', 'wavelet', 'multitaper', ...
         'wavelet_signal', 'multitaper_signal'})} = 'signal'
    opts.showact   (1,1) logical = true
    opts.xwindow         double = []
    opts.maxpoints (1,1) double {mustBePositive} = 20000
    opts.smoothsec (1,1) double {mustBePositive} = 5
    opts.tffreqlim (1,2) double {mustBePositive} = [1 40]
    opts.tfnfreq   (1,1) double {mustBePositive} = 40
    opts.tfcycles  (1,2) double {mustBePositive} = [3 10]
    opts.tfnorm    (1,:) char {mustBeMember(opts.tfnorm, ...
        {'db', 'percent', 'raw'})} = 'db'
    opts.tfsmooth  (1,1) double {mustBeNonnegative} = 1
    opts.tfwindow  (1,1) double {mustBePositive} = 4
    opts.tftw      (1,1) double {mustBePositive} = 3
    opts.tftapers  (1,1) double {mustBeNonnegative} = 0
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

%%% The '_signal' suffix is a display choice, not a different transform, so it is
%%% split off here and everything downstream sees the plain method name.
issig  = endsWith(opts.acttype, '_signal') && showact;
method = erase(opts.acttype, '_signal');
istf   = any(strcmp(method, {'wavelet', 'multitaper'})) && showact;

%%% Two columns per component: that makes the grid divisible both by the number
%%% of maps in row 2 and by the two panels sharing row 1, whatever ncomps is.
%%% The spectrograms need one more for their colorbar, and the columns are
%%% counted in thirds so that carving it out costs the traces a fiftieth of the
%%% width rather than a ninth.
unit     = 1 + 2 * double(istf);
ncols    = 2 * ncomps * unit;                   % columns the panels span
gridcols = ncols + double(istf);                % columns the grid actually has
halfrow  = ncols / 2;                           % the two panels sharing row 1
topoWide = 2 * unit;                            % one component map
topoSpan = 2;                                   % maps get double height
nrows    = 1 + topoSpan + ncomps * double(showact);

%%% A signal strip takes a fifth of its component's block, so the figure grows by
%%% a quarter of the block to make room for it - otherwise the strips would be
%%% carved out of the spectrograms and every one of them would come out shorter
%%% than it is without them.
sizerows = nrows + 0.25 * ncomps * double(issig);

fig = figure('Color', 'w', 'Name', 'GED diagnostics', 'NumberTitle', 'off', ...
    'Position', [80 80 1250 min(1250, 220 + 105 * sizerows)]);
tl = tiledlayout(nrows, gridcols, 'TileSpacing', 'compact', 'Padding', 'compact');

palette  = gedpalette(ncomps);
cmapAxes = gobjects(0);

%% ----------------------------------------- eigenspectrum and power spectra
%%% Eigenspectrum. The elbow says how many directions actually separate S from R;
%%% the dashed line is the permutation threshold, when one was computed.
nexttile(1, [1 halfrow]);
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
axS = nexttile(1 + halfrow, [1 halfrow]);
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
    ax = nexttile(gridcols + (i - 1) * topoWide + 1, [topoSpan topoWide]);
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
    %%% Scale the axis to what is actually on screen: with a 30 s window over an
    %%% 8 h recording, labelling in hours would put every tick at 0.00.
    tref = t(end);
    if ~isempty(opts.xwindow), tref = opts.xwindow; end
    if tref > 7200,    tscale = 3600; tlabel = 'Time (h)';
    elseif tref > 120, tscale = 60;   tlabel = 'Time (min)';
    else,              tscale = 1;    tlabel = 'Time (s)';
    end
    tplot  = t / tscale;
    xlims  = [0 max(tplot(end), eps)];
    axTime = gobjects(1, ncomps);

    %%% The window the axes will actually open on, worked out here rather than
    %%% left to addtimescroll at the end. Setting the full range first and
    %%% narrowing it afterwards costs a whole-recording spectrogram that is
    %%% thrown away before anyone sees it - which for a night of data is a
    %%% minute of waiting for a view that was never asked for.
    winlen = opts.xwindow / tscale;
    xview  = xlims;
    if ~isempty(winlen) && isscalar(winlen) && isfinite(winlen) && ...
            winlen > 0 && winlen < diff(xlims)
        xview = xlims(1) + [0 winlen];
    end

    %%% One spectrogram column per pixel or so. maxpoints is a budget for min/max
    %%% pairs along a line and is far more than an image can show, so it only
    %%% acts as a ceiling here.
    tf   = struct('method', method, ...
        'freqlim', opts.tffreqlim, 'nfreq', opts.tfnfreq, ...
        'cycles', opts.tfcycles, 'norm', opts.tfnorm, 'smooth', opts.tfsmooth, ...
        'window', opts.tfwindow, 'tw', opts.tftw, 'ntapers', opts.tftapers, ...
        'ncols', min(opts.maxpoints, 2000), 'xinit', xview);
    tfclim = zeros(0, 2);

    axSig = gobjects(1, ncomps * double(issig));

    for i = 1:ncomps
        row = 1 + topoSpan + i;
        if issig
            %%% A layout of its own for the pair. TileSpacing is a property of a
            %%% whole layout, so the gap between a spectrogram and its trace can
            %%% only be closed by giving the two a layout to themselves - doing
            %%% it on the parent would close up every other gap in the figure
            %%% along with it. Five rows split four to one is the 80/20.
            blk = tiledlayout(tl, 5, 1, 'TileSpacing', 'none', 'Padding', 'none');
            blk.Layout.Tile     = (row - 1) * gridcols + 1;
            blk.Layout.TileSpan = [1 ncols];
            ax   = nexttile(blk, 1, [4 1]);
            asig = nexttile(blk, 5, [1 1]);
        else
            ax   = nexttile(tl, (row - 1) * gridcols + 1, [1 ncols]);
            asig = gobjects(0);
        end

        if istf
            %%% single: the series is kept alive behind every one of these axes
            %%% for the life of the figure, since the transform re-runs on each
            %%% change of the x limits, and a night per component adds up.
            tfclim(end + 1, :) = drawtfmap(ax, single(act(i, :)), srate, tscale, tf); %#ok<AGROW>
            cmapAxes(end + 1)  = ax; %#ok<AGROW>
            ylabel(ax, sprintf('#%d (Hz)', comps(i)), 'Color', palette(i, :), ...
                'FontWeight', 'bold');
        else
            if drawactivation(ax, tplot, act(i, :), srate, method, palette(i, :), ...
                    opts.maxpoints, opts.smoothsec)
                cmapAxes(end + 1) = ax; %#ok<AGROW>
            end
            ylabel(ax, sprintf('#%d', comps(i)), 'Color', palette(i, :), 'FontWeight', 'bold');
        end
        xlim(ax, xview);
        axTime(i) = ax;

        if issig
            %%% The trace the power above was computed from, decimated the same
            %%% way the 'signal' mode does it - so it stays honest at screen
            %%% resolution and redraws at full detail as the window moves.
            drawactivation(asig, tplot, act(i, :), srate, 'signal', palette(i, :), ...
                opts.maxpoints, opts.smoothsec);
            %%% No y ticks: a fifth of a panel has no room to read a scale off,
            %%% and the strip is here to show the shape of the trace - artefacts,
            %%% clipping, what the spectrogram was actually made from.
            set(asig, 'YTick', []);
            %%% Its own limits, and robust ones. drawactivation fits the axis
            %%% tightly to the whole recording, so one artefact sets the scale
            %%% and every ordinary stretch collapses to a flat line - which at a
            %%% fifth of a panel's height leaves nothing to see at all. Trimming
            %%% the outer fifth of a percent fixes that, and keeping the limits
            %%% fixed rather than refitting them per window keeps the strip
            %%% honest while it scrolls. Sorting a subsample, not the whole row:
            %%% a percentile needs nothing like a night of data to be steady.
            sub = sort(act(i, 1:max(1, floor(size(act, 2) / 20000)):end));
            k   = max(1, round(0.002 * numel(sub)));
            lo  = sub(k);
            hi  = sub(end - k + 1);
            ylim(asig, [lo hi] + max(0.05 * (hi - lo), eps) * [-1 1]);
            xlim(asig, xview);
            axSig(i) = asig;
            %%% The trace carries the time axis for the pair, so the spectrogram
            %%% never labels one - its ticks would collide with the panel below.
            set(ax, 'XTickLabel', []);
        end

        bottom = axTime(i);
        if issig, bottom = axSig(i); end
        if i < ncomps
            set(bottom, 'XTickLabel', []);
        else
            xlabel(bottom, tlabel);
        end
    end

    %%% Link the time axes: scrolling or zooming one moves them all. For the
    %%% trace modes the amplitude axis is linked too, which is what makes their
    %%% heights comparable; a spectrogram's y axis is frequency, so linking it
    %%% would say nothing and its own limits are already identical.
    %%% linkprop rather than linkaxes, because plotgednight has to link two
    %%% different sets of axes on two different properties and linkaxes would
    %%% overwrite one with the other; the handle has to be kept alive, hence the
    %%% appdata.
    %%% One x link over everything that carries time - the spectrograms and, when
    %%% they are there, the traces underneath them - so a scroll moves the whole
    %%% stack as one.
    %%%
    %%% The amplitude link is a different matter. It goes on the trace modes,
    %%% where a shared range is what makes one component's height mean the same
    %%% as another's. It deliberately does not go on the signal strips: a
    %%% spectrogram's y axis is frequency, so there is nothing to link there, and
    %%% each strip is normalised against its own recording anyway - exactly as
    %%% the spectrogram above it is. Forcing one range across all of them would
    %%% only let the loudest component flatten the rest.
    xaxes = [axTime axSig];
    yaxes = gobjects(0);
    if ~istf, yaxes = axTime; end

    links = {};
    if numel(xaxes) > 1
        links{end + 1} = linkprop(xaxes, 'XLim');
    end
    if numel(yaxes) > 1
        links{end + 1} = linkprop(yaxes, 'YLim');
        %%% A common amplitude range, set once: each axis was scaled to its own
        %%% component, so linking alone would hand every trace whichever range
        %%% happened to be set last.
        yl = get(yaxes, 'YLim');
        yl = vertcat(yl{:});
        ylim(yaxes(1), [min(yl(:, 1)) max(yl(:, 2))]);
    end
    setappdata(fig, 'gedAxisLinks', links);

    addtimescroll(fig, tl, xaxes, opts.xwindow / tscale, xlims);

    %%% One colour scale for every component, and one bar to read it by. Each
    %%% panel proposed limits from its own recording; the median of those keeps a
    %%% single odd component from stretching the scale flat for the rest.
    %%%
    %%% The bar is given its own column, spanning the trace rows and no others,
    %%% so it starts at the top of the first spectrogram and stops at the bottom
    %%% of the last - rather than running the height of the figure past the maps
    %%% and the eigenspectrum, which it says nothing about.
    %%%
    %%% Placed by tile rather than by setting Position from the measured extent
    %%% of the axes: a tiled axes does not report a Position that follows the
    %%% layout, so anything measured off it is a stale number and the bar lands
    %%% off the edge of the figure.
    %%%
    %%% It hangs off an invisible axes of its own rather than off one of the
    %%% spectrograms, because with the signal strips those sit inside nested
    %%% layouts - and a colorbar reads Layout.Tile against the layout its own
    %%% axes belongs to, which would be the nested one, five tiles tall. A parent
    %%% tile index means nothing there and the bar silently disappears.
    if istf
        clim   = [median(tfclim(:, 1)) median(tfclim(:, 2))];
        cbTile = (1 + topoSpan) * gridcols + gridcols;
        set(axTime, 'CLim', clim);

        axCB = nexttile(tl, cbTile, [ncomps 1]);
        axis(axCB, 'off');
        set(axCB, 'CLim', clim);
        cmapAxes(end + 1) = axCB;

        cb = colorbar(axCB);
        cb.Layout.Tile     = cbTile;
        cb.Layout.TileSpan = [ncomps 1];
        switch opts.tfnorm
            case 'db',      ylabel(cb, 'Power (dB vs. median)');
            case 'percent', ylabel(cb, 'Power (% of median)');
            case 'raw',     ylabel(cb, 'Power (dB)');
        end
    end
end

%%% Set after every topoplot call, not inside the loop: EEGLAB's topoplot applies
%%% its own colormap and background from icadefs to the figure it lands in, and
%%% would otherwise overwrite whatever the earlier tiles were given.
for ax = [topoAxes cmapAxes]
    colormap(ax, cmap);
end
set(fig, 'Color', 'w');
end

