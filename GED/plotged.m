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
%                'signal_wavelet', 'signal_multitaper'
%                            the same pair with the vertical ratio switched:
%                            the trace takes four fifths of the block and the
%                            spectrogram shrinks to a strip underneath it.
%                            Reach for this the other way round from above -
%                            when the trace is what needs reading closely and
%                            the spectrogram is there only for context.
%              The first three are the same modes as plotgednight.
%   showact    Draw the activation traces at all. Default: true.
%   events     Logical matrix [components x samples] marking events with true.
%              Same number of columns as the activations, and indexed by
%              component number down the rows - row 3 marks component #3
%              whether or not it is the third one shown. Drawn over that
%              component's activation trace or spectrogram in its own colour: a
%              run of true samples as a shaded band, an isolated true sample as
%              a thin line. For the '..._signal' and 'signal_...' pairings,
%              drawn on the signal strip as well. Default: none.
%   scoring    One sleep-stage digit per epoch, as scoreloader returns them:
%              -3 N3, -2 N2, -1 N1, 0 Wake, 1 REM. Given it, a hypnogram is
%              drawn full width above the first activation trace, spanning the
%              whole recording. It does not scroll with xwindow; instead a
%              vertical line sweeps across it to show where the visible window
%              sits in the night. Default: none.
%   epochlength  Scoring epoch length in seconds, for the hypnogram. Default: 30.
%   keptepochs   Epoch indices that entered the GED. Shaded as a band in the
%              hypnogram, so it is visible which part of the night the filters
%              were built on. If GED.comp is itself the concatenation of those
%              epochs (as after bidsfun_subcomp), the scoring is subset with
%              them so the hypnogram still lines up with the traces. Default: none.
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
%     plotged(GEDs(1).ged, 'acttype', 'signal_wavelet', 'xwindow', 30)
%                                                  % same pair, trace given the room
%     ev = false(size(GEDs(1).ged.comp));
%     ev(1, spindlepeaks) = true;                  % samples to mark on #1
%     plotged(GEDs(1).ged, 'events', ev, 'xwindow', 30)
%     plotged(GEDs(1).ged, 'scoring', scoring, 'keptepochs', GEDs(1).keptepochs, ...
%         'xwindow', 30)                           % hypnogram on top, cursor tracks
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
         'wavelet_signal', 'multitaper_signal', ...
         'signal_wavelet', 'signal_multitaper'})} = 'signal'
    opts.showact   (1,1) logical = true
    opts.events          logical = logical([])
    opts.scoring         double = []
    opts.epochlength (1,1) double {mustBePositive} = 30
    opts.keptepochs      double = []
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

%%% The 'signal' pairing is a display choice, not a different transform, so it
%%% is split off here and everything downstream sees the plain method name.
%%% It can be named either way round - '<method>_signal' or 'signal_<method>' -
%%% and swapsig records which, since that also decides which of the pair gets
%%% the four-fifths share of the block below.
swapsig = startsWith(opts.acttype, 'signal_') && showact;
issig   = (swapsig || endsWith(opts.acttype, '_signal')) && showact;
if swapsig
    method = erase(opts.acttype, 'signal_');
else
    method = erase(opts.acttype, '_signal');
end
istf   = any(strcmp(method, {'wavelet', 'multitaper'})) && showact;

%%% Event overlay (optional): a logical [components x samples] mask, indexed by
%%% component number. Checked here so a size mistake is caught before anything
%%% is drawn, and dropped when there are no activations to draw it on.
events = opts.events;
if ~isempty(events) && showact
    if size(events, 2) ~= size(GED.comp, 2)
        error('plotged:eventsSize', ...
            ['events must have the same number of samples as the activations ' ...
             '(%d); got %d.'], size(GED.comp, 2), size(events, 2));
    end
    if size(events, 1) < max(comps)
        error('plotged:eventsRows', ...
            ['events has %d row(s), but component #%d is being plotted; give it ' ...
             'one row per component.'], size(events, 1), max(comps));
    end
elseif ~showact
    events = logical([]);
end

%%% Hypnogram (optional): drawn only when there are traces for it to sit above.
%%% The scoring and the activations have to describe the same stretch of
%%% recording; when GED.comp is a concatenation of the kept epochs, keptepochs
%%% repairs the mismatch by subsetting the scoring the same way - the same
%%% reconciliation plotgednight does.
scoring    = opts.scoring(:)';
keptepochs = opts.keptepochs;
hasHypno   = ~isempty(scoring) && showact;
if hasHypno
    nEpochsData = floor(size(GED.comp, 2) / (opts.epochlength * srate));
    if abs(numel(scoring) - nEpochsData) > 2
        if ~isempty(keptepochs) && numel(keptepochs) == nEpochsData
            scoring    = scoring(keptepochs);
            keptepochs = [];        % already applied; nothing left to shade
            fprintf(['plotged: scoring subset to the %d kept epoch(s) to match the ' ...
                     'stage-selected recording.\n'], numel(scoring));
        else
            warning('plotged:scoringMismatch', ...
                ['the scoring covers %d epoch(s) but the data covers %d; the ' ...
                 'hypnogram and the traces will not line up. Pass keptepochs, or ' ...
                 'the scoring for exactly the data ged() was given.'], ...
                numel(scoring), nEpochsData);
        end
    end
end

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
hypRows  = double(hasHypno);                    % the whole-night hypnogram
nrows    = 1 + topoSpan + hypRows + ncomps * double(showact);

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

    %%% The hypnogram sits directly above the first trace and shows the whole
    %%% night at once. It is deliberately kept out of the scroll link below, so
    %%% xwindow pages the traces while the hypnogram holds still; a cursor line,
    %%% added once the traces exist, then tracks the visible window across it.
    %%%
    %%% Its own time base, not the traces': with a 30 s window the traces label
    %%% in seconds, which over a whole night would read 0, 500, 1000, ... Pick
    %%% the unit from the night's length instead. hyp2plot rescales a position
    %%% on the trace axis to this one, and is 1 whenever no window is asked for.
    axHyp = gobjects(0);
    if hasHypno
        tnight = t(end);
        if tnight > 7200,    hypscale = 3600; hyplabel = 'Time (h)';
        elseif tnight > 120, hypscale = 60;   hyplabel = 'Time (min)';
        else,                hypscale = 1;    hyplabel = 'Time (s)';
        end
        hyp2plot = tscale / hypscale;

        hrow  = 2 + topoSpan;
        axHyp = nexttile(tl, (hrow - 1) * gridcols + 1, [1 ncols]);
        drawhypnogram(axHyp, scoring, opts.epochlength / hypscale, keptepochs, 1);
        xlim(axHyp, [0 max(tnight / hypscale, eps)]);
        set(axHyp, 'XTickLabel', []);
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
        row = 1 + topoSpan + hypRows + i;
        if issig
            %%% A layout of its own for the pair. TileSpacing is a property of a
            %%% whole layout, so the gap between a spectrogram and its trace can
            %%% only be closed by giving the two a layout to themselves - doing
            %%% it on the parent would close up every other gap in the figure
            %%% along with it. Five rows split four to one is the 80/20.
            blk = tiledlayout(tl, 5, 1, 'TileSpacing', 'none', 'Padding', 'none');
            blk.Layout.Tile     = (row - 1) * gridcols + 1;
            blk.Layout.TileSpan = [1 ncols];
            %%% 'signal_wavelet'/'signal_multitaper' switch this ratio round
            %%% from 'wavelet_signal'/'multitaper_signal': the trace becomes
            %%% the four-fifths panel and the spectrogram the strip beneath it.
            if swapsig
                asig = nexttile(blk, 1, [4 1]);
                ax   = nexttile(blk, 5, [1 1]);
            else
                ax   = nexttile(blk, 1, [4 1]);
                asig = nexttile(blk, 5, [1 1]);
            end
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
            %%% clipping, what the spectrogram was actually made from. swapsig
            %%% hands the trace four fifths instead, precisely so its scale can
            %%% be read, so there the ticks stay.
            if ~swapsig
                set(asig, 'YTick', []);
            end
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
            %%% Whichever of the pair sits at the bottom carries the time axis;
            %%% the one on top never labels one - its ticks would collide with
            %%% the panel below. swapsig puts the trace on top instead of the
            %%% spectrogram, so it is the one cleared in that case.
            if swapsig
                set(asig, 'XTickLabel', []);
            else
                set(ax, 'XTickLabel', []);
            end
        end

        bottom = axTime(i);
        if issig && ~swapsig, bottom = axSig(i); end
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

    %%% The "you are here" cursor on the hypnogram. Only drawn when the traces
    %%% show a window rather than the whole night; it hangs off the linked time
    %%% axis, so the scrollbar, a zoom and a pan all carry it alike.
    if hasHypno && ~isequal(xview, xlims)
        cur = xline(axHyp, mean(xview) * hyp2plot, '-', 'Color', [0.10 0.45 0.95], ...
            'LineWidth', 2);
        cur.HandleVisibility = 'off';
        addlistener(axTime(1), 'XLim', 'PostSet', ...
            @(~, ~) movecursor(cur, axTime(1), hyp2plot));

        %%% Click anywhere on the hypnogram to jump the window there, centred on
        %%% the point clicked. The stairs and bands are made click-through so the
        %%% event always reaches the axes; the linked time axis and its listeners
        %%% (this cursor, the scrollbar) then follow from the one xlim change.
        seeklen = diff(xview);
        axHyp.ButtonDownFcn = @(src, ~) seekhypno(src, axTime(1), hyp2plot, seeklen, xlims);
        set(allchild(axHyp), 'HitTest', 'off');
        %%% The hypnogram runs on its own time base now, so it keeps its tick
        %%% labels - they are how the cursor's position is read off. Without a
        %%% window the two axes match and the labels would only be a duplicate
        %%% row, so there they stay hidden.
        set(axHyp, 'XTickLabelMode', 'auto');
        xlabel(axHyp, hyplabel);
    end

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
        cbTile = (1 + topoSpan + hypRows) * gridcols + gridcols;
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

    %%% Event overlay, last of all - once every axis has its final y-limits, so
    %%% the bands can be built to span exactly those and nothing they draw can
    %%% drag the amplitude scale the way a placeholder range added mid-loop
    %%% would. Drawn on the trace/spectrogram, and on the signal strip too when
    %%% there is one - a pairing has the marks read against both views.
    if ~isempty(events)
        s2p = 1 / (srate * tscale);
        for i = 1:ncomps
            drawevents(axTime(i), events(comps(i), :), s2p, palette(i, :));
            if issig
                drawevents(axSig(i), events(comps(i), :), s2p, palette(i, :));
            end
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


function drawevents(ax, mask, s2p, col)
%%% One line for the point events and one patch for the spans - never an object
%%% per marked sample. A night of detections is thousands of xlines, each with
%%% its own limit listeners, and every scroll then redraws the lot and the
%%% figure crawls. Here contiguous true runs collapse to a single patch face
%%% each and isolated true samples share one NaN-separated line, so a component
%%% costs two graphics objects however many events it has.
%%%
%%% s2p converts a sample index to the axis' plot-time units: t = (k-1)*s2p.
%%% Lines and span edges are a darkened opaque form of the component colour so
%%% they read over a saturated spectrogram as well as over a pale trace; the
%%% span interior is the plain colour at low alpha. The markers are built to
%%% span the axis' current y-limits and re-spanned by one listener if those are
%%% later zoomed; they never write YLim themselves, so they cannot pull the
%%% amplitude scale.
mask = logical(mask(:)');
if ~any(mask), return; end

edges  = diff([false, mask, false]);
starts = find(edges == 1);
stops  = find(edges == -1) - 1;
ispt   = starts == stops;

yl      = ax.YLim;
cdark   = 0.55 * col;
handles = gobjects(0);

if any(ispt)
    x = (starts(ispt) - 1) * s2p;
    n = numel(x);
    handles(end + 1) = line(ax, reshape([x; x; nan(1, n)], [], 1), ...
        repmat([yl(1); yl(2); NaN], n, 1), ...
        'Color', cdark, 'LineWidth', 0.75, 'HandleVisibility', 'off');
end
if any(~ispt)
    xa = (starts(~ispt) - 1) * s2p;
    xb = (stops(~ispt)  - 1) * s2p;
    m  = numel(xa);
    handles(end + 1) = patch(ax, 'XData', [xa; xb; xb; xa], ...
        'YData', repmat([yl(1); yl(1); yl(2); yl(2)], 1, m), ...
        'FaceColor', col, 'FaceAlpha', 0.30, 'EdgeColor', cdark, ...
        'LineWidth', 0.75, 'HandleVisibility', 'off');
end

addlistener(ax, 'YLim', 'PostSet', @(~, ~) respanevents(ax, handles));
end

% -------------------------------------------------------------------------
function respanevents(ax, handles)
% Stretch the event markers to the axis' current y-limits, without touching
% YLim itself - so a zoom on one linked trace re-spans every component's marks.
yl = ax.YLim;
for k = 1:numel(handles)
    h = handles(k);
    if ~isgraphics(h), continue, end
    y = h.YData;
    if strcmp(h.Type, 'line')
        y(1:3:end) = yl(1);
        y(2:3:end) = yl(2);
    else
        y([1 2], :) = yl(1);
        y([3 4], :) = yl(2);
    end
    h.YData = y;
end
end

% -------------------------------------------------------------------------
function movecursor(cur, ax, k)
% Keep the hypnogram's "you are here" line on the centre of the visible window,
% rescaled from the trace axis' time unit to the hypnogram's own.
if isgraphics(cur) && isgraphics(ax)
    cur.Value = mean(ax.XLim) * k;
end
end

% -------------------------------------------------------------------------
function seekhypno(axHyp, axTrace, k, winlen, xfull)
% Jump the scroll window to the point clicked on the hypnogram, centred there
% and clamped to the recording. k rescales the hypnogram's time unit back to
% the trace axis'; the linked axes and their listeners follow from the xlim.
if ~isgraphics(axHyp) || ~isgraphics(axTrace), return, end
cx = axHyp.CurrentPoint(1, 1) / k;
lo = max(xfull(1), min(cx - winlen / 2, xfull(2) - winlen));
xlim(axTrace, [lo lo + winlen]);
end

