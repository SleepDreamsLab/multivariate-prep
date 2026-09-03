function clim = drawtfmap(ax, x, srate, tscale, tf)
% DRAWTFMAP  One component's time-frequency power into one axes, computed lazily.
%
%   Used by plotged for acttype 'wavelet' and 'multitaper'. Returns the colour
%   limits it would like, so the caller can settle on one scale for every
%   component - which is what makes the panels comparable to each other.
%
%   Which of the two transforms runs is the tf.method field, and it is the only
%   thing that differs between them here: both put a column at the same instant,
%   both are normalised per frequency the same way, and both are recomputed on
%   the same events.
%
%   Nothing is precomputed over the whole recording. The transform runs over the
%   visible window only, at whatever resolution that window can show, and again
%   whenever the x limits change - from the scrollbar, from zooming, from panning
%   alike. That is the only way both ends of the range stay usable: a whole night
%   held at the resolution a 30 s window wants would be hundreds of megabytes per
%   component, and one held at the resolution a whole night wants would be a
%   hundred blocky columns once you scrolled into it. Computing 30 s of wavelet
%   convolution takes a few milliseconds, so paying per view is cheaper than
%   storing either.
%
%   ax       axes to draw into
%   x        the component time series (single is fine, and is what the caller
%            should hand over - it is kept alive for the life of the figure)
%   srate    sampling rate, Hz
%   tscale   seconds per x-axis unit, so the axes can be labelled in minutes or
%            hours while the transform still works in samples
%   tf       struct with fields method ('wavelet' or 'multitaper'), freqlim,
%            nfreq, norm, ncols, xinit (the x limits the figure will open on, in
%            axis units), cycles and smooth for the wavelet, window, tw and
%            ntapers for the multitaper

arguments
    ax
    x      (1,:) {mustBeNumeric, mustBeReal}
    srate  (1,1) double {mustBePositive}
    tscale (1,1) double {mustBePositive}
    tf     (1,1) struct
end

%%% Built here, not left to tfmorlet, because the axis has to be labelled with
%%% the frequencies before any of them have been computed.
hi    = min(tf.freqlim(2), srate / 2);
lo    = min(tf.freqlim(1), hi);
freqs = linspace(lo, hi, tf.nfreq);
nfreq = numel(freqs);

%%% Everything a transform needs, gathered once, so the two places that run one
%%% - the reference just below and every later redraw - cannot drift apart.
p       = tf;
p.x     = x;
p.srate = srate;
p.freqs = freqs;

[ref, clim] = referencepower(p);

%%% The axis is Hz, on a linear scale, and the rows are evenly spaced in Hz to
%%% match - which is what imagesc assumes. The two go together: a log-spaced bank
%%% drawn against a linear axis would put every row at the wrong height, so the
%%% spacing of the bank is not a free choice once the scale is fixed.
step = 0;
if nfreq > 1, step = (freqs(end) - freqs(1)) / (nfreq - 1); end
ydata = [freqs(1) freqs(end)];

h = imagesc(ax, [0 1], ydata, nan(nfreq, 2));
set(ax, 'YDir', 'normal', 'CLim', clim);
%%% Half a row of padding at each end: imagesc places a row at its centre, so
%%% without it the top and bottom rows are cut in half.
ylim(ax, [freqs(1) - step/2, freqs(end) + step/2]);
%%% Before the listener exists, so the first transform runs over the window the
%%% figure is about to open on and not over whatever imagesc left behind.
xlim(ax, tf.xinit);
box(ax, 'off');

d = struct();
d.p      = p;
d.tscale = tscale;
d.ydata  = ydata;
d.norm   = tf.norm;
d.ncols  = tf.ncols;
d.ref    = ref;
d.h      = h;
d.n      = numel(x);
d.i1     = 0;
d.i2     = -1;
d.dec    = 0;
setappdata(ax, 'gedTF', d);
%%% The listener's lifetime is tied to the axes, so it needs no other owner -
%%% same arrangement as the traces in drawactivation.
addlistener(ax, 'XLim', 'PostSet', @(~, ~) refreshtf(ax));
refreshtf(ax);
end

% -------------------------------------------------------------------------
function refreshtf(ax)
% Recompute the visible window, coalescing the burst of events a scrollbar drag
% produces. Each one of those would otherwise queue a transform that is already
% stale by the time it runs, and the figure would lag a second behind the thumb.

if getappdata(ax, 'gedTFbusy')
    setappdata(ax, 'gedTFpending', true);
    return
end
setappdata(ax, 'gedTFbusy', true);
done = onCleanup(@() clearbusy(ax));

again = true;
while again
    setappdata(ax, 'gedTFpending', false);
    updatetf(ax);
    again = isgraphics(ax) && isequal(getappdata(ax, 'gedTFpending'), true);
end
end

% -------------------------------------------------------------------------
function clearbusy(ax)
% The axes can be gone by now, if the figure was closed mid-transform.

if isgraphics(ax), setappdata(ax, 'gedTFbusy', false); end
end

% -------------------------------------------------------------------------
function updatetf(ax)

d = getappdata(ax, 'gedTF');
if isempty(d) || ~isgraphics(d.h), return, end

%%% The visible stretch, in samples. Uniform sampling makes this arithmetic - no
%%% search over a multi-million-sample time vector on every pan event.
xl = xlim(ax);
i1 = max(1,   floor(xl(1) * d.tscale * d.p.srate) + 1);
i2 = min(d.n, ceil( xl(2) * d.tscale * d.p.srate) + 1);
if i2 - i1 < 2, return, end

%%% Roughly one column per pixel: more is invisible, less is blocky. The bins
%%% double as the anti-alias filter for the decimation, which is why tfmorlet
%%% averages the power into them rather than sampling it.
dec = max(1, floor((i2 - i1 + 1) / d.ncols));

%%% What is already on the axes may cover the new window: the image is drawn
%%% wider than the view, so a small scroll needs no transform at all.
if i1 >= d.i1 && i2 <= d.i2 && dec == d.dec
    return
end

%%% Half a window of margin on either side, which is what makes that check pay
%%% off while still only transforming a fraction more data than is shown.
margin = round(0.25 * (i2 - i1 + 1));
j1 = max(1,   i1 - margin);
j2 = min(d.n, i2 + margin);

[pow, t] = transform(d.p, [j1 j2], d.p.srate / dec);
if numel(t) < 2, return, end

set(d.h, 'XData', [t(1) t(end)] / d.tscale, 'YData', d.ydata, ...
    'CData', normalisepower(pow, d.ref, d.norm));

d.i1 = j1; d.i2 = j2; d.dec = dec;
setappdata(ax, 'gedTF', d);
end

% -------------------------------------------------------------------------
function [pow, t] = transform(p, range, outrate)
% Run whichever transform this axes was built for. The single place the two
% methods are told apart, so everything downstream - normalisation, colour
% limits, caching, the time axis - is shared by construction.

switch p.method
    case 'multitaper'
        [pow, t] = tfmulti(p.x, p.srate, 'freqs', p.freqs, 'window', p.window, ...
            'tw', p.tw, 'ntapers', p.ntapers, 'range', range, 'outrate', outrate);
    otherwise
        [pow, t] = tfmorlet(p.x, p.srate, 'freqs', p.freqs, 'cycles', p.cycles, ...
            'range', range, 'outrate', outrate, 'smoothsigma', p.smooth);
end
end

% -------------------------------------------------------------------------
function C = normalisepower(pow, ref, mode)
% Power on a readable scale. Raw wavelet power is unreadable across a broad
% band: 1/f alone leaves the bottom rows a solid bright bar and flattens
% everything above them, whatever is actually happening up there. Dividing each
% row by its own reference is what makes a spindle burst at 13 Hz and a slow
% wave at 1 Hz show up as the same kind of feature.

pow = double(pow);
switch mode
    case 'db',      C = 10 * log10(pow ./ ref);
    case 'percent', C = 100 * (pow ./ ref - 1);
    otherwise,      C = 10 * log10(pow);          % 'raw'
end
end

% -------------------------------------------------------------------------
function [ref, clim] = referencepower(p)
% The per-frequency reference and the colour limits, fixed once, up front.
%
% Both have to describe the whole recording rather than whichever window happens
% to be on screen first. A colour scale that refitted itself on every scroll
% would make no two windows comparable, and a reference taken from one window
% would call that window average by construction.
%
% Sampled rather than transformed in full: half a dozen minutes spread across
% the night, which costs a fraction of a second, and a median over them, which
% is a steadier estimate than a mean over everything would have been anyway -
% arousals and movement artefacts move a mean and leave a median alone.

n     = numel(p.x);
snip  = min(n, round(30 * p.srate));
nsnip = max(1, min(12, floor(n / snip)));
start = round(linspace(1, max(1, n - snip + 1), nsnip));

%%% 5 Hz is ample for an estimate pooled over minutes, and cuts the work. Run
%%% through the same dispatcher as the display, so the reference is made by the
%%% method it will be dividing.
pow = transform(p, [start(:) start(:) + snip - 1], 5);

ref = median(double(pow), 2);
ref(~isfinite(ref) | ref <= 0) = eps;

C = normalisepower(pow, ref, p.norm);
C = C(isfinite(C));
if strcmp(p.norm, 'raw')
    %%% Unreferenced power has no meaningful centre, so the limits just bracket
    %%% the data.
    clim = quantiles(C, [0.02 0.98]);
else
    %%% Symmetric about zero, because zero is where the reference sits and half
    %%% the recording lies on either side of it. A scale that were not centred
    %%% would put the midpoint of a diverging colormap - the white of a bwr, the
    %%% pale middle of any of them - somewhere other than "typical for this
    %%% frequency", and every panel would read as biased one way.
    clim = quantiles(abs(C), 0.98) * [-1 1];
end
if ~(clim(2) > clim(1)), clim = clim(1) + [0 1]; end
end

% -------------------------------------------------------------------------
function v = quantiles(x, p)
% Quantiles without the Statistics toolbox - trimming the tails keeps a handful
% of artefact samples from washing out the whole colour scale.

x = sort(x(:));
if isempty(x), v = [0 1]; return, end
v = x(min(numel(x), max(1, round(p * numel(x)))))';
end

%%% No custom frequency ticks any more: on a linear axis MATLAB's own choice is
%%% already round numbers, and it re-picks them when the panel is resized.
