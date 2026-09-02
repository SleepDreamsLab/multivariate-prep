function isheat = drawactivation(ax, tplot, x, srate, acttype, color, maxpoints, smoothsec)
% DRAWACTIVATION  One component's time course into one axes.
%
%   Shared by plotged and plotgednight. Returns true when the axes was drawn as
%   a heat strip, so the caller knows it still needs the component colormap
%   applied (EEGLAB's topoplot would otherwise overwrite it - see the note where
%   the callers re-colour at the end).
%
%   acttype 'signal'    the component time series itself, min/max decimated
%           'envelope'  its smoothed absolute amplitude, as a line
%           'heat'      the same envelope, as a colour strip
%
%   Decimation follows the view, not the recording. The full-resolution trace is
%   kept on the axes and redrawn from it whenever the x limits change, so a
%   window always gets the full maxpoints budget. Decimating once over the whole
%   recording instead would leave a 30 s window holding only the handful of
%   points that happened to fall inside it, and the min/max pairs of the
%   decimator would then read as a spiky oscillation that is not in the data.

arguments
    ax
    tplot     double
    x         double
    srate     (1,1) double
    acttype   (1,:) char
    color     (1,3) double
    maxpoints (1,1) double = 20000
    smoothsec (1,1) double = 5
end

%%% Both non-raw modes show the same quantity - the smoothed absolute amplitude.
%%% Across long recordings the raw trace collapses into a solid band, while its
%%% envelope shows when the component actually switches on.
isheat = strcmp(acttype, 'heat');
if isheat || strcmp(acttype, 'envelope')
    y = movmean(abs(x), max(1, round(smoothsec * srate)));
else
    y = x;
end

if isheat
    [td, cd] = binaverage(tplot, y, maxpoints);
    h = imagesc(ax, td, [0 1], cd);
    set(ax, 'YTick', []);
    %%% Colour limits are fixed to the whole recording, so scrolling does not
    %%% silently rescale the colours under the reader.
    clim(ax, [0 prctileish(y, 99)]);
else
    [td, yd] = minmaxdecimate(tplot, y, maxpoints);
    h = plot(ax, td, yd, 'Color', color, 'LineWidth', 0.7 + 0.4 * strcmp(acttype, 'envelope'));
    axis(ax, 'tight');
    if strcmp(acttype, 'envelope')
        %%% An amplitude is measured from zero, so the axis starts there: with a
        %%% tight lower limit, ordinary ripple in a flat envelope would look like
        %%% a component switching on and off.
        ylim(ax, [0 max(max(yd) * 1.05, eps)]);
    end
end
box(ax, 'off');

%%% Keep the full-resolution trace with the axes and re-decimate on every limit
%%% change - from the scrollbar, from zooming, from panning alike.
dt = 0;
if numel(tplot) > 1, dt = tplot(2) - tplot(1); end
setappdata(ax, 'gedTrace', struct('t0', tplot(1), 'dt', dt, 'n', numel(y), ...
    'y', y, 'h', h, 'maxpoints', maxpoints, 'isheat', isheat));
%%% The listener's lifetime is tied to the axes, so it needs no other owner.
addlistener(ax, 'XLim', 'PostSet', @(~, ~) refreshtrace(ax));
end

% -------------------------------------------------------------------------
function refreshtrace(ax)
% Redraw the visible stretch at full detail.

d = getappdata(ax, 'gedTrace');
if isempty(d) || ~isgraphics(d.h) || d.dt <= 0
    return
end

%%% Samples are uniformly spaced, so the visible range is arithmetic - no search
%%% over a multi-million-sample vector on every pan event.
xl = xlim(ax);
i1 = max(1,   floor((xl(1) - d.t0) / d.dt) + 1);
i2 = min(d.n, ceil( (xl(2) - d.t0) / d.dt) + 1);
if i2 - i1 < 2
    return
end

idx = i1:i2;
tv  = d.t0 + (idx - 1) * d.dt;
if d.isheat
    [td, cd] = binaverage(tv, d.y(idx), d.maxpoints);
    set(d.h, 'XData', td, 'CData', cd);
else
    [td, yd] = minmaxdecimate(tv, d.y(idx), d.maxpoints);
    set(d.h, 'XData', td, 'YData', yd);
end
end

% -------------------------------------------------------------------------
function [td, yd] = minmaxdecimate(t, y, maxpoints)
% Decimate for display by keeping the min and max of each bin, in time order.
% Plain striding would alias hours of oscillation into whatever the stride
% happened to land on; this keeps the envelope truthful.

n = numel(y);
if n <= maxpoints
    td = t; yd = y;
    return
end
nbin = max(1, floor(maxpoints / 2));
per  = floor(n / nbin);
n2   = nbin * per;
Y    = reshape(y(1:n2), per, nbin);
T    = reshape(t(1:n2), per, nbin);

[mn, imn] = min(Y, [], 1);
[mx, imx] = max(Y, [], 1);
tmn = T(sub2ind(size(T), imn, 1:nbin));
tmx = T(sub2ind(size(T), imx, 1:nbin));

first = tmn <= tmx;                       % whichever extreme comes first in time
td = zeros(1, 2 * nbin);
yd = zeros(1, 2 * nbin);
td(1:2:end) = min(tmn, tmx);  td(2:2:end) = max(tmn, tmx);
yd(1:2:end) = mn .* first + mx .* ~first;
yd(2:2:end) = mx .* first + mn .* ~first;
end

% -------------------------------------------------------------------------
function [td, yd] = binaverage(t, y, nbin)
% Bin-average a trace down to nbin points, for the heatmap view.

n = numel(y);
if n <= nbin
    td = t; yd = y;
    return
end
per = floor(n / nbin);
n2  = nbin * per;
yd  = mean(reshape(y(1:n2), per, nbin), 1);
td  = mean(reshape(t(1:n2), per, nbin), 1);
end

% -------------------------------------------------------------------------
function v = prctileish(x, p)
% Percentile without the Statistics toolbox - keeps a few huge samples from
% washing out the heatmap's colour scale.

x = sort(x(:));
if isempty(x), v = 1; return, end
v = x(min(numel(x), max(1, round(p / 100 * numel(x)))));
if v <= 0, v = max(x(end), eps); end
end
