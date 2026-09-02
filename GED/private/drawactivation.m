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
%   'line' is accepted as a deprecated alias for 'signal'.

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

if strcmp(acttype, 'line'), acttype = 'signal'; end     % pre-rename callers

%%% Both non-raw modes show the same quantity - the smoothed absolute amplitude.
%%% Across long recordings the raw trace collapses into a solid band, while its
%%% envelope shows when the component actually switches on.
isheat = strcmp(acttype, 'heat');
switch acttype
    case 'heat'
        env      = envelope_(x, srate, smoothsec);
        [td, yd] = binaverage(tplot, env, maxpoints);
        imagesc(ax, td, [0 1], yd);
        set(ax, 'YTick', []);
        clim(ax, [0 prctileish(yd, 99)]);

    case 'envelope'
        env      = envelope_(x, srate, smoothsec);
        [td, yd] = minmaxdecimate(tplot, env, maxpoints);
        plot(ax, td, yd, 'Color', color, 'LineWidth', 1.1);
        axis(ax, 'tight');
        %%% An amplitude is measured from zero, so the axis starts there: with a
        %%% tight lower limit, ordinary ripple in a flat envelope would look like
        %%% a component switching on and off.
        ylim(ax, [0 max(max(yd) * 1.05, eps)]);

    otherwise   % 'signal'
        [td, yd] = minmaxdecimate(tplot, x, maxpoints);
        plot(ax, td, yd, 'Color', color, 'LineWidth', 0.7);
        axis(ax, 'tight');
end
box(ax, 'off');
end

% -------------------------------------------------------------------------
function env = envelope_(x, srate, smoothsec)
% Smoothed absolute amplitude. A moving average of |x| rather than a Hilbert
% envelope, so no Signal Processing Toolbox is needed.

env = movmean(abs(x), max(1, round(smoothsec * srate)));
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
