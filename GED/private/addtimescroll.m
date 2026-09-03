function addtimescroll(fig, tl, axTime, winlen, xfull)
% ADDTIMESCROLL  Show a fixed-length time window, with a scrollbar to move it.
%
%   Used by plotged and plotgednight when the caller asks for a window shorter
%   than the recording (xwindow). Without it a whole night is squeezed into one
%   axis, where individual spindles or slow waves are a few pixels wide; with it
%   the axes show winlen at a time and the slider pages through the recording.
%
%   The axes are already linked by the caller, so moving one moves all - the
%   slider therefore only has to set the limits of the first.
%
%   fig      figure to put the slider in
%   tl       the TiledChartLayout, shrunk to make room at the bottom
%   axTime   the linked time axes
%   winlen   window length, in the same units as the axes (not seconds)
%   xfull    [start stop] of the full recording, in those same units
%
%   Does nothing when the window covers the whole recording anyway.

if isempty(winlen) || ~isfinite(winlen) || winlen <= 0 || winlen >= diff(xfull)
    return
end

%%% Free a strip at the bottom of the figure. The layout owns the whole figure
%%% otherwise, and the slider would sit on top of the bottom axis labels.
tl.OuterPosition = [0 0.06 1 0.94];

vmax = max(xfull(2) - winlen, xfull(1) + eps);
sld  = uicontrol(fig, 'Style', 'slider', 'Units', 'normalized', ...
    'Position', [0.07 0.008 0.88 0.032], ...
    'Min', xfull(1), 'Max', vmax, 'Value', xfull(1), ...
    'TooltipString', 'Scroll through the recording');

%%% Arrow click nudges by a tenth of a window, trough click pages by a whole one.
span = vmax - xfull(1);
sld.SliderStep = [min(0.1 * winlen / span, 1) min(winlen / span, 1)];

setwindow(axTime, xfull(1), winlen);
sld.Callback = @(src, ~) setwindow(axTime, src.Value, winlen);

%%% ContinuousValueChange makes it scroll while the thumb is dragged, instead of
%%% only jumping when it is released.
addlistener(sld, 'ContinuousValueChange', @(src, ~) setwindow(axTime, src.Value, winlen));
end

% -------------------------------------------------------------------------
function setwindow(axs, v, winlen)
% Axes get deleted when the user closes the figure mid-drag, so check first.

axs = axs(isgraphics(axs));
if isempty(axs), return, end
xlim(axs(1), [v v + winlen]);      % linked: the rest follow
end
