function drawhypnogram(ax, scoring, epochdur, keptepochs, lw)
% DRAWHYPNOGRAM  Standard hypnogram into one axes.
%
%   Wake on top, then REM, N1, N2, N3 going down. REM is redrawn over the trace
%   in its own colour, the way sleep labs read it. Shared by plotged and
%   plotgednight.
%
%   ax          target axes
%   scoring     one sleep-stage digit per epoch: -3 N3, -2 N2, -1 N1, 0 Wake,
%               1 REM (as scoreloader returns them)
%   epochdur    epoch length on the x axis, in whatever units the axis uses
%   keptepochs  epoch indices that fed the GED, shaded as a background band;
%               [] to draw none
%   lw          line width for the stage trace

if nargin < 4, keptepochs = []; end
if nargin < 5 || isempty(lw), lw = 1; end

digits = [ 1  0 -1 -2 -3];          % W  REM  N1  N2  N3
ypos   = [ 5  4  3  2  1];
names  = {'N3', 'N2', 'N1', 'REM', 'W'};

scoring = scoring(:)';
y = nan(size(scoring));
for k = 1:numel(digits)
    y(scoring == digits(k)) = ypos(k);
end

x = (0:numel(scoring)) * epochdur;              % epoch edges
hold(ax, 'on');

%%% Which epochs fed the GED. Drawn first, as a background band, and run-length
%%% encoded so a night of alternating stages does not become 900 patches.
if ~isempty(keptepochs)
    mask = false(1, numel(scoring));
    mask(keptepochs(keptepochs >= 1 & keptepochs <= numel(scoring))) = true;
    for run = findruns(mask)
        patch(ax, x([run(1) run(2) + 1 run(2) + 1 run(1)]), [0.4 0.4 5.6 5.6], ...
            [0.85 0.90 0.97], 'EdgeColor', 'none');
    end
end

stairs(ax, x, [y y(end)], 'Color', [0.20 0.20 0.20], 'LineWidth', lw);

%%% REM on top of the trace, in red.
for run = findruns(scoring == 1)
    idx = run(1):run(2);
    stairs(ax, x([idx idx(end) + 1]), [y(idx) y(idx(end))], ...
        'Color', [0.84 0.19 0.15], 'LineWidth', lw + 0.6);
end

set(ax, 'YTick', ypos(end:-1:1), 'YTickLabel', names, 'YLim', [0.4 5.6]);
ylabel(ax, 'Stage'); box(ax, 'off');
end
