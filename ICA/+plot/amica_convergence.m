function fig = amica_convergence(LL, opts)
% ICA.PLOT.AMICA_CONVERGENCE  AMICA log-likelihood convergence plot.
%
%   ica.plot.amica_convergence(LL)
%   ica.plot.amica_convergence(LL, Label=root, SavePath=fullfile(outDir, id))
%
%   LL       – mod.LL(:) from loadmodout15; trailing unused (zero) iterations
%              are trimmed automatically.
%   MinDLL   – mindll convergence threshold to mark on the improvement panel
%              (default 5e-7).
%   Label    – text (e.g. the AMICA output folder) annotated on the figure.
%   SavePath – base path for saving; '_convergence.png' appended.
%              '' (default) = don't save.

arguments
    LL (:,1) double
    opts.MinDLL   (1,1) double = 5e-7
    opts.Label    char         = ''
    opts.SavePath char         = ''
end

last = find(LL ~= 0, 1, 'last');
if isempty(last), error('ica:plot:amica_convergence:emptyLL', 'LL is empty or all zeros'); end
LL = LL(1:last);
if numel(LL) < 2, error('ica:plot:amica_convergence:tooShort', 'LL has %d usable iteration(s)', numel(LL)); end

it   = (1:numel(LL))';
dLL  = diff(LL);
col  = [0.20 0.36 0.60];
colA = [0.85 0.33 0.20];

fig = figure('Color', 'w', 'Position', [100 100 900 620]);
ax1 = axes(fig, 'Position', [0.09 0.56 0.86 0.34]);
ax2 = axes(fig, 'Position', [0.09 0.10 0.86 0.34]);

% --- trajectory ---
hold(ax1, 'on'); box(ax1, 'off'); grid(ax1, 'on');
plot(ax1, it, LL, 'LineWidth', 1.6, 'Color', col);
yline(ax1, LL(end), ':', sprintf('final %.4f', LL(end)), 'Color', [.45 .45 .45], ...
      'LabelHorizontalAlignment', 'left', 'FontSize', 9);
[~, iBest] = max(LL);
if iBest < numel(LL)                        % LL can dip after a late overshoot
    plot(ax1, iBest, LL(iBest), 'o', 'MarkerSize', 6, ...
         'MarkerFaceColor', colA, 'MarkerEdgeColor', 'none');
    text(ax1, iBest, LL(iBest), sprintf('  peak @ %d', iBest), 'FontSize', 9, 'Color', colA);
end
ylabel(ax1, 'log-likelihood');
ax1.GridAlpha = 0.12;
srt = sort(LL); lo = srt(max(1, ceil(0.02*numel(srt))));   % clip the early cliff
if max(LL) > lo, ylim(ax1, [lo max(LL) + 0.02*(max(LL)-lo)]); end
title(ax1, 'AMICA convergence', 'FontWeight', 'normal', 'FontSize', 12);

% --- per-iteration improvement ---
hold(ax2, 'on'); box(ax2, 'off'); grid(ax2, 'on');
pos = dLL > 0;
h = gobjects(0); lbl = {};
if any(pos)
    h(end+1) = semilogy(ax2, it(find(pos)+1), dLL(pos), '-', 'LineWidth', 1.1, 'Color', col);
    lbl{end+1} = 'increase';
end
if any(~pos)
    h(end+1) = semilogy(ax2, it(find(~pos)+1), abs(dLL(~pos)), '.', 'MarkerSize', 7, 'Color', colA);
    lbl{end+1} = 'decrease (overshoot)';
end
set(ax2, 'YScale', 'log');
yline(ax2, opts.MinDLL, '--', sprintf('mindll = %.10f', opts.MinDLL), 'Color', [.45 .45 .45], 'FontSize', 9);
if numel(h) > 1
    legend(ax2, h, lbl, 'Box', 'off', 'Location', 'northeast', 'FontSize', 9);
end
xlabel(ax2, 'iteration'); ylabel(ax2, '\DeltaLL per iteration');
ax2.GridAlpha = 0.12;

linkaxes([ax1 ax2], 'x');
xlim(ax1, [1 numel(LL)]);

% --- inset: last 10% (created last so it sits on top) ---
n0  = max(1, round(0.9*numel(LL)));
axI = axes(fig, 'Position', [0.60 0.62 0.30 0.15]);
hold(axI, 'on'); box(axI, 'on'); grid(axI, 'on');
plot(axI, it(n0:end), LL(n0:end), 'LineWidth', 1.2, 'Color', col);
axI.FontSize = 8; axI.GridAlpha = 0.12;
xlim(axI, [it(n0) it(end)]);
title(axI, 'last 10%', 'FontWeight', 'normal', 'FontSize', 9);

% --- verdict ---
tail = mean(dLL(max(1, end-99):end));
if tail > 1e-9, verdict = 'STILL CLIMBING (hit the iteration ceiling)';
else,           verdict = 'converged'; end
fprintf('\nAMICA LL: %d iters, final %.5f, peak %.5f @ %d\n', ...
        numel(LL), LL(end), max(LL), iBest);
fprintf('mean dLL over last 100 iters: %.3e  -->  %s\n', tail, verdict);

if ~isempty(opts.Label)
    annotation(fig, 'textbox', [0.01 0.005 0.3 0.02], 'String', opts.Label, ...
        'FontSize', 6, 'EdgeColor', 'none', 'Interpreter', 'none');
end

save_fig(fig, opts.SavePath, 'convergence');
end
