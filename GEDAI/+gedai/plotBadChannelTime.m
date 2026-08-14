function plotBadChannelTime(corr, removed, savefile, opts)
% PLOTBADCHANNELTIME  Reconstruction correlation over the recording, per channel.
%
%   gedai.plotBadChannelTime(corr, removed, savefile)
%
%   corr and removed are the per-window correlations and the bad channel mask returned
%   by clean_channels (EEG.etc.badchans.corr / .mask). Flagged channels are drawn on
%   top, a sample of retained channels below a red divider for reference.
%
%   The topoplot in gedai.plotBadChannels answers "which channels", this answers "when
%   and how badly": whether a channel is broken all night or only in bouts, and how far
%   past threshold it sits. A contact that dries out shows up as a channel that is clean
%   early and degrades monotonically; one that settles in does the reverse. Both are
%   removed for the whole recording, so it is worth seeing which you have.

arguments
    corr      double
    removed   logical
    savefile  char
    opts.windowseconds (1,1) double = 5    % clean_channels default window length
    opts.nrefchans     (1,1) double = 6
    opts.nbins         (1,1) double = 300
    opts.clim          (1,2) double = [0.3 1]
    opts.title         char         = ''
end

W    = size(corr, 2);
flag = find(removed(:))';
good = find(~removed(:))';
if isempty(good), warning('gedai:plotBadChannelTime:noGoodChans', 'No retained channels.'); return, end
good = good(round(linspace(1, numel(good), min(opts.nrefchans, numel(good)))));
rows = [flag(:); good(:)];

%%% Average the windows into bins so a full night stays legible
nb    = max(1, min(opts.nbins, W));
edges = round(linspace(0, W, nb+1));
M     = zeros(numel(rows), nb);
for iRow = 1:numel(rows)
    for iBin = 1:nb
        M(iRow,iBin) = mean(corr(rows(iRow), edges(iBin)+1:edges(iBin+1)));
    end
end
binHours = (edges(1:end-1) + edges(2:end)) / 2 * opts.windowseconds / 3600;

figure('Color', 'w', 'Position', [100 100 1100 620]);
imagesc(binHours, 1:numel(rows), M); clim(opts.clim)
colormap(flipud(gray));
cb = colorbar; cb.Label.String = sprintf('RANSAC correlation (%g-s windows)', opts.windowseconds);
set(gca, 'YTick', 1:numel(rows), 'FontSize', 8, 'YTickLabel', ...
    [arrayfun(@(c) sprintf('E%d *', c), flag(:), 'uni', 0); ...
     arrayfun(@(c) sprintf('E%d',   c), good(:), 'uni', 0)]);
yline(numel(flag) + 0.5, 'r-', 'LineWidth', 2);
xlabel('hours into recording')
title({opts.title, sprintf(['%d flagged (*), %d retained channels shown for reference; ' ...
    'brighter = worse'], numel(flag), numel(good))}, 'Interpreter', 'none')

print(gcf, savefile, '-dpng', '-r100');
close
end
