function plotBadChannelTime(corr, removed, savefile, opts)
% PLOTBADCHANNELTIME  Reconstruction correlation over the recording, per channel.
%
%   gedai.plotBadChannelTime(corr, removed, savefile)
%
%   corr and removed are the per-window correlations and the bad channel mask returned
%   by clean_channels (EEG.etc.badchans.corr / .mask). Every channel is drawn, ordered
%   top-to-bottom by proportion of low-correlation windows (opts.corrThreshold) - the same
%   quantity the correlation criterion itself decides on - so the worst channels sit at the
%   top regardless of whether they crossed maxBrokenTime and were actually removed. Flagged
%   (removed) channels are marked with a "*" in their row label.
%
%   The topoplot in gedai.plotBadChannels answers "which channels", this answers "when
%   and how badly": whether a channel is broken all night or only in bouts, and how far
%   past threshold it sits. A contact that dries out shows up as a channel that is clean
%   early and degrades monotonically; one that settles in does the reverse.

arguments
    corr      double
    removed   logical
    savefile  char
    opts.windowseconds (1,1) double = 5    % clean_channels default window length
    opts.corrThreshold (1,1) double = 0.7  % clean_channels default correlation threshold
    opts.nbins         (1,1) double = 300
    opts.clim          (1,2) double = [0.3 1]
    opts.title         char         = ''
end

[nChan, W]  = size(corr);
lowcorrprop = sum(corr < opts.corrThreshold, 2) / W;
[~, order]  = sort(lowcorrprop, 'descend');

%%% Average the windows into bins so a full night stays legible
nb    = max(1, min(opts.nbins, W));
edges = round(linspace(0, W, nb+1));
M     = zeros(nChan, nb);
for iBin = 1:nb
    M(:, iBin) = mean(corr(order, edges(iBin)+1:edges(iBin+1)), 2);
end
binHours = (edges(1:end-1) + edges(2:end)) / 2 * opts.windowseconds / 3600;

figure('Color', 'w', 'Position', [100 100 1100 max(720, 8*nChan)]);
imagesc(binHours, 1:nChan, M); clim(opts.clim)
% colormap(flipud(gray));
colormap(gray);
cb = colorbar; cb.Label.String = sprintf('RANSAC correlation (%g-s windows)', opts.windowseconds);
labels = compose('E%d', order(:));
labels(removed(order)) = strcat(labels(removed(order)), ' *');
set(gca, 'YTick', 1:nChan, 'FontSize', 8, 'YTickLabel', cellstr(labels))
xlabel('hours into recording')
title({opts.title, sprintf(['%d/%d channels flagged (*); ordered by proportion of low-' ...
    'correlation windows, worst at top; darker = worse'], nnz(removed), nChan)}, ...
    'Interpreter', 'none')

gedai.printFigure(gcf, savefile);
close
end
