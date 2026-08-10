function [fig, EEG] = amica_iclabel_bars(EEG, opts)
% ICA.PLOT.AMICA_ICLABEL_BARS  Vertical bar summary of ICLabel certainty per component.
%
%   ica.plot.amica_iclabel_bars(EEG)
%   ica.plot.amica_iclabel_bars(EEG, SavePath=fullfile(outDir, id))
%   [fig, EEG] = ica.plot.amica_iclabel_bars(EEG)   % reuse the classified EEG,
%                                                    % e.g. for ica.plot.amica_topographies
%
%   Runs eeg_checkset(EEG,'ica') + iclabel(EEG) — recomputing EEG.icaact (the
%   slow step, see ica.plot.amica_topographies for a quick unlabelled version)
%   — then plots one bar per component: bar height is the probability of its
%   top ICLabel class, bar colour/label is that class, x-tick is the
%   component index. Static axes only, so it prints/saves cleanly (unlike
%   pop_selectcomps's interactive topography grid). Skipped if EEG is already
%   classified (EEG.etc.ic_classification.ICLabel present), so it's safe to
%   call after another function has already run ICLabel.
%
%   EEG               – EEGLAB struct with icaweights/icasphere/icawinv/icachansind
%                       and EEG.data populated (EEG.data is required to
%                       recompute icaact for ICLabel).
%   Components        – component indices to include (default: all).
%   ArtefactThreshold – ICLabel probability threshold for pre-marking components
%                       as rejected in EEG.reject.gcompreject (default 0.5).
%                       That's what makes pop_selectcomps (via
%                       ica.plot.amica_topographies) show a component's button
%                       red — ICLabel alone only drives the topoplot title.
%   ICLabelClasses    – ICLabel classes to treat as artefacts (default [2 3 4 5 6]).
%                       1=Brain, 2=Muscle, 3=Eye, 4=Heart, 5=Line Noise, 6=Channel Noise, 7=Other
%   SavePath          – base path for saving; '_iclabel_bars.png' appended.
%                       '' (default) = don't save.

arguments
    EEG
    opts.Components        (1,:) double = []
    opts.ArtefactThreshold (1,1) double = 0.5
    opts.ICLabelClasses    (1,:) double = [2, 3, 4, 5, 6]
    opts.SavePath          char         = ''
end

if isempty(opts.Components)
    opts.Components = 1:size(EEG.icaweights, 1);
end

if ~isfield(EEG, 'etc') || ~isfield(EEG.etc, 'ic_classification') || ~isfield(EEG.etc.ic_classification, 'ICLabel')
    EEG = eeg_checkset(EEG, 'ica');
    EEG = iclabel(EEG);
end

%%% Pre-mark artefact components so pop_selectcomps (via amica_topographies)
%%% highlights them red, mirroring ica.selectcomps's pre-marking logic.
nAllComps = size(EEG.icaweights, 1);
[allProbs, allTopClass] = max(EEG.etc.ic_classification.ICLabel.classifications, [], 2);
artifactComps = ismember(allTopClass, opts.ICLabelClasses) & allProbs >= opts.ArtefactThreshold;
EEG.reject.gcompreject = artifactComps';
fprintf('ICLabel pre-marked %d/%d components for rejection\n', sum(artifactComps), nAllComps);

classNames = EEG.etc.ic_classification.ICLabel.classes;
classProbs = EEG.etc.ic_classification.ICLabel.classifications(opts.Components, :);
[certainty, topClass] = max(classProbs, [], 2);

classColors = [
    0.20 0.55 0.30;   % Brain
    0.85 0.55 0.10;   % Muscle
    0.75 0.20 0.65;   % Eye
    0.80 0.20 0.20;   % Heart
    0.30 0.55 0.85;   % Line Noise
    0.55 0.40 0.25;   % Channel Noise
    0.55 0.55 0.55];  % Other

nComps = numel(opts.Components);
fig = figure('Color', 'w', 'Position', [100 100 max(500, 18*nComps + 150) 500]);
ax = axes(fig); hold(ax, 'on'); box(ax, 'off');

bar(ax, 1:nComps, certainty, 0.7, 'FaceColor', 'flat', ...
    'CData', classColors(topClass, :), 'EdgeColor', 'none');

for i = 1:nComps
    text(ax, i, certainty(i) + 0.01, classNames{topClass(i)}, ...
        'FontSize', 7, 'Rotation', 90, 'HorizontalAlignment', 'left', 'Color', classColors(topClass(i), :));
end

set(ax, 'XTick', 1:nComps, 'XTickLabel', arrayfun(@(c) sprintf('IC%d', c), opts.Components, 'uni', 0), ...
    'FontSize', 8, 'TickDir', 'out');
xtickangle(ax, 90);
ylim(ax, [0 1]);
xlim(ax, [0.5 nComps + 0.5]);
ylabel(ax, 'ICLabel certainty (top class probability)');
title(ax, 'AMICA components — ICLabel classification', 'FontWeight', 'normal');
grid(ax, 'on'); ax.GridAlpha = 0.12; ax.XGrid = 'off';

legH = gobjects(1, numel(classNames));
for c = 1:numel(classNames)
    legH(c) = patch(ax, nan, nan, classColors(c,:), 'EdgeColor', 'none', 'DisplayName', classNames{c});
end
legend(ax, legH, 'Location', 'eastoutside', 'Box', 'off', 'FontSize', 8);

save_fig(fig, opts.SavePath, 'iclabel_bars');
end
