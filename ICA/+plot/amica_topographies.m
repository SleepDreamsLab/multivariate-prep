function figs = amica_topographies(EEG, opts)
% ICA.PLOT.AMICA_TOPOGRAPHIES  Component topography grids (pop_selectcomps).
%
%   ica.plot.amica_topographies(EEG)
%   ica.plot.amica_topographies(EEG, SavePath=fullfile(outDir, id))
%
%   Just forwards EEG to pop_selectcomps and saves the resulting page(s).
%   If EEG is already ICLabel-classified (EEG.etc.ic_classification.ICLabel
%   present — e.g. EEG returned by ica.plot.amica_iclabel_bars), pop_selectcomps
%   titles each topoplot with its class + probability automatically. Passing
%   an unclassified EEG skips the EEG.icaact recompute and ICLabel run
%   entirely, for a quick unlabelled version.
%
%   EEG        – EEGLAB struct with icaweights/icasphere/icawinv/icachansind
%               and chanlocs populated.
%   Components – component indices to display (default: all).
%   SavePath   – base path for saving; '_topographies_pageN.png' appended
%               per figure. '' (default) = don't save.
%   Close      – close the figure(s) after saving (default true).
%
%   pop_selectcomps opens one figure per page (~35 components each).

arguments
    EEG
    opts.Components (1,:) double  = []
    opts.SavePath    char         = ''
    opts.Close       (1,1) logical = true
end

if isempty(opts.Components)
    opts.Components = 1:size(EEG.icaweights, 1);
end

if ~isfield(EEG, 'reject') || ~isfield(EEG.reject, 'gcompreject') || isempty(EEG.reject.gcompreject)
    EEG.reject.gcompreject = false(1, size(EEG.icawinv, 2));
end

% pop_selectcomps tags each figure 'selcomp<rand>' (not the literal
% 'selcomp'), so an exact-Tag findobj never matches. Diff all figures
% instead — a >35-component request opens several windows at once.
preExisting = findobj(0, 'Type', 'figure');
pop_selectcomps(EEG, opts.Components);
figs = setdiff(findobj(0, 'Type', 'figure'), preExisting);

hasLabels = isfield(EEG, 'etc') && isfield(EEG.etc, 'ic_classification') && isfield(EEG.etc.ic_classification, 'ICLabel');
suffix    = 'topographies_page%d';
if hasLabels, suffix = 'topographies_iclabel_page%d'; end

for iFig = 1:numel(figs)
    save_fig(figs(iFig), opts.SavePath, sprintf(suffix, iFig));
end

if opts.Close
    close(figs);
end
end
