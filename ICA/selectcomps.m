function ALLICA = selectcomps(ALLICA, opts)
% SELECTCOMPS  Launch pop_selectcomps on an ALLICA struct.
%
%   ica.selectcomps(ALLICA)
%   ica.selectcomps(ALLICA, components=1:32)
%   ica.selectcomps(ALLICA, ArtefactThreshold=0.7, ICLabelClasses=[2 3 4])
%
%   ALLICA            – array of ICA structs from wrapper.ICA_whole_night.
%                       ALLICA(1) must have .data populated (call ica.reinsertData first).
%   components        – component indices to display (default: all).
%   ArtefactThreshold – ICLabel probability threshold for pre-marking components (default: 0.5).
%   ICLabelClasses    – ICLabel classes to treat as artefacts (default: [2 3 4 5 6]).
%                       1=Brain, 2=Muscle, 3=Eye, 4=Heart, 5=Line Noise, 6=Channel Noise, 7=Other
%
%   pop_selectcomps requires the variable to be named "EEG", so ALLICA(1)
%   is copied into a local variable with that name before the call.
%   Components flagged as artefacts by ICLabel are pre-marked in
%   EEG.reject.gcompreject so they appear highlighted in the GUI.

arguments
    ALLICA
    opts.components        {mustBeInteger, mustBePositive} = []
    opts.ArtefactThreshold (1,1) double                   = 0.5
    opts.ICLabelClasses    (1,:) double                   = [2, 3, 4, 5, 6]
    opts.ManualQC          (1,1) logical                  = true
end

%%% Copy first ICA struct into a variable named EEG (required by pop_selectcomps)
EEG = ALLICA(1);

%%% Default: show all components
if isempty(opts.components)
    opts.components = 1:size(EEG.icaweights, 1);
end

%%% Pre-mark artefact components from ICLabel (mirrors ICA_remove_whole_night)
nComps = size(EEG.icaweights, 1);
EEG.reject.gcompreject = false(1, nComps);

if isfield(EEG, 'etc') && isfield(EEG.etc, 'ic_classification')
    classifications        = EEG.etc.ic_classification.ICLabel.classifications;
    [topProbs, topClasses] = max(classifications, [], 2);
    artifactComps          = ismember(topClasses, opts.ICLabelClasses) & topProbs >= opts.ArtefactThreshold;
    EEG.reject.gcompreject = artifactComps';
    fprintf('ICLabel pre-marked %d/%d components for rejection\n', sum(artifactComps), nComps);
else
    warning('ica:selectcomps:noICLabel', 'No ICLabel classifications found in ALLICA(1).etc — no components pre-marked.');
end

if opts.ManualQC
    pop_selectcomps(EEG, opts.components);
    
    %%% Wait for user to finish selecting components
    fprintf('Select/deselect components in the figure, then press Space to continue...\n');
    pause;
    
    %%% Read updated rejection flags from the figure's UserData (set by GUI callbacks)
    fig = findobj(0, 'Tag', 'selcomp');
    if ~isempty(fig)
        EEG = get(fig(1), 'UserData');
    else
        warning('ica:selectcomps:figNotFound', 'Could not find selcomp figure — rejection labels not updated.');
    end
end

%%% Write updated rejection labels back into ALLICA
ALLICA(1).reject.gcompreject = EEG.reject.gcompreject;

end
