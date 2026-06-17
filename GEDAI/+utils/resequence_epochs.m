function EEG = resequence_epochs(EEG, ndxepochs)
% GEDAI.UTILS.RESEQUENCE_EPOCHS  Reorder 30-s epochs in a continuous EEG structure.
%
%   EEG = gedai.utils.resequence_epochs(EEG, ndxepochs)
%
%   Inputs
%   ------
%   EEG       : Continuous EEGLAB EEG structure (EEG.trials == 1).
%               The data are expected to consist of back-to-back 30-second
%               epochs (any trailing samples that do not fill a full epoch
%               are left untouched).
%   ndxepochs : Integer permutation vector of length nEpochs.
%               Epoch i of the output corresponds to epoch ndxepochs(i) of
%               the input.  Obtain this via [~, ndxepochs] = sort(epochIdx)
%               to restore chronological order after stage-grouped processing.
%
%   Output
%   ------
%   EEG : Same structure with EEG.data reordered epoch-wise.
%         Boundary events are removed because their latency values are no
%         longer valid after reordering; use epochIdx to map epochs back to
%         the original Scoring vector.

nEpochs      = numel(ndxepochs);
epochSamples = floor(EEG.pnts / nEpochs);
nSamples     = epochSamples * nEpochs;   % samples that span full epochs

dataEpoched             = reshape(EEG.data(:, 1:nSamples), [EEG.nbchan, epochSamples, nEpochs]);
dataEpoched             = dataEpoched(:, :, ndxepochs);
EEG.data(:, 1:nSamples) = reshape(dataEpoched, [EEG.nbchan, nSamples]);

% Remove boundary events – latencies point into the old ordering
if ~isempty(EEG.event)
    EEG.event(strcmpi({EEG.event.type}, 'boundary')) = [];
    EEG = eeg_checkset(EEG);
end

end
