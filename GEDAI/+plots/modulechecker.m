function [] = modulechecker(EEG, Artefacts, EpochLength)
arguments
    EEG
    Artefacts
    EpochLength = 30;
end

ArtefactSignal = EEG.data;
ArtefactSignal(~Artefacts) = nan;
sprep.plot.eeglab_scroll(EEG, ArtefactSignal, EpochLength)