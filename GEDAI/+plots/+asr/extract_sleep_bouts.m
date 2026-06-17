function [X] = extract_sleep_bouts(EEG, Stage, opts)
arguments (Input)
    EEG
    Stage
    opts.ScoringDuration = 30
    opts.BoutDuration = []
    opts.ThisBout = 1
end

% Extract EEG of sleep stage
EEG_SleepStage  = pop_epoch(EEG, Stage, [0 opts.ScoringDuration]);
EEG_SleepStage  = eeg_epoch2continuous(EEG_SleepStage); 

% Bout Length
if isempty(opts.BoutDuration)
    BoutLength = EEG_SleepStage.pnts ;
else
    BoutLength = opts.BoutDuration * 60 * EEG.srate;    
end
BoutStarts = 1 : BoutLength : (EEG_SleepStage.pnts);

    
% Extract the current bout
idxwindow = BoutStarts(opts.ThisBout) : BoutStarts(opts.ThisBout) + BoutLength - 1;

% Extract respective data
X = EEG_SleepStage;
X.data = EEG_SleepStage.data(:, idxwindow);
X = eeg_checkset(X);
