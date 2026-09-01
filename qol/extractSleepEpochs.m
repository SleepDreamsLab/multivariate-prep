function [EEG, keepIdx] = extractSleepEpochs(EEG, scoringDigits, opts)
% EXTRACTSLEEPEPOCHS  Restrict a continuous recording to the requested sleep stages.
%
%   [EEG, keepIdx] = extractSleepEpochs(EEG, scoringDigits)
%   [EEG, keepIdx] = extractSleepEpochs(EEG, scoringDigits, 'flatten', false)
%
%   Cuts the recording into epochlength-second epochs, matches them against the
%   sleep scoring, keeps only the epochs whose stage is in keepTheseStages, and
%   (by default) flattens the survivors back into one continuous recording.
%
%   This is the epoching step of GEDAI_StageSpecific, but nothing here dilates
%   stage runs into neighbouring epochs (gedai.dilateStages) or reassigns N1 to
%   a neighbour (gedai.killN1): the scoring digits are used exactly as read.
%
%   scoringDigits   one stage digit per epoch, as returned by scoreloader:
%                   -3 N3, -2 N2, -1 N1, 0 Wake, 1 REM
%
%   Name-value, all optional:
%     epochlength      sleep-epoch length in seconds                       (30)
%     keepTheseStages  stage digits to keep. Default keeps everything: [-3 -2 -1 0 1]
%     flatten          concatenate the surviving epochs back into one continuous
%                      recording                                         (true)
%
%                      Leave it true for anything that needs a contiguous time
%                      series -- ged(), for one, segments the recording into its
%                      own covariance windows and cannot take a stack of disjoint
%                      epochs. Set it false to get the epoched EEG back (trials =
%                      numel(keepIdx)), e.g. for per-epoch statistics.
%
%   Outputs:
%     EEG        the stage-selected recording, continuous or epoched per flatten
%     keepIdx    indices, into the pre-selection epoching, of the epochs kept
%
%   Errors if the scoring and the recording disagree on epoch count, or if no
%   epoch matches keepTheseStages.
arguments
    EEG (1,1) struct
    scoringDigits double
    opts.epochlength (1,1) double {mustBePositive} = 30
    opts.keepTheseStages (1,:) double = [-3 -2 -1 0 1]
    opts.flatten (1,1) logical = true
end

%%% Correct scoring length if needed - the scoring can run a few epochs
%%% longer than the recording it was scored from.
nEpochs = floor(EEG.pnts / (opts.epochlength * EEG.srate));
while numel(scoringDigits) > nEpochs, scoringDigits(end) = []; end

%%% Epoch into fixed-length sleep epochs, same as GEDAI_StageSpecific: strip
%%% boundary/epoch events first, since eeg_regepochs silently drops epochs
%%% that overlap one, which would desync the epochs from scoringDigits.
if ~isempty(EEG.event)
    EEG.event(strcmpi({EEG.event.type}, 'boundary')) = [];
    EEG.event(contains({EEG.event.type}, 'Epoch')) = [];
    EEG = eeg_checkset(EEG);
end
EEG = eeg_regepochs(EEG, 'recurrence', opts.epochlength, ...
    'limits', [0 opts.epochlength], 'eventtype', sprintf('Epoch%ds', opts.epochlength));

if numel(scoringDigits) ~= EEG.trials
    error('extractSleepEpochs:scoringMismatch', ...
        'Scoring has %d epoch(s) but the recording has %d %ds-epoch(s).', ...
        numel(scoringDigits), EEG.trials, opts.epochlength);
end

%%% Keep only the requested stages - no dilation, no N1 reassignment: the
%%% scoring digits are used exactly as read.
keepIdx = find(ismember(scoringDigits, opts.keepTheseStages));
if isempty(keepIdx)
    error('extractSleepEpochs:noStageEpochs', ...
        'No epochs match keepTheseStages = [%s].', num2str(opts.keepTheseStages));
end
fprintf('Keeping %d/%d epochs (stages [%s])\n', ...
    numel(keepIdx), numel(scoringDigits), num2str(opts.keepTheseStages));
EEG = pop_select(EEG, 'trial', keepIdx);

if opts.flatten
    EEG = eeg_epoch2continuous(EEG);
end
end
