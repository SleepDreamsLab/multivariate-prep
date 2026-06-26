function GEDAI_BIDS(BIDS, opts)
%GEDAI_BIDS  Run GEDAI artefact-rejection on a BIDS sleep EEG dataset.
%
%   GEDAI_BIDS(BIDS)
%   GEDAI_BIDS(BIDS, Name, Value, ...)
%
%   Required
%   --------
%   BIDS              bids.layout object.
%
%   Non-standard BIDS paths  (all optional)
%   ----------------------------------------
%   scoringpath       Directory containing sleep-scoring files (.json or .csv).
%                     Default: <BIDS root>/derivatives/scoring/scores/Manual_Checked
%   sfppath           BIDS root path, passed to the study-specific SFP
%                     resolver (dispatched via fileID on "PM" or "DROP").
%                     Default: <BIDS root>
%   leadfielddir      Root directory for brainstorm leadfields, organised as
%                     <sub>/<ses>/headmodel_surf_openmeeg.mat.
%                     Default: <BIDS root>/../Data_Analysis/Brainstorm_db/Leadfield_PM/data
%
%   Saving
%   ------
%   savepath          Output directory. Default: <BIDS root>/derivatives/GEDAI
%   refresh           Force re-run even if a cache file exists. Default: false.
%
%   EEG
%   ---
%   tasklabel         BIDS task label(s) to query. Default: {'Sleep','sleep'}.
%   recordinglabel    BIDS recording label to query. Default: '125Hz'.
%   noteegchannels    Channel indices to drop before processing.
%                     Default: [257:264].
%   targetsrate       Resample EEG to this rate (Hz). 0 = no resampling.
%                     Default: 125.
%   net               EEG net identifier passed to chans1020. Default: 'EGI256'.
%
%   Subject filter
%   --------------
%   subjectfilter     Cell array of subject ID strings; process only files
%                     whose fileID contains at least one entry.
%                     Default: {'hpmam006'}.
%
%   GEDAI
%   -----
%   runmode           'StageSpecific' (separate GEDAI per stage) or
%                     'WholeNight' (single run over all stages).
%                     Default: 'StageSpecific'.
%   epochlength       Sleep-epoch length in seconds. Default: 30.
%   runs              Cell array of GEDAI run-config structs. Default: one
%                     run with ICAtype='none' and canonical settings.
%   epochstoplot      Epoch indices included in diagnostic figures.
%                     Default: auto-selected (first + middle of each stage).
%
%   See dependancies.m.

arguments
    BIDS   % bids.layout object

    %--- Non-standard BIDS paths ---
    opts.scoringpath      char = fullfile(BIDS.pth, 'derivatives\scoring\scores\Manual_Checked')
    opts.sfppath          char = BIDS.pth
    opts.leadfielddir     char = fullfile(BIDS.pth, '..', 'Data_Analysis\Brainstorm_db\Leadfield_PM\data')

    %--- Saving ---
    opts.savepath         char = fullfile(BIDS.pth, 'derivatives', 'GEDAI')
    opts.refresh (1,1) logical = false

    %--- EEG ---
    opts.tasklabel               = {'Sleep', 'sleep'}
    opts.recordinglabel     char = '125Hz'
    opts.noteegchannels          = [257:264]
    opts.targetsrate (1,1) double = 125
    opts.net                char = 'EGI256'
    opts.removeDC  (1,1) logical = true;
    opts.removeLN  (1,1) logical = true;

    %--- Subject filter ---
    opts.subjectfilter      cell = {}

    %--- GEDAI ---
    opts.runmode {mustBeMember(opts.runmode, {'WholeNight', 'StageSpecific', 'StateWise'})} = 'StageSpecific'
    opts.epochlength (1,1) double = 30
    opts.runs                    = []
    opts.epochstoplot            = []
    opts.prefix             char = '';
end


%%% Initiate variables
KeepTime = [];
if isempty(opts.runs), opts.runs = gedai.defaultRuns(); end

%%% Query EEG files
filesEEG = bids.query(BIDS, 'data', 'extension', '.vhdr', ...
    'task', opts.tasklabel, 'recording', opts.recordinglabel);
if isempty(filesEEG); error('GEDAI_BIDS:noFiles', 'No matching EEG files found in BIDS layout.'); end

%%% Scoring files
scoringfiles = gedai.collectScoringFiles(opts.scoringpath);

%%% Loop over EEG files
for ifile = 1:numel(filesEEG)
    eegFile  = filesEEG{ifile};
    p        = bids.internal.parse_filename(eegFile);
    fileID   = strjoin(cellfun(@(k) [k '-' p.entities.(k)], fieldnames(p.entities), 'uni', 0), '_');

    %%% Subject filter
    if ~isempty(opts.subjectfilter) && ~contains(fileID, opts.subjectfilter)
        continue
    end
    fprintf('\n=== %s ===\n', fileID)
    fprintf('Save path → %s\n', opts.savepath)

    %%% Import EEG
    EEG = eeg_import(eegFile);

    %%% Drop non-EEG channels
    EEG = pop_select(EEG, 'nochannel', intersect(1:EEG.nbchan, opts.noteegchannels));

    %%% Optional downsampling
    if opts.targetsrate > 0 && EEG.srate ~= opts.targetsrate
        fprintf('Resampling %d → %d Hz ...\n', EEG.srate, opts.targetsrate)
        EEG = pop_resample(EEG, opts.targetsrate);
    end

    %%% Build filters
    fprintf('Building filters (srate = %d Hz) ...\n', EEG.srate)
    EEG_DCFilter_NumDen = filterbank(EEG.srate, 'DC_RCSquareFilt');

    %%% DC removal
    if opts.removeDC
        D = tic; fprintf('\nDC removal ...\n')
        EEG.data = filtfilt(EEG_DCFilter_NumDen(1,:), EEG_DCFilter_NumDen(2,:), double(EEG.data'))';
        KeepTime.DCRemoval = toc(D);
    end

    %%% Zapline
    if opts.removeLN
        D = tic; fprintf('\nZapline plus ...\n')
        [EEG.data, plotHandles, analyticsResults] = clean_data_with_zapline_plus( ...
            EEG.data, EEG.srate, ...
            'noisefreqs', 'line', ...   % auto-detect line freq, or specify [50]
            'plotResults', 0);
        KeepTime.Zapline = toc(D);
        fprintf('ZapLine-plus: %.2f min\n', KeepTime.Zapline/60);
    end
end