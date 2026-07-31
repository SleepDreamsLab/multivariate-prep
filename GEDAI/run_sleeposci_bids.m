function failures = run_sleeposci_bids(BIDS, opts)
% RUN_GEDAI_BIDS  Run GEDAI artefact-rejection on pre-filtered BIDS EEG data.
%
%   run_gedai_bids(BIDS)
%   run_gedai_bids(BIDS, Name, Value, ...)
%
%   Required
%   --------
%   BIDS              bids.layout object.
%
%   Input paths
%   -----------
%   filteredpath      Root of the filtered derivatives, organised as
%                     <sub>/<ses>/<fileID>_desc-<filtdesc>_eeg.dat.
%                     Default: <BIDS root>/derivatives/prep-ged/A_filtered
%   filtdesc          desc label used when building the filtered filename.
%                     Default: 'filt'
%   scoringpath       Directory containing sleep-scoring files (.json or .csv).
%                     Default: <BIDS root>/derivatives/scoring/scores/Manual_Checked
%   sfppath           Path passed to the SFP resolver.
%                     Default: <BIDS root>
%   leadfielddir      Root for Brainstorm leadfields (<sub>/<ses>/headmodel_surf_openmeeg.mat).
%                     Default: <BIDS root>/../Data_Analysis/Brainstorm_db/Leadfield_PM/data
%
%   Output paths
%   ------------
%   savepath          Root for GEDAI outputs.
%                     Default: <BIDS root>/derivatives/prep-ged/GEDAI
%   figpath           Root for all figures.
%                     Default: <BIDS root>/derivatives/prep-ged/figures
%   geddesc           desc label for GEDAI BrainVision output files.
%                     Default: 'filtGEDAI'
%   refresh           Force re-run even if a cache file exists. Default: false.
%
%   EEG
%   ---
%   tasklabel         BIDS task label(s) to query. Default: {'Sleep','sleep'}.
%   acqlabel    BIDS recording label to query. Default: '125Hz'.
%   noteegchannels    Channel indices to drop. Default: 257:264.
%   net               EEG net identifier passed to chans1020. Default: 'EGI256'.
%
%   Subject filter
%   --------------
%   subjectfilter     Cell array of subject ID strings; {} = all subjects.
%
%   GEDAI
%   -----
%   runmode           'StageSpecific', 'StateWise', or 'WholeNight'. Default: 'StageSpecific'.
%   epochlength       Sleep-epoch length in seconds. Default: 30.
%   runs              Cell array of GEDAI run-config structs. Default: gedai.defaultRuns().
%   epochstoplot      Epoch indices for diagnostic figures. Default: auto.
%   prefix            Prefix prepended to the run savename. Default: ''.
%
%   See dependancies.m.

arguments
    BIDS

    %--- Derivative folder ---
    opts.derivfolder      char = 'prep-ged'

    %--- Input paths ---
    opts.inputpath        char = ''
    opts.inputdesc        char = 'filt2ged'
    opts.inputfileext     char = '.set'    
    opts.scoringpath      char = fullfile(BIDS.pth, 'derivatives', 'scoring', 'scores', 'Manual_Checked')
    opts.sfppath          char = BIDS.pth
    opts.leadfielddir     char = fullfile(BIDS.pth, '..', 'Data_Analysis', 'Brainstorm_db', 'Leadfield_PM', 'data')

    %--- Output paths ---
    opts.savepath         char = ''
    opts.figpath          char = ''
    opts.geddesc          char = 'filt2ged2sleeposci'
    opts.refresh (1,1) logical = false
    opts.savefileext      char = '.set'

    %--- EEG ---
    opts.tasklabel                      = {'Sleep', 'sleep'}
    opts.acqlabel    char               = '125Hz'
    opts.noteegchannels    (1,:) double = 257:300
    opts.net               char         = 'EGI256'

    %--- Subject filter ---
    opts.subjectfilter     cell          = {}

    %--- GEDAI ---
    opts.runmode {mustBeMember(opts.runmode, {'WholeNight', 'StageSpecific', 'StateWise'})} = 'StageSpecific'
    opts.epochlength (1,1) double = 30
    opts.runs                     = []
    opts.epochstoplot             = []
    opts.prefix            char   = ''
end

if isempty(opts.inputpath), opts.inputpath = fullfile(BIDS.pth, 'derivatives', opts.derivfolder); end
if isempty(opts.savepath),  opts.savepath  = fullfile(BIDS.pth, 'derivatives', opts.derivfolder); end
if isempty(opts.figpath),   opts.figpath   = fullfile(BIDS.pth, 'derivatives', opts.derivfolder, 'figures'); end

KeepTime = struct();
if isempty(opts.runs), opts.runs = gedai.defaultRuns();
end

%%% Query EEG files from BIDS (used for entity extraction and subject iteration)
filesEEG = bids.query(BIDS, 'data', 'extension', '.vhdr', ...
    'task', opts.tasklabel, 'acq', opts.acqlabel);
if isempty(filesEEG); error('run_gedai_bids:noFiles', 'No matching EEG files found in BIDS layout.');
end

%%% Scoring files
if ~isempty(opts.scoringpath); scoringfiles = gedai.collectScoringFiles(opts.scoringpath);
end

%%% Loop over EEG files
failures = {};
for ifile = 1:numel(filesEEG)
    eegFile = filesEEG{ifile};
    p       = bids.internal.parse_filename(eegFile);
    fileID  = strjoin(cellfun(@(k) [k '-' p.entities.(k)], fieldnames(p.entities), 'uni', 0), '_');
    subDir  = fullfile(['sub-' p.entities.sub], ['ses-' p.entities.ses]);

    %%% Subject filter
    if ~isempty(opts.subjectfilter) && ~contains(fileID, opts.subjectfilter)
        continue
    end
    fprintf('\n=== %s ===\n', fileID)
    
    %%% Try block
%     try

    %%% Resolve filtered input file
    filtFile = fullfile(opts.inputpath, subDir, [fileID '_desc-' opts.inputdesc '_eeg' opts.inputfileext]);
    if ~isfile(filtFile)
        fprintf('[skip] %s: filtered file not found (%s)\n', fileID, filtFile)
        continue
    end
    fprintf('Input  → %s\n', filtFile)
    fprintf('Output → %s\n', opts.savepath)

    %%% Find matching scoring file
    if isempty(opts.scoringpath)
        scoringpath = fullfile(BIDS.pth, subDir);
        scoringfiles = gedai.collectScoringFiles(scoringpath);
    end      
    scoringFile = gedai.matchScoringFile(p.entities, scoringfiles);
    fprintf('Scoring → %s\n', scoringFile)
    if isempty(scoringFile)
        error('run_gedai_bids:noScoring', 'No scoring file matched for %s.', fileID);
    end

    %%% Load sleep scoring
    fprintf('\nReading %s ...\n', scoringFile)
    scoringDigits = scoreloader(scoringFile);    

    %%% Import filtered EEG
    D = tic; fprintf('\nImporting filtered EEG ...\n')
    EEG = eeg_import(filtFile);
    KeepTime.EEGimport = toc(D);
end