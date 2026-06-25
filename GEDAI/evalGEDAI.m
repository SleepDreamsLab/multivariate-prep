function evalGEDAI(BIDS, savenames, opts)
arguments
    BIDS
    savenames cell

    %--- Non-standard BIDS paths ---
    opts.scoringpath      char = fullfile(BIDS.pth, 'derivatives\scoring\scores\Manual_Checked')    

    %--- Subject filter ---
    opts.subjectfilter      cell = {}

    %--- EEG ---
    opts.tasklabel               = {'Sleep', 'sleep'}
    opts.recordinglabel     char = '125Hz'
    opts.net                char = 'EGI256'   

    %--- GEDAI ---
    opts.runmode {mustBeMember(opts.runmode, {'WholeNight', 'StageSpecific', 'StateWise'})} = 'StageSpecific'    
    opts.epochlength (1,1) double = 30    
    opts.savepath           char = fullfile(BIDS.pth, 'derivatives', 'GEDAI')
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

    %%% Find matching scoring file (match on sub+ses only; recording label may differ)
    scoringFile = gedai.matchScoringFile(p.entities, scoringfiles);
    if isempty(scoringFile); error('GEDAI_BIDS:noScoring', 'No scoring file matched for %s.', fileID); end

    %%% Load sleep scoring
    fprintf('\nReading %s ...\n', scoringFile)
    scoringDigits = scoreloader(scoringFile); 

    %%% Epochs to plot
    epochsToPlot = gedai.resolveEpochsToPlot(opts.epochstoplot, scoringDigits);

    %%% Replace isolated N1 epochs at stage boundaries with neighbour stage
    scoringDigits_NoN1 = gedai.killN1(scoringDigits);


    %%% GEDAI runs
%     for iRun = 1:numel(opts.runs)
    for savename = savenames
        savename = savename{1};
%         r        = opts.runs{iRun}; % GEDAI parameters
%         savename = gedai.buildSaveName(r, EEG.srate);
%         fprintf('Run %d/%d: %s\n', iRun, numel(opts.runs), savename)
    
        % Define GEDAI logic
%         clear EEGgedai
%         savename = [opts.prefix opts.runmode '_' savename];
    
        % load GEDAI run
        gedaifile = fullfile(opts.savepath, 'EEG', savename, [fileID '.mat']);
        load(gedaifile)

        % temporary stuff
        EEG             = pop_interp(EEG, EEG.urchanlocs, 'spherical');
        EEGclean        = pop_interp(EEGclean, EEG.urchanlocs, 'spherical');
        EEG.chanlocs    = readlocs('C:\Postdoc\Code\exploratory-prep\locfiles\electrodes.tsv')


        %%% Evaluation plot
        run.eval_clean(EEG, EEGgedai, scoringDigits_NoN1(ndxepochs), ...
            'EpochLength', opts.epochlength, 'WelchWindow', 4, ...
            'EpochsToPlot', epochsToPlot, 'refresh', opts.refresh, ...
            'SavePath', fullfile(opts.savepath, 'Figures', ['EGI_' savename], fileID, fileID))
        close all;        
    
    end
end