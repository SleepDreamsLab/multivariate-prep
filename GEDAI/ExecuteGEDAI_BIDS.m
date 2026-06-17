%%% Start clean
clearvars
clc; close all

%%% Variables
refresh = false; 
KeepTime = [];

run('dependancies.m')
run('Parameters/ParametersPM.m')

%%% Projects to run
[BIDS_PM] = qol.smartcache( ...
    @() qol.bidswizard({'PM'}, 'W:', 'Data_BIDS'), ...
    fullfile(pwd, 'TestFiles', 'BIDS_PM.mat'), false, {'BIDS_PM'});

%%% Loop over database entries
Database = organon.databasePM(BIDS_PM, 'recording', '125Hz', ...
    'relative_save_path', 'derivatives\GEDAI', ...
    'relative_scoring_path', 'derivatives\GEDAI\KCs');

%%% Add chanlocs
chanlocs = {fullfile(BIDS_PM{1}.pth, 'channels.tsv'), fullfile(BIDS_PM{1}.pth, 'electrodes.tsv')};
Database.Chanlocs = repmat(chanlocs, height(Database), 1);

%%% Epochs to plot
EpochsToPlot = [84, 1005, 42, 565, 321, 995];


%%% Loop over database entries
% parfor (ifile = 1:size(Database, 1), 3)
for ifile = find(contains(Database.FileID, 'sub-hpmam006_ses-S1'))
    fileID = Database.FileID{ifile};
    savepath = Database.Save{ifile};
    fprintf('Processing %s ...\n', fileID)
    fprintf('Saving files to %s ...\n', savepath)

    %%% Import sleep scoring
    scoringHero     = jsondecode(fileread(Database.Scoring{ifile}));
    scoringDigits   = [scoringHero{1}.digit];

    %%% Import EEG
    EEG = eeg_import(Database.EEG{ifile});

    %%% Import channel locations
    [EEG, channelData, elecData] = bids_importchanlocs(EEG, Database.Chanlocs{ifile, 1}, Database.Chanlocs{ifile, 2});
    EEG = qol.bids_fixchanlocs(EEG);

    %%% Select channels of interest
    EEG = pop_select(EEG, 'nochannel', Parameters.Channels.NotEEG);  
    
    %%% Build filters (original sampling rate)
    if ~exist("EEG_DCFilter_NumDen")
        fprintf('Building filters...\n')
        EEG_DCFilter_NumDen = filterbank(EEG.srate, 'DC_RCSquareFilt');
        EEG_HiPassFilt_IIR  = filterbank(EEG.srate, 'EEG_HiPassFilt_IIR'); % EEG_HiPassFilt_IIR, ICA_HiPassFilt_IIR
        EEG_NotchFilt_IIR   = filterbank(EEG.srate, 'EEG_NotchFilt_IIR2'); % 'EEG_NotchFilt_IIR2', 'EMG_NotchFilt_IIR'
    end

    %%% DC Removal
    D = tic;            
    fprintf('Removing DC ...\n')   
    EEG_DCFilter_NumDen = filterbank(EEG.srate, 'DC_RCSquareFilt');
    EEG.data = filtfilt(EEG_DCFilter_NumDen(1,:), EEG_DCFilter_NumDen(2,:), double(EEG.data'))'; 
    KeepTime.DCRemoval = toc(D);

%     %%% High-pass filter
%     D = tic;
%     fprintf('High-pass filtering ...\n')
%     EEG.data = filtfilt(EEG_HiPassFilt_IIR, EEG.data')';
%     KeepTime.HPFilter = toc(D);    

    %%% Remove line noise
    D = tic; EEG.data = double(EEG.data);    
    EEG_NotchFilt_IIR = filterbank(EEG.srate, 'EEG_NotchFilt_IIR2');
    for ifilt = 1:numel(EEG_NotchFilt_IIR)
        fprintf('%dHz notch filtering ...\n', 50*ifilt) 
        EEG.data = filtfilt(EEG_NotchFilt_IIR{ifilt}, EEG.data')';  % 50 Hz
    end
    KeepTime.NotchFilter = toc(D); 

    %%% Average re-reference
    EEG = pop_reref(EEG, []);   

    %%% Real channel locations
    coordinates = [extractBefore(Database.EEG{ifile}, '_task'), '_coordinates.sfp'];
    chanlocs = readlocs(coordinates);
    EEG.chanlocs = chanlocs(4:EEG.nbchan+3); 

    %%% Remove reference channel if present
    EEG = pop_select(EEG, 'nochannel', 257);    

    %%% Kill N1
    N1 = -1; scoringDigits_NoN1 = scoringDigits;
    while any(scoringDigits_NoN1 == N1)
        firstN1 = find([0 diff(scoringDigits_NoN1 == N1)] == 1);
        lastN1  = find([diff(scoringDigits_NoN1 == N1), 0] == -1);
    
        scoringDigits_NoN1(firstN1) = scoringDigits_NoN1(firstN1 - 1);
        scoringDigits_NoN1(lastN1)  = scoringDigits_NoN1(lastN1 + 1);
    end

    %%% Assign 10-20 labels
    EEG = gedai.extract_1020(EEG, false);


    %%% ------------------------------------------------------------------- %%%
    %%%                      Bad channel detection                          %%%

    % Bad channel detection
    chancorr_crit = 0.7;
    line_crit = Inf;
    channel_crit_maxbad_time = 0.5;
    num_samples = 25;
    D = tic;
    [EEG, removed_channels, corrs] = clean_channels(EEG, chancorr_crit, line_crit, [], channel_crit_maxbad_time, num_samples);
    KeepTime.BadChannelDetection = toc(D);

    % ### leadfield from brainstorm
    bstorm  = load('\\vs03.herseninstituut.knaw.nl\VS03-SandD-4\PM\Data_Analysis\Brainstorm_db\Leadfield_PM\data\sub-hpmam006\ses-S1\headmodel_surf_openmeeg.mat');
%     bstorm  = load('\\vs03.herseninstituut.knaw.nl\VS03-SandD-4\PM\Data_Analysis\Brainstorm_db\GEDAI_Leadfield\data\sub-hpmam006_ses-S1\shift_warp\headmodel_surf_openmeeg.mat');
    B = bstorm.Gain(setdiff(1:256, find(removed_channels)), :);            % [nchan x 3*nsources], unconstrained
    B = B - sum(B,1)/(size(B,1)+1);      % GEDAI's non-rank-deficient avg ref
    lfCOV = B*B';                        % drop in place of interp_mont_GEDAI's lfCOV
     
%     %%% Leadfield covariance matrix
%     L = load('fsavLEADFIELD_4_GEDAI.mat');
%     leadfield_EEG      = L.leadfield4GEDAI.EEG;
%     leadfield_EEG.data = L.leadfield4GEDAI.Gain - sum(L.leadfield4GEDAI.Gain, 1) / (size(L.leadfield4GEDAI.Gain, 1) + 1); 
%     interpolated_EEG   = interp_mont_GEDAI(leadfield_EEG, EEG.chanlocs);
%     lfCOV              = interpolated_EEG.data * interpolated_EEG.data';    
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%                         GEDAI + ICA                                 %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    %%% Epochs to plot
    EpochsToPlot = [84, 1005, 42, 565, 321, 995];
    for score = -3:1
        sepochs = find(scoringDigits==score);
        EpochsToPlot   = [EpochsToPlot, sepochs( round(sum(scoringDigits==score)/2) )];
        EpochsToPlot   = [EpochsToPlot, sepochs( 1 )];
        if score == -2
            EpochsToPlot   = [EpochsToPlot, sepochs( [5:5:20] )];    
        end    
    end    

    %%
    %%% RUNS
    refresh = false; clear EEGclean;
    runs = { ...
        struct('GEDAIMode','auto','GEDAILowCutOffFreq',0.1,'GEDAIEpochSize',12,'GEDAIBroadbandEpochSize',10,'boost1',1,'boost2',1,'broadbandOnly',false,'percentileThreshold',98, 'WeightKC', 0, 'BBMinThreshold', -2, 'computeSENSAI', false, 'ICAtype', 'none', 'GEDAIEnovaChannelThreshold', Inf); ...
        struct('GEDAIMode','auto','GEDAILowCutOffFreq',0.1,'GEDAIEpochSize',12,'GEDAIBroadbandEpochSize',10,'boost1',1,'boost2',1,'broadbandOnly',false,'percentileThreshold',98, 'WeightKC', 0, 'BBMinThreshold', -2, 'computeSENSAI', false, 'ICAtype', 'picard', 'GEDAIEnovaChannelThreshold', Inf); ...
%         struct('GEDAIMode','auto+','GEDAILowCutOffFreq',0.1,'GEDAIEpochSize',12,'GEDAIBroadbandEpochSize',10,'boost1',1,'boost2',1,'broadbandOnly',false,'percentileThreshold',98, 'WeightKC', 0, 'BBMinThreshold', -2, 'computeSENSAI', false, 'ICAtype', 'none', 'GEDAIEnovaChannelThreshold', Inf); ...
%         struct('GEDAIMode','auto-','GEDAILowCutOffFreq',0.1,'GEDAIEpochSize',12,'GEDAIBroadbandEpochSize',10,'boost1',1,'boost2',1,'broadbandOnly',false,'percentileThreshold',98, 'WeightKC', 0, 'BBMinThreshold', -2, 'computeSENSAI', false, 'ICAtype', 'none', 'GEDAIEnovaChannelThreshold', Inf); ...
    };
    
    for iRun = 1:numel(runs)
        r = runs{iRun};
    
        %%% Savename
        savename = sprintf('BadChan_BstormAUTOC_ICA%s_BBonly%d_LowCutOff%d_BBEpochSize%d_%s_PrcThresh%d_BBMinThresh%d_GEDAIEnovaChannelThreshold%d_SENSAI%d_WeightKC%dAllStages_b1x%d_b2x%d_Srate%d', ...
            r.ICAtype, r.broadbandOnly, r.GEDAILowCutOffFreq*10, r.GEDAIBroadbandEpochSize, r.GEDAIMode, r.percentileThreshold, r.BBMinThreshold, r.GEDAIEnovaChannelThreshold, r.computeSENSAI, r.WeightKC*100, ...
            round(r.boost1*10), round(r.boost2*10), EEG.srate);  
    
        %%% GEDAI
        clear EEGgedai EEGraw
        refresh = false; close all;
        [EEGgedai, ndxepochs, KeepTime] = qol.smartcache( ...
            @() gedai.GEDAI_PerStage(EEG, scoringDigits_NoN1, ...
                {[-2], [-3], [0], [1]}, KeepTime, ...
                'EpochLength', 30, ...
                'GEDAIMode', r.GEDAIMode, ...
                'GEDAIEpochSize', r.GEDAIEpochSize, ...
                'GEDAILowCutOffFreq', r.GEDAILowCutOffFreq, ...
                'BBEpochSize', r.GEDAIBroadbandEpochSize, ...
                'BroadbandOnly', r.broadbandOnly, ...
                'GEDAIEnovaChannelThreshold', r.GEDAIEnovaChannelThreshold, ...
                'PercentileThreshold', r.percentileThreshold, ...
                'BBMinThreshold', r.BBMinThreshold, ...
                'ComputeSENSAI', r.computeSENSAI, ...
                'ICAtype', r.ICAtype, ...
                'RefCOV', {lfCOV, lfCOV, lfCOV, lfCOV}), ...
            fullfile(savepath, 'EEG', [savename '_' fileID '.mat']), ...
            refresh, {'EEGgedai', '', 'ndxepochs', 'KeepTime'});
        fprintf('GEDAI took %.2f min\n', KeepTime.GEDAI/60)
    %     ndxepochsAdj = unique(ndxepochs(2:2:end)/2);

        %%% Performance plots
        gedai.eval_clean(EEG, EEGgedai, scoringDigits_NoN1(ndxepochs), ...
          'TopoBandLims', [1.2, 2.1; .7, 1.5; -.7, -.1; -.3, .5; .0, 0.8], ...                    
          'EpochLength', 30, 'WelchWindow', 4, 'EpochsToPlot', EpochsToPlot, 'refresh', refresh, ...
          'SavePath', fullfile(savepath, 'Figures', ['GEDAIonly_' savename], fileID, [fileID]))
        close all;    
    
        %%% Remove ICA
        if isfield(EEGgedai.etc, 'ic_classification')
            EEGgedai = ica.selectcomps(EEGgedai, 'ArtefactThreshold', 0.5, 'ManualQC', false);
            EEGgedai = pop_subcomp(EEGgedai, find(EEG.reject.gcompreject), 0);
    
            %%% Performance plots
            gedai.eval_clean(EEG, EEGgedai, scoringDigits_NoN1(ndxepochs), ...
              'TopoBandLims', [1.2, 2.1; .7, 1.5; -.7, -.1; -.3, .5; .0, 0.8], ...                   
              'EpochLength', 30, 'WelchWindow', 4, 'EpochsToPlot', EpochsToPlot, 'refresh', refresh, ...
              'SavePath', fullfile(savepath, 'Figures', savename, fileID, [fileID]))
            close all;
        end
    end
end    