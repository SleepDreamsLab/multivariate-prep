%%% Start clean
clearvars
clc; close all

%%% Dependencies
run('..\dependancies.m')

%% ===================== CONFIG =====================
%--- Subjects/sessions to process — the whole pipeline below loops over every
%   subject x session x AMICA-run-variant combination.
subjects = {'sub-drop0001'};
sessions = {'ses-t1', 'ses-t2', 'ses-t3', 'ses-t4', 'ses-t5'};

%--- AMICA run variant(s) per subject/session: the "desc-..." suffix after
%   <sub>_<ses>_task-sleep_run-01_ in the AMICA output folder name.
runDescs = { ...
    'zc2gedWakeBBAutoPlusFSAutoPlus2amica-Rej1Nmodel1F32Dll5e7MaxIter600Stride4_ica', ...
    'zc2gedWakeBBAutoPlusFSAutoPlus2amica-Rej0Nmodel1F32Dll5e7MaxIter600Stride4_ica'};
runDescs = {'zc2gedWakeBBAutoPlusFSAutoPlus2amica-Rej0Nmodel1F32Dll5e7MaxIter600Stride4_ica'};

%--- "desc" of the pre-AMICA preprocessed EEG file used for the component-topography plot
eegDesc = 'zc2gedWakeBBAutoPlusFSAutoPlus';

%--- Base directories (subject/session are inserted per iteration below)
pamicaRoot = '\\vs03.herseninstituut.knaw.nl\VS03-SandC-1\data\nin\data-drop\rawdata\derivatives\pamica';
prepRoot   = '\\vs03.herseninstituut.knaw.nl\VS03-SandC-1\data\nin\data-drop\rawdata\derivatives\prep-ged';
scoreRoot  = '\\vs03.herseninstituut.knaw.nl\VS03-SandC-1\data\nin\data-drop\derivatives\scores\final';
figRoot    = fullfile(pamicaRoot, 'figures');

useICLabel = true;   % true = ICLabel certainty bar chart (slow: recomputes icaact); false = topographies only (fast)

for iSub = 1:numel(subjects)
    subj = subjects{iSub};

    for iSes = 1:numel(sessions)
        sess = sessions{iSes};
        fprintf('\n=== %s %s ===\n', subj, sess);

        fold1       = fullfile(pamicaRoot, subj, sess);
        foldOut     = fullfile(figRoot, subj, sess);
        scoringFile = fullfile(scoreRoot, sprintf('%s_%s_task-sleep_run-01_eeg.csv', subj, sess));
        eegFile     = fullfile(prepRoot, subj, sess, sprintf('%s_%s_task-sleep_run-01_desc-%s_eeg.set', subj, sess, eegDesc));

        scoringDigits = scoreloader(scoringFile);
        EEGraw        = fast_eeg_import(eegFile);

        for iRun = 1:numel(runDescs)
            close all
            fold2LLcur = sprintf('%s_%s_task-sleep_run-01_desc-%s', subj, sess, runDescs{iRun});
            rootLL     = fullfile(fold1, fold2LLcur);
            outDir     = fullfile(foldOut, fold2LLcur);

            %% ===================== CONVERGENCE + LOG-LIKELIHOOD-BY-STAGE =====================
            mod = loadmodout15(rootLL);

            plot.amica_convergence(mod.LL(:), 'Label', rootLL, ...
                'SavePath', fullfile(outDir, fold2LLcur));

            plot.amica_lt_by_stage(mod.Lht, mod.v, scoringDigits, 'Label', rootLL, ...
                'SavePath', fullfile(outDir, fold2LLcur));

            %% ===================== COMPONENT TOPOGRAPHIES =====================
            pam = load(rootLL);            % icasphere + chanlabels
            mod = loadmodout15(rootLL);    % W, S — byte-identical to Fortran

            P = double(pam.icasphere);                        % k x nchan, the rank projection
            S = double(mod.S(1:mod.num_pcs, 1:mod.data_dim));  % S is num_pcs x data_dim

            labels = string(pam.chanlabels(:));
            assert(size(P,2) == numel(labels), 'P columns vs chanlabels mismatch');
            assert(size(S,2) == size(P,1),     'amicaout rank vs P mismatch');

            [tf, idx] = ismember(labels, string({EEGraw.chanlocs.labels}'));
            assert(all(tf), 'channels missing from EEG: %s', strjoin(labels(~tf), ', '));

            % Each AMICA model has its own unmixing matrix (and so its own components) —
            % loop over them so every model gets its own topography/ICLabel figures.
            for h = 1:mod.num_models
                EEG = EEGraw;
                EEG.icaweights  = double(mod.W(:,:,h));    % amicaout lives in the P-space, so fold P into the sphere
                EEG.icasphere   = S * P;
                EEG.icawinv     = pinv(EEG.icaweights * EEG.icasphere);
                EEG.icachansind = idx(:)';

                saveBase = fold2LLcur;
                if mod.num_models > 1
                    saveBase = sprintf('%s_model%d', fold2LLcur, h);
                end

                if useICLabel
                    [~, EEG] = plot.amica_iclabel_bars(EEG, 'SavePath', fullfile(outDir, saveBase), 'ICLabelClasses', [2:7], 'ArtefactThreshold', .0);
                end
                plot.amica_topographies(EEG, 'SavePath', fullfile(outDir, saveBase));
            end
        end
    end
end
