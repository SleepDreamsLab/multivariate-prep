%% EEGCompareViewer0
% Overlay two EEG recordings channel-by-channel for direct visual
% comparison, with epoch/stage navigation driven by a sleep scoring file.
% GUI lives in +viewer/eegCompareViewer.m (shared with EEGCompareViewer.m).
% Requires EEGLAB (pop_interp, pop_select, ...) already on the path.
clc;
run(fullfile(fileparts(mfilename('fullpath')), '..', 'dependancies.m'))
addpath(fileparts(mfilename('fullpath')))

%% ===================== USER SETTINGS =====================
% file1       = '\\vs03.herseninstituut.knaw.nl\VS03-SandC-1\data\nin\data-drop\rawdata\derivatives\prep-ged\sub-drop0001\ses-t1\sub-drop0001_ses-t1_task-sleep_run-01_desc-zc2gedWakeBBAuto_eeg.set';   % first EEG file
file1       = '\\vs03.herseninstituut.knaw.nl\VS03-SandC-1\data\nin\data-drop\rawdata\derivatives\prep-ged\sub-drop0001\ses-t1\sub-drop0001_ses-t1_task-sleep_run-01_desc-zc_eeg.vhdr';   % first EEG file
file2       = '\\vs03.herseninstituut.knaw.nl\VS03-SandC-1\data\nin\data-drop\rawdata\derivatives\prep-ged\sub-drop0001\ses-t1\sub-drop0001_ses-t1_task-sleep_run-01_desc-zc2gedWakeBBAutoFSAutoPlus_eeg.set';   % second EEG file
% file2       = '\\vs03.herseninstituut.knaw.nl\VS03-SandC-1\data\nin\data-drop\rawdata\derivatives\prep-ged\sub-drop0001\ses-t1\sub-drop0001_ses-t1_task-sleep_run-01_desc-zc2gedWakeBBAutoPlusFSAutoPlus_eeg.set';   % second EEG file
% file1       = '\\vs03.herseninstituut.knaw.nl\VS03-SandC-1\data\nin\data-drop\rawdata\derivatives\prep-ged\sub-drop0001\ses-t1\sub-drop0001_ses-t1_task-sleep_run-01_desc-zc2gedWakeBBAutoPlusFSAutoPlus_eeg.set';   % second EEG file
% file2       = '\\vs03.herseninstituut.knaw.nl\VS03-SandC-1\data\nin\data-drop\rawdata\derivatives\prep-ged\sub-drop0001\ses-t1\sub-drop0001_ses-t1_task-sleep_run-01_desc-zc2gedWakeBBAutoFSAutoPlusN2N3REMFSAutoMinus_eeg.set';   % second EEG file
scoringFile = '\\vs03.herseninstituut.knaw.nl\VS03-SandC-1\data\nin\data-drop\derivatives\scores\final\sub-drop0001_ses-t1_task-sleep_run-01_eeg.csv';                         % sleep scoring file, '' if none

net = 'EGI256'; % net type, passed to chans1020
plotDecimation = 5; % only plot every Nth sample (display only; does not affect underlying data)

%% ===================== IMPORT + PREP =====================
EEG1 = fast_eeg_import(file1);
EEG2 = fast_eeg_import(file2);

scoringDigits = [];
if ~isempty(scoringFile)
    scoringDigits = scoreloader(scoringFile);
end


EEG1.data = EEG1.data - sum(EEG1.data, 1) / (size(EEG1.data, 1) + 1);

[EEG1, chanmap] = chans1020(EEG1, false, 'add_eog', true, 'net', net, 'chanprefix', 'E');
[EEG2, chanmap] = chans1020(EEG2, false, 'add_eog', true, 'net', net, 'chanprefix', 'E');

viewer.eegCompareViewer(EEG1, EEG2, scoringDigits, file1, file2, plotDecimation, fieldnames(chanmap));
