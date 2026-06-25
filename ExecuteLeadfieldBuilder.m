%%% Start clean
clearvars
clc; close all

%%% Dependancies
run('dependancies.m')

%%% Projects to run
% [BIDS_PM] = smartcache( ...
%     @() bidswizard({'PM'}, '\\vs03.herseninstituut.knaw.nl\VS03-SandD-4', 'Data_BIDS'), ...
%     fullfile(pwd, 'BidsFiles', 'BIDS_PM.mat'), false, {'BIDS_PM'});
% BIDS = BIDS_PM{1};

[BIDS_DROP] = smartcache( ...
    @() bidswizard({'data-drop'}, '\\vs03.herseninstituut.knaw.nl\VS03-SandC-1\data\nin', 'rawdata'), ...
    fullfile(pwd, 'BidsFiles', 'BIDS_DROP.mat'), false, {'BIDS_DROP'});
BIDS = BIDS_DROP{1};

%%% Build leadfield matrix
brainstorm;
build_leadfield_bids(BIDS, 'ProtocolName', 'DROP_Leadfields', 'SubjectFilter', {'drop0001'}, ...
    'QCDir', 'R:\data\nin\data-drop\derivatives\GEDAI\Leadfield\Figures');