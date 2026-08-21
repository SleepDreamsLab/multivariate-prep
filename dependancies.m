
% External dependancies ---------------------------------------------------
% *************************************************************************
% GEDAI
addpath(genpath(['..' filesep 'GEDAI-master'])) % https://github.com/SvennoNito/GEDAI-master
% --> switch from main to sleep-fast branch!

% BIDS
addpath(genpath(['..' filesep 'bids-matlab'])) % https://github.com/bids-standard/bids-matlab

% EEGLAB extensions
addpath(genpath(['..' filesep 'cleanline'])) % https://github.com/sccn/cleanline
addpath(genpath(['..' filesep 'zapline-plus'])) % https://github.com/MariusKlug/zapline-plus
addpath(genpath(['..' filesep 'bva-io'])) % https://github.com/sccn/bva-io
addpath(genpath(['..' filesep 'pAMICA'])) % https://github.com/sccn/pAMICA

% Brainstorm
addpath(genpath(['..' filesep 'brainstorm3'])) % https://github.com/brainstorm-tools/brainstorm3

% EEGLAB
eeglabpath = ['..' filesep 'eeglab']; % https://github.com/sccn/eeglab
addpath(eeglabpath) 

% EEGLAB startup
if exist('pop_loadset', 'file') ~= 2
    % Checks if EEGLAB has previously been called in this session. If not,
    % call EEGLAB and close the pop-up window.
    eeglab; close;
end

% Sophias preprocessing packages
% addpath('..\chART')
addpath(genpath(['..' filesep 'eeg-oscillations'])) % https://github.com/SvennoNito/eeg-oscillations
% -> oscip.fit_fooof

% This folder
addpath('qol')
addpath('patches')
addpath('GEDAI')
addpath('BidsFiles')
addpath('chanlocs')
addpath('ICA')

