
% External dependancies ---------------------------------------------------
% *************************************************************************
% GEDAI
addpath(genpath(['..' filesep 'GEDAI-master'])) % https://github.com/SvennoNito/dusk2dawn

% BIDS
addpath(['..' filesep 'bids-matlab']) % https://github.com/bids-standard/bids-matlab

% Sophias preprocessing packages
% addpath('..\chART')
addpath(['..' filesep 'eeg-oscillations']) % https://github.com/SvennoNito/eeg-oscillations
% -> oscip.fit_fooof

% EEGLAB
eeglabpath = ['..' filesep 'eeglab']; % https://sccn.ucsd.edu/eeglab/download.php
addpath(eeglabpath) 

% EEGLAB extensions
addpath(['..' filesep 'cleanline'])
addpath(['..' filesep 'zapline-plus'])
addpath(['..' filesep 'bva-io'])

% Brainstorm
addpath(['..' filesep 'brainstorm3'])

% EEGLAB startup
if exist('pop_loadset', 'file') ~= 2
    % Checks if EEGLAB has previously been called in this session. If not,
    % call EEGLAB and close the pop-up window.
    eeglab; close;
end

% This folder
addpath('qol')
addpath('patches')
addpath('GEDAI')
addpath('BidsFiles')
addpath('chanlocs')
addpath('ICA')

