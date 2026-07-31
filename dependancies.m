
% External dependancies ---------------------------------------------------
% *************************************************************************
% GEDAI
addpath(genpath('..\GEDAI-master')) % https://github.com/SvennoNito/dusk2dawn

% BIDS
addpath('..\bids-matlab') % https://github.com/bids-standard/bids-matlab

% Sophias preprocessing packages
% addpath('..\chART')
addpath('..\eeg-oscillations') % https://github.com/SvennoNito/eeg-oscillations
% -> oscip.fit_fooof

% EEGLAB
eeglabpath = '..\EEGLAB_2025.0.0'; % https://sccn.ucsd.edu/eeglab/download.php
addpath(eeglabpath) 

% Brainstorm
addpath('..\brainstorm3')

% EEGLAB requires extensions for importing BrainVision files.
% Run `eeglab`, then install needed plugins via: File -> Manage EEGLAB extensions
% Extensions: 
% -> EEG-BIDS 
% -> bva-io1.73

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

