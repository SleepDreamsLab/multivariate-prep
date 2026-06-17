
% External dependancies ---------------------------------------------------
% *************************************************************************
addpath(genpath(pwd))

% BIDS
addpath('..\bids-matlab') % https://github.com/bids-standard/bids-matlab

%%% Sophias preprocessing packages
% addpath('..\chART')
addpath('..\episl-preprocessing\functions') % https://github.com/snipeso/episl-preprocessing
addpath('..\sleep-prep') % https://github.com/SvennoNito/sleep-prep
addpath('..\eeg-oscillations') % https://github.com/SvennoNito/eeg-oscillations
addpath('..\dusk2dawn') % https://github.com/SvennoNito/dusk2dawn
addpath(genpath('..\GEDAI-master')) % https://github.com/SvennoNito/dusk2dawn

%%% MATLAB package namespaces (parent dir must be on path, not the +pkg dir itself)
% +qol, +organon: live in this repo (multivariate-prep) — already on path via addpath(pwd)
% +gedai: C:\Postdoc\Code\+gedai
addpath('C:\Postdoc\Code')
% +ica: C:\Postdoc\Code\exploratory-prep\+ica
addpath('C:\Postdoc\Code\exploratory-prep')

% %%% Fieldtrip (only needed for brainvision files saved as INT32
% addpath('..\sleeptrip')

%%% EEGLAB
eeglabpath = '..\EEGLAB_2025.0.0'; % https://sccn.ucsd.edu/eeglab/download.php
addpath(eeglabpath) 

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
