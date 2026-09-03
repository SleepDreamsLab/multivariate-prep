% External dependancies ---------------------------------------------------
% *************************************************************************
% All paths below are anchored to THIS file's folder, not to the current
% folder, so dependancies.m can be run from anywhere (e.g. from Pipelines/).
here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end % e.g. when pasted into the command window
sibling = @(name) fullfile(here, '..', name); % repos checked out next to this one

% GEDAI
addpath(genpath(sibling('GEDAI-master'))) % https://github.com/SvennoNito/GEDAI-master
% --> switch from main to sleep-fast branch!

% BIDS
addpath(genpath(sibling('bids-matlab'))) % https://github.com/bids-standard/bids-matlab

% EEGLAB extensions
addpath(genpath(sibling('cleanline'))) % https://github.com/sccn/cleanline
addpath(genpath(sibling('zapline-plus'))) % https://github.com/MariusKlug/zapline-plus
addpath(genpath(sibling('bva-io'))) % https://github.com/sccn/bva-io
addpath(genpath(sibling('pAMICA'))) % https://github.com/sccn/pAMICA
addpath(genpath(sibling('ICLabel'))) % https://github.com/sccn/ICLabel
addpath(genpath(sibling('firfilt'))) % https://github.com/sccn/firfilt
% --> both as sibling checkouts, like the other sccn plugins above. The copies inside
%     the eeglab clone are git submodules and stay empty unless you run
%     `git submodule update --init` in there, so bidsfun_iclabel relies on these.
%     firfilt is not optional: ICLabel's eeg_rpsd calls its windows() on every run.

% Brainstorm
addpath(genpath(sibling('brainstorm3'))) % https://github.com/brainstorm-tools/brainstorm3

% EEGLAB
eeglabpath = sibling('eeglab'); % https://github.com/sccn/eeglab
addpath(eeglabpath)

% EEGLAB startup
if exist('pop_loadset', 'file') ~= 2
    % Checks if EEGLAB has previously been called in this session. If not,
    % call EEGLAB and close the pop-up window.
    eeglab; close;
end

% Sophias preprocessing packages
% addpath(sibling('chART'))
addpath(genpath(sibling('eeg-oscillations'))) % https://github.com/SvennoNito/eeg-oscillations
% -> oscip.fit_fooof

% This folder
addpath(here)                        % -> filterbank, style
addpath(fullfile(here, 'qol'))
addpath(fullfile(here, 'patches'))
addpath(fullfile(here, 'GEDAI'))
addpath(fullfile(here, 'BidsFiles'))
addpath(fullfile(here, 'chanlocs'))
addpath(fullfile(here, 'ICA'))
addpath(fullfile(here, 'GED'))       % -> ged (generalized eigendecomposition)
addpath(fullfile(here, 'Leadfield')) % -> ExecuteLeadfieldBuilder, build_leadfield_bids
addpath(fullfile(here, 'SleepOsci')) % -> run_sleeposci_bids
addpath(fullfile(here, 'colormaps')) % -> slanCM, 200 colormaps (FEX #120088)
clear here sibling
