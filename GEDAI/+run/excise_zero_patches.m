function [EEG, nPatches] = excise_zero_patches(EEG, opts)
% EXCISE_ZERO_PATCHES  Cut all-zero patches out of a continuous EEG struct.
%
% Amplifiers that crash and restart sometimes have the lost time padded with a patch of
% exact zeros so the recording keeps its wall-clock length. Such a patch carries no
% signal, but it does damage every filtering step applied to it: an IIR high-pass run
% through filtfilt smears it into a large exponential transient, and Zapline-plus errors
% out on any chunk that falls entirely inside it.
%
% This function removes those patches so the whole filter chain only ever sees real data.
% The removed samples are recorded in EEG.etc.zeropatches, and run.restore_zero_patches
% puts them back (as zeros, at their original sample positions) once filtering is done.
%
% IMPORTANT: run this on raw data, BEFORE any filtering. The patch is only exactly zero
% until a high-pass filter touches it.
%
% USAGE:
%   EEG = run.excise_zero_patches(EEG)
%   [EEG, nPatches] = run.excise_zero_patches(EEG, minseconds=5)
%
% OPTIONAL NAME-VALUE:
%   minseconds   only patches longer than this are removed   (default 5)
%
% OUTPUTS:
%   EEG        — EEG struct with the patches cut out of .data and .pnts/.times/.xmax
%                updated to match; .etc.zeropatches holds what is needed to restore it
%   nPatches   — number of patches removed
%
% See also: run.restore_zero_patches, run.run_filter

arguments
    EEG                  struct
    opts.minseconds (1,1) double = 5
end

nPatches = 0;

if EEG.trials > 1
    warning('run:excise_zero_patches:notContinuous', ...
        'Data is epoched (%d trials) - skipping zero-patch removal.', EEG.trials)
    return
end

%%% Find samples that are exactly zero on every channel
zeroSample = all(EEG.data == 0, 1);

%%% Group them into runs, keep only runs longer than the minimum
d         = diff([false, zeroSample, false]);
runStarts = find(d ==  1);
runEnds   = find(d == -1) - 1;
runLength = runEnds - runStarts + 1;

isPatch   = runLength > round(opts.minseconds * EEG.srate);
runStarts = runStarts(isPatch);
runEnds   = runEnds(isPatch);
runLength = runLength(isPatch);

nPatches = numel(runStarts);
if nPatches == 0
    return
end

%%% Build the keep mask
keepMask = true(1, EEG.pnts);
for iPatch = 1:nPatches
    keepMask(runStarts(iPatch):runEnds(iPatch)) = false;
end

fprintf('Removing %d all-zero patch(es) totalling %.1f s before filtering (restored before saving):\n', ...
    nPatches, sum(runLength) / EEG.srate);
for iPatch = 1:nPatches
    fprintf('  patch %d: %.1f - %.1f s (%.1f s)\n', iPatch, ...
        (runStarts(iPatch)-1)/EEG.srate, runEnds(iPatch)/EEG.srate, runLength(iPatch)/EEG.srate);
end

%%% Record what is needed to put it back, then cut
EEG.etc.zeropatches = struct( ...
    'keepmask',  keepMask, ...
    'origpnts',  EEG.pnts, ...
    'origxmax',  EEG.xmax, ...
    'starts',    runStarts, ...
    'lengths',   runLength, ...
    'minseconds', opts.minseconds, ...
    'restored',  false);

EEG.data  = EEG.data(:, keepMask);
EEG.pnts  = size(EEG.data, 2);
EEG.xmax  = EEG.xmin + (EEG.pnts - 1) / EEG.srate;
EEG.times = EEG.xmin*1000 + (0:EEG.pnts-1) / EEG.srate * 1000;

end
