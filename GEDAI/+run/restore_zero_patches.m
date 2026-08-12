function EEG = restore_zero_patches(EEG)
% RESTORE_ZERO_PATCHES  Put back the all-zero patches removed by run.excise_zero_patches.
%
% Reinserts zeros at the original sample positions so the saved recording keeps its
% original length and timing. Call this once filtering is done, immediately before
% saving.
%
% Safe to call unconditionally: it is a no-op if nothing was ever excised, or if the
% patches have already been restored.
%
% USAGE:
%   EEG = run.restore_zero_patches(EEG);
%
% See also: run.excise_zero_patches, run.run_filter

arguments
    EEG struct
end

if ~isfield(EEG, 'etc') || ~isfield(EEG.etc, 'zeropatches') || isempty(EEG.etc.zeropatches)
    return
end

zp = EEG.etc.zeropatches;
if zp.restored
    return
end

if size(EEG.data, 2) ~= sum(zp.keepmask)
    error('run:restore_zero_patches:lengthMismatch', ...
        ['Data has %d samples but %d were expected from the excision record. ' ...
         'The number of samples must not change between excising and restoring.'], ...
        size(EEG.data, 2), sum(zp.keepmask));
end

fprintf('Restoring %d all-zero patch(es) totalling %.1f s.\n', ...
    numel(zp.starts), sum(zp.lengths) / EEG.srate);

fullData = zeros(size(EEG.data, 1), zp.origpnts, 'like', EEG.data);
fullData(:, zp.keepmask) = EEG.data;

EEG.data  = fullData;
EEG.pnts  = zp.origpnts;
EEG.xmax  = zp.origxmax;
EEG.times = EEG.xmin*1000 + (0:EEG.pnts-1) / EEG.srate * 1000;

EEG.etc.zeropatches.restored = true;

end
