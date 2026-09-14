function x = readOneEpoch(fid, nbchan, iEpoch, epochSamples)
% READONEEPOCH  Read one epoch [nbchan x epochSamples] from an open EEGLAB .fdt file.
% The .fdt is float32, multiplexed (all channels of sample 1, then sample 2, ...), which is
% how pop_saveset writes it and how fast_eeg_import reads it back.
offset = (iEpoch - 1) * epochSamples * nbchan * 4;
if fseek(fid, offset, 'bof') ~= 0
    error('bench:readEpoch:seek', 'Cannot seek to epoch %d: %s', iEpoch, ferror(fid));
end
x = fread(fid, [nbchan, epochSamples], 'float32=>single');
if numel(x) ~= nbchan * epochSamples
    error('bench:readEpoch:short', 'Epoch %d is truncated in the .fdt file.', iEpoch);
end
end
