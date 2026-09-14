function data = readEpochs(fdtFile, nbchan, epochs, epochSamples, opts)
% READEPOCHS  Read selected epochs from an EEGLAB .fdt without loading the whole night.
%
%   data = bench.readEpochs(fdtFile, nbchan, epochs, epochSamples)
%   data = bench.readEpochs(..., 'demean', false)
%
%   epochs    epoch indices (1-based) in the recording; read in the given order and
%             concatenated
%   demean    remove each channel's mean per epoch. The epochs come from different parts
%             of the night, so a slow offset would otherwise become a step at every epoch
%             boundary of the concatenated benchmark recording.            (true)
%
%   data      [nbchan x numel(epochs)*epochSamples] single

arguments
    fdtFile char
    nbchan (1,1) double
    epochs (1,:) double
    epochSamples (1,1) double
    opts.demean (1,1) logical = true
end

fid = fopen(fdtFile, 'r');
if fid < 0
    error('bench:readEpochs:open', 'Cannot open %s.', fdtFile);
end
closer = onCleanup(@() fclose(fid));

data = zeros(nbchan, numel(epochs) * epochSamples, 'single');
for k = 1:numel(epochs)
    x = readOneEpoch(fid, nbchan, epochs(k), epochSamples);
    if opts.demean
        x = x - mean(x, 2);
    end
    data(:, (k - 1) * epochSamples + 1:k * epochSamples) = x;
end
end
