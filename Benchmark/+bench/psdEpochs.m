function [P, f] = psdEpochs(dataList, weights, stageEpoch, stages, srate, epochLength, winSec)
% PSDEPOCHS  Mean Welch PSD per stage of a weighted sum of data matrices, epoch by epoch.
%
%   [P, f] = bench.psdEpochs({B, A}, [1 1], stageEpoch, stages, srate, epochLength, 4)
%
%   Each 30-s epoch is summed from the matrices in dataList (so B + A never has to exist as a
%   whole) and passed through pwelch on its own - windows never straddle the discontinuity
%   between two concatenated epochs.
%
%   dataList   cell of [nChan x nSamples] matrices of equal size
%   weights    one weight per matrix
%   P          [numel(stages) x nFreq x nChan], NaN for a stage without epochs
%   f          [nFreq x 1] Hz

arguments
    dataList cell
    weights (1,:) double
    stageEpoch (1,:) double
    stages (1,:) double
    srate (1,1) double
    epochLength (1,1) double
    winSec (1,1) double = 4
end

epochSamples = round(epochLength * srate);
nChan  = size(dataList{1}, 1);
window = round(winSec * srate);
nFreq  = floor(window / 2) + 1;
P = nan(numel(stages), nFreq, nChan);
f = (0:nFreq - 1)' * srate / window;

for iStage = 1:numel(stages)
    epochs = find(stageEpoch == stages(iStage));
    if isempty(epochs), continue; end
    total = zeros(nFreq, nChan);
    for iEpoch = epochs
        cols = (iEpoch - 1) * epochSamples + 1:iEpoch * epochSamples;
        x = zeros(nChan, epochSamples);
        for d = 1:numel(dataList)
            x = x + weights(d) * double(dataList{d}(:, cols));
        end
        [pxx, f] = pwelch(x', hann(window), round(window / 2), window, srate);
        total = total + pxx;
    end
    P(iStage, :, :) = total / numel(epochs);
end
end
