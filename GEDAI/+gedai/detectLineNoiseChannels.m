function [mask, znoise] = detectLineNoiseChannels(data, srate, opts)
% DETECTLINENOISECHANNELS  Flag channels whose line noise SURVIVED cleaning.
%
%   [mask, znoise] = gedai.detectLineNoiseChannels(data, srate)
%   [mask, znoise] = gedai.detectLineNoiseChannels(data, srate, threshold=4)
%
%   data   channels x samples, AFTER Zapline/CleanLine
%   srate  sampling rate in Hz
%
%   Same statistic clean_channels uses - the ratio of >45 Hz to <45 Hz robust amplitude,
%   robust-z-scored across channels - but applied after line-noise removal rather than
%   before, which asks a different and more useful question. Before cleaning it measures
%   how much mains a channel picks up, and on a HydroCel net that ranks the electrodes
%   under the lead bundle at the vertex highest: real, but Zapline's job, and those
%   channels have intact EEG (reconstruction correlation ~0.99). After cleaning it
%   measures what Zapline could not fix, which is the thing that actually warrants
%   removal. Measured on sub-drop0001: pre-Zapline this flags the vertex ring, at
%   z = 41/37/32 against a p95 of 2; post-Zapline it flags peripheral channels instead.
%
%   Returns a mask over the channels of `data` as given - if bad channels were already
%   removed upstream, it is in that reduced frame, not the full montage.
%
% See also: clean_channels, gedai.detectFlatChannels

arguments
    data                     double
    srate              (1,1) double
    opts.threshold     (1,1) double = 4
end

nCh = size(data, 1);
if srate <= 100
    %% clean_channels only computes this above 100 Hz - below it there is no band left
    %% above the 45-50 Hz lowpass to form the numerator from.
    mask = false(nCh, 1); znoise = zeros(1, nCh);
    return
end

fprintf('Scanning for residual line noise ...\n');

%%% Split at 45-50 Hz exactly as clean_channels does, so the two are comparable
B = firls(100, [2*[0 45 50]/srate 1], [1 1 0 0]);
X = zeros(size(data, 2), nCh);
for c = nCh:-1:1
    X(:,c) = filtfilt(B, 1, data(c,:)');
end

md        = @(A) median(abs(A - median(A)));      % clean_channels' local mad()
noisiness = md(data' - X) ./ md(X);
znoise    = (noisiness - median(noisiness)) ./ (md(noisiness) * 1.4826);
mask      = znoise(:) > opts.threshold;

if any(mask)
    fprintf('Residual line noise above z = %g on channel(s): %s\n', ...
        opts.threshold, mat2str(find(mask)'));
else
    fprintf('No channel above z = %g after cleaning (max %.1f).\n', ...
        opts.threshold, max(znoise));
end
end
