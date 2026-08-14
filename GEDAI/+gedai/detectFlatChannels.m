function [flatmask, flatprop] = detectFlatChannels(data, srate, opts)
% DETECTFLATCHANNELS  Flag channels that stop recording for too much of the session.
%
%   [flatmask, flatprop] = gedai.detectFlatChannels(data, srate)
%   [flatmask, flatprop] = gedai.detectFlatChannels(data, srate, threshold=0.5)
%
%   data    channels x samples
%   srate   sampling rate in Hz
%
%   Why this exists: a dead electrode is invisible to both criteria in clean_channels.
%   Its correlation with the RANSAC reconstruction is 0/0 and its noise-to-signal ratio
%   is 0/0, both NaN, and NaN fails every comparison - so `corrs < corr_threshold` and
%   `znoise > noise_threshold` are false and the channel is silently retained. The one
%   dead channel in sub-drop0001/ses-t1 (E133, flat from roughly the middle of the
%   night) was only caught because it was noisy while it was dying.
%
%   A window counts as flat when the signal does not move by more than `threshold`
%   peak-to-peak across it. Real channels clear that easily: thermal noise alone puts a
%   live electrode well above it, while a dead one sits at exactly one value.
%
%   The window grid and the duration rule match the correlation criterion in
%   clean_channels, so the two are directly comparable: a channel is flagged once it has
%   been flat for longer than maxbrokentime, given as a fraction of the recording
%   (0 < x < 1) or in seconds (x >= 1).
%
%   Run this BEFORE any re-referencing. Average referencing turns a dead channel into
%   minus the common average, which has plenty of variance and is not flat at all.
%
% See also: clean_channels, gedai.plotBadChannels

arguments
    data                  double
    srate           (1,1) double
    opts.windowseconds (1,1) double = 5      % clean_channels default
    opts.maxbrokentime (1,1) double = 0.5    % as passed to clean_channels
    opts.threshold     (1,1) double = 0.5    % peak-to-peak, in data units (uV)
end

[C, S]     = size(data);
window_len = opts.windowseconds * round(srate);
wnd        = 0:window_len-1;
offsets    = 1:window_len:S-window_len;      % same grid as clean_channels
W          = numel(offsets);

if W < 1
    flatmask = false(C,1); flatprop = zeros(C,1);
    return
end

if opts.maxbrokentime > 0 && opts.maxbrokentime < 1
    max_broken_time = S * opts.maxbrokentime;
else
    max_broken_time = round(srate) * opts.maxbrokentime;
end

fprintf('Scanning for flat channels...\n');
flatwin = false(C, W);
for o = 1:W
    seg = data(:, offsets(o)+wnd);
    flatwin(:,o) = (max(seg, [], 2) - min(seg, [], 2)) <= opts.threshold;
end

flatprop = sum(flatwin, 2) ./ W;
flatmask = sum(flatwin, 2) * window_len > max_broken_time;

if any(flatmask)
    fprintf('Flat channels: %s\n', strjoin(arrayfun(@(c) ...
        sprintf('%d (%.0f%% of recording)', c, 100*flatprop(c)), find(flatmask)', 'uni', 0), ', '));
elseif any(flatprop > 0)
    fprintf('Flat stretches seen but under the duration limit (worst: channel %d, %.0f%%).\n', ...
        find(flatprop == max(flatprop), 1), 100*max(flatprop));
end
end
