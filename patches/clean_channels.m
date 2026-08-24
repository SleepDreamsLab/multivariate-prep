function [signal,removed_channels, corrs, znoise] = clean_channels(signal,corr_threshold,noise_threshold,window_len,max_broken_time,num_samples,subset_size,window_stride,ramsaver)
% Remove channels with abnormal data from a continuous data set.
% Signal = clean_channels(Signal,CorrelationThreshold,LineNoiseThreshold,WindowLength,MaxBrokenTime,NumSamples,SubsetSize)
%
% This is an automated artifact rejection function which ensures that the data contains no channels
% that record only noise for extended periods of time. If channels with control signals are
% contained in the data these are usually also removed. The criterion is based on correlation: if a
% channel has lower correlation to its robust estimate (based on other channels) than a given threshold
% for a minimum period of time (or percentage of the recording), it will be removed.
%
% In:
%   Signal          : Continuous data set, assumed to be appropriately high-passed (e.g. >0.5Hz or
%                     with a 0.5Hz - 2.0Hz transition band).
%
%   CorrelationThreshold : Correlation threshold. If a channel is correlated at less than this value
%                          to its robust estimate (based on other channels), it is considered abnormal in
%                          the given time window. Default: 0.85.
%                     
%   LineNoiseThreshold : If a channel has more line noise relative to its signal than this value, in
%                        standard deviations from the channel population mean, it is considered abnormal.
%                        Default: 4.
%
%
%   The following are detail parameters that usually do not have to be tuned. If you cannot get
%   the function to do what you want, you might consider adapting these to your data.
%   
%   WindowLength    : Length of the windows (in seconds) for which correlation is computed; ideally
%                     short enough to reasonably capture periods of global artifacts or intermittent 
%                     sensor dropouts, but not shorter (for statistical reasons). Default: 5.
% 
%   MaxBrokenTime : Maximum time (either in seconds or as fraction of the recording) during which a 
%                   retained channel may be broken. Reasonable range: 0.1 (very aggressive) to 0.6
%                   (very lax). The default is 0.4.
%
%   NumSamples : Number of RANSAC samples. This is the number of samples to generate in the random
%                sampling consensus process. The larger this value, the more robust but also slower 
%                the processing will be. Default: 50.
%
%   SubsetSize : Subset size. This is the size of the channel subsets to use for robust reconstruction,
%                as a fraction of the total number of channels. Default: 0.25.
%
%   RamSaver : (local addition) if true, the >50Hz removal filter is applied one channel
%              at a time instead of to the whole data matrix at once, trading speed for a
%              lower peak memory footprint. Default: false.
%
% Out:
%   Signal : data set with bad channels removed
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2014-05-12

% Copyright (C) Christian Kothe, SCCN, 2014, christian@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify it under the terms of the GNU
% General Public License as published by the Free Software Foundation; either version 2 of the
% License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
% even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License along with this program; if not,
% write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
% USA

if ~exist('corr_threshold','var') || isempty(corr_threshold) corr_threshold = 0.8; end
if ~exist('noise_threshold','var') || isempty(noise_threshold) noise_threshold = 4; end
if ~exist('window_len','var') || isempty(window_len) window_len = 5; end
if ~exist('max_broken_time','var') || isempty(max_broken_time) max_broken_time = 0.4; end
if ~exist('num_samples','var') || isempty(num_samples) num_samples = 50; end
if ~exist('subset_size','var') || isempty(subset_size) subset_size = 0.25; end
if ~exist('reset_rng','var') || isempty(reset_rng) reset_rng = true; end
% SPEED (local addition): evaluate every window_stride-th window instead of every one.
% The correlation criterion is a *proportion* of windows below threshold compared against
% max_broken_time, so subsampling costs precision, not correctness: at stride 2 on a 9-h
% night the proportion is still estimated from ~3300 windows, s.e. ~0.9%. Only channels
% sitting within about a percentage point of the threshold can flip. 1 = original.
if ~exist('window_stride','var') || isempty(window_stride) window_stride = 1; end
% RAM SAVER (local addition): filter one channel at a time instead of the whole
% signal.data matrix at once. filtfilt on the full transposed matrix briefly holds
% several full-size copies (transpose, double-cast, filtered output) at once; looping
% keeps only one channel's vector live at a time, at the cost of speed. Off by default.
if ~exist('ramsaver','var') || isempty(ramsaver) ramsaver = false; end

subset_size = round(subset_size*size(signal.data,1)); 

% flag channels
if max_broken_time > 0 && max_broken_time < 1  %#ok<*NODEF>
    max_broken_time = size(signal.data,2)*max_broken_time;
else
    max_broken_time = round(signal.srate)*max_broken_time;
end

% RANSAC (below) needs X in double: the projector P is double, and a single x double
% matmul falls off the BLAS fast path and converts per window; the noisiness ratio also
% needs double, because mad(data-X) is a small difference of large numbers that cancels
% badly in single. Tested, rejected - don't retry.
%
% signal.data itself, however, is deliberately NOT cast to double here. It used to be
% (`signal.data = double(signal.data)`), unconditionally, ramsaver or not - a full extra
% double-sized copy of the whole recording (~13.7 GB at 256 ch x 8 h x 250 Hz) held
% alive purely so the noisiness ratio below could subtract it from X. That copy is never
% used for anything else: the caller (detectBadChannels/smartcache) discards this
% function's `signal` output entirely and keeps only removed_channels/corrs/znoise. So
% under ramsaver, X and the noisiness ratio are now computed one channel at a time
% instead, and signal.data is cast to double only per-channel, transiently. Peak memory
% drops by roughly two full-recording double copies (the old signal.data cast, and its
% transpose at the old `signal.data'-X`) - the only full-size double array still held is
% X itself, which the RANSAC step below needs intact regardless.
[C,S] = size(signal.data);
window_len = window_len*round(signal.srate);
wnd = 0:window_len-1;
offsets = 1:window_len*window_stride:S-window_len;
W = length(offsets);

fprintf('Scanning for bad channels...\n');

if signal.srate > 100
    % remove signal content above 50Hz
    B = firls(100,[2*[0 45 50]/signal.srate 1],[1 1 0 0]);
    % determine z-scored level of EM noise-to-signal ratio for each channel
    % TRIED AND REJECTED: dividing by median(mad(X,1)) - the typical channel's
    % low-frequency amplitude - instead of by each channel's own, on the theory that the
    % vertex electrodes were being flagged for having little EEG rather than much mains.
    % They are not: with the denominator removed entirely the same three channels still
    % topped the montage on sub-drop0001, and by more (100/81/70 vs 41/37/32), while the
    % flagged count rose from 13 to 19. Those electrodes genuinely carry far more >45 Hz
    % energy than the rest of the head - plausibly because the net's lead bundle exits at
    % the vertex - so no denominator and no threshold can separate them. Kept upstream.
    if ramsaver
        % Looping from C down to 1 preallocates X on the first (largest) iteration -
        % same trick as the original commented-out loop above. Each channel's double
        % cast, filter output, and noisiness ratio are computed together and the raw
        % double copy dropped immediately, so no full-recording double copy of
        % signal.data (or of signal.data'-X) is ever held alongside X.
        noisiness = zeros(1,C);
        for c=C:-1:1
            fprintf('Filtering channel %d/%d ...\n', c, C);
            raw_c = double(signal.data(c,:))';
            X(:,c) = filtfilt(B,1,raw_c);
            noisiness(c) = mad(raw_c-X(:,c))./mad(X(:,c),1);
        end
    else
        X = filtfilt(B, 1, double(signal.data'));
        noisiness = mad(double(signal.data)'-X)./mad(X,1);
    end
    znoise = (noisiness - median(noisiness)) ./ (mad(noisiness,1)*1.4826);
    % trim channels based on that
    noise_mask = znoise > noise_threshold;
else
    X = double(signal.data)';
    noise_mask = false(C,1)'; % transpose added. Otherwise gives an error below at removed_channels = removed_channels | noise_mask';  (by Ozgur Balkan)
end

if ~(isfield(signal.chanlocs,'X') && isfield(signal.chanlocs,'Y') && isfield(signal.chanlocs,'Z') && all([length([signal.chanlocs.X]),length([signal.chanlocs.Y]),length([signal.chanlocs.Z])] > length(signal.chanlocs)*0.5))
    error('clean_channels:bad_chanlocs','To use this function most of your channels should have X,Y,Z location measurements.'); end

% get the matrix of all channel locations [3xN]
[x,y,z] = deal({signal.chanlocs.X},{signal.chanlocs.Y},{signal.chanlocs.Z});
usable_channels = find(~cellfun('isempty',x) & ~cellfun('isempty',y) & ~cellfun('isempty',z));
locs = [cell2mat(x(usable_channels));cell2mat(y(usable_channels));cell2mat(z(usable_channels))];
% Skip the reindex when it would be a no-op (the normal case once chanlocs-less
% channels are already filtered out upstream): X is full-recording double, so
% X(:,usable_channels) otherwise allocates a second ~13.7 GB copy for nothing.
if numel(usable_channels) ~= C || ~isequal(usable_channels(:)', 1:C)
    X = X(:,usable_channels);
end
  
% caculate all-channel reconstruction matrices from random channel subsets   
if reset_rng
    rng('default')
end
P = calc_projector(locs,num_samples,subset_size);
corrs = zeros(length(usable_channels),W);
        
% calculate each channel's correlation to its RANSAC reconstruction for each window
timePassedList = zeros(W,1);
for o=1:W
    tic; % makoto
    XX = X(offsets(o)+wnd,:);
    YY = sort(reshape(XX*P,length(wnd),length(usable_channels),num_samples),3);
    YY = YY(:,:,round(end/2));
	corrs(:,o) = sum(XX.*YY)./(sqrt(sum(XX.^2)).*sqrt(sum(YY.^2)));
    timePassedList(o) = toc; % makoto
    medianTimePassed = median(timePassedList(1:o));
    fprintf('clean_channel: %3.0d/%d blocks, %.1f minutes remaining.\n', o, W, medianTimePassed*(W-o)/60); % makoto
end
        
flagged = corrs < corr_threshold;
        
% mark all channels for removal which have more flagged samples than the maximum number of
% ignored samples
removed_channels = false(C,1);
% window_stride scales the flagged count back onto the full recording, so max_broken_time
% keeps meaning the same fraction of the night whatever the stride.
removed_channels(usable_channels) = sum(flagged,2)*window_len*window_stride > max_broken_time;
removed_channels = removed_channels | noise_mask';

% apply removal
if mean(removed_channels) > 0.75
    error('clean_channels:bad_chanlocs','More than 75%% of your channels were removed -- this is probably caused by incorrect channel location measurements (e.g., wrong cap design).');
elseif any(removed_channels)
    try
        signal = pop_select(signal,'nochannel',find(removed_channels));
    catch e
        if ~exist('pop_select','file')
            disp('Apparently you do not have EEGLAB''s pop_select() on the path.');
        else
            disp('Could not select channels using EEGLAB''s pop_select(); details: ');
            hlp_handleerror(e,1);
        end
        fprintf('Removing %i channels and dropping signal meta-data.\n',nnz(removed_channels));
        if length(signal.chanlocs) == size(signal.data,1)
            signal.chanlocs = signal.chanlocs(~removed_channels); end
        signal.data = signal.data(~removed_channels,:);
        signal.nbchan = size(signal.data,1);
        [signal.icawinv,signal.icasphere,signal.icaweights,signal.icaact,signal.stats,signal.specdata,signal.specicaact] = deal([]);
    end
    if isfield(signal.etc,'clean_channel_mask')
        signal.etc.clean_channel_mask(signal.etc.clean_channel_mask) = ~removed_channels;
    else
        signal.etc.clean_channel_mask = ~removed_channels;
    end
end



% calculate a bag of reconstruction matrices from random channel subsets
function P = calc_projector(locs,num_samples,subset_size)
%stream = RandStream('mt19937ar','Seed',435656);
rand_samples = {};
for k=num_samples:-1:1
    tmp = zeros(size(locs,2));
    subset = randsample(1:size(locs,2),subset_size);
%    subset = randsample(1:size(locs,2),subset_size,stream);
    tmp(subset,:) = real(sphericalSplineInterpolate(locs(:,subset),locs))';
    rand_samples{k} = tmp;
end
P = horzcat(rand_samples{:});


function Y = randsample(X,num)
Y = [];
while length(Y)<num
    pick = round(1 + (length(X)-1).*rand());
    Y(end+1) = X(pick);
    X(pick) = [];
end
% 
% function Y = randsample(X,num,stream)
% Y = [];
% while length(Y)<num
%     pick = round(1 + (length(X)-1).*rand(stream));
%     Y(end+1) = X(pick);
%     X(pick) = [];
% end

function Y = mad(X,flag) %#ok<INUSD>
Y = median(abs(bsxfun(@minus,X,median(X))));
