function [spindleMask, spindleSamples] = detect_spindle(eegData, sfreq, options)
    % DETECT_SPINDLE Detect sleep spindles in EEG signals using multitaper spectral analysis.
    %
    % Directly adapts the MT-KCD approach (Oliveira et al. 2020) for sigma-band
    % spindle detection. Candidate regions are identified from median-normalized
    % sigma-band spectrogram power. Spindle boundaries are refined using the
    % Hilbert envelope of the sigma-band filtered signal.
    %
    % USAGE:
    %   spindleMask = detect_spindle(eegData, sfreq)
    %   spindleMask = detect_spindle(eegData, sfreq, 'fmin', 11, 'fmax', 16, 'amin', 10)
    %   spindleMask = detect_spindle(eegData, sfreq, 'stages', stageVector, 'stageFilter', [2 3])
    %   [spindleMask, spindleSamples] = detect_spindle(eegData, sfreq)
    %
    % INPUT:
    %   eegData      : [n_channels x n_samples] EEG signal matrix (in µV)
    %   sfreq        : Sampling frequency (Hz)
    %   options      : Name-value pairs (see below)
    %
    % NAME-VALUE PAIRS:
    %   'fmin'       : Lower sigma band bound (Hz). Default: 11.0
    %   'fmax'       : Upper sigma band bound (Hz). Default: 16.0
    %   'amin'       : Minimum peak envelope amplitude of sigma-band signal (µV). Default: 10.0
    %   'dmin_s'     : Minimum spindle duration (seconds). Default: 0.5
    %   'dmax_s'     : Maximum spindle duration (seconds). Default: 2.0
    %   'q'          : Percentile rank for candidate-region threshold (0–100). Default: 95.0
    %   'stages'     : Optional sleep stage vector [1 x n_epochs]. Each value is a stage number.
    %   'stageFilter': Stages to include (e.g., [2 3] for N2 and N3). If empty, all stages included.
    %   'epochLen'   : Epoch length in seconds (required if using stages). Default: 30 seconds.
    %
    % OUTPUT:
    %   spindleMask    : [n_channels x n_samples] logical matrix, same orientation as
    %                    (possibly transposed) eegData. True where a sample falls
    %                    inside a detected spindle event on that channel.
    %   spindleSamples : [n_channels x 1] cell array, one cell per channel
    %                    (row into eegData, after the transpose above).
    %                    Each cell holds a further [n_events_ch x 1] cell
    %                    array of that channel's detected (and
    %                    artifact-passed) spindle events, each a row
    %                    vector of sample indices (n1:n2, into the
    %                    original eegData columns) spanned by that event.
    %
    % EXAMPLE:
    %   % Load EEG data
    %   load('eeg_data.mat');  % Contains eegData (channels x samples), sfreq
    %
    %   % Detect spindles with default parameters
    %   spindleMask = detect_spindle(eegData, 200);
    %
    %   % Detect spindles in N2/N3 only
    %   stages = [1 1 2 2 2 3 3 1 1];  % Sleep stages for 9 epochs
    %   spindleMask = detect_spindle(eegData, 200, 'stages', stages, ...
    %                           'stageFilter', [2 3], 'epochLen', 30);

    arguments
        eegData (:,:) {mustBeNumeric}
        sfreq (1,1) {mustBePositive}
        options.fmin (1,1) {mustBePositive} = 11.0
        options.fmax (1,1) {mustBePositive} = 16.0
        options.amin (1,1) {mustBePositive} = 10.0
        options.dmin_s (1,1) {mustBePositive} = 0.5
        options.dmax_s (1,1) {mustBePositive} = 2.0
        options.q (1,1) {mustBeInRange(options.q, 0, 100)} = 95.0
        options.stages (:,1) {mustBeInteger} = []
        options.stageFilter (:,1) {mustBeInteger} = []
        options.epochLen (1,1) {mustBePositive} = 30
        options.plotFlag (1,1) logical = false
        options.ampUpper (1,1) double = 150
        options.ampLower (1,1) double = -150
    end

    % Transpose if needed (ensure channels x samples)
    if size(eegData, 1) > size(eegData, 2)
        eegData = eegData';
    end

    n_channels = size(eegData, 1);
    n_samples  = size(eegData, 2);

    % Initialize accumulators. allChannels(i) records which channel allEvents(i,:) came from.
    allEvents   = [];
    allChannels = [];
    allSignals  = {};

    fprintf('detect_spindle: starting detection on %d channel(s), %d samples @ %.1f Hz\n', ...
            n_channels, n_samples, sfreq);

    % Process each channel
    for ch = 1:n_channels
        ch_tic = tic;
        fprintf('detect_spindle: channel %d/%d - starting\n', ch, n_channels);
        signal = double(eegData(ch, :));

        % Apply stage filtering if provided
        if ~isempty(options.stages)
            signal = filterByStages(signal, options.stages, options.stageFilter, ...
                                   options.epochLen, sfreq);
            if isempty(signal)
                fprintf('detect_spindle: channel %d/%d - no samples remain after stage filtering, skipping\n', ...
                        ch, n_channels);
                continue;
            end
        end

        % Detect spindles in this channel
        events_ch = detectSpindle_singleChannel(signal, sfreq, ...
                                               options.fmin, options.fmax, ...
                                               options.amin, options.dmin_s, options.dmax_s, ...
                                               options.q);

        if ~isempty(events_ch)
            for i = 1:size(events_ch, 1)
                n1 = max(1, round(events_ch(i,1) * sfreq));
                n2 = min(length(signal), round(events_ch(i,2) * sfreq));
                allSignals{end+1} = signal(n1:n2);
            end
            allEvents   = [allEvents; events_ch];
            allChannels = [allChannels; repmat(ch, size(events_ch, 1), 1)];
        end

        fprintf('detect_spindle: channel %d/%d - done, %d candidate event(s) found (%.1f s)\n', ...
                ch, n_channels, size(events_ch, 1), toc(ch_tic));
    end

    fprintf('detect_spindle: all channels processed, %d candidate event(s) total before artifact rejection\n', ...
            size(allEvents, 1));

    % Artifact rejection: remove events with any sample outside [ampLower, ampUpper]
    rejEvents = [];
    if ~isempty(allEvents)
        keep = cellfun(@(s) max(s) <= options.ampUpper && min(s) >= options.ampLower, allSignals)';
        n_removed = sum(~keep);
        fprintf('Artifact rejection: removed %d / %d spindle events (amp outside [%.0f, %.0f] µV)\n', ...
                n_removed, numel(keep), options.ampLower, options.ampUpper);
        rejEvents   = allEvents(~keep, :);
        allEvents   = allEvents(keep,  :);
        allChannels = allChannels(keep);
    end

    % Plot overlaid waveforms
    if options.plotFlag
        plotSpindleEvents(eegData, allEvents, rejEvents, sfreq, options.amin, options.ampUpper, options.ampLower);
    end

    % Paint kept events into a channels x samples logical mask
    fprintf('detect_spindle: building %d x %d spindle mask from %d kept event(s)...\n', ...
            n_channels, n_samples, size(allEvents, 1));
    spindleMask = false(n_channels, n_samples);
    spindleSamples = cell(n_channels, 1);
    for ch = 1:n_channels
        spindleSamples{ch} = {};
    end
    for i = 1:size(allEvents, 1)
        n1 = max(1, round(allEvents(i,1) * sfreq));
        n2 = min(n_samples, round(allEvents(i,2) * sfreq));
        ch = allChannels(i);
        spindleMask(ch, n1:n2) = true;
        spindleSamples{ch}{end+1, 1} = n1:n2;
    end

    fprintf('detect_spindle: done, %d / %d sample(s) marked as spindle across %d channel(s)\n', ...
            nnz(spindleMask), numel(spindleMask), n_channels);
end

%% =========================================================================
% Helper function: Filter signal by sleep stages
% =========================================================================

function signalFiltered = filterByStages(signal, stages, stageFilter, epochLen, sfreq)
    if isempty(stageFilter)
        signalFiltered = signal;
        return;
    end

    n_epochs = length(stages);
    n_samples_epoch = round(epochLen * sfreq);

    mask = false(1, length(signal));
    for e = 1:n_epochs
        if ismember(stages(e), stageFilter)
            idx_start = (e - 1) * n_samples_epoch + 1;
            idx_end = min(e * n_samples_epoch, length(signal));
            if idx_start <= length(signal)
                mask(idx_start:idx_end) = true;
            end
        end
    end

    mask = mask(1:length(signal));
    signalFiltered = signal;
    signalFiltered(~mask) = 0;
end

%% =========================================================================
% Main detection routine for a single channel
% =========================================================================

function events = detectSpindle_singleChannel(x, sfreq, fmin, fmax, amin, dmin_s, dmax_s, q)
    x = double(x(:));  % Ensure column vector
    N = length(x);

    % Scale sample-based parameters to actual sfreq (same as MT-KCD)
    L       = max(2, round(sfreq));                  % 1 s window
    delta_j = max(1, round(0.05 * sfreq));           % 0.05 s step
    delta_f = 4.0;                                   % Hz, spectral resolution (fixed)
    Ishort  = max(1, round(0.5 * sfreq / delta_j));  % ≈ 0.5 s of spectrogram cols
    Ibackg  = 10 * Ishort;                           % ≈ 5 s
    Lsmth   = max(1, round(0.15 * sfreq));           % 0.15 s smoothing
    Lbackg  = delta_j * Ibackg;
    Dmin    = round(dmin_s * sfreq);
    Dmax    = round(dmax_s * sfreq);

    % 1. Broadband bandpass (same as MT-KCD)
    fprintf('  [1/5] broadband bandpass filtering (0.3-35 Hz)...\n');
    step_tic = tic;
    x_broad = bandpassFilter(x, sfreq, 0.3, 35.0);
    fprintf('  [1/5] done (%.1f s)\n', toc(step_tic));

    % 2. Multitaper spectrogram
    fprintf('  [2/5] computing multitaper spectrogram (%d column(s) expected)...\n', ceil(N / delta_j));
    step_tic = tic;
    [SG, J, R] = computeSpectrogram(x_broad, sfreq, L, delta_j, delta_f);
    fprintf('  [2/5] done (%.1f s)\n', toc(step_tic));

    % 3. Candidate regions from median-normalized sigma-band power
    fprintf('  [3/5] identifying candidate regions (sigma band %.1f-%.1f Hz)...\n', fmin, fmax);
    step_tic = tic;
    regions = identifyCandidateRegions(SG, J, R, sfreq, fmin, fmax, Ishort, Ibackg, q);
    fprintf('  [3/5] done, %d candidate region(s) found (%.1f s)\n', size(regions, 1), toc(step_tic));
    if isempty(regions)
        events = [];
        return;
    end

    % 4. Sigma-band envelope for boundary detection
    fprintf('  [4/5] computing sigma-band envelope for boundary refinement...\n');
    step_tic = tic;
    x_sigma = bandpassFilter(x, sfreq, fmin, fmax);
    envelope = abs(hilbert(x_sigma));
    env_smooth = centralMovingAverage(envelope, Lsmth);
    env_backg  = centralMovingAverage(env_smooth, Lbackg);
    fprintf('  [4/5] done (%.1f s)\n', toc(step_tic));

    % 5. Find spindle boundaries within candidate regions
    fprintf('  [5/5] refining spindle boundaries within candidate regions...\n');
    step_tic = tic;
    spindles = [];
    for r = 1:size(regions, 1)
        j1 = regions(r, 1);
        j2 = regions(r, 2);

        % Convert spectrogram column indices (1-based) to sample indices (1-based)
        n_start = max(1, (j1 - 1) * delta_j + 1);
        n_end   = min(N, (j2 - 1) * delta_j + Dmax);

        above_bg = env_smooth(n_start:n_end) > env_backg(n_start:n_end);

        % Find rising and falling edges via diff trick
        diff_above = diff([0; above_bg(:); 0]);
        starts_local = find(diff_above == 1);
        ends_local   = find(diff_above == -1);

        for i = 1:min(length(starts_local), length(ends_local))
            % starts_local(i) = index in above_bg where it first becomes True
            % ends_local(i)   = index in above_bg that is the first False after True
            n1  = n_start + starts_local(i) - 1;
            n2  = n_start + ends_local(i)   - 1;
            dur = n2 - n1;

            if n2 > n1
                peak = max(env_smooth(n1:n2-1));
            else
                peak = 0.0;
            end

            if dur >= Dmin && dur < Dmax && peak >= amin
                spindles = [spindles; n1, n2];
            end
        end
    end

    % Convert sample indices to seconds
    if isempty(spindles)
        events = [];
    else
        events = spindles / sfreq;
    end
    fprintf('  [5/5] done, %d spindle(s) passed duration/amplitude criteria (%.1f s)\n', ...
            size(spindles, 1), toc(step_tic));
end

%% =========================================================================
% Internal helpers
% =========================================================================

function y = bandpassFilter(x, sfreq, low, high)
    nyq = sfreq / 2.0;
    Wn = [low/nyq, min(high/nyq, 0.999)];
    [b, a] = butter(4, Wn, 'bandpass');
    y = filtfilt(b, a, x);
end

function y = centralMovingAverage(x, window)
    window = max(1, round(window));
    y = movmean(x, window, 'Endpoints', 'shrink');
end

function p = nextPow2(n)
    p = 2^ceil(log2(n));
end

function [SG, J, R] = computeSpectrogram(x, sfreq, L, delta_j, delta_f)
    % Vectorized multitaper spectrogram.
    %
    % Mathematically identical to the original per-column, per-taper loop
    % (same zero-padding of the last, partial window; same tapers; same
    % FFT length and scaling), but instead of calling fft() once per
    % column per taper (hundreds of thousands of tiny calls for a
    % whole-night recording), it builds a windowed L x nCols matrix and
    % calls fft() once per taper across a whole block of columns at once.
    % This lets MATLAB/FFTW batch the FFTs, which is dramatically faster
    % than the scalar loop. Columns are processed in blocks to bound peak
    % memory use for very long recordings.
    N = length(x);
    x = x(:);
    TW = (L * delta_f) / (2.0 * sfreq);
    K = max(1, round(2 * TW) - 1);

    [tapers, ~] = dpss(L, TW, K);

    R = nextPow2(L);
    J = ceil(N / delta_j);
    half = R / 2 + 1;

    SG = zeros(half, J);

    % Zero-pad the tail so every window (including the last, partial one)
    % can be sliced out directly, matching the original seg(1:avail)=... behavior.
    nPad = max(0, (J - 1) * delta_j + L - N);
    x_padded = [x; zeros(nPad, 1)];

    % Bound the windowed-matrix size to roughly 5e7 elements (~400 MB as double)
    % per block, regardless of window length L.
    blockCols  = max(1, min(J, floor(5e7 / max(L, 1))));
    nBlocks    = ceil(J / blockCols);

    progress_tic = tic;
    report_every_block = max(1, round(nBlocks / 10));  % ~10 progress updates
    blockNum = 0;

    for jStart = 1:blockCols:J
        blockNum = blockNum + 1;
        jEnd  = min(J, jStart + blockCols - 1);
        nCols = jEnd - jStart + 1;

        colStarts = (jStart - 1 : jEnd - 1) * delta_j;   % 0-based window starts
        idxMat    = (0:L-1)' + colStarts + 1;            % L x nCols (1-based indices)
        Xblock    = x_padded(idxMat);                    % L x nCols windowed segments

        SWblock = zeros(half, nCols);
        for k = 1:K
            tapered = tapers(:, k) .* Xblock;            % L x nCols
            fft_out = fft(tapered, R, 1);                % R x nCols
            SWblock = SWblock + (1.0 / sfreq) * abs(fft_out(1:half, :)).^2 / K;
        end

        SG(:, jStart:jEnd) = SWblock;

        if mod(blockNum, report_every_block) == 0 || blockNum == nBlocks
            fprintf('    spectrogram: columns %d-%d/%d (%.0f%%, %.1f s elapsed)\n', ...
                    jStart, jEnd, J, 100 * jEnd / J, toc(progress_tic));
        end
    end

    SG = 10.0 * log10(SG + 1.0);
end

function regions = identifyCandidateRegions(SG, J, R, sfreq, fmin, fmax, Ishort, Ibackg, q)
    freqs = (0:(size(SG, 1)-1))' * (sfreq / R);

    % Sigma-band mask
    mask = freqs >= fmin & freqs <= fmax;
    if ~any(mask)
        regions = [];
        return;
    end

    % dB above per-frequency median baseline
    SG_sigma = SG(mask, :);
    SG_norm  = SG_sigma - median(SG_sigma, 2);

    % Sum normalized sigma-band power per spectrogram column
    C = sum(SG_norm, 1)';

    % Short and background moving averages (same as MT-KCD)
    Cshort = centralMovingAverage(C, Ishort);
    Cbackg = centralMovingAverage(C, Ibackg);
    Cdiff  = Cshort - Cbackg;

    % Percentile threshold
    thresh = prctile(Cdiff, q);

    above = Cdiff >= thresh;
    regions = [];
    in_region = false;
    j1 = 0;

    for j = 1:J
        if above(j) && ~in_region
            j1 = j;
            in_region = true;
        elseif ~above(j) && in_region
            regions = [regions; j1, j - 1];
            in_region = false;
        end
    end

    if in_region
        regions = [regions; j1, J];
    end
end

%% =========================================================================
% Stacked plot of detected spindle events
% =========================================================================

function plotSpindleEvents(eegData, events_s, rejEvents_s, sfreq, amin, ampUpper, ampLower)
    n_samples = size(eegData, 2);
    signal    = double(eegData(1, :));
    pad_s     = 0.5;

    figure('Name', 'Detected Sleep Spindles');
    ax = axes;
    hold(ax, 'on');

    for i = 1:size(rejEvents_s, 1)
        mid_s   = (rejEvents_s(i,1) + rejEvents_s(i,2)) / 2;
        half_s  = (rejEvents_s(i,2) - rejEvents_s(i,1)) / 2;
        n_start = max(1, round((mid_s - half_s - pad_s) * sfreq));
        n_end   = min(n_samples, round((mid_s + half_s + pad_s) * sfreq));
        seg     = signal(n_start:n_end);
        t       = (0:length(seg)-1) / sfreq - (half_s + pad_s);
        h = plot(ax, t, seg, 'Color', [0.85 0.15 0.15], 'LineWidth', 0.8);
        h.Color(4) = 0.1;
    end

    for i = 1:size(events_s, 1)
        mid_s   = (events_s(i,1) + events_s(i,2)) / 2;
        half_s  = (events_s(i,2) - events_s(i,1)) / 2;
        n_start = max(1, round((mid_s - half_s - pad_s) * sfreq));
        n_end   = min(n_samples, round((mid_s + half_s + pad_s) * sfreq));
        seg     = signal(n_start:n_end);
        t       = (0:length(seg)-1) / sfreq - (half_s + pad_s);
        h = plot(ax, t, seg, 'Color', [0.15 0.15 0.15], 'LineWidth', 0.8);
        h.Color(4) = 0.1;
    end

    all_events = [events_s; rejEvents_s];
    max_dur = 0;
    if ~isempty(all_events)
        max_dur = max(all_events(:,2) - all_events(:,1));
    end

    xline(ax, 0, 'b--', 'Spindle centre');
    yline(ax,  amin, 'b:', sprintf('%.0f µV (amin)', amin));
    yline(ax, -amin, 'b:', sprintf('-%.0f µV', amin));
    if isfinite(ampUpper)
        yline(ax, ampUpper, 'r:', sprintf('%.0f µV', ampUpper));
    end
    if isfinite(ampLower)
        yline(ax, ampLower, 'r:', sprintf('%.0f µV', ampLower));
    end
    xlabel(ax, 'Time relative to spindle onset (s)');
    ylabel(ax, 'Amplitude (\muV)');
    title(ax, sprintf('Sleep spindles: %d kept (dark), %d rejected (red) — ch 1', ...
                      size(events_s, 1), size(rejEvents_s, 1)));
    if max_dur > 0
        xlim(ax, [-(max_dur/2 + pad_s), max_dur/2 + pad_s]);
    end
    box(ax, 'off');
    hold(ax, 'off');
end
