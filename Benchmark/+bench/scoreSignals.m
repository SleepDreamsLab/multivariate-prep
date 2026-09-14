function [metrics, detail] = scoreSignals(S, Xhat, Bhat, B, events, classes, stageEpoch, srate, epochLength, opts)
% SCORESIGNALS  Over-cleaning metrics for injected graphoelements: SER and detector recovery.
%
%   [metrics, detail] = bench.scoreSignals(S, Xhat, Bhat, B, events, classes, stageEpoch, srate, epochLength)
%
%   All data are [nChan x nSamples], average-referenced:
%     S     injected ground truth
%     B     background (the benchmark input without injections)
%     Bhat  GEDAI(B)
%     Xhat  GEDAI(B + S)
%   The uncleaned contaminated input X = B + S is formed channel by channel where needed.
%
%   Paired design
%   -------------
%   GEDAI is nonlinear and the background is real EEG with its own artifacts, so the part of
%   the output that belongs to the injection is estimated as Shat = Xhat - Bhat: identical
%   background, identical GEDAI configuration, only the injection differs. Whatever GEDAI
%   does to the background alone cancels; what remains is how it treated the injected
%   events plus any change in how the background was cleaned because they were there. Both
%   count as error.
%
%   Per class and stage (and pooled as stage 'all'); windows are each event +/- eventpad
%   ----------------------------------------------------------------------------------
%   SER_dB             10 log10( sum S^2 / sum (Shat - S)^2 ), pooled over channels and
%                      windows. Signal-to-error ratio in the sense of Somers et al. (2018,
%                      J Neural Eng 15:036007), here on injected events instead of clean
%                      segments. Higher is better; 0 dB = error as large as the signal.
%   SER_raw_dB         the same with Shat = Xhat - B (unpaired, Somers' form). Only
%                      meaningful on a synthetic background, where B itself is clean.
%   gain               <Shat, S> / <S, S>: fraction of the injected field preserved
%                      (1 = intact, 0 = removed); insensitive to uncorrelated error.
%   SER_event_median, SER_event_p10   distribution over events (p10 = the worst decile)
%   gain_peakchan_median   gain at each event's peak channel
%   ampratio_median, ampratio_p10     band-limited amplitude of Shat over S at the peak
%                      channel: envelope maximum for spindles, peak-to-peak for slow waves
%   morphcorr_median   correlation of Shat and S at the peak channel within the event
%   nEvents
%
%   Detector recovery - per class, on each distinct peak channel of that class
%   --------------------------------------------------------------------------
%   Thresholds are calibrated once on the uncleaned input X and frozen for Xhat, B and Bhat
%   (see bench.detectSpindles / bench.detectSlowWaves). An injected spindle counts as
%   recovered when a detection on its peak channel overlaps it (+/- matchtolerance); an
%   injected slow wave when a detected trough lies within matchtolerance of its own trough
%   (window overlap alone is met by chance in dense N3 background).
%   recall_in, recall_out      fraction of injected events detected in X / in Xhat
%   recall_retention           recall_out / recall_in
%   density_in, density_out    detections per minute in X / Xhat
%   density_bg_in, density_bg_out   the same in B / Bhat (real background events)
%   excess_density_in, excess_density_out   density minus the matching background density:
%                              the injected events the detector finds, per minute
%   density_retention          excess_density_out / excess_density_in
%   injected_density           injected events per minute
%   detamp_ratio_median        detector amplitude Xhat / X for events recovered in both
%   detfreq_shift_median       Hz, detected spindle frequency Xhat - X (spindles)
%   detduration_ratio_median   detected duration Xhat / X
%
%   Name-value
%   ----------
%   eventpad        s around each event                                       (0.25)
%   validmask       [1 x nSamples] logical, samples usable for detection      (all)
%   matchtolerance  s                                                         (0.2)
%   erplength       s, length of the stored event-locked averages             (3)
%
%   Outputs
%   -------
%   metrics   table: class, stage, metric, value
%   detail    .events  per-event table (id, class, SER_dB, gain, gain_peakchan, ampratio,
%                      morphcorr, hit_in, hit_out)
%             .erp.<class>  event-locked means of S and Shat at the peak channel (t, S, Shat)
%             .example.<class>  the median-SER event: t, S, Shat, X, Xhat at its peak channel

arguments
    S single
    Xhat single
    Bhat single
    B single
    events table
    classes struct
    stageEpoch (1,:) double
    srate (1,1) double
    epochLength (1,1) double
    opts.eventpad (1,1) double = 0.25
    opts.validmask (1,:) logical = logical([])
    opts.matchtolerance (1,1) double = 0.2
    opts.erplength (1,1) double = 3
end

nSamples      = size(S, 2);
epochSamples  = round(epochLength * srate);
stageAtSample = repelem(stageEpoch, epochSamples);
valid         = opts.validmask;
if isempty(valid), valid = true(1, nSamples); end
pad    = round(opts.eventpad * srate);
nErp   = round(opts.erplength * srate);
rows   = {};
detail = struct('events', table(), 'erp', struct(), 'example', struct());
eventTables = {};

for k = 1:numel(classes)
    c = classes(k);
    ev = events(strcmp(events.trial_type, c.name), :);
    nEv = height(ev);
    if nEv == 0, continue; end

    %%% ---- Per-event SER, gain, amplitude and morphology ----
    [s2, e2, eRaw2, cross] = deal(zeros(nEv, 1));
    [gainPc, ampRatio, morphCorr] = deal(nan(nEv, 1));
    erpS = nan(nEv, nErp); erpShat = nan(nEv, nErp);
    for i = 1:nEv
        a = max(1, ev.sample_start(i) - pad);
        b = min(nSamples, ev.sample_end(i) + pad);
        truth = double(S(:, a:b));
        estimate = double(Xhat(:, a:b)) - double(Bhat(:, a:b));
        s2(i)    = sum(truth.^2, 'all');
        e2(i)    = sum((estimate - truth).^2, 'all');
        eRaw2(i) = sum((double(Xhat(:, a:b)) - double(B(:, a:b)) - truth).^2, 'all');
        cross(i) = sum(estimate .* truth, 'all');

        pc = ev.peakchan_idx(i);
        tPc = truth(pc, :);
        ePc = estimate(pc, :);
        gainPc(i) = (ePc * tPc') / (tPc * tPc');
        r = corrcoef(tPc, ePc);
        morphCorr(i) = r(1, 2);
        ampRatio(i) = amplitudeRatio(S, Xhat, Bhat, pc, a, b, c, srate);

        cols = ev.sample_start(i):min(nSamples, ev.sample_start(i) + nErp - 1);
        erpS(i, 1:numel(cols))    = S(pc, cols);
        erpShat(i, 1:numel(cols)) = Xhat(pc, cols) - Bhat(pc, cols);
    end
    serEvent = 10 * log10(s2 ./ e2);
    gainEvent = cross ./ s2;

    %%% ---- Detector recovery on each peak channel ----
    [hitIn, hitOut, rec] = detectorRecovery(S, Xhat, Bhat, B, ev, c, stageAtSample, valid, srate, ...
        round(opts.matchtolerance * srate));

    %%% ---- Aggregate per stage and pooled ----
    stageList = unique(ev.stage)';
    for iStage = 0:numel(stageList)
        if iStage == 0
            sel = true(nEv, 1); label = 'all'; recStage = rec.all;
        else
            sel = ev.stage == stageList(iStage); label = bench.stageLabel(stageList(iStage));
            recStage = rec.(label);
        end
        add = @(name, value) {c.name, label, name, value};
        rows = [rows; ...
            add('nEvents',              nnz(sel)); ...
            add('SER_dB',               10 * log10(sum(s2(sel)) / sum(e2(sel)))); ...
            add('SER_raw_dB',           10 * log10(sum(s2(sel)) / sum(eRaw2(sel)))); ...
            add('gain',                 sum(cross(sel)) / sum(s2(sel))); ...
            add('SER_event_median',     median(serEvent(sel))); ...
            add('SER_event_p10',        prctile(serEvent(sel), 10)); ...
            add('gain_peakchan_median', median(gainPc(sel), 'omitnan')); ...
            add('ampratio_median',      median(ampRatio(sel), 'omitnan')); ...
            add('ampratio_p10',         prctile(ampRatio(sel), 10)); ...
            add('morphcorr_median',     median(morphCorr(sel), 'omitnan')); ...
            add('recall_in',            mean(hitIn(sel))); ...
            add('recall_out',           mean(hitOut(sel))); ...
            add('recall_retention',     mean(hitOut(sel)) / mean(hitIn(sel)))]; %#ok<AGROW>
        for f = fieldnames(recStage)'
            rows(end+1, :) = add(f{1}, recStage.(f{1})); %#ok<AGROW>
        end
    end

    %%% ---- Detail ----
    eventTables{end+1} = table(ev.id, repmat({c.name}, nEv, 1), ev.stage, serEvent, gainEvent, ...
        gainPc, ampRatio, morphCorr, hitIn, hitOut, 'VariableNames', {'id', 'trial_type', 'stage', ...
        'SER_dB', 'gain', 'gain_peakchan', 'ampratio', 'morphcorr', 'hit_in', 'hit_out'}); %#ok<AGROW>
    detail.erp.(c.name) = struct('t', (0:nErp - 1) / srate, 'S', mean(erpS, 1, 'omitnan'), ...
        'Shat', mean(erpShat, 1, 'omitnan'), 'n', nEv);
    [~, order] = sort(serEvent);
    iEx = order(ceil(nEv / 2));
    pc = ev.peakchan_idx(iEx);
    a = max(1, ev.sample_start(iEx) - srate);
    b = min(nSamples, ev.sample_end(iEx) + srate);
    detail.example.(c.name) = struct('t', (0:b - a) / srate, 'channel', ev.peakchan{iEx}, ...
        'S', S(pc, a:b), 'Shat', Xhat(pc, a:b) - Bhat(pc, a:b), 'X', B(pc, a:b) + S(pc, a:b), ...
        'Xhat', Xhat(pc, a:b), 'SER_dB', serEvent(iEx));
end

metrics = cell2table(rows, 'VariableNames', {'class', 'stage', 'metric', 'value'});
if ~isempty(eventTables)
    detail.events = vertcat(eventTables{:});
end
end

% -------------------------------------------------------------------------
function ratio = amplitudeRatio(S, Xhat, Bhat, pc, a, b, c, srate)
% Band-limited amplitude of the recovered event over the injected one, at the peak channel.
% One second of padding keeps the filter's edge transients out of the event window.
nSamples = size(S, 2);
a2 = max(1, a - round(srate));
b2 = min(nSamples, b + round(srate));
filtered = bandpassZeroPhase([double(S(pc, a2:b2)); double(Xhat(pc, a2:b2)) - double(Bhat(pc, a2:b2))], ...
    srate, c.detectband, 4);
core = a - a2 + 1:b - a2 + 1;
switch c.kind
    case 'spindle'
        envelope = abs(hilbert(filtered'))';
        ratio = max(envelope(2, core)) / max(envelope(1, core));
    case 'slowwave'
        ratio = (max(filtered(2, core)) - min(filtered(2, core))) / (max(filtered(1, core)) - min(filtered(1, core)));
    otherwise
        error('bench:scoreSignals:kind', 'Unknown kind %s.', c.kind);
end
end

% -------------------------------------------------------------------------
function [hitIn, hitOut, rec] = detectorRecovery(S, Xhat, Bhat, B, ev, c, stageAtSample, valid, srate, tolerance)
nEv = height(ev);
[hitIn, hitOut] = deal(false(nEv, 1));
[ampIn, ampOut, freqIn, freqOut, durIn, durOut] = deal(nan(nEv, 1));
stageMask = ismember(stageAtSample, c.stages) & valid;
stageList = unique(ev.stage)';
counts = zeros(numel(stageList), 4);                 % X, Xhat, B, Bhat
channels = unique(ev.peakchan_idx)';

for ch = channels
    signals = {double(B(ch, :)) + double(S(ch, :)), double(Xhat(ch, :)), double(B(ch, :)), double(Bhat(ch, :))};
    detections = cell(1, 4);
    switch c.kind
        case 'spindle'
            [detections{1}, threshold] = bench.detectSpindles(signals{1}, srate, stageMask, 'band', c.detectband);
            for v = 2:4
                detections{v} = bench.detectSpindles(signals{v}, srate, stageMask, 'band', c.detectband, ...
                    'threshold', threshold);
            end
        case 'slowwave'
            [detections{1}, thresholds] = bench.detectSlowWaves(signals{1}, srate, stageMask, 'band', c.detectband);
            for v = 2:4
                detections{v} = bench.detectSlowWaves(signals{v}, srate, stageMask, 'band', c.detectband, ...
                    'thresholds', thresholds);
            end
        otherwise
            error('bench:scoreSignals:kind', 'Unknown kind %s.', c.kind);
    end

    here = find(ev.peakchan_idx == ch);
    if strcmp(c.kind, 'slowwave')
        %%% Trough of the injected wave at this channel, matched to detected troughs
        troughs = zeros(numel(here), 1);
        for i = 1:numel(here)
            [~, iMin] = min(S(ch, ev.sample_start(here(i)):ev.sample_end(here(i))));
            troughs(i) = ev.sample_start(here(i)) + iMin - 1;
        end
        [hIn, jIn]   = matchPeaks(troughs, detections{1}.negpeaksample, tolerance);
        [hOut, jOut] = matchPeaks(troughs, detections{2}.negpeaksample, tolerance);
    else
        [hIn, jIn]   = matchEvents(ev.sample_start(here), ev.sample_end(here), detections{1}.onset, detections{1}.offset, tolerance);
        [hOut, jOut] = matchEvents(ev.sample_start(here), ev.sample_end(here), detections{2}.onset, detections{2}.offset, tolerance);
    end
    hitIn(here) = hIn;
    hitOut(here) = hOut;
    both = hIn & hOut;
    ampIn(here(both))  = detections{1}.amplitude(jIn(both));
    ampOut(here(both)) = detections{2}.amplitude(jOut(both));
    durIn(here(both))  = detections{1}.duration(jIn(both));
    durOut(here(both)) = detections{2}.duration(jOut(both));
    if strcmp(c.kind, 'spindle')
        freqIn(here(both))  = detections{1}.freq(jIn(both));
        freqOut(here(both)) = detections{2}.freq(jOut(both));
    end
    for v = 1:4
        onsetStage = stageAtSample(detections{v}.onset);
        counts(:, v) = counts(:, v) + sum(onsetStage(:) == stageList, 1)';
    end
end

%%% Densities per minute of eligible data, averaged over the channels searched
nChannels = numel(channels);
minutes = arrayfun(@(st) nnz(stageMask & stageAtSample == st) / srate / 60, stageList);
for iStage = 0:numel(stageList)
    if iStage == 0
        sel = true(nEv, 1); use = true(1, numel(stageList)); label = 'all';
    else
        sel = ev.stage == stageList(iStage); use = (1:numel(stageList)) == iStage;
        label = bench.stageLabel(stageList(iStage));
    end
    density = sum(counts(use, :), 1) / (nChannels * sum(minutes(use)));
    both = sel & hitIn & hitOut;
    rec.(label) = struct( ...
        'density_in',               density(1), ...
        'density_out',              density(2), ...
        'density_bg_in',            density(3), ...
        'density_bg_out',           density(4), ...
        'excess_density_in',        density(1) - density(3), ...
        'excess_density_out',       density(2) - density(4), ...
        'density_retention',        (density(2) - density(4)) / (density(1) - density(3)), ...
        'injected_density',         nnz(sel) / sum(minutes(use)), ...
        'detamp_ratio_median',      median(ampOut(both) ./ ampIn(both)), ...
        'detfreq_shift_median',     median(freqOut(both) - freqIn(both)), ...
        'detduration_ratio_median', median(durOut(both) ./ durIn(both)));
end
end
