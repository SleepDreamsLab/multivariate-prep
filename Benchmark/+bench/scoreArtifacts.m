function [metrics, detail] = scoreArtifacts(A, Xhat, Bhat, B, events, type, cfg, stageEpoch, srate, epochLength, psdB, psdBhat, opts)
% SCOREARTIFACTS  Under-cleaning metrics for injected artifacts: ARR, RRMSE, CC, PSD ratios.
%
%   [metrics, detail] = bench.scoreArtifacts(A, Xhat, Bhat, B, events, type, cfg, stageEpoch, ...
%       srate, epochLength, psdB, psdBhat)
%
%   All data are [nChan x nSamples], average-referenced:
%     A     injected artifact (ground truth), already multiplied by the condition's scale
%     B     background;  Bhat = GEDAI(B)
%     Xhat  GEDAI(B + A)
%   psdB, psdBhat: struct with .P [nStages x nFreq x nChan], .f, .stages (bench.psdEpochs),
%   computed once per GEDAI run and shared across artifact conditions.
%
%   The residual is Xhat - Bhat (paired design, see bench.scoreSignals): what the injected
%   artifact left behind, plus any change in how the background was cleaned because of it.
%
%   Windows: EOG and EMG use each event +/- eventpad; ECG uses every sample of cfg.stages.
%
%   Per artifact type and stage (and pooled as stage 'all')
%   -------------------------------------------------------
%   ARR_dB             10 log10( sum A^2 / sum (Xhat - Bhat)^2 ). Artifact-to-residue ratio
%                      of Somers et al. (2018, J Neural Eng 15:036007), whose A - Ahat
%                      equals Xhat - B; here B is replaced by Bhat. Higher is better; 0 dB =
%                      nothing removed (or as much damage done as removed).
%   ARR_raw_dB         Somers' unpaired form with B. Meaningful on a synthetic background.
%   residual_gain      <Xhat - Bhat, A> / <A, A>: fraction of the artifact left in the output
%   RRMSE_t_in         ||A|| / ||B||: temporal relative RMSE before cleaning
%   RRMSE_t_out        ||Xhat - Bhat|| / ||Bhat||: after cleaning (Zhang et al. 2021,
%                      J Neural Eng 18:056057, with GEDAI(B) as the clean reference)
%   CC_in, CC_out      median over channels of corr(B + A, B) and corr(Xhat, Bhat)
%
%   Per stage only (PSDs are per stage; 'all' averages stages weighted by epoch count)
%   RRMSE_s_in, RRMSE_s_out   spectral RRMSE over 0.5 Hz to 0.4 x srate, median over
%                      channels: ||P(B + A) - P(B)|| / ||P(B)|| and ||P(Xhat) - P(Bhat)|| / ||P(Bhat)||
%   psd_contamination_dB   10 log10( band power of B + A / band power of B ): how much the
%                      artifact added in cfg.psdbands, power pooled over channels
%   psd_residual_dB    10 log10( band power of Xhat / band power of Bhat ): what is left
%   psd_residual_excess    (P(Xhat) - P(Bhat)) / (P(B + A) - P(B)) in-band: 0 = removed,
%                      1 = untouched, < 0 = more removed than the artifact added
%
%   Name-value
%   ----------
%   eventpad     s                                                          (0.25)
%   validmask    [1 x nSamples] logical, samples to score                   (all)
%   welchwindow  s, must match psdB/psdBhat                                  (4)
%
%   Outputs
%   -------
%   metrics   table: class, stage, metric, value (class = artifact type)
%   detail    .spectra(k): stage, f, ratioIn_dB, ratioOut_dB (median over channels)
%             .topo: channel vectors for stage 'all' - ARR_dB, psd_residual_dB,
%                    psd_contamination_dB

arguments
    A single
    Xhat single
    Bhat single
    B single
    events table
    type {mustBeMember(type, {'eog', 'emg', 'ecg'})}
    cfg struct
    stageEpoch (1,:) double
    srate (1,1) double
    epochLength (1,1) double
    psdB struct
    psdBhat struct
    opts.eventpad (1,1) double = 0.25
    opts.validmask (1,:) logical = logical([])
    opts.welchwindow (1,1) double = 4
end

nSamples      = size(A, 2);
epochSamples  = round(epochLength * srate);
stageAtSample = repelem(stageEpoch, epochSamples);
valid         = opts.validmask;
if isempty(valid), valid = true(1, nSamples); end

if strcmp(type, 'ecg')
    window = ismember(stageAtSample, cfg.stages);
else
    window = eventWindowMask(nSamples, events.sample_start, events.sample_end, round(opts.eventpad * srate));
end
window = window & valid;

stages = cfg.stages(ismember(cfg.stages, stageEpoch));
nEpochs = arrayfun(@(st) nnz(stageEpoch == st), stages);

%%% PSDs of the contaminated input and the cleaned output
PX    = bench.psdEpochs({B, A}, [1 1], stageEpoch, stages, srate, epochLength, opts.welchwindow);
[PXhat, f] = bench.psdEpochs({Xhat}, 1, stageEpoch, stages, srate, epochLength, opts.welchwindow);
[~, iB] = ismember(stages, psdB.stages);
[~, iBhat] = ismember(stages, psdBhat.stages);
PB    = psdB.P(iB, :, :);
PBhat = psdBhat.P(iBhat, :, :);

bands  = cfg.psdbands;
inBand = any(f >= bands(:, 1)' & f <= bands(:, 2)', 2);
broad  = f >= 0.5 & f <= 0.4 * srate;

rows = {};
detail = struct('spectra', struct('stage', {}, 'f', {}, 'ratioIn_dB', {}, 'ratioOut_dB', {}), 'topo', struct());
for iStage = 0:numel(stages)
    if iStage == 0
        label = 'all';
        mask = window;
        weights = nEpochs / sum(nEpochs);
    else
        label = bench.stageLabel(stages(iStage));
        mask = window & stageAtSample == stages(iStage);
        weights = double((1:numel(stages)) == iStage);
    end
    if ~any(mask), continue; end
    add = @(name, value) {type, label, name, value};

    sums = windowSums(A, Xhat, Bhat, B, mask, epochSamples);
    rows = [rows; ...
        add('ARR_dB',        10 * log10(sum(sums.a2) / sum(sums.d2))); ...
        add('ARR_raw_dB',    10 * log10(sum(sums.a2) / sum(sums.dRaw2))); ...
        add('residual_gain', sum(sums.da) / sum(sums.a2)); ...
        add('RRMSE_t_in',    sqrt(sum(sums.a2) / sum(sums.b2))); ...
        add('RRMSE_t_out',   sqrt(sum(sums.d2) / sum(sums.bhat2))); ...
        add('CC_in',         median(sums.ccIn, 'omitnan')); ...
        add('CC_out',        median(sums.ccOut, 'omitnan'))]; %#ok<AGROW>

    %%% Spectral metrics: stage PSDs combined with epoch-count weights
    pX    = squeeze(sum(weights(:) .* PX, 1));      % [nFreq x nChan]
    pXhat = squeeze(sum(weights(:) .* PXhat, 1));
    pB    = squeeze(sum(weights(:) .* PB, 1));
    pBhat = squeeze(sum(weights(:) .* PBhat, 1));
    rrmseIn  = sqrt(sum((pX(broad, :) - pB(broad, :)).^2, 1)) ./ sqrt(sum(pB(broad, :).^2, 1));
    rrmseOut = sqrt(sum((pXhat(broad, :) - pBhat(broad, :)).^2, 1)) ./ sqrt(sum(pBhat(broad, :).^2, 1));
    bandX = sum(pX(inBand, :), 'all');  bandB = sum(pB(inBand, :), 'all');
    bandXhat = sum(pXhat(inBand, :), 'all');  bandBhat = sum(pBhat(inBand, :), 'all');
    rows = [rows; ...
        add('RRMSE_s_in',           median(rrmseIn)); ...
        add('RRMSE_s_out',          median(rrmseOut)); ...
        add('psd_contamination_dB', 10 * log10(bandX / bandB)); ...
        add('psd_residual_dB',      10 * log10(bandXhat / bandBhat)); ...
        add('psd_residual_excess',  (bandXhat - bandBhat) / (bandX - bandB))]; %#ok<AGROW>

    detail.spectra(end+1) = struct('stage', label, 'f', f, ...
        'ratioIn_dB', median(10 * log10(pX ./ pB), 2), 'ratioOut_dB', median(10 * log10(pXhat ./ pBhat), 2));
    if iStage == 0
        detail.topo = struct('ARR_dB', 10 * log10(sums.a2 ./ sums.d2), ...
            'psd_residual_dB', 10 * log10(sum(pXhat(inBand, :), 1) ./ sum(pBhat(inBand, :), 1))', ...
            'psd_contamination_dB', 10 * log10(sum(pX(inBand, :), 1) ./ sum(pB(inBand, :), 1))');
    end
end
metrics = cell2table(rows, 'VariableNames', {'class', 'stage', 'metric', 'value'});
end

% -------------------------------------------------------------------------
function sums = windowSums(A, Xhat, Bhat, B, mask, epochSamples)
% Per-channel sums over the masked samples, accumulated epoch by epoch so no copy of the
% full recording is made. Correlations come from the running first and second moments.
nChan = size(A, 1);
z = zeros(nChan, 1);
sums = struct('a2', z, 'd2', z, 'dRaw2', z, 'da', z, 'b2', z, 'bhat2', z);
moments = struct('n', 0, 'xIn', z, 'yIn', z, 'xxIn', z, 'yyIn', z, 'xyIn', z, ...
    'xOut', z, 'yOut', z, 'xxOut', z, 'yyOut', z, 'xyOut', z);
for iEpoch = 1:size(A, 2) / epochSamples
    cols = (iEpoch - 1) * epochSamples + 1:iEpoch * epochSamples;
    cols = cols(mask(cols));
    if isempty(cols), continue; end
    a = double(A(:, cols));  b = double(B(:, cols));
    xh = double(Xhat(:, cols));  bh = double(Bhat(:, cols));
    d = xh - bh;
    sums.a2    = sums.a2 + sum(a.^2, 2);
    sums.d2    = sums.d2 + sum(d.^2, 2);
    sums.dRaw2 = sums.dRaw2 + sum((xh - b).^2, 2);
    sums.da    = sums.da + sum(d .* a, 2);
    sums.b2    = sums.b2 + sum(b.^2, 2);
    sums.bhat2 = sums.bhat2 + sum(bh.^2, 2);

    x = b + a;
    moments.n     = moments.n + numel(cols);
    moments.xIn   = moments.xIn + sum(x, 2);    moments.yIn   = moments.yIn + sum(b, 2);
    moments.xxIn  = moments.xxIn + sum(x.^2, 2); moments.yyIn  = moments.yyIn + sum(b.^2, 2);
    moments.xyIn  = moments.xyIn + sum(x .* b, 2);
    moments.xOut  = moments.xOut + sum(xh, 2);   moments.yOut  = moments.yOut + sum(bh, 2);
    moments.xxOut = moments.xxOut + sum(xh.^2, 2); moments.yyOut = moments.yyOut + sum(bh.^2, 2);
    moments.xyOut = moments.xyOut + sum(xh .* bh, 2);
end
sums.ccIn  = pearson(moments.n, moments.xIn, moments.yIn, moments.xxIn, moments.yyIn, moments.xyIn);
sums.ccOut = pearson(moments.n, moments.xOut, moments.yOut, moments.xxOut, moments.yyOut, moments.xyOut);
end

% -------------------------------------------------------------------------
function r = pearson(n, sx, sy, sxx, syy, sxy)
r = (n * sxy - sx .* sy) ./ sqrt((n * sxx - sx.^2) .* (n * syy - sy.^2));
end
