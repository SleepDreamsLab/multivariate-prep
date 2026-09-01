function GED = ged(data, opts)
% GED  Generalized eigendecomposition: spatial filters that maximise S over R.
%
%   GED = ged(EEG, Name, Value, ...)
%   GED = ged(X, 'srate', 500, Name, Value, ...)
%
%   Builds two channel covariance matrices - S from the data feature to enhance,
%   R from the data feature that serves as reference - and solves S*W = R*W*L.
%   Each column of W is a spatial filter, each diagonal element of L is the
%   multivariate S:R ratio along that filter, and the filter with the largest
%   eigenvalue is the one that maximises the contrast the caller specified.
%
%   Everything here follows Cohen (2022), NeuroImage 247:118809, "A tutorial on
%   generalized eigendecomposition for denoising, contrast enhancement, and
%   dimension reduction in multichannel electrophysiology". Section numbers in the
%   comments below refer to that paper.
%
%   The data are never modified and nothing is written to disk; the caller decides
%   what to do with the filters.
%
%   Required
%   --------
%   data      Either an EEGLAB struct (uses .data, .srate, .times, .chanlocs) or a
%             numeric array, channels x time or channels x time x trials, in which
%             case 'srate' is required. The data must already be cleaned (3.1): GED
%             cannot tell brain from artefact, it only separates two covariance
%             matrices. Reduced rank is fine (3.4) and interpolation is not needed.
%             The reference scheme does not matter, as long as the same one is used
%             to build the filters and to apply them.
%
%   Choosing the contrast (3.2) - the one decision that matters most
%   ----------------------------------------------------------------
%   contrast  Which feature pair to separate. One of:
%
%     'spectral'   (default) S = narrowband filtered data, R = the same segments
%                  broadband ('refmode' 'broadband') or filtered at the flanking
%                  frequencies ('refmode' 'neighbour'). The classic way to pull out
%                  an oscillation - alpha, a spindle band, midfrontal theta.
%                  Options: peakfreq, fwhm, refmode, neighbourdist.
%     'window'     S = a within-trial time window, R = a baseline window (task vs.
%                  pre-trial). Options: swindow, rwindow (seconds, need 'times').
%     'condition'  S = one set of trials, R = another. Options: strials, rtrials
%                  (indices into the 3rd dimension), optionally swindow/rwindow.
%     'erp'        S = covariance of the trial average, R = average of the
%                  single-trial covariances. Enhances phase-locked activity.
%                  Optionally restricted to swindow. Averaging shrinks the
%                  variance in S, so the eigenvalues land well below 1 here -
%                  read them relative to each other, or set covnorm 'trace'.
%     'data'       S and R computed from two arrays passed in directly.
%                  Options: sdata, rdata. For anything not covered above.
%     'cov'        S and R passed in ready-made. Options: S, R. Covariance
%                  cleaning, permutation testing and cross-validation are then
%                  unavailable, since those all need the individual segments.
%
%   Covariance construction (3.3)
%   -----------------------------
%   segdur        Segment length in seconds for continuous (2-D) data. Covariances
%                 are computed per segment and averaged, each segment mean-centred
%                 on its own. Default: 2. Epoched data use one segment per trial.
%                 For narrowband contrasts keep this well above one cycle of
%                 peakfreq - the check below warns if it is not.
%   covoutlierz   Drop segment covariances whose Euclidean distance to the average
%                 covariance is more than this many standard deviations away; those
%                 are multivariate outliers that would bias the average. Default:
%                 2.31 (p < .01, the value used in the paper). [] disables.
%   channorm      'none' (default) keeps the data in their own units, which is what
%                 you want for single-modality EEG. 'pooled' mean-centres each
%                 channel and divides everything by the pooled standard deviation,
%                 preserving relative channel variance. 'zscore' normalises each
%                 channel separately - only for combining modalities (6.1), as it
%                 distorts the between-channel covariances and thus the maps.
%   covnorm       'none' (default), 'trace' (scale each covariance to a mean
%                 eigenvalue of 1) or 'norm' (unit Frobenius norm). Needed when S
%                 and R live on different scales and you want to permutation test
%                 or read the eigenvalues against the null value of 1 (2.4).
%
%   Conditioning the solution (3.4, 3.10)
%   -------------------------------------
%   shrinkage     Shrinkage regularisation gamma applied to R:
%                 Rreg = R*(1-gamma) + gamma*alpha*I, alpha = mean eigenvalue of R.
%                 Default: 0.01. As little as possible, as much as necessary: at
%                 gamma = 1 the GED degenerates into a PCA on S.
%   pcacompress   Compress the data with a PCA before the GED (two-stage GED), for
%                 many channels, severely reduced-rank covariances or complex
%                 solutions. Default: false. Filters and maps come back in the
%                 original channel space either way.
%   pcadims       'rank' (default, lossless), a number >= 1 read as an explicit
%                 number of components, or a number < 1 read as a variance
%                 threshold in percent (e.g. 0.1).
%   complexaction What to do when the solution is complex-valued (3.8, a sign that
%                 the covariances are ill-conditioned or hard to separate):
%                 'warn' (default, keeps the real part), 'error' or 'ignore'.
%
%   Filters, maps and component time series (3.5, 3.6)
%   --------------------------------------------------
%   unitnorm      Scale each eigenvector to unit length. Default: true.
%   signfix       Flip each eigenvector so that the largest-magnitude channel in
%                 its map is positive, removing the arbitrary eigenvector sign.
%                 Default: true.
%   ncomps        How many component time series to compute. Default: 10.
%                 Filters, eigenvalues and maps are always returned for all.
%   applyto       Data to project through the filters, same channels in the same
%                 order. Default: the input data - note that with a spectral
%                 contrast this means the filter is built on narrowband data and
%                 applied to the broadband signal, which is the intended use.
%   zscorecomp    z-score each component time series. Default: false.
%
%   Evaluating the solution (2.3, 2.4)
%   ----------------------------------
%   nperm         Permutation iterations for the null distribution of the largest
%                 eigenvalue (maxT correction). Segment covariances are randomly
%                 reassigned to S and R, so the spatiotemporal structure of the
%                 data is preserved. Default: 0 (off). 1000 is a typical value.
%   permseed      rng seed for the permutation, for reproducibility. Default: 1.
%   cvfolds       k-fold cross-validation of the eigenvalues: fit the filters on
%                 k-1 folds of the segments and evaluate the S:R ratio on the
%                 held-out fold. Default: 0 (off).
%   plot          Draw the diagnostic figure - eigenspectrum, maps and component
%                 spectra for the first components (3.7, Fig. 4). Always look at
%                 this before trusting the top component. Default: false. Call
%                 plotged(GED) to draw it again later, e.g. for a GED struct
%                 stored by bidsfun_subcomp, which never plots on its own.
%   verbose       Print a summary. Default: true.
%
%   Output
%   ------
%   GED.filters      channels x components, the eigenvectors (spatial filters).
%   GED.evals        components x 1, the eigenvalues (the S:R ratios), sorted
%                    descending. The null-hypothesis value is 1 when S and R are on
%                    the same scale.
%   GED.maps         channels x components, the component maps (S*w). These, not
%                    the eigenvectors, are what you plot and interpret.
%   GED.comp         ncomps x time (x trials), the component time series.
%   GED.S, GED.R     The two covariance matrices, after normalisation.
%   GED.Rreg         R after shrinkage - the matrix actually handed to eig.
%   GED.covS, covR   The per-segment covariances that survived the outlier check.
%   GED.perm         .maxnull, .p (per component), .crit95, .nperm
%   GED.cv           .lambda (components x folds), .mapcorr, .nfolds
%   GED.diagnostics  ranks, condition numbers, trace ratio, dropped segments,
%                    whether the solution was complex.
%   GED.info         Everything that was asked for, plus srate, labels, chanlocs.
%   GED.apply        Function handle projecting new data: comp = GED.apply(X).
%
%   Examples
%   --------
%     % An alpha component from resting-state data, with a look at the result
%     G = ged(EEG, 'contrast', 'spectral', 'peakfreq', 10, 'fwhm', 3, 'plot', true);
%     alpha = G.comp(1, :);
%
%     % A sleep spindle filter against the neighbouring frequencies, tested
%     G = ged(EEG, 'contrast', 'spectral', 'peakfreq', 13.5, 'fwhm', 2, ...
%             'refmode', 'neighbour', 'covnorm', 'trace', 'nperm', 1000);
%     fprintf('top component p = %.3f\n', G.perm.p(1));
%
%     % Task vs. baseline on epoched data
%     G = ged(EEG, 'contrast', 'window', 'swindow', [0 0.8], 'rwindow', [-0.5 0]);
%
%     % Build the filter on one recording, apply it to another (2.3, point 3)
%     G    = ged(EEGa, 'contrast', 'spectral', 'peakfreq', 11);
%     comp = G.apply(EEGb.data);
%
% Reference:
%   Cohen MX (2022). A tutorial on generalized eigendecomposition for denoising,
%   contrast enhancement, and dimension reduction in multichannel
%   electrophysiology. NeuroImage 247:118809.

arguments
    data
    opts.srate         (1,1) double = NaN
    opts.times               double = []
    opts.chanlocs                   = []

    opts.contrast      (1,:) char {mustBeMember(opts.contrast, ...
        {'spectral', 'window', 'condition', 'erp', 'data', 'cov'})} = 'spectral'

    opts.peakfreq      (1,1) double = 10
    opts.fwhm          (1,1) double = 3
    opts.refmode       (1,:) char {mustBeMember(opts.refmode, {'broadband', 'neighbour'})} = 'broadband'
    opts.neighbourdist       double = []

    opts.swindow             double = []
    opts.rwindow             double = []
    opts.strials             double = []
    opts.rtrials             double = []
    opts.sdata                      = []
    opts.rdata                      = []
    opts.S                   double = []
    opts.R                   double = []

    opts.segdur        (1,1) double = 2
    opts.covoutlierz         double = 2.31
    opts.channorm      (1,:) char {mustBeMember(opts.channorm, {'none', 'pooled', 'zscore'})} = 'none'
    opts.covnorm       (1,:) char {mustBeMember(opts.covnorm,  {'none', 'trace', 'norm'})}    = 'none'

    opts.shrinkage     (1,1) double {mustBeInRange(opts.shrinkage, 0, 1)} = 0.01
    opts.pcacompress   (1,1) logical = false
    opts.pcadims                     = 'rank'
    opts.complexaction (1,:) char {mustBeMember(opts.complexaction, {'warn', 'error', 'ignore'})} = 'warn'

    opts.unitnorm      (1,1) logical = true
    opts.signfix       (1,1) logical = true
    opts.ncomps        (1,1) double  = 10
    opts.applyto                     = []
    opts.zscorecomp    (1,1) logical = false

    opts.nperm         (1,1) double  = 0
    opts.permseed      (1,1) double  = 1
    opts.cvfolds       (1,1) double  = 0
    opts.plot          (1,1) logical = false
    opts.verbose       (1,1) logical = true
end

%% ------------------------------------------------------------------ the data
[X, srate, times, chanlocs, labels] = unpack(data, opts);
X = normalisechannels(X, opts.channorm);
[nchan, npnts, ntrials] = size(X, 1, 2, 3);

if opts.verbose
    fprintf('\n=== GED (%s contrast) ===\n', opts.contrast);
    fprintf('%d channels, %d time points, %d trial(s), %g Hz\n', nchan, npnts, ntrials, srate);
end

%% ------------------------------------------- covariance matrices, per segment
%%% Every contrast ends up here: a stack of S covariances and a stack of R
%%% covariances, one per data segment. Keeping the segments instead of a single
%%% pooled matrix is what makes outlier cleaning, permutation testing and
%%% cross-validation possible further down (3.3, 2.4).
segsamples = max(2, round(opts.segdur * srate));
switch opts.contrast
    case 'spectral'
        checkcyclesperseg(opts, srate, ntrials, npnts, segsamples);
        covS = covstack(gaussfilter(X, srate, opts.peakfreq, opts.fwhm), segsamples);
        if strcmp(opts.refmode, 'broadband')
            covR = covstack(X, segsamples);
        else
            d = opts.neighbourdist;
            if isempty(d), d = 2 * opts.fwhm; end
            %%% Flanking bands, pooled: the reference is "the same topography at
            %%% neighbouring frequencies", so the two sides are averaged rather
            %%% than treated as separate segments.
            covR = (covstack(gaussfilter(X, srate, opts.peakfreq - d, opts.fwhm), segsamples) + ...
                    covstack(gaussfilter(X, srate, opts.peakfreq + d, opts.fwhm), segsamples)) / 2;
        end

    case 'window'
        covS = covstack(X(:, windowsamples(opts.swindow, times, 'swindow'), :), segsamples);
        covR = covstack(X(:, windowsamples(opts.rwindow, times, 'rwindow'), :), segsamples);

    case 'condition'
        if ntrials < 2
            error('ged:needTrials', 'The condition contrast needs epoched (3-D) data.');
        end
        if isempty(opts.strials) || isempty(opts.rtrials)
            error('ged:needTrialIndices', 'Give both strials and rtrials for the condition contrast.');
        end
        sIdx = ':'; rIdx = ':';
        if ~isempty(opts.swindow), sIdx = windowsamples(opts.swindow, times, 'swindow'); end
        if ~isempty(opts.rwindow), rIdx = windowsamples(opts.rwindow, times, 'rwindow'); end
        covS = covstack(X(:, sIdx, opts.strials), segsamples);
        covR = covstack(X(:, rIdx, opts.rtrials), segsamples);
        %%% Unequal trial counts let the larger condition dominate the eigenvalue,
        %%% which is the trial-count confound the paper warns about (2.3).
        if opts.verbose && numel(opts.strials) ~= numel(opts.rtrials)
            fprintf(['[note] %d S trials vs %d R trials - unequal amounts of data can bias\n' ...
                     '       the contrast towards the condition that contributes more.\n'], ...
                numel(opts.strials), numel(opts.rtrials));
        end

    case 'erp'
        if ntrials < 2
            error('ged:needTrials', 'The erp contrast needs epoched (3-D) data.');
        end
        sIdx = ':';
        if ~isempty(opts.swindow), sIdx = windowsamples(opts.swindow, times, 'swindow'); end
        Xwin = X(:, sIdx, :);
        covS = covstack(mean(Xwin, 3), size(Xwin, 2));   % covariance of the trial average
        covR = covstack(Xwin, segsamples);               % average of the single-trial covariances

    case 'data'
        if isempty(opts.sdata) || isempty(opts.rdata)
            error('ged:needData', 'Give both sdata and rdata for the data contrast.');
        end
        covS = covstack(normalisechannels(double(opts.sdata), opts.channorm), segsamples);
        covR = covstack(normalisechannels(double(opts.rdata), opts.channorm), segsamples);

    case 'cov'
        if isempty(opts.S) || isempty(opts.R)
            error('ged:needCov', 'Give both S and R for the cov contrast.');
        end
        covS = opts.S;
        covR = opts.R;
end

if size(covS, 1) ~= size(covR, 1)
    error('ged:sizeMismatch', ...
        ['S spans %d channels and R spans %d - the two covariance matrices must cover the ' ...
         'same channels, in the same order.'], size(covS, 1), size(covR, 1));
end

%%% Multivariate outlier segments (3.3). One stretch of muscle activity or a lost
%%% electrode connection can dominate an averaged covariance matrix.
[covS, droppedS] = dropoutliers(covS, opts.covoutlierz);
[covR, droppedR] = dropoutliers(covR, opts.covoutlierz);
if opts.verbose && (droppedS + droppedR) > 0
    fprintf('Dropped %d/%d S and %d/%d R segment covariances as outliers (z > %.2f)\n', ...
        droppedS, droppedS + size(covS, 3), droppedR, droppedR + size(covR, 3), opts.covoutlierz);
end

S = normalisecov(mean(covS, 3), opts.covnorm);
R = normalisecov(mean(covR, 3), opts.covnorm);

%% -------------------------------------------------------------- conditioning
diagnostics = struct( ...
    'rankS',      rank(S), ...
    'rankR',      rank(R), ...
    'condS',      cond(S), ...
    'condR',      cond(R), ...
    'traceratio', trace(S) / trace(R), ...
    'nsegS',      size(covS, 3), ...
    'nsegR',      size(covR, 3), ...
    'droppedS',   droppedS, ...
    'droppedR',   droppedR, ...
    'complex',    false, ...
    'npcadims',   NaN);

%%% Two-stage GED (3.10): compress with a PCA, separate with a GED. V carries the
%%% filters back into the channel space afterwards, so the caller never has to
%%% think about the compressed space.
V = eye(nchan);
if opts.pcacompress
    V = pcabasis((mean(covS, 3) + mean(covR, 3)) / 2, opts.pcadims);
    S = V' * S * V;
    R = V' * R * V;
    diagnostics.npcadims = size(V, 2);
    if opts.verbose
        fprintf('PCA compression: %d channels -> %d dimensions\n', nchan, size(V, 2));
    end
end

%%% Shrinkage regularisation (3.4, Eq. 15). Scaling R down by (1-gamma) keeps its
%%% trace, and with it the total energy of the eigenspectrum, unchanged.
Rreg = shrink(R, opts.shrinkage);

%% ----------------------------------------------------------------- the solve
[W, L] = eig(S, Rreg);
evals  = diag(L);

%%% Complex solutions (3.8). R\S is not symmetric, so complex conjugate pairs can
%%% appear - nearly always a sign that the covariances are ill-conditioned, or that
%%% S and R are simply not separable.
imagshare = max(abs(imag(evals))) / max(abs(real(evals)));
if imagshare > 1e-8
    diagnostics.complex = true;
    msg = sprintf(['the GED returned complex solutions (the largest imaginary part is %.2g of the ' ...
        'largest real eigenvalue). Use more data, a more separable contrast, PCA compression ' ...
        '(pcacompress) or more shrinkage.'], imagshare);
    switch opts.complexaction
        case 'error', error('ged:complexSolution', '%s', msg);
        case 'warn',  warning('ged:complexSolution', '%s', msg);
    end
end
evals = real(evals);
W     = real(W);

[evals, order] = sort(evals, 'descend');
W = W(:, order);
W = V * W;                       % back into the channel space (identity unless compressed)

if opts.unitnorm
    W = W ./ vecnorm(W);         % eigenvectors carry an arbitrary scale (3.6)
end

%%% Component maps (3.6). The eigenvectors themselves are not physiologically
%%% interpretable - they also suppress irrelevant channels - so what gets plotted
%%% is w'S: the source projecting outwards onto the electrodes.
Sfull = normalisecov(mean(covS, 3), opts.covnorm);
maps  = Sfull * W;

%%% The eigenvector sign is arbitrary (3.5); fix it on the strongest channel of the
%%% map so that maps and ERPs can be averaged across recordings.
if opts.signfix
    [~, peak] = max(abs(maps), [], 1);
    flip = sign(maps(sub2ind(size(maps), peak, 1:size(maps, 2))));
    flip(flip == 0) = 1;
    W    = W    .* flip;
    maps = maps .* flip;
end

%% ---------------------------------------------------- component time series
ncomps = min(opts.ncomps, size(W, 2));
if isempty(opts.applyto)
    projdata = X;
else
    projdata = normalisechannels(double(opts.applyto), opts.channorm);
end
if size(projdata, 1) ~= nchan
    error('ged:applyChannels', ...
        'applyto has %d channels but the filters were built on %d.', size(projdata, 1), nchan);
end
comp = project(projdata, W(:, 1:ncomps));
if opts.zscorecomp
    comp = (comp - mean(comp, 2)) ./ std(comp, 0, 2);
end

%% ------------------------------------------------------ inferential statistics
perm = struct('maxnull', [], 'p', [], 'crit95', NaN, 'nperm', 0);
if opts.nperm > 0
    if strcmp(opts.contrast, 'cov')
        warning('ged:noPermutation', ...
            'Permutation testing needs the segment covariances, which the cov contrast does not provide.');
    else
        if strcmp(opts.covnorm, 'none') && (diagnostics.traceratio < 0.5 || diagnostics.traceratio > 2)
            warning('ged:scaleMismatch', ...
                ['S and R differ in scale by a factor of %.2g, so the observed eigenvalues and the ' ...
                 'permuted ones are not centred on the same value. Set covnorm to ''trace'' (2.4).'], ...
                diagnostics.traceratio);
        end
        perm = permutetest(covS, covR, V, opts);
        if opts.verbose
            fprintf('Permutation (%d iterations): critical lambda = %.3f, %d/%d components at p < .05\n', ...
                opts.nperm, perm.crit95, nnz(perm.p < 0.05), numel(perm.p));
        end
    end
end

cv = struct('lambda', [], 'mapcorr', [], 'nfolds', 0);
if opts.cvfolds > 1
    if strcmp(opts.contrast, 'cov')
        warning('ged:noCrossvalidation', ...
            'Cross-validation needs the segment covariances, which the cov contrast does not provide.');
    else
        cv = crossvalidate(covS, covR, V, maps, opts, ncomps);
        if opts.verbose
            fprintf('Cross-validation (%d folds): held-out lambda of component 1 = %.3f (SD %.3f)\n', ...
                cv.nfolds, mean(cv.lambda(1, :)), std(cv.lambda(1, :)));
        end
    end
end

%% ---------------------------------------------------------------- the output
GED             = struct();
GED.filters     = W;
GED.evals       = evals;
GED.maps        = maps;
GED.comp        = comp;
GED.S           = S;
GED.R           = R;
GED.Rreg        = Rreg;
GED.covS        = covS;
GED.covR        = covR;
GED.perm        = perm;
GED.cv          = cv;
GED.diagnostics = diagnostics;
GED.info          = opts;
GED.info.srate    = srate;
GED.info.labels   = labels;
GED.info.chanlocs = chanlocs;
GED.apply       = @(Y) project(normalisechannels(double(Y), opts.channorm), W(:, 1:ncomps));

if opts.verbose
    fprintf('Rank S/R: %d/%d of %d. Condition number of R: %.3g. trace(S)/trace(R) = %.3g\n', ...
        diagnostics.rankS, diagnostics.rankR, nchan, diagnostics.condR, diagnostics.traceratio);
    fprintf('Top eigenvalues: %s\n', num2str(evals(1:min(5, end))', '%.3f  '));
    fprintf(['[remember] the largest eigenvalue separates S from R best mathematically, not\n' ...
             '           necessarily physiologically - look at the maps before committing (3.7).\n']);
end

if opts.plot
    plotged(GED, 'ncomps', ncomps);
end
end

%% ========================================================================= %%
%  Local functions
%% ========================================================================= %%

function [X, srate, times, chanlocs, labels] = unpack(data, opts)
% Accept an EEGLAB struct or a raw array and return everything in one shape.

if isstruct(data) && isfield(data, 'data')
    X        = double(data.data);
    srate    = double(data.srate);
    chanlocs = data.chanlocs;
    times    = [];
    if isfield(data, 'times') && ~isempty(data.times)
        times = double(data.times(:))' / 1000;   % EEGLAB keeps milliseconds
    end
elseif isnumeric(data)
    X        = double(data);
    srate    = opts.srate;
    times    = opts.times(:)';
    chanlocs = opts.chanlocs;
    if isnan(srate)
        error('ged:noSrate', 'Give the sampling rate with ''srate'' when passing a numeric array.');
    end
else
    error('ged:badInput', 'data must be an EEGLAB struct or a numeric array.');
end

%%% Explicit name-value arguments win over whatever the struct carried.
if ~isempty(opts.times),    times    = opts.times(:)'; end
if ~isempty(opts.chanlocs), chanlocs = opts.chanlocs;  end
if ~isnan(opts.srate),      srate    = opts.srate;     end

if isempty(times)
    times = (0:size(X, 2) - 1) / srate;
end
labels = {};
if ~isempty(chanlocs) && isfield(chanlocs, 'labels')
    labels = {chanlocs.labels};
end
if any(~isfinite(X(:)))
    error('ged:nonFinite', 'The data contain NaN or Inf, which would make the covariances meaningless.');
end
end

% -------------------------------------------------------------------------
function X = normalisechannels(X, mode)
% Channel scaling (3.3). 'pooled' keeps the relative channel variances, which the
% component maps depend on; 'zscore' does not, and is for multimodal data only.

switch mode
    case 'none'
        % covariance matrices keep the units of the data, which is usually right
    case 'pooled'
        X = X - mean(X, 2);
        X = X / std(X(:));
    case 'zscore'
        X = (X - mean(X, 2)) ./ std(X, 0, 2);
end
end

% -------------------------------------------------------------------------
function C = covstack(X, segsamples)
% One covariance matrix per segment, each mean-centred on its own (3.3). Epoched
% data give one covariance per trial; continuous data are cut into segments of
% segsamples points, with a trailing remainder shorter than one segment ignored.

[nchan, npnts, ntrials] = size(X, 1, 2, 3);

if ntrials > 1
    C = zeros(nchan, nchan, ntrials);
    for t = 1:ntrials
        C(:, :, t) = onecov(X(:, :, t));
    end
    return
end

if segsamples >= npnts
    segsamples = npnts;
end
nseg = max(1, floor(npnts / segsamples));
C    = zeros(nchan, nchan, nseg);
for s = 1:nseg
    C(:, :, s) = onecov(X(:, (s - 1) * segsamples + (1:segsamples)));
end
end

function C = onecov(Y)
Y = Y - mean(Y, 2);                  % mean offsets would steer the solution (2.1)
C = (Y * Y') / (size(Y, 2) - 1);
C = (C + C') / 2;                    % symmetric to machine precision
end

% -------------------------------------------------------------------------
function [C, ndropped] = dropoutliers(C, zthresh)
% Reject segment covariances that sit far away from the average covariance (3.3).

ndropped = 0;
if isempty(zthresh) || size(C, 3) < 4, return, end

mu = mean(C, 3);
d  = zeros(size(C, 3), 1);
for s = 1:size(C, 3)
    d(s) = norm(C(:, :, s) - mu, 'fro');
end
if std(d) == 0, return, end

keep = (d - mean(d)) / std(d) <= zthresh;
if ~any(keep)
    warning('ged:allOutliers', 'Every segment covariance looked like an outlier; keeping all of them.');
    return
end
ndropped = nnz(~keep);
C = C(:, :, keep);
end

% -------------------------------------------------------------------------
function C = normalisecov(C, mode)
% Covariance normalisation (3.3). Only needed when S and R are on different scales
% and the absolute eigenvalues have to mean something.

switch mode
    case 'none'
    case 'trace', C = C / (trace(C) / size(C, 1));   % mean eigenvalue becomes 1
    case 'norm',  C = C / norm(C, 'fro');
end
end

% -------------------------------------------------------------------------
function Rreg = shrink(R, gamma)
% Shrinkage regularisation, Eq. 15: Rreg = R(1-gamma) + gamma*alpha*I, with alpha
% the average eigenvalue of R (= trace(R)/M). The trace is preserved.

if gamma == 0
    Rreg = R;
    return
end
alpha = trace(R) / size(R, 1);
Rreg  = R * (1 - gamma) + gamma * alpha * eye(size(R, 1));
end

% -------------------------------------------------------------------------
function V = pcabasis(C, dims)
% First stage of a two-stage GED (3.10): the PCA basis of the data covariance.

[V, D] = eig((C + C') / 2);
[d, o] = sort(diag(D), 'descend');
V      = V(:, o);
d      = max(d, 0);

if ischar(dims) || isstring(dims)
    n = rank(C);                                  % lossless: keep the whole rank
elseif dims >= 1
    n = min(round(dims), size(V, 2));             % an explicit number of components
else
    n = nnz(100 * d / sum(d) > dims);             % a variance threshold, in percent
end
V = V(:, 1:max(n, 1));
end

% -------------------------------------------------------------------------
function Y = gaussfilter(X, srate, peakfreq, fwhm)
% Narrowband filter by a frequency-domain Gaussian, specified by its full width at
% half maximum in Hz: no toolbox, no phase distortion, no filter order to pick.
% Applied to the whole time series rather than to the covariance window alone, so
% that edge artefacts stay out of the covariance matrices (3.3).

if peakfreq <= 0
    error('ged:badFrequency', 'peakfreq must be positive (asked for %g Hz).', peakfreq);
end
if peakfreq > srate / 2
    error('ged:aboveNyquist', 'peakfreq (%g Hz) is above the Nyquist frequency (%g Hz).', ...
        peakfreq, srate / 2);
end
npnts = size(X, 2);
hz    = linspace(0, srate, npnts);
s     = fwhm * (2 * pi - 1) / (4 * pi);        % Gaussian width from the FWHM
gauss = exp(-0.5 * ((hz - peakfreq) / s).^2);
gauss = gauss / max(gauss);                    % gain-normalised: no amplitude bias
Y     = 2 * real(ifft(fft(X, [], 2) .* gauss, [], 2));
end

% -------------------------------------------------------------------------
function idx = windowsamples(window, times, name)
% Turn a [start stop] window in seconds into sample indices.

if isempty(window)
    error('ged:noWindow', 'This contrast needs %s, in seconds.', name);
end
if isempty(times)
    error('ged:noTimes', 'No time vector available - pass ''times'' (in seconds) to use %s.', name);
end
idx = find(times >= window(1) & times <= window(2));
if numel(idx) < 2
    error('ged:emptyWindow', '%s = [%g %g] s selects %d sample(s).', ...
        name, window(1), window(2), numel(idx));
end
end

% -------------------------------------------------------------------------
function checkcyclesperseg(opts, srate, ntrials, npnts, segsamples)
% A covariance of narrowband data needs at least one cycle per segment (3.3).

if ntrials > 1
    seglen = npnts / srate;
else
    seglen = segsamples / srate;
end
if seglen * opts.peakfreq < 1
    warning('ged:shortSegments', ...
        ['covariance segments are %.0f ms long, less than one cycle of %g Hz. Use longer segments ' ...
         '(segdur) or the covariance matrices will be unstable.'], seglen * 1000, opts.peakfreq);
end
end

% -------------------------------------------------------------------------
function comp = project(X, W)
% Apply the spatial filters: comp = w'X, for 2-D and 3-D data alike (3.6).

X = X - mean(X, 2);
if ismatrix(X)
    comp = W' * X;
else
    [~, npnts, ntrials] = size(X);
    comp = reshape(W' * reshape(X, size(X, 1), []), size(W, 2), npnts, ntrials);
end
end

% -------------------------------------------------------------------------
function perm = permutetest(covS, covR, V, opts)
% Null distribution of the largest eigenvalue (2.4). The segment covariances are
% randomly reassigned to S and R, which keeps the spatiotemporal structure of real
% data in the null - covariances of random numbers would give a far too liberal
% threshold. Taking the largest eigenvalue of every iteration corrects for the
% comparisons across components (maxT / extreme-value correction).

rng(opts.permseed);
pool = cat(3, covS, covR);
nS   = size(covS, 3);
nAll = size(pool, 3);

maxnull = zeros(opts.nperm, 1);
for p = 1:opts.nperm
    shuffled = randperm(nAll);
    Sp = V' * normalisecov(mean(pool(:, :, shuffled(1:nS)),       3), opts.covnorm) * V;
    Rp = V' * normalisecov(mean(pool(:, :, shuffled(nS + 1:end)), 3), opts.covnorm) * V;
    maxnull(p) = max(real(eig(Sp, shrink(Rp, opts.shrinkage))));
end

%%% The observed eigenvalues are recomputed here in exactly the same way, so that
%%% the comparison is like for like even when the caller asked for sign fixing,
%%% unit norming or anything else that touches the returned solution.
observed = sort(real(eig( ...
    V' * normalisecov(mean(covS, 3), opts.covnorm) * V, ...
    shrink(V' * normalisecov(mean(covR, 3), opts.covnorm) * V, opts.shrinkage))), 'descend');

perm         = struct();
perm.maxnull = maxnull;
perm.p       = arrayfun(@(lambda) (1 + nnz(maxnull >= lambda)) / (1 + opts.nperm), observed);
perm.crit95  = percentile(maxnull, 95);
perm.nperm   = opts.nperm;
end

% -------------------------------------------------------------------------
function cv = crossvalidate(covS, covR, V, maps, opts, ncomps)
% k-fold cross-validation (2.3). The filters are fitted on the training segments
% and their S:R ratio is evaluated on the held-out segments, which is the ratio
% the filter would reach on data it has never seen. Components are matched across
% folds by rank, so read the map correlations alongside the eigenvalues: a low
% correlation means that fold found a different component, not a worse one.

k       = min(opts.cvfolds, min(size(covS, 3), size(covR, 3)));
foldS   = foldassignment(size(covS, 3), k);
foldR   = foldassignment(size(covR, 3), k);
lambda  = nan(ncomps, k);
mapcorr = nan(ncomps, k);

for f = 1:k
    Str = normalisecov(mean(covS(:, :, foldS ~= f), 3), opts.covnorm);
    Rtr = normalisecov(mean(covR(:, :, foldR ~= f), 3), opts.covnorm);
    Ste = normalisecov(mean(covS(:, :, foldS == f), 3), opts.covnorm);
    Rte = normalisecov(mean(covR(:, :, foldR == f), 3), opts.covnorm);

    [Wf, Lf] = eig(V' * Str * V, shrink(V' * Rtr * V, opts.shrinkage));
    [~, o]   = sort(real(diag(Lf)), 'descend');
    Wf       = V * real(Wf(:, o(1:ncomps)));
    Wf       = Wf ./ vecnorm(Wf);

    lambda(:, f) = diag(Wf' * Ste * Wf) ./ diag(Wf' * Rte * Wf);
    for c = 1:ncomps
        mapcorr(c, f) = abs(spatialcorr(Str * Wf(:, c), maps(:, c)));
    end
end

cv = struct('lambda', lambda, 'mapcorr', mapcorr, 'nfolds', k);
end

function fold = foldassignment(n, k)
% Spread the segments over k folds in random order.

fold = mod(randperm(n), k) + 1;
end

function r = spatialcorr(a, b)
% Pearson correlation of two component maps, without the Statistics toolbox.

a = a(:) - mean(a(:));
b = b(:) - mean(b(:));
r = (a' * b) / (norm(a) * norm(b));
end

% -------------------------------------------------------------------------
function v = percentile(x, p)
% Linear-interpolation percentile, so that no toolbox is needed.

x = sort(x(:));
n = numel(x);
if n == 1, v = x; return, end
pos = (p / 100) * n + 0.5;
lo  = min(max(floor(pos), 1), n);
hi  = min(max(ceil(pos),  1), n);
v   = x(lo) + (pos - floor(pos)) * (x(hi) - x(lo));
end

%%% The diagnostic figure itself lives in plotged.m, so it can be called again
%%% later on any stored GED struct (ged.m calls it above when 'plot' is true).
