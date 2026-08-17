function [EEG, KeepTime] = run_filter(EEG, opts)
% RUN_FILTER  Resample, remove DC offset, and suppress line noise on an EEG struct.
%
% USAGE:
%   EEG = run.run_filter(EEG)
%   [EEG, KeepTime] = run.run_filter(EEG, removeDC=true, zapline=true)
%
% INPUTS:
%   EEG   — EEGLAB EEG struct (already imported by the caller)
%
% OPTIONAL NAME-VALUE:
%   noteegchannels  channel indices to drop before processing  (default 257:300)
%   targetsrate     resample to this rate in Hz; 0 = skip      (default 125)
%   removeDC        apply DC-removal filter                    (default true)
%   zapline        apply Zapline-plus line-noise removal      (default true)
%   noiseCompDetectSigma  sigma of the per-chunk outlier detector that picks how many
%                   DSS components to remove. Higher = fewer components = gentler.
%                   5 is the upper bound of Zapline-plus' own adaptive range (default 5)
%   adaptiveSigma   let Zapline-plus tune noiseCompDetectSigma by reprocessing the whole
%                   recording once per 0.25 step. Off: on these data it always walks to
%                   the 5.0 ceiling, so 9 passes are computed and 8 discarded (default false)
%   KeepTime        struct of prior timings to merge into the output (default [])
%   cleanline       apply CleanLine after Zapline               (default true)
%   badchannels     detect and remove bad channels with clean_channels, after DC
%                   removal and before Zapline. EEG.chanlocs must carry X/Y/Z. The
%                   full montage stays available in EEG.urchanlocs, so the removed
%                   channels can be interpolated back downstream  (default false)
%   badchannelsonly stop after bad channel detection: no channel is removed, and Zapline
%                   and CleanLine are skipped. Lets bidsfun_detect_badchans reach the
%                   detection step through this same code path, so the data detection
%                   sees is defined in exactly one place              (default false)
%   badchanstride   evaluate every Nth window in the correlation criterion. The criterion
%                   is a proportion of windows, so this trades precision, not
%                   correctness: at 2 the proportion still rests on ~3300 windows of a
%                   9-h night. Set 1 for the original behaviour   (default 2)
%   flatthreshold   peak-to-peak, in data units (uV), below which a window counts as
%                   flat for gedai.detectFlatChannels. A channel flat for more than
%                   half the recording is removed, the same duration rule the
%                   correlation criterion uses                    (default 0.5)
%   badchanavgref   average-reference the data for the bad channel detection only, and
%                   undo it afterwards. Off, a single-electrode reference makes both
%                   clean_channels criteria misfire on the ring of channels next to
%                   that electrode - see the note at the detection block (default true)
%   badchanfile     .mat cache path for the bad channel results; required when
%                   badchannels is true
%   refresh         recompute the bad channel cache even if it exists (default false)
%   JsonFile        full path for a JSON sidecar with processing timings;
%                   '' = skip                                  (default '')
%   zeropatchseconds   cut out all-zero patches (amplifier crash padding) longer
%                   than this many seconds before filtering; 0 = skip (default 5)
%   restorezeropatches  put those patches back before returning. Set false if the
%                   caller does more filtering and restores itself with
%                   run.restore_zero_patches before saving        (default true)
%
% Methods secton:
%
% Continuous EEG (256 channels, 250 Hz) was high-pass filtered to remove DC offset and then 
% cleaned of power-line artifacts using Zapline-plus (Klug & Kloosterman, 2022), an extension 
% of the Zapline algorithm (de Cheveigné, 2020). Zapline separates the data into a 
% line-frequency-dominated and a residual component using a comb filter, isolates the artifactual 
% subspace via DSS applied to the line frequency and its harmonics, and removes it by spatial 
% projection, thereby avoiding the spectral distortion introduced by notch filtering. The target 
% frequency was fixed at 50 Hz; at a 250 Hz sampling rate the fundamental and its second harmonic 
% (100 Hz) were removed jointly. To accommodate non-stationarity of the line artifact across the 
% [8-h] recordings, data were processed in fixed 300-s chunks, with the individual noise peak 
% re-estimated within each chunk (search window 50 ± 3 Hz) and the number of removed components 
% determined adaptively per chunk by iterative outlier detection on the DSS component scores
% (minimum 1 component). The outlier-detection threshold was held fixed at sigma = 5, the upper
% bound of the algorithm's adaptive range and therefore its most conservative setting; the
% iterative threshold adaptation was disabled because it converged to this bound on every pilot
% recording, so a fixed threshold yields identical output while also equalising cleaning
% strength across recordings. Residual artifacts at 50 and 100 Hz were subsequently attenuated using
% a custom parallelised reimplementation CleanLine (Mullen, 2012; EEGLAB), leveraging 
% multi-taper sinusoidal regression, .
% Within 4-s windows advanced in 2-s steps, the Thomson F-test for a deterministic sinusoid 
% was evaluated at 50 and 100 Hz using 7 Slepian tapers (time–bandwidth product 4, 
% resolution bandwidth 2 Hz) and an FFT zero-padded to 4096 points (0.061 Hz bin spacing). 
% Sinusoids reaching p < 0.01 were fitted by least squares and subtracted; target frequencies 
% not reaching significance in a given window and channel were left unmodified. Fitted 
% components were crossfaded sigmoidally across the 2-s overlap between adjacent windows.
% Unlike the default CleanLine implementation, which subtracts an estimated sinusoid at each
% target frequency regardless of significance, subtraction here was conditional on the F-test.
%
% Bad channels:
%
% Bad channels were identified and removed before line-noise correction. For detection only, the
% recording was re-referenced to the common average; the data carried forward retained the original
% vertex (Cz) reference. This is necessary because under a single-electrode reference the electrodes
% immediately surrounding the reference carry a near-zero potential difference, and both criteria
% below are ratios in which that difference forms the denominator, so those channels are
% systematically misclassified (here, the median reconstruction correlation of the ring of
% electrodes adjacent to the reference fell from 0.99 to 0.73 when detection was run on
% vertex-referenced data, while every other channel was unaffected). Three criteria were evaluated
% over consecutive 5-s windows, and a channel was rejected when any of them held for more than 50%
% of the recording. (i) Flat: the signal varied by no more than 0.5 uV peak-to-peak within the
% window. This criterion is required because a dead electrode is invisible to the two that follow -
% both are ratios that evaluate to 0/0 for a constant signal, and the resulting NaN fails every
% threshold comparison. (ii) Line noise: the ratio of the >45 Hz to the <45 Hz robust amplitude
% (median absolute deviation), the two bands separated by a 100th-order least-squares FIR lowpass
% with a 45-50 Hz transition, expressed as a robust z-score across channels and thresholded at
% z = 4. Because this criterion keys on power-line contamination, detection preceded Zapline-plus,
% while the 50 and 100 Hz components were still present. (iii) Reconstruction correlation: 25 random
% subsets of 25% of the channels were each used to predict all channels by spherical-spline
% interpolation (Perrin et al., 1989); the element-wise median of the 25 predictions was taken as
% the consensus estimate (RANSAC; Fischler & Bolles, 1981), and the Pearson correlation between each
% channel's <45 Hz signal and its consensus prediction computed per window, with a threshold of
% r = 0.7. Windows were evaluated at a stride of two, leaving the proportion estimated from
% approximately 3,300 windows per recording. Channels flagged by any criterion were removed before
% Zapline-plus so that they could not contribute to the estimation of its spatial filters; the full
% montage was retained alongside the data so that removed channels could be restored by
% spherical-spline interpolation wherever whole-head output was required. Detection used
% clean_channels from the clean_rawdata plugin for EEGLAB (Kothe & Makeig, 2013; Delorme & Makeig,
% 2004), locally modified to return the per-window correlations and to accept a window stride; the
% flat-line criterion was implemented separately.

arguments
    EEG                  struct
    opts.noteegchannels  (1,:) double  = 257:300
    opts.targetsrate     (1,1) double  = 0
    opts.removeDC        (1,1) logical = true
    opts.zapline         (1,1) logical = true
    opts.KeepTime        struct        = struct()
    opts.cleanline       (1,1) logical = true
    opts.JsonFile        (1,1) string  = ""
    opts.noisefreqs                    = 50
    opts.adaptiveNremove (1,1) logical = true
    opts.fixedNremove    (1,1) double  = 1
    opts.chunkLength     (1,1) double  = 300
    opts.noiseCompDetectSigma (1,1) double  = 5
    opts.adaptiveSigma   (1,1) logical = false
    opts.plotResults     (1,1) logical = true
    opts.zeropatchseconds   (1,1) double  = 5
    opts.restorezeropatches (1,1) logical = true
    opts.badchannels     (1,1) logical = false
    opts.badchannelsonly (1,1) logical = false
    opts.badchanavgref   (1,1) logical = true
    opts.badchanstride   (1,1) double  = 2
    opts.flatthreshold   (1,1) double  = 0.5
    opts.badchanfile     char          = ''
    opts.refresh         (1,1) logical = false
end

KeepTime = opts.KeepTime;

%%% Drop non-EEG channels
EEG = pop_select(EEG, 'nochannel', intersect(1:EEG.nbchan, opts.noteegchannels));

%%% Optional downsampling
if opts.targetsrate > 0 && EEG.srate ~= opts.targetsrate
    D = tic; fprintf('Resampling %d → %d Hz ...\n', EEG.srate, opts.targetsrate)
    EEG = pop_resample(EEG, opts.targetsrate);
    KeepTime.Downsample = toc(D);
end

%%% Cut out all-zero patches (amplifier crash padding) so no filter ever sees them.
%%% Must happen before the DC filter: filtfilt would smear them into large transients.
if opts.zeropatchseconds > 0
    EEG = run.excise_zero_patches(EEG, 'minseconds', opts.zeropatchseconds);
end

%%% DC removal
if opts.removeDC
    
    %%% Build filters
    fprintf('Building filters (srate = %d Hz) ...\n', EEG.srate)
    EEG_DCFilter_NumDen = filterbank(EEG.srate, 'DC_RCSquareFilt');

    %%% DC removal
    D = tic; fprintf('\nDC removal ...\n')
    EEG.data = filtfilt(EEG_DCFilter_NumDen(1,:), EEG_DCFilter_NumDen(2,:), double(EEG.data'))';
    KeepTime.DCRemoval = toc(D);
end

%%% Bad channel detection and removal
%%% Sits here, between the DC filter and Zapline, for two reasons. Detection needs the
%%% line noise: the second criterion in clean_channels is the ratio of >50 Hz to <45 Hz
%%% amplitude, which Zapline and CleanLine are about to flatten. Removal belongs here
%%% too: Zapline estimates its DSS spatial filter from the channel covariance, so a
%%% channel dominated by line noise both biases how many components are removed per
%%% chunk and gets its artefact smeared over the montage by the projection. The full
%%% montage is preserved in EEG.urchanlocs for interpolation.
%%%
%%% Detection runs on an average-referenced copy (see badchanavgref). Both criteria in
%%% clean_channels are ratios of unreconstructible noise to reconstructible signal, and
%%% referencing to a single electrode drives the denominator towards zero for the ring
%%% of channels next to it: subtracting Cz cancels the spatially smooth field (brain
%%% signal) but not per-electrode line pickup, which is impedance- and lead-driven and
%%% therefore rough. Measured on sub-drop0001/ses-t1, that ring went from a median
%%% RANSAC correlation of 0.99 post-Zapline to 0.73 pre-Zapline while every other
%%% channel stayed put - i.e. the vertex electrodes, the most valuable ones for sleep,
%%% were being flagged for the reference geometry rather than for being broken. The
%%% 45-50 Hz lowpass clean_channels applies before correlating only rejects ~31 dB at
%%% 50 Hz, so ~3% of the line amplitude reaches the correlation band; that is nothing
%%% against a normal channel's EEG and everything against a near-reference channel's.
if opts.badchannels
    if isempty(opts.badchanfile)
        error('run_filter:noBadChanFile', ...
            'Bad channel detection needs a cache path; pass it via the badchanfile option.')
    end
    if isempty(EEG.chanlocs) || ~isfield(EEG.chanlocs, 'X') || isempty([EEG.chanlocs.X])
        error('run_filter:noChanlocs', ...
            'Bad channel detection needs EEG.chanlocs with X/Y/Z coordinates.')
    end

    %%% Keep a record of the full montage before anything is dropped
    if ~isfield(EEG, 'urchanlocs') || isempty(EEG.urchanlocs)
        EEG.urchanlocs = EEG.chanlocs;
        for iCh = 1:numel(EEG.chanlocs)
            EEG.chanlocs(iCh).urchan = iCh;
        end
    end

    %%% Parameters for both criteria, kept in one place: the calls below read them from
    %%% here, and bidsfun_hp_zap_cleanline writes them to the JSON sidecar, so what is recorded is
    %%% necessarily what was run.
    %%% noiseThreshold is Inf on purpose: the line-noise criterion is computed and
    %%% reported, but no longer removes anything. On these recordings the electrodes it
    %%% ranks highest are consistently the ones around the vertex, and they are not
    %%% artefacts of the statistic - with the per-channel denominator removed they came
    %%% out further ahead still (100/81/70 sigma above the population). They genuinely
    %%% carry the most >45 Hz energy on the head, plausibly because the net's lead bundle
    %%% exits there. But their reconstruction correlation is ~0.99, i.e. the EEG in them
    %%% is intact, they sit where sleep activity matters most, and removing mains is what
    %%% the very next stage is for. znoise stays in the .mat, the channels.tsv and the
    %%% topoplot, so a pathological channel is still visible; it just no longer decides.
    %%% Set noiseThreshold back to a finite value to restore the old behaviour - the
    %%% topoplot marker follows the same threshold, so there is only one number to set.
    bcp = struct( ...
        'corrThreshold',       0.7, ...
        'noiseThreshold',      Inf, ...
        'windowSeconds',       5, ...
        'maxBrokenTime',       0.5, ...
        'numSamples',          25, ...
        'subsetSizeFraction',  0.25, ...
        'windowStride',        opts.badchanstride, ...
        'averageReferenced',   opts.badchanavgref, ...
        'flatThresholdMicroV', opts.flatthreshold, ...
        'flatMaxBrokenTime',   0.5);

    %%% Flat channels, on the data as recorded. Must come before the average reference:
    %%% a dead channel reads as one constant value here, but subtracting the common
    %%% average turns it into minus that average, which has real variance and is not
    %%% flat at all. Same window grid and duration rule as the correlation criterion.
    D = tic; fprintf('\nFlat channel detection ...\n')
    [flatmask, flatprop] = gedai.detectFlatChannels(EEG.data, EEG.srate, ...
        'windowseconds',  bcp.windowSeconds, ...
        'maxbrokentime',  bcp.flatMaxBrokenTime, ...
        'threshold',      bcp.flatThresholdMicroV);
    KeepTime.FlatChannelDetection = toc(D);

    D = tic; fprintf('\nBad channel detection ...\n')

    %%% Average-reference for detection only, then undo it: the data handed to Zapline
    %%% and written out keeps the reference it came in with.
    avgRef = [];
    if opts.badchanavgref
        fprintf('Average-referencing for detection ...\n')
        avgRef   = sum(EEG.data, 1) / (size(EEG.data, 1) + 1);   % +1: implicit reference channel
        EEG.data = EEG.data - avgRef;
    end

    %%% The flat mask is folded in inside the cached call, so the mask on disk is the
    %%% one actually applied - bidsfun_gedai reads it back to index the leadfield.
    [removed_channels, corrs, znoise, flatprop] = smartcache( ...
        @() detectBadChannels(EEG, flatmask, flatprop, bcp), ...
        opts.badchanfile, opts.refresh, ...
        {'', 'removed_channels', 'corr', 'znoise', 'flatprop'});

    if ~isempty(avgRef)
        EEG.data = EEG.data + avgRef;
        clear avgRef
    end

    EEG.etc.badchans = struct('mask', removed_channels, 'corr', corrs, ...
        'znoise', znoise, 'flatprop', flatprop, 'params', bcp);
    KeepTime.BadChannelDetection = toc(D);

    if opts.badchannelsonly
        %%% Detection stage: leave the montage intact so the caller can label every
        %%% channel, and skip the filtering below entirely.
        fprintf('%d/%d channels flagged as bad (detection only, data left intact).\n', ...
            nnz(removed_channels), EEG.nbchan)
    else
        fprintf('Removing %d/%d channels flagged as bad.\n', nnz(removed_channels), EEG.nbchan)
        EEG = pop_select(EEG, 'nochannel', find(removed_channels));
    end

    EEG.etc.filterparams.BadChannels = bcp;
    EEG.etc.filterparams.BadChannels.nRemoved = nnz(removed_channels);
end

%%% Zapline
%%% Skipped under badchannelsonly: that mode exists so bidsfun_detect_badchans can reach
%%% the detection step through exactly this code path, rather than reimplementing the
%%% drop / resample / excise / DC chain that decides what detection actually sees.
if opts.zapline && ~opts.badchannelsonly
    D = tic; fprintf('\nZapline plus ...\n')
    [EEG.data, zaplineConfig, analyticsResults] = clean_data_with_zapline_plus( ...
        double(EEG.data), EEG.srate, ...
        'noisefreqs',      opts.noisefreqs, ...
        'adaptiveNremove', opts.adaptiveNremove, ...
        'fixedNremove',    opts.fixedNremove, ...
        'chunkLength',     opts.chunkLength, ...
        'noiseCompDetectSigma', opts.noiseCompDetectSigma, ...
        'adaptiveSigma',   opts.adaptiveSigma, ...
        'plotResults',     opts.plotResults);
    EEG.etc.zapline.config    = zaplineConfig;
    EEG.etc.zapline.analytics = analyticsResults;
    EEG.etc.filterparams.Zapline = struct( ...
        'noisefreqs',           opts.noisefreqs, ...
        'adaptiveNremove',      opts.adaptiveNremove, ...
        'fixedNremove',         opts.fixedNremove, ...
        'chunkLengthSeconds',   opts.chunkLength, ...
        'noiseCompDetectSigma', opts.noiseCompDetectSigma, ...
        'adaptiveSigma',        opts.adaptiveSigma);
    KeepTime.Zapline = toc(D);
    fprintf('ZapLine-plus: %.2f min\n', KeepTime.Zapline / 60);
end

%%% Cleanline
if opts.cleanline && ~opts.badchannelsonly

%     % base CleanLine directly on what ZapLine actually found and treated
%     linefreqs = zaplineConfig.noisefreqs;
%     fprintf('ZapLine-plus detected and cleaned: %s Hz\n', mat2str(round(linefreqs,2)));
% 
%     % only pass freqs that had a real noise ratio to CleanLine
%     keepIdx = analyticsResults.ratioNoiseRaw > 1.5;   % adjust threshold based on what you see
%     linefreqs = zaplineConfig.noisefreqs(keepIdx);
%     fprintf('Only keep those with noise ratio >1.5: %s Hz\n', mat2str(round(linefreqs,2)));

    D = tic; fprintf('\nClean line ...\n')
    % EEG.data = double(EEG.data);
    clp = struct('linefreqs', [opts.noisefreqs opts.noisefreqs*2], ...
        'winsize', 4, 'winstep', 2, 'sigtest', true, 'pad', 2);
    EEG = cleanline_fast(EEG, 'linefreqs', clp.linefreqs, ...
        'winsize', clp.winsize, 'winstep', clp.winstep, ...
        'sigtest', clp.sigtest, 'pad', clp.pad);
    EEG.etc.filterparams.CleanLine = clp;

    % EEG = pop_cleanline(EEG, ...
    %     'chanlist', 1:EEG.nbchan, ...
    %     'linefreqs', 50, ...
    %     'winsize', 4, ...
    %     'winstep', 2, ...
    %     'computepower', false);
    KeepTime.Cleanline = toc(D);

end

%%% Put the all-zero patches back so the data regains its original length and timing.
%%% bidsfun_hp_zap_cleanline sets this to false and restores later, after its own second Zapline
%%% pass, so that pass also runs on patch-free data.
if opts.restorezeropatches
    EEG = run.restore_zero_patches(EEG);
end

%%% Write JSON sidecar
if opts.JsonFile ~= ""
    sidecarjson(KeepTime, char(opts.JsonFile));
end
end

% -------------------------------------------------------------------------
function [signal, removed_channels, corrs, znoise, flatprop] = detectBadChannels(EEG, flatmask, flatprop, bcp)
% Union of the clean_channels criteria and the flat-line criterion, as one cached unit.
% Kept together so the mask written to the cache is the mask actually applied to the
% data: bidsfun_gedai reads it back to pick the matching rows of the leadfield, and a
% cache holding only part of the criteria would silently desync from the saved file.
    [signal, removed_channels, corrs, znoise] = clean_channels(EEG, ...
        bcp.corrThreshold, bcp.noiseThreshold, bcp.windowSeconds, bcp.maxBrokenTime, ...
        bcp.numSamples, bcp.subsetSizeFraction, bcp.windowStride);
    removed_channels = removed_channels(:) | flatmask(:);
end
