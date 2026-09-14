# GEDAI ground-truth benchmark

`bidsfun_gedai_benchmark` measures GEDAI cleaning quality on real recordings by injecting
signal and artifacts whose ground truth is known, running GEDAI through the production code path
(`run.GEDAI_StageSpecific`), and scoring what comes out.

| Question | Injected | Headline metrics |
|---|---|---|
| Does GEDAI remove brain signal it should keep? (over-cleaning) | fast/slow spindles, K-complexes, slow oscillations | SER, preserved gain, amplitude ratio, morphology correlation, detector recall/density with frozen thresholds |
| Does GEDAI leave artifacts in? (under-cleaning) | EOG (blinks, saccades, slow eye movements), EMG, ECG | ARR, residual gain, temporal/spectral RRMSE, CC, in-band PSD ratios |

## 1. Core design

**Semi-synthetic, paired.** Per recording, `nepochs` epochs per stage are read from the GEDAI
input (`desc-hpzc`) and concatenated: this is the background *B*. GEDAI runs once on *B* alone
(baseline, *B̂*) and once per condition on *B + S* or *B + A* (*X̂*). Everything is scored on
the paired difference *X̂ − B̂*:

- Real sleep EEG is never artifact-free and already contains real spindles and slow waves.
  Comparing *X̂* with *B* would count GEDAI's legitimate cleaning of *B* as error. Comparing
  with *B̂* cancels it: the background and the GEDAI configuration are identical, and only the
  injection differs.
- GEDAI is nonlinear, with stage-group-level thresholds. *X̂ − B̂* therefore contains what
  survived of the injection plus any change in how the background was cleaned *because* of the
  injection. Both count as error. `gain` / `residual_gain` (projections onto the truth) separate
  the two: they ignore error uncorrelated with the injection.
- The unpaired forms (`SER_raw_dB`, `ARR_raw_dB`, exactly Somers et al. 2018) are also reported.
  They are the right numbers on `'background', 'synthetic'`, where *B* is clean by construction.

**Null control.** The `null` condition runs GEDAI a second time on the unchanged background.
Two runs on identical input should agree. Whatever difference remains is the floor under every
paired metric (`null_RRMSE`, `null_SNR_dB`, `null_maxabs_uV`). An SER or ARR is only meaningful
where the injection's energy is well above that floor.

**Same code path.** Stage groups, modes per stage (`GEDAIMode_dict`), and leadfield reference
covariance (`gedai.loadrefcov`) are identical to `bidsfun_gedai`. A group is re-run whenever any
of its epochs was touched by an injection, because thresholds are group-level. Untouched groups
are copied from the baseline run.

**Separate conditions.** Graphoelements and each artifact type are injected in separate GEDAI
runs, so every error can be attributed to one cause. `signalscales` / `artifactscales` repeat a
condition at other amplitudes on the same realisation (dose-response). Every entry of `runs` is
scored on the same injections, so GEDAI configurations can be compared head to head.

## 2. Signal source: parametric templates vs neural mass models

The benchmark needs events whose **timing, extent and amplitude are known exactly**, with
realistic **spatial topography** (what a spatial filter acts on) and plausible
**spectral/temporal morphology**. The options:

| Option | What it is | Feasibility here | Ground truth |
|---|---|---|---|
| **Parametric templates** (implemented, default) | Hann-envelope chirped sinusoids (spindles); area-matched biphasic Gaussian pairs (KCs, SOs) | MATLAB-native, no dependencies, milliseconds per event | Exact: onset, duration, amplitude, frequency, source patch |
| **Schellenberger Costa et al. 2016** thalamocortical NMM ([PLoS Comput Biol 12:e1005022](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005022); [ModelDB 226474](https://modeldb.science/226474)) | One cortical + one thalamic population model. Produces N2/N3 EEG with spindles, K-complexes and slow oscillations, and their timing relations. Cortex after Weigenand et al. 2014 (PLoS Comput Biol 10:e1003923). | C++ with a MATLAB interface on ModelDB. Needs a mex build; not tried here. | Model output is a single time course: events must be *labelled* by thresholding model variables, so the truth is model-derived, not designed |
| **neurolib** ALN cortex + thalamus (Cakan et al. 2023, Cogn Comput; Jajcay et al. 2022, [Front Comput Neurosci 16:769860](https://pmc.ncbi.nlm.nih.gov/articles/PMC9120371/); [code](https://github.com/jajcayn/thalamocortical_model_study)) | Adaptive cortical node (SOs) coupled to a thalamic node (fast spindles). Reproduces SO–spindle phase coupling. | Python, pip-installable, fits the repo's `uv` setup. Hours of single-node simulation are cheap. Not tried here. | Same labelling issue. Also one node, so no spatial information. |
| **Conductance-based networks** (Bazhenov et al. 2002, J Neurosci 22:8691; Krishnan et al. 2016, eLife 5:e18607) | Thousands of Hodgkin–Huxley-type neurons | Hours of CPU for minutes of activity, and a non-trivial LFP-to-EEG mapping | Labelling issue, and far more cost |

**Recommendation.** For a benchmark of a *spatial* cleaning method, biophysical realism inside
the event barely changes what is tested. What matters is realistic spatial extent, topography,
propagation, amplitude, density and spectrum, and all of those come from the leadfield
projection and the class parameters, whatever produces the time course. Full biophysical
simulation is not worth it. NMMs do add realistic *variability* and *coupling* of waveforms, so
the stage exposes a hook: `templatebank` takes waveforms per class (e.g. KCs, spindles and SOs
cut from a Costa-2016 or neurolib run, labelled offline), and placement, projection, scaling and
scoring stay unchanged. Generate a bank offline; do not embed an ODE solver in the pipeline.

Bank format (`.mat` with variable `bank`):

```matlab
bank.kcomplex.srate     = 250;
bank.kcomplex.waveforms = {w1, w2, ...};   % vectors, any length, arbitrary units
bank.spindle_fast       = struct('srate', 1000, 'waveforms', {{...}});
```

Waveforms are resampled to the recording rate and normalised: envelope peak 1 for spindles,
peak-to-peak 1 for slow waves.

## 3. Spatial models

**Graphoelements: cortical leadfield.** Each event picks an anchor channel (`anchors`). The
"source under the channel" is found from the leadfield alone, among sources whose strongest gain
falls on that channel (`bench.anchorSource`), so no electrode/MRI coregistration is needed.
Sources within `patchradius` fire with Gaussian weights through normal-constrained BEM gain.
Slow waves and KCs can travel front to back (`speed`; Massimini et al. 2004, 1.2–7 m/s), so their
field is not rank-1. Each event is scaled on its average-referenced field at the peak channel and
sign-corrected to be surface-negative. KCs and SOs can carry a nested fast spindle on their
up-state (`coupleclass`).

**Artifacts: deliberately not the cortical leadfield.** GEDAI keeps what the brain leadfield
explains. An artifact generated by cortical sources would be kept by construction, and the
benchmark would report under-cleaning that says nothing about real artifacts.

- **EOG**: a corneo-retinal dipole at each eye, placed from the nasion fiducial (`eyeoffset`),
  in a homogeneous sphere (`bench.sphereLeadfield`, the Zhang 1995 closed form, checked against
  Brainstorm's `bst_eeg_sph`). The sphere is fitted to the electrodes within 10 cm of the eyes.
  A whole-head fit is centred well above ear level on EGI nets and puts the eyes on its surface
  (1.01 R on sub-drop0007); the local fit puts them at 0.8 R with a 6 mm residual. Conjugate
  gaze rotations give blinks (vertical transient), saccades (step, hold, return) and slow eye
  movements (drift).
- **EMG**: band-limited (20–100 Hz) noise with a shared plus an independent part per channel,
  weighted by a Gaussian spread around muscle sites (temporalis, frontalis, occipital/neck,
  masseter), optionally bilateral. Sites are fixed positions in the fiducial head frame
  (`sitepositions`), mapped to the nearest channel. Modelled in electrode space because EMG
  comes from muscle directly beneath the electrodes. Locally high-rank.
- **ECG**: a far-field cardiac dipole about 32 cm below the head, infinite homogeneous medium,
  driven by a three-lead vectorcardiogram. Waves are Gaussian in cardiac phase (McSharry et al.
  2003; 3-lead form after Sameni et al. 2007), with heart-rate variability and respiratory sinus
  arrhythmia. Continuous in every contaminated epoch.

All fields are computed at the data channels and the reference electrode (Cz), re-referenced to
Cz as the data on disk are, then average-referenced exactly as `bidsfun_gedai` does
(`bench.avgRef`).

## 4. Detectors

Existing detectors found nearby:

- `ScoringHero/detect_spindle.m`: multitaper, after MT-KCD
- `GEDAI-master/tests/validate_whole_night_sleep.m`: single-channel Fz proxies
- `dusk2dawn/d2d_detectSlowWaves.m`

None allows **frozen thresholds**, which this benchmark needs. A detector calibrated on each
signal separately is relative: if GEDAI shrank every event by 30%, it would find the same events
and hide the loss. The benchmark therefore ships two small, standard detectors:

- `bench.detectSpindles`: Mölle et al. (2002) band-pass RMS, mean + 1.5 SD, 0.5–3 s
- `bench.detectSlowWaves`: Massimini et al. (2004) zero crossings, with thresholds relative to
  the candidate distribution

Both are calibrated **once on the uncleaned input *B + S*** and applied unchanged to *X̂*, *B*
and *B̂*. Detection runs on each event's peak channel. Slow waves are matched by trough time
(±0.2 s), because window overlap is met by chance in N3. On synthetic data, recall with frozen
thresholds falls from 1.00 to about 0.1–0.2 when slow waves are attenuated to 10%; with
recalibrated thresholds it only falls to about 0.5 (`tests/test_benchmark_units.m`). For publication, cross-check
with an established detector such as YASA (Vallat & Walker 2021, eLife 10:e70092) or A7 (Lacourse
et al. 2019, J Neurosci Methods 316:3), fed the four signals with fixed thresholds.

## 5. Metrics

Written to `<fileID>_desc-<desc>_metrics.tsv` as a long table (`run, runname, condition, scale,
class, stage, metric, value`). Stage `all` pools the stages. Full definitions are in the headers
of `bench.scoreSignals` and `bench.scoreArtifacts` and in the JSON sidecar.

- Over-cleaning:
  - `SER_dB`, `SER_raw_dB`, `gain`
  - `SER_event_median/p10`, `gain_peakchan_median`, `ampratio_median/p10`, `morphcorr_median`
  - `recall_in/out/retention`
  - `density_in/out`, `density_bg_in/out`, `excess_density_in/out`, `density_retention`
  - `detamp_ratio_median`, `detfreq_shift_median`, `detduration_ratio_median`
  - `SER_background_dB` (GEDAI on the background alone)
- Under-cleaning:
  - `ARR_dB`, `ARR_raw_dB`, `residual_gain`
  - `RRMSE_t_in/out`, `RRMSE_s_in/out` (Zhang et al. 2021, EEGdenoiseNet), `CC_in/out`
  - `psd_contamination_dB`, `psd_residual_dB`, `psd_residual_excess` in the artifact bands
    (EOG 0.5–4 Hz, EMG 20–45 and 55–95 Hz, ECG 1–20 Hz)

## 6. Limitations

- **Subset, not whole night.** GEDAI thresholds come from the data they clean; 40 epochs per
  stage are fewer than a night's stage group. Paired comparisons stay valid, but absolute
  cleaning strength can differ somewhat from production. Raise `nepochs` if RAM allows.
- **Injected densities shift GEDAI's thresholds.** That is intended, since real nights contain
  these events, but results depend on `density`. Keep densities physiological, or sweep them.
- **Head models are approximations.** The sphere eye model ignores the orbit and skull
  conductivity; the ECG infinite medium ignores the neck path; eye positions come from
  anthropometric offsets. Topographies are right in shape (frontal EOG, lateral-inferior ECG
  gradient), not subject-exact.
- **Amplitudes are average-referenced**, so they read lower than mastoid-referenced norms.
  The defaults are plausible, not fitted to DROP.
- **ICA is never run** in the benchmark, whatever `r.ICAtype` says.

## 7. Running

```matlab
fails.benchmark = bidsfun_gedai_benchmark(BIDS, defaults{:}, ...
    'inputdesc', filtdesc, 'desc', 'hpzcgedbench', ...
    'leadfielddir', leadfieldpath, 'scoringpath', scoringpath, ...
    'runs', gedai.defaultRuns(), 'nepochs', 40);
```

Memory: about 8 × (4 stages × `nepochs` × 7500 × 256 × 4 bytes). That is roughly 10 GB at
`nepochs` 40, plus GEDAI's own working memory. Runtime is one baseline plus up to four
conditions of GEDAI per run, each on the touched stage groups only.

Unit checks (no data, no GEDAI): `run('Benchmark/tests/test_benchmark_units.m')`.
