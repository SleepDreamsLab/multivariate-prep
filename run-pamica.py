r"""
AMICA (pamica) on BIDS derivative recordings, either EEGLAB .set or BrainVision
.vhdr files.

Per recording: import -> Chebyshev-II high-pass -> rank projection (GEDAI may
have left the data rank deficient) -> AMICA -> unmixing matrices written in the
form EEGLAB expects in EEG.icaweights / EEG.icasphere / EEG.icawinv -> optionally
ICLabel (mne-icalabel) over those components.

The ICLabel pass is deliberately independent of the AMICA fit: it keys off the
_ica.mat, so a recording whose .mat already exists gets labelled without refitting.
That matters because a fit is hours on a GPU and a relabel is a minute -- see
run_iclabel() and RUN_ICLABEL below. The MATLAB equivalent, ICA/bidsfun_iclabel.m,
writes the same _iclabels.tsv schema from the canonical EEGLAB implementation; use
either, not both, per recording.

.venv\Scripts\Activate.ps1
python run-pamica.py
"""

import json
import time
from datetime import datetime
from pathlib import Path

import mne
import numpy as np
import scipy.io as sio
import torch
from pamica import AMICA
from scipy.signal import cheb2ord, cheby2, sosfiltfilt

# time.sleep(1 * 30 * 60)

BIDS_MAT      = Path(__file__).parent / "BidsFiles" / "BIDS_DROP.mat"
DERIV_IN_DIR  = "prep-zc-ged"  # derivatives subfolder to read the desc-* .set files from
DERIV_OUT_DIR = "prep-zc-ged"    # derivatives subfolder to write AMICA output under
DESC          = "hpzcged" # zc
OUT_DESC      = "pamica"# True # "zc2gedWakeBBAutoPlusFSAutoPlus2AmicaF32DllAutoStride4Rej0Nmodel1"  # desc entity for the AMICA output filename; None = same as DESC (the input's own desc)
SUBJECTS      = None # ["drop0001"]  # None = all subjects
SESSIONS      = None  # None # ["t1"]  # None = all sessions
TASKS         = ["Sleep", "sleep"]  # None = all tasks

DATA_EXTENSIONS = (".set", ".vhdr")  # preference order when both exist for a recording

# Input files can still be getting written by another machine as this script starts.
# The __main__ loop scans all requested recordings each pass, processing whichever
# already exist immediately and deferring the rest -- a missing file never blocks
# checking the next one. Once a full pass leaves recordings still missing, it waits
# WAIT_LOOP_MINUTES before rescanning just those. If no scan pass turns up a new
# file for WAIT_MAX_MINUTES straight, it gives up and raises for whatever's left.
WAIT_MAX_MINUTES = 180  # give up if no new file appears across scans for this long
WAIT_LOOP_MINUTES = 5   # how long to wait between rescans of the still-missing set


MAX_ITER = 700  # EEGLAB-AMICA's usual budget; pamica's fit default is lower
DO_REJECT = False  # Fortran-style outlier rejection; costs ~2x GPU memory (pamica clones the full array
                   # every iteration under do_reject, core.py:1710) -- OOMs on this recording size/GPU
DO_NEWTON = True   # Fortran-parity Newton preconditioner (tune newt_start/newtrate via fit() kwargs)
PASS_FRQ, STOP_FRQ, PASS_RIPPLE, STOP_ATTEN = 1.6, 0.8, 0.05, 30
RANK_TOL = 1e-7  # rank_projection: eigenvalue-ratio cutoff for the kept subspace

# AMICA knobs below are usually left alone -- named here (rather than left as call-site
# literals) so their values get recorded in amica_runtime.json for provenance.
N_MODELS   = 1              # AMICA: number of models
N_MIX      = 3              # AMICA: mixture components per source
BLOCK_SIZE = 8192           # AMICA: E-step accumulation chunk size; pure chunking, doesn't change results
DTYPE      = torch.float32  # AMICA: compute dtype
DEVICE     = "cuda"         # AMICA: compute device
MIN_DLL     = 5e-7    # 1e-9 is pamica's own default (Fortran parity). Note it's below float32's LL-
                      # accumulation noise floor (~1e-7-1e-6) and pamica logs a warning at construction
                      # if dtype=torch.float32 with min_dll < 1e-6 -- raise this or switch DTYPE to
                      # torch.float64 if that matters for your run.
MAXINCS     = 5       # ...for this many consecutive iterations
USE_MIN_DLL = True    # gate for the above; False disables the check entirely

# Second, independent early-stop pAMICA added alongside min_dll (Fortran
# use_grad_norm/min_nd): stops once the RMS weight-update norm falls to/below
# min_nd, regardless of the LL-gain check above. This is specifically what fixes
# the do_newton=True-on-CUDA case where lrate sits at newtrate and oscillates
# instead of annealing -- previously nothing caught that, so max_iter was the
# only real stop. Defaults match pamica's own (both True/1e-7); named here only
# so they land in amica_runtime.json like every other AMICA knob.
USE_GRAD_NORM = True  # gate for the min_nd stop; independent of USE_MIN_DLL
MIN_ND = 1e-7          # RMS weight-update-norm threshold for the above

# Decimation: ICA wants roughly MIN_SAMPLES_FACTOR x k^2 samples (k = rank kept by
# rank_projection) to be well-determined; take every STRIDEth sample instead of
# every sample when there's more data than that needs, up to MAX_STRIDE. Kept as
# an "every Xth sample" stride rather than a plain sample cap so it still spans
# the full recording instead of just its first portion.
MIN_SAMPLES_FACTOR = 30  # target sample count = k^2 * this
MAX_STRIDE = 4            # largest allowed stride; never take every 5th sample or coarser


# ICLabel (mne-icalabel), run over the components AMICA just fitted -- or over the
# ones an earlier run already wrote, since this keys off the _ica.mat rather than the
# fit. Writes <mat_stem>_iclabels.tsv (+ .json data dictionary) beside the .mat and
# folds the probability matrix back into the .mat itself.
RUN_ICLABEL = True
ICLABEL_BACKEND = None   # None = mne-icalabel picks: torch if installed, else onnx

# ICLabel sees this many minutes, as evenly spaced chunks rather than the head of the
# recording, so every sleep stage contributes in roughly the proportion it occupies.
# 0 = hand over the whole recording (the default here).
#
# Mind the memory if raising this above 0: _eeg_autocorr_welch stacks the activations
# through a 50%-overlapping 3-s window index in float64, so it transiently needs well
# over an order of magnitude more than the recording itself -- a 236-component 8-h
# night at 250 Hz needs on the order of 40 GB there. Set this to e.g. 60 if that OOMs
# on your machine; both data-derived features (median 1-s PSD, mean 3-s
# autocorrelation) are window averages and converge well inside an hour.
ICLABEL_MINUTES = 0
ICLABEL_CHUNK_SECONDS = 60

# Marking components bad: top class among these, at or above this probability.
# Indices are into ICLABEL_CLASSES below (0-based), i.e. everything except brain and
# "other" -- the same set ICA/bidsfun_iclabel.m defaults to.
ICLABEL_ARTEFACT_CLASSES = (1, 2, 3, 4, 5)
ICLABEL_THRESHOLD = 0.5

# GEDAI average-references before it runs, so the derivative this reads already is
# common-average even though nothing in the file says so. Declaring it keeps
# mne-icalabel from warning once per recording about a reference that is in fact
# correct. Set False if you ever feed this something that is NOT average referenced.
ICLABEL_ASSUME_AVGREF = True

# mne-icalabel reports 'muscle artifact' / 'eye blink' / 'heart beat'; the MATLAB
# original reports 'Muscle' / 'Eye' / 'Heart'. The network and the column order are
# the same, so the labels are renamed to the MATLAB spelling and the .tsv written
# here is interchangeable with bidsfun_iclabel.m's.
ICLABEL_CLASSES = ("Brain", "Muscle", "Eye", "Heart", "Line Noise", "Channel Noise", "Other")
ICLABEL_COLUMNS = ("p_brain", "p_muscle", "p_eye", "p_heart",
                   "p_line_noise", "p_channel_noise", "p_other")


def _unwrap(x):
    """Peel the 1x1 cell / array wrappers left by loadmat."""
    while not isinstance(x, dict):
        x = np.atleast_1d(x)[0]
    return x


def load_bids(mat_path):
    """Load the BIDS layout struct dumped from MATLAB (bids-matlab's BIDS object)."""
    mat = sio.loadmat(mat_path, simplify_cells=True)
    return _unwrap(mat[next(k for k in mat if not k.startswith("__"))])

def deriv_entries(bids, desc, deriv_in):
    """Build (base, candidates) pairs for every raw .vhdr recording in the BIDS
    struct that matches the SUBJECTS/SESSIONS/TASKS filters, mapped onto its
    prep-ged derivative file -- WITHOUT checking whether that file actually
    exists yet (see the scan/retry loop in __main__ for that).

    Mirrors the MATLAB filtering (bids.query(..., 'extension', '.vhdr')) then
    bids.internal.parse_filename entity join, using the ext/entities fields
    that BIDS_DROP.mat already carries per recording instead of re-deriving
    them from the filename. The derivative can be either an EEGLAB .set or a
    BrainVision .vhdr file; both extensions are tried per recording (in
    DATA_EXTENSIONS order).
    """
    entries = []
    for sub in np.atleast_1d(bids["subjects"]):
        for rec in np.atleast_1d(sub["eeg"]):
            if rec["ext"] != ".vhdr":
                continue
            ent = rec["entities"]
            if SUBJECTS is not None and ent["sub"] not in SUBJECTS:
                continue
            if SESSIONS is not None and ent.get("ses") not in SESSIONS:
                continue
            if TASKS is not None and ent.get("task") not in TASKS:
                continue
            file_id = "_".join(f"{k}-{v}" for k, v in ent.items())
            folders = [f"sub-{ent['sub']}"] + ([f"ses-{ent['ses']}"] if "ses" in ent else [])
            base = deriv_in.joinpath(*folders, f"{file_id}_desc-{desc}_eeg")
            candidates = [base.with_suffix(ext) for ext in DATA_EXTENSIONS]
            entries.append((base, candidates))
    return entries


def _scan_pass(entries):
    """One pass over `entries` (list of (base, candidates)): returns (found,
    still_missing) -- found is the list of resolved Paths (existing now),
    still_missing is the (base, candidates) entries that don't exist yet, to
    retry on a later pass. A missing entry never blocks checking the next one.
    """
    found = []
    still_missing = []
    for base, candidates in entries:
        hits = [c for c in candidates if c.is_file()]
        if not hits:
            still_missing.append((base, candidates))
            continue
        if len(hits) > 1:
            print(f"both {' and '.join(c.suffix for c in hits)} present for {base.name}, using {hits[0].suffix}")
        found.append(hits[0])
    return found, still_missing


def _fmt_small(x):
    """Compact filename-safe tag for a small positive float: 1e-09 -> '1e9'."""
    mantissa, exponent = f"{x:.0e}".split("e")
    return f"{mantissa}e{abs(int(exponent))}"


def build_auto_desc(stride):
    """Auto-generate an output desc by appending "2amica-" plus this run's key
    parameters onto DESC: RejX (DO_REJECT), NmodelX (N_MODELS), F32/F64 (DTYPE),
    DllX (MIN_DLL), MaxIterX (MAX_ITER), StrideX (the decimation stride used)."""
    dtype_tag = "F32" if DTYPE == torch.float32 else "F64"
    return (
        f"{DESC}2amica-Rej{int(DO_REJECT)}Nmodel{N_MODELS}{dtype_tag}"
        f"Dll{_fmt_small(MIN_DLL)}MaxIter{MAX_ITER}Stride{stride}"
    )


def output_paths(data_file, deriv_out, out_desc, stride=None):
    """Derive this recording's AMICA output paths: deriv_out/sub-X/[ses-Y/]<mat_stem>,
    where <mat_stem> is data_file's stem with its desc swapped to out_desc and its
    trailing _eeg replaced with _ica.

    out_desc: None = same desc as the input; True = auto-generate from this run's
    parameters via build_auto_desc() (requires stride, since that's only known
    after decimation -- call this after rank_projection/decimation, not before);
    a string = use it literally.

    Returns (mat_path, bin_dir, json_path): the .mat file, the same-named folder for
    write_amica_output()'s binaries (no extension), and a same-named json.
    """
    stem = data_file.stem
    if not stem.endswith("_eeg"):
        raise ValueError(f"expected a derivative filename ending in '_eeg', got: {stem}")
    parts = stem.split("_")
    sub = next(p for p in parts if p.startswith("sub-"))
    ses = next((p for p in parts if p.startswith("ses-")), None)

    if out_desc is None:
        effective_desc = DESC
    elif out_desc is True:
        if stride is None:
            raise ValueError("OUT_DESC=True requires stride -- call output_paths() after decimation")
        effective_desc = build_auto_desc(stride)
    else:
        effective_desc = out_desc
    mat_stem = stem[: -len("_eeg")].replace(f"desc-{DESC}", f"desc-{effective_desc}") + "_ica"

    session_dir = deriv_out / sub / ses if ses else deriv_out / sub
    return (
        session_dir / f"{mat_stem}.mat",
        session_dir / mat_stem,
        session_dir / f"{mat_stem}.json",
    )


def highpass(X, sfreq):
    """Minimum-order Chebyshev-II high-pass, zero phase (= designfilt + filtfilt)."""
    order, wn = cheb2ord(PASS_FRQ, STOP_FRQ, PASS_RIPPLE, STOP_ATTEN, fs=sfreq)
    print(f"  high-pass: order {order}, Wn {wn:.4g} (pass {PASS_FRQ} Hz / stop {STOP_FRQ} Hz @ {sfreq} Hz)")
    sos = cheby2(order, STOP_ATTEN, wn, btype="highpass", output="sos", fs=sfreq)
    t0 = time.perf_counter()
    Xf = sosfiltfilt(sos, X, axis=1)
    print(f"  high-pass: filtered {X.shape[0]} channels x {X.shape[1]} samples in {time.perf_counter() - t0:.1f} s")
    return Xf


def rank_projection(X, tol=RANK_TOL):
    """Orthonormal k x n projection onto the non-degenerate subspace."""
    print(f"  rank projection: eigendecomposing {X.shape[0]}x{X.shape[0]} covariance...")
    t0 = time.perf_counter()
    d, V = np.linalg.eigh(X @ X.T / X.shape[1])
    d, V = d[::-1], V[:, ::-1]
    ratio = d / d[0]
    k = int(np.sum(ratio > tol))
    tail = " ".join(f"{r:.1e}" for r in ratio[max(k - 3, 0):k + 3])
    print(f"  rank {k}/{len(d)} at tol {tol:.0e}; eigenvalue ratios across the cut: {tail}"
          f" ({time.perf_counter() - t0:.1f} s)")
    return V[:, :k].T


def read_raw(data_file, preload=True):
    """Load a derivative recording with the MNE reader matching its extension.

    preload=False is what the ICLabel path uses: it only ever wants a handful of
    chunks, and preloading a full night to throw most of it away is the one thing
    that makes that step expensive.
    """
    if data_file.suffix == ".set":
        return mne.io.read_raw_eeglab(data_file, preload=preload)
    if data_file.suffix == ".vhdr":
        return mne.io.read_raw_brainvision(data_file, preload=preload)
    raise ValueError(f"unsupported derivative extension: {data_file.suffix} ({data_file})")


def run_amica(data_file, deriv_out):
    # Early out. With a literal OUT_DESC the output name is fully determined by the
    # input name, so an already-fitted night costs one stat() here instead of a full
    # preload plus an nchan x nchan eigendecomposition before the late check below --
    # which matters now that the point of revisiting a finished night is usually just
    # to give it the ICLabel pass it never got.
    #
    # OUT_DESC=True cannot use this: it encodes the decimation stride in the filename,
    # and the stride is only known after preprocessing. That case still falls through.
    if OUT_DESC is not True:
        mat_path, _, _ = output_paths(data_file, deriv_out, OUT_DESC)
        if mat_path.is_file():
            print(f"skip {data_file.stem} (already done: {mat_path.name})")
            return None, mat_path

    raw = read_raw(data_file)
    picks = mne.pick_types(raw.info, eeg=True, exclude="bads")
    if raw.info["bads"]:
        print(f"  {len(raw.info['bads'])} channel(s) flagged bad in the file, excluded: {raw.info['bads']}")
    labels = [raw.ch_names[i] for i in picks]
    sfreq = raw.info["sfreq"]
    X = raw.get_data(picks=picks) * 1e6  # MNE volts -> EEGLAB microvolts
    del raw

    print(f"{data_file.stem}: {len(labels)} channels")
    # X = highpass(X, sfreq)
    # X -= X.mean(axis=1, keepdims=True)
    P = rank_projection(X)  # kept in float64: cheap, and wants a clean rank cut
    Xr = (P @ X).astype(np.float32)
    del X

    # Decimate to every STRIDEth sample if there's more data than ICA needs (see
    # MIN_SAMPLES_FACTOR/MAX_STRIDE above). Picks the largest stride in
    # 1..MAX_STRIDE that still leaves >= k^2 * MIN_SAMPLES_FACTOR samples; falls
    # back to stride=1 (no decimation) if even that isn't enough headroom.
    k, n_full = Xr.shape
    min_samples = k**2 * MIN_SAMPLES_FACTOR
    stride = next((s for s in range(MAX_STRIDE, 0, -1) if n_full // s >= min_samples), 1)
    if stride > 1:
        Xr = Xr[:, ::stride]
    print(f"  decimation: k={k}, target >= {min_samples} samples (k^2 x {MIN_SAMPLES_FACTOR}); "
          f"stride={stride} -> {Xr.shape[1]}/{n_full} samples")

    # Resolved only now: with OUT_DESC=True the filename depends on stride, which
    # is only known post-decimation -- so the "already done" check has to live
    # here too, after preprocessing, rather than before it like a plain-desc run.
    mat_path, bin_dir, json_path = output_paths(data_file, deriv_out, OUT_DESC, stride)
    if mat_path.is_file():
        print(f"skip {data_file.stem} (already done: {mat_path.name})")
        return None, mat_path

    model = AMICA(n_models=N_MODELS, n_mix=N_MIX, device=DEVICE)
    t0 = time.perf_counter()
    model.fit(Xr, max_iter=MAX_ITER, block_size=BLOCK_SIZE, dtype=DTYPE, do_reject=DO_REJECT, do_newton=DO_NEWTON,
              min_dll=MIN_DLL, maxincs=MAXINCS, use_min_dll=USE_MIN_DLL,
              use_grad_norm=USE_GRAD_NORM, min_nd=MIN_ND)
    elapsed = time.perf_counter() - t0
    print(f"  AMICA fit took {elapsed / 60:.1f} min ({elapsed:.1f} s)")

    mat_path.parent.mkdir(parents=True, exist_ok=True)
    # Always written first: byte-identical to the Fortran reference, unaffected by
    # anything below, and the fallback if the .mat convention check fails.
    model.write_amica_output(str(bin_dir))  # NB: in the rank-reduced space

    # get_unmixing_matrix() returns W^T *alone* -- only PART of the unmixing.
    # rank_projection() only rotates/decorrelates (no variance normalization), so
    # with do_sphere defaulting True, pamica does real additional whitening of Xr
    # internally; the actual unmixing that operates on Xr is W^T @ sphere (this is
    # exactly how pamica's own mir() composes it, core.py:2086). Using W^T alone
    # here was the bug behind the earlier 0.012 (and a later 236/236) convention-
    # check failures -- not float32 pinv conditioning as first suspected.
    sphere = model.model_.sphere.detach().cpu().numpy()
    W_full = model.get_unmixing_matrix() @ sphere

    # Sanity check: W_full @ data should reproduce model.transform(data). These are
    # two independent pamica code paths (transform() re-runs the forward pass;
    # W_full reads the stored W/sphere tensors directly), so agreement here is a
    # real check on what's about to be saved, not a tautology. transform() also
    # subtracts data-space mean/center offsets that W_full @ data doesn't -- but
    # those are per-source constants, invisible to a correlation check, so this
    # still validates the part that matters (the linear/sphere composition). A
    # failure here must stop the .mat from being written rather than just warn.
    ref, chk = model.transform(Xr[:, :5000]), W_full @ Xr[:, :5000]
    r = np.abs(np.array([np.corrcoef(a, b)[0, 1] for a, b in zip(ref, chk)]))
    n_bad = int(np.sum(r < 0.999))
    if n_bad:
        raise RuntimeError(
            f"{data_file.stem}: mixing-matrix convention check failed -- "
            f"{n_bad}/{len(r)} components below 0.999 |corr| "
            f"(median {np.median(r):.6f}, min {r.min():.6f}). "
            f"amicaout was written; amica_eeglab.mat was not."
        )

    # Diagnostic only: A and W are maintained as mutual (approximate) inverses
    # within pamica's own internal EM updates, so this checks that internal
    # consistency -- it does NOT exercise the sphere-composition above (both sides
    # omit it identically), so it will not catch that class of bug on its own.
    W_pinv = np.linalg.pinv(model.get_mixing_matrix())
    r_pinv = np.abs(np.array([np.corrcoef(a, b)[0, 1] for a, b in zip(model.get_unmixing_matrix(), W_pinv)]))
    print(f"  diagnostic: get_unmixing_matrix() vs pinv(get_mixing_matrix()) -- "
          f"median |corr|={np.median(r_pinv):.6f}, min={r_pinv.min():.6f}, "
          f"{int(np.sum(r_pinv < 0.999))}/{len(r_pinv)} below 0.999")

    order = model.variance_order()
    icaweights = W_full[order, :]  # EEGLAB order: IC1 = highest variance
    icasphere = P  # kept in float64 throughout
    # loadmodout15's formula: A = pinv(W @ S), cast to float64 before inverting.
    icawinv = np.linalg.pinv(icaweights.astype(np.float64) @ icasphere)

    sio.savemat(
        mat_path,
        {
            "icaweights": icaweights,
            "icasphere": icasphere,
            "icawinv": icawinv,
            "chanlabels": np.array(labels, dtype=object),
            "setfile": str(data_file),
            "srate": sfreq,
            "final_ll": model.final_ll_,
            "stop_reason": model.stop_reason_,
            "n_iter": len(model.ll_history_),
            "rank_kept": Xr.shape[0],
            "n_chan": len(labels),
            "dtype": str(Xr.dtype),
            "LL": np.asarray(model.ll_history_, dtype=np.float64),
        },
    )
    with open(json_path, "w") as f:
        json.dump(
            {
                "PamicaParameters": {
                    "max_iter": MAX_ITER,
                    "do_reject": DO_REJECT,
                    "do_newton": DO_NEWTON,
                    "pass_frq": PASS_FRQ,
                    "stop_frq": STOP_FRQ,
                    "pass_ripple": PASS_RIPPLE,
                    "stop_atten": STOP_ATTEN,
                    "rank_tol": RANK_TOL,
                    "n_models": N_MODELS,
                    "n_mix": N_MIX,
                    "block_size": BLOCK_SIZE,
                    "dtype": str(DTYPE),
                    "device": DEVICE,
                    "min_samples_factor": MIN_SAMPLES_FACTOR,
                    "max_stride": MAX_STRIDE,
                    "stride": stride,
                    "n_samples": Xr.shape[1],
                    "min_dll": MIN_DLL,
                    "maxincs": MAXINCS,
                    "use_min_dll": USE_MIN_DLL,
                    "use_grad_norm": USE_GRAD_NORM,
                    "min_nd": MIN_ND,
                },
                "ProcessingDurationsMinutes": elapsed / 60,
                "GeneratedDate": datetime.now().isoformat(),
            },
            f,
            indent=2,
        )
    return model.final_ll_, mat_path


def _eeglab_style_ica(icaweights, icasphere, info):
    """Wrap an externally computed unmixing (icaweights @ icasphere) in an mne ICA.

    Mirrors mne.preprocessing.read_ica_eeglab, which is MNE's own supported route for
    bringing an EEGLAB-shaped decomposition in -- the same route mne-icalabel's docs
    point at. The only difference is the source: the matrices come from memory here
    instead of out of a .set file.

    Why the SVD. mne-icalabel reads the decomposition back as
    unmixing_matrix_ @ pca_components_, so those two have to multiply out to the full
    channel-space unmixing. AMICA's icaweights (ncomp x k) and icasphere (k x nchan)
    do not fit those slots directly when k < nchan -- the rank projection makes
    icaweights non-square in the wrong way -- so the product is re-split by SVD into
    a square (u*s) and an orthonormal v, which do.
    """
    from mne.preprocessing import ICA

    n_components = icaweights.shape[0]
    unmixing = icaweights @ icasphere                      # ncomp x nchan
    u, s, v = np.linalg.svd(unmixing, full_matrices=False)

    ica = ICA(method="imported_eeglab", n_components=n_components)
    ica.current_fit = "eeglab"
    ica.ch_names = info["ch_names"]
    ica.n_pca_components = None
    ica.n_components_ = n_components
    ica.pre_whitener_ = np.ones((len(info["ch_names"]), 1))
    ica.pca_mean_ = np.zeros(len(info["ch_names"]))
    ica.unmixing_matrix_ = u * s
    ica.pca_components_ = v
    ica.pca_explained_variance_ = s * s
    ica.info = info
    ica._update_mixing_matrix()
    ica._update_ica_names()
    ica.reject_ = None

    # What mne-icalabel will actually reconstruct and apply to the data
    # (features._retrieve_eeglab_icawinv), recomputed here from public attributes.
    # If the split above were wrong, the components would still get labels -- just
    # labels for the wrong sources -- so this has to be an error, not a warning.
    roundtrip = ica.unmixing_matrix_ @ ica.pca_components_[:n_components, :]
    err = np.abs(roundtrip - unmixing).max() / np.abs(unmixing).max()
    if err > 1e-6:
        raise RuntimeError(
            f"ICA round-trip mismatch ({err:.2e}): unmixing_matrix_ @ pca_components_ "
            f"does not reproduce icaweights @ icasphere."
        )
    return ica


def iclabel_subset(data_file, ch_names, minutes, chunk_seconds):
    """Load `minutes` of `data_file` as evenly spaced chunks, as an mne Raw.

    Returns (raw, minutes_used). minutes <= 0, or a recording shorter than that,
    loads everything. Channels are picked by name, in the order the decomposition
    was fitted with, so raw.ch_names lines up with the columns of icasphere.
    """
    import mne
    from mne.io import RawArray

    raw = read_raw(data_file, preload=False)
    picks = mne.pick_channels(raw.ch_names, include=list(ch_names), ordered=True)
    info = mne.pick_info(raw.info, picks)
    sfreq = raw.info["sfreq"]
    n_total = raw.n_times

    n_want = int(round(minutes * 60 * sfreq))
    if minutes <= 0 or n_want >= n_total:
        data = raw.get_data(picks=picks)
    else:
        chunk = min(int(round(chunk_seconds * sfreq)), n_total)
        n_chunk = max(1, n_want // chunk)
        starts = np.unique(np.linspace(0, n_total - chunk, n_chunk).round().astype(int))
        data = np.hstack([raw.get_data(picks=picks, start=s, stop=s + chunk) for s in starts])
    del raw

    if ICLABEL_ASSUME_AVGREF:
        with info._unlock():
            info["custom_ref_applied"] = True
    out = RawArray(data, info, verbose="ERROR")
    return out, out.n_times / sfreq / 60


def run_iclabel(mat_path):
    """Label the components in `mat_path` with ICLabel and write the results beside it.

    Reads the decomposition back out of the .mat rather than taking it from the fit,
    so this works the same whether AMICA just ran or ran last week.

    Writes three things:
      <stem>_iclabels.tsv   BIDS table, one row per component (the file to read or
                            hand-edit when deciding what to subtract)
      <stem>_iclabels.json  its data dictionary, plus which implementation ran
      <stem>_ica.mat        updated in place with 'ic_classification' (EEGLAB's own
                            nested shape), 'gcompreject' and 'varfrac', so MATLAB can
                            pick the whole decomposition up with ICA/loadica.m

    Returns the (n_components, 7) probability matrix, or None if the labels were
    already there.
    """
    from mne_icalabel.iclabel import iclabel_label_components

    # <fileID>_desc-X_ica.mat -> <fileID>_desc-X_iclabels.tsv, i.e. the _ica suffix is
    # replaced rather than appended. That is the name ICA/bidsfun_iclabel.m writes too,
    # so whichever implementation produced the labels, downstream code opens one path.
    stem = mat_path.stem[: -len("_ica")] if mat_path.stem.endswith("_ica") else mat_path.stem
    tsv_path = mat_path.with_name(stem + "_iclabels.tsv")
    json_path = mat_path.with_name(stem + "_iclabels.json")

    mat = sio.loadmat(mat_path, simplify_cells=True)
    if "ic_classification" in mat and tsv_path.is_file():
        print(f"  ICLabel: already labelled ({tsv_path.name})")
        return None

    ch_names = [str(c) for c in np.atleast_1d(mat["chanlabels"])]
    data_file = Path(str(mat["setfile"]))
    if not data_file.is_file():
        raise FileNotFoundError(
            f"{mat_path.name} was fitted on {data_file}, which is not there any more -- "
            f"ICLabel needs the data, not just the matrices."
        )

    # Checked before the data is touched: loading a subset of a full night only to
    # find the matrices contradict each other wastes minutes for no reason.
    unmixing = np.asarray(mat["icaweights"], dtype=float) @ np.asarray(mat["icasphere"], dtype=float)
    if not np.allclose(unmixing, np.linalg.pinv(np.asarray(mat["icawinv"], dtype=float)),
                       rtol=1e-5, atol=1e-8):
        raise RuntimeError(
            f"{mat_path.name}: icawinv is not the pseudo-inverse of icaweights @ icasphere "
            f"-- the three matrices disagree, so any labels derived from them would "
            f"describe a decomposition that does not exist."
        )

    t0 = time.perf_counter()
    raw, minutes_used = iclabel_subset(data_file, ch_names, ICLABEL_MINUTES, ICLABEL_CHUNK_SECONDS)
    print(f"  ICLabel: {minutes_used:.1f} min over {len(ch_names)} channels, "
          f"{mat['icaweights'].shape[0]} components")
    ica = _eeglab_style_ica(
        np.asarray(mat["icaweights"], dtype=float),
        np.asarray(mat["icasphere"], dtype=float),
        raw.info,
    )
    probs = iclabel_label_components(raw, ica, inplace=False, backend=ICLABEL_BACKEND)
    elapsed = time.perf_counter() - t0

    top = probs.argmax(axis=1)
    top_prob = probs.max(axis=1)
    is_artefact = np.isin(top, ICLABEL_ARTEFACT_CLASSES) & (top_prob >= ICLABEL_THRESHOLD)

    # Share of back-projected variance: var(activation) scaled by the squared norm of
    # the component's scalp projection. Computed on the ICLabel subset, which is what
    # is loaded here -- a variance ordering, not a claim about the whole night.
    act = ica.unmixing_matrix_ @ ica.pca_components_ @ (raw.get_data() * 1e6)
    backvar = act.var(axis=1) * (np.asarray(mat["icawinv"], dtype=float) ** 2).sum(axis=0)
    varfrac = backvar / backvar.sum()

    names = [ICLABEL_CLASSES[i] for i in top]
    header = ["component", "label", "probability", "status", "varfrac", *ICLABEL_COLUMNS]
    with open(tsv_path, "w", newline="\n") as f:
        f.write("\t".join(header) + "\n")
        for i in range(probs.shape[0]):
            row = [f"IC{i + 1}", names[i], f"{top_prob[i]:.10g}",
                   "bad" if is_artefact[i] else "good", f"{varfrac[i]:.10g}"]
            row += [f"{p:.10g}" for p in probs[i]]
            f.write("\t".join(row) + "\n")
    print(f"  ICLabel: {int(is_artefact.sum())}/{probs.shape[0]} flagged bad -> {tsv_path.name}")
    for cls in range(len(ICLABEL_CLASSES)):
        n = int((top == cls).sum())
        if n:
            print(f"    {ICLABEL_CLASSES[cls]:<14} {n:3d}")

    artefact_names = [ICLABEL_CLASSES[i] for i in ICLABEL_ARTEFACT_CLASSES]
    with open(json_path, "w") as f:
        json.dump(
            {
                "component":   {"Description": "Component identifier, IC<n>, in the row order of icaweights (AMICA components ordered by decreasing back-projected variance)."},
                "label":       {"Description": "Most probable ICLabel class, one of: " + ", ".join(ICLABEL_CLASSES) + "."},
                "probability": {"Description": "Posterior probability of the class named in 'label'."},
                "status":      {"Description": "Whether the component is considered artefactual.",
                                "Levels": {"good": "Retained.",
                                           "bad": f"Top class in [{', '.join(artefact_names)}] with probability >= {ICLABEL_THRESHOLD}."}},
                "varfrac":     {"Description": "Share of total back-projected variance: var(activation) * ||icawinv column||^2, normalised across components, computed on the ICLabel subset. A proxy for pvaf, not pvaf itself.",
                                "Units": "fraction"},
                **{col: {"Description": f'ICLabel posterior probability of class "{cls}".'}
                   for col, cls in zip(ICLABEL_COLUMNS, ICLABEL_CLASSES)},
                "GeneratedBy": {
                    "Name": "run-pamica.py / mne-icalabel",
                    "Description": "ICLabel via mne-icalabel, not the MATLAB EEGLAB plugin. Same network and class order; the two implementations agree closely but are not bit-identical.",
                    "MinutesRequested": ICLABEL_MINUTES,
                    "MinutesClassified": round(minutes_used, 2),
                    "AssumedAverageReference": ICLABEL_ASSUME_AVGREF,
                    "Backend": ICLABEL_BACKEND or "auto",
                    "DurationMinutes": round(elapsed / 60, 3),
                    "GeneratedDate": datetime.now().isoformat(),
                },
            },
            f,
            indent=2,
        )

    # Folded back into the .mat in the exact shape EEGLAB keeps it in, so MATLAB can
    # assign it across without rearranging anything:
    #   EEG.etc.ic_classification = S.ic_classification
    #   EEG.reject.gcompreject    = S.gcompreject
    # scipy writes a nested dict as a nested struct and an object array as a cell
    # array, which is what EEG.etc.ic_classification.ICLabel.classes has to be.
    # ICA/loadica.m does this together with the channel matching; see it for the rest.
    mat.pop("__header__", None); mat.pop("__version__", None); mat.pop("__globals__", None)
    mat["ic_classification"] = {
        "ICLabel": {
            "classes": np.array(ICLABEL_CLASSES, dtype=object),
            "classifications": probs.astype(np.float64),
            # EEGLAB's iclabel() writes 'default' here. Saying which implementation
            # produced it matters: the two are not bit-identical, and this field is
            # the only thing that travels with the numbers into MATLAB.
            "version": "default (mne-icalabel)",
        }
    }
    mat["iclabel_labels"] = np.array(names, dtype=object)
    mat["gcompreject"] = is_artefact.reshape(1, -1)     # EEGLAB wants 1 x ncomp
    mat["varfrac"] = varfrac.reshape(-1, 1)
    sio.savemat(mat_path, mat)
    return probs


if __name__ == "__main__":
    bids = load_bids(BIDS_MAT)
    rawdata = Path(bids["pth"])  # BIDS.pth is <project_root>\rawdata; derivatives live inside it
    deriv_in = rawdata / "derivatives" / DERIV_IN_DIR
    deriv_out = rawdata / "derivatives" / DERIV_OUT_DIR

    # Scan/retry loop for input files another machine may still be writing: each
    # pass checks every still-missing recording once (a miss never blocks checking
    # the next), processes whatever's newly found immediately, then waits
    # WAIT_LOOP_MINUTES before rescanning the remainder. Gives up (raises) once
    # WAIT_MAX_MINUTES pass with no scan turning up anything new.
    pending = deriv_entries(bids, DESC, deriv_in)
    last_progress = time.monotonic()
    while pending:
        found, pending = _scan_pass(pending)
        if found:
            last_progress = time.monotonic()
            for data_file in found:
                final_ll, mat_path = run_amica(data_file, deriv_out)
                if final_ll is not None:
                    print(f"{data_file.stem}: final LL = {final_ll:.5f}")
                # Runs for skipped recordings too: the fit is what is expensive and
                # already-done, the labels may simply never have been made.
                if RUN_ICLABEL and mat_path.is_file():
                    try:
                        run_iclabel(mat_path)
                    except Exception as exc:  # noqa: BLE001
                        # The AMICA output is written and valid at this point; losing
                        # the whole batch over a labelling problem would be worse than
                        # carrying on and relabelling later.
                        print(f"  ICLabel FAILED for {mat_path.name}: {exc!r}")
        if not pending:
            break
        stalled_for = time.monotonic() - last_progress
        if stalled_for >= WAIT_MAX_MINUTES * 60:
            names = ", ".join(base.name for base, _ in pending)
            raise FileNotFoundError(
                f"giving up after {WAIT_MAX_MINUTES} min with no new input files -- "
                f"still missing: {names}"
            )
        print(f"  {len(pending)} file(s) still missing, rechecking in {WAIT_LOOP_MINUTES} min "
              f"(giving up after {WAIT_MAX_MINUTES} min total without progress)")
        time.sleep(WAIT_LOOP_MINUTES * 60)