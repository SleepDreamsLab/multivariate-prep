r"""
AMICA (pamica) on BIDS derivative recordings, either EEGLAB .set or BrainVision
.vhdr files.

Per recording: import -> Chebyshev-II high-pass -> rank projection (GEDAI may
have left the data rank deficient) -> AMICA -> unmixing matrices written in the
form EEGLAB expects in EEG.icaweights / EEG.icasphere / EEG.icawinv.

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


MAX_ITER = 600  # EEGLAB-AMICA's usual budget; pamica's fit default is lower
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


def read_raw(data_file):
    """Load a derivative recording with the MNE reader matching its extension."""
    if data_file.suffix == ".set":
        return mne.io.read_raw_eeglab(data_file, preload=True)
    if data_file.suffix == ".vhdr":
        return mne.io.read_raw_brainvision(data_file, preload=True)
    raise ValueError(f"unsupported derivative extension: {data_file.suffix} ({data_file})")


def run_amica(data_file, deriv_out):
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
        return None

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
    return model.final_ll_


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
                final_ll = run_amica(data_file, deriv_out)
                if final_ll is not None:
                    print(f"{data_file.stem}: final LL = {final_ll:.5f}")
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