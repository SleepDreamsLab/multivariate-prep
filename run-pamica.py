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
import threading
import time
from pathlib import Path

import mne
import numpy as np
import scipy.io as sio
import torch
from pamica import AMICA
from scipy.signal import cheb2ord, cheby2, sosfiltfilt

BIDS_MAT      = Path(__file__).parent / "BidsFiles" / "BIDS_DROP.mat"
DERIV_IN_DIR  = "prep-ged"  # derivatives subfolder to read the desc-* .set files from
DERIV_OUT_DIR = "pamica"    # derivatives subfolder to write AMICA output under
DESC          = "zc2gedWakeBBAutoplusFSAutoPlus" # zc
SUBJECTS      = ["drop0001"]  # None = all subjects
SESSIONS      = ["t1"]  # None = all sessions
TASKS         = ["Sleep", "sleep"]  # None = all tasks


MAX_ITER = 1000  # EEGLAB-AMICA's usual budget; pamica's fit default is lower
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
DTYPE      = torch.float64  # AMICA: compute dtype
DEVICE     = "cuda"         # AMICA: compute device
MEM_MONITOR_INTERVAL = 60   # seconds between GPU memory polls during fit(); 0/None disables


def _unwrap(x):
    """Peel the 1x1 cell / array wrappers left by loadmat."""
    while not isinstance(x, dict):
        x = np.atleast_1d(x)[0]
    return x


def load_bids(mat_path):
    """Load the BIDS layout struct dumped from MATLAB (bids-matlab's BIDS object)."""
    mat = sio.loadmat(mat_path, simplify_cells=True)
    return _unwrap(mat[next(k for k in mat if not k.startswith("__"))])


DATA_EXTENSIONS = (".set", ".vhdr")  # preference order when both exist for a recording


def deriv_paths(bids, desc, deriv_in):
    """Map each raw .vhdr recording in the BIDS struct onto its prep-ged derivative file.

    Mirrors the MATLAB filtering (bids.query(..., 'extension', '.vhdr')) then
    bids.internal.parse_filename entity join, using the ext/entities fields
    that BIDS_DROP.mat already carries per recording instead of re-deriving
    them from the filename.

    The derivative can be either an EEGLAB .set or a BrainVision .vhdr file;
    both extensions are tried per recording (in DATA_EXTENSIONS order), and
    run_amica() picks the matching MNE reader off the extension it gets back.
    """
    paths = []
    missing = []
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
            found = [c for c in candidates if c.is_file()]
            if not found:
                missing.append(base)
                continue
            if len(found) > 1:
                print(f"both {' and '.join(c.suffix for c in found)} present for {base.name}, using {found[0].suffix}")
            paths.append(found[0])

    for base in missing:
        print(f"missing: {base} ({'/'.join(DATA_EXTENSIONS)})")
    return paths


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


class GPUMemoryMonitor:
    """Background thread that periodically prints torch's own CUDA memory
    counters while a blocking call (like AMICA's fit()) runs on the main
    thread -- fit() has no per-iteration hook to poll from, so this is the
    only way to see peak allocation over the course of a fit rather than
    only its final value. Reports this process's own allocation, unlike
    nvidia-smi which shows total GPU usage across all processes.
    """

    def __init__(self, interval=MEM_MONITOR_INTERVAL, label=""):
        self.interval = interval
        self.label = label
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._run, daemon=True)

    def _run(self):
        while not self._stop.wait(self.interval):
            allocated = torch.cuda.memory_allocated() / 1e9
            peak = torch.cuda.max_memory_allocated() / 1e9
            reserved = torch.cuda.memory_reserved() / 1e9
            print(f"  [mem{self.label}] allocated={allocated:.2f} GB peak={peak:.2f} GB reserved={reserved:.2f} GB")

    def __enter__(self):
        if self.interval:
            self._thread.start()
        return self

    def __exit__(self, *exc_info):
        self._stop.set()
        if self.interval:
            self._thread.join()


def read_raw(data_file):
    """Load a derivative recording with the MNE reader matching its extension."""
    if data_file.suffix == ".set":
        return mne.io.read_raw_eeglab(data_file, preload=True)
    if data_file.suffix == ".vhdr":
        return mne.io.read_raw_brainvision(data_file, preload=True)
    raise ValueError(f"unsupported derivative extension: {data_file.suffix} ({data_file})")


def run_amica(data_file, out_dir):
    raw = read_raw(data_file)
    picks = mne.pick_types(raw.info, eeg=True, exclude="bads")
    if raw.info["bads"]:
        print(f"  {len(raw.info['bads'])} channel(s) flagged bad in the file, excluded: {raw.info['bads']}")
    labels = [raw.ch_names[i] for i in picks]
    sfreq = raw.info["sfreq"]
    X = raw.get_data(picks=picks) * 1e6  # MNE volts -> EEGLAB microvolts
    del raw

    print(f"{data_file.stem}: {len(labels)} channels")
    X = highpass(X, sfreq)
    X -= X.mean(axis=1, keepdims=True)
    P = rank_projection(X)  # kept in float64: cheap, and wants a clean rank cut
    Xr = (P @ X).astype(np.float32)
    del X

    model = AMICA(n_models=N_MODELS, n_mix=N_MIX, device=DEVICE)
    t0 = time.perf_counter()
    with GPUMemoryMonitor(label=f" {data_file.stem}"):
        model.fit(Xr, max_iter=MAX_ITER, block_size=BLOCK_SIZE, dtype=DTYPE, do_reject=DO_REJECT, do_newton=DO_NEWTON)
    elapsed = time.perf_counter() - t0
    print(f"  AMICA fit took {elapsed / 60:.1f} min ({elapsed:.1f} s)")

    out_dir.mkdir(parents=True, exist_ok=True)
    # Always written first: byte-identical to the Fortran reference, unaffected by
    # anything below, and the fallback if the .mat convention check fails.
    model.write_amica_output(str(out_dir / "amicaout"))  # NB: in the rank-reduced space

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
        out_dir / "amica_eeglab.mat",
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
    with open(out_dir / "amica_runtime.json", "w") as f:
        json.dump(
            {
                "elapsed_seconds": elapsed,
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

    for data_file in deriv_paths(bids, DESC, deriv_in):
        out_dir = deriv_out / data_file.stem
        if (out_dir / "amica_eeglab.mat").is_file():
            print(f"skip {data_file.stem} (already done)")
            continue
        print(f"{data_file.stem}: final LL = {run_amica(data_file, out_dir):.5f}")