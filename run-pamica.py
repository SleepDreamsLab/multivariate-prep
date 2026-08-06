r"""
AMICA (pamica) on BIDS derivative EEGLAB .set files.

Per recording: import -> Chebyshev-II high-pass -> rank projection (GEDAI may
have left the data rank deficient) -> AMICA -> unmixing matrices written in the
form EEGLAB expects in EEG.icaweights / EEG.icasphere / EEG.icawinv.

.venv\Scripts\Activate.ps1
python run-pamica.py
"""

import json
import time
from pathlib import Path

import mne
import numpy as np
import scipy.io as sio
from pamica import AMICA
from scipy.signal import cheb2ord, cheby2, sosfiltfilt

BIDS_MAT      = Path(__file__).parent / "BidsFiles" / "BIDS_DROP.mat"
DERIV_IN_DIR  = "prep-ged"  # derivatives subfolder to read the desc-* .set files from
DERIV_OUT_DIR = "pamica"    # derivatives subfolder to write AMICA output under
DESC          = "zc2gedWakeBBAuto"
SUBJECTS      = ["drop0001"]  # None = all subjects
SESSIONS      = ["t1"]  # None = all sessions
TASKS         = ["Sleep", "sleep"]  # None = all tasks

MAX_ITER = 2000  # EEGLAB-AMICA's usual budget; pamica's fit default is lower
PASS_FRQ, STOP_FRQ, PASS_RIPPLE, STOP_ATTEN = 1.6, 0.8, 0.05, 30


def _unwrap(x):
    """Peel the 1x1 cell / array wrappers left by loadmat."""
    while not isinstance(x, dict):
        x = np.atleast_1d(x)[0]
    return x


def load_bids(mat_path):
    """Load the BIDS layout struct dumped from MATLAB (bids-matlab's BIDS object)."""
    mat = sio.loadmat(mat_path, simplify_cells=True)
    return _unwrap(mat[next(k for k in mat if not k.startswith("__"))])


def deriv_paths(bids, desc, deriv_in):
    """Map each raw .vhdr recording in the BIDS struct onto its prep-ged .set file.

    Mirrors the MATLAB filtering (bids.query(..., 'extension', '.vhdr')) then
    bids.internal.parse_filename entity join, using the ext/entities fields
    that BIDS_DROP.mat already carries per recording instead of re-deriving
    them from the filename.
    """
    paths = []
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
            paths.append(deriv_in.joinpath(*folders, f"{file_id}_desc-{desc}_eeg.set"))

    missing = [p for p in paths if not p.is_file()]
    for p in missing:
        print(f"missing: {p}")
    return [p for p in paths if p.is_file()]


def highpass(X, sfreq):
    """Minimum-order Chebyshev-II high-pass, zero phase (= designfilt + filtfilt)."""
    order, wn = cheb2ord(PASS_FRQ, STOP_FRQ, PASS_RIPPLE, STOP_ATTEN, fs=sfreq)
    print(f"  high-pass: order {order}, Wn {wn:.4g} (pass {PASS_FRQ} Hz / stop {STOP_FRQ} Hz @ {sfreq} Hz)")
    sos = cheby2(order, STOP_ATTEN, wn, btype="highpass", output="sos", fs=sfreq)
    t0 = time.perf_counter()
    Xf = sosfiltfilt(sos, X, axis=1)
    print(f"  high-pass: filtered {X.shape[0]} channels x {X.shape[1]} samples in {time.perf_counter() - t0:.1f} s")
    return Xf


def rank_projection(X, tol=1e-7):
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


def run_amica(set_file, out_dir):
    raw = mne.io.read_raw_eeglab(set_file, preload=True)
    picks = mne.pick_types(raw.info, eeg=True, exclude="bads")
    if raw.info["bads"]:
        print(f"  {len(raw.info['bads'])} channel(s) flagged bad in the file, excluded: {raw.info['bads']}")
    labels = [raw.ch_names[i] for i in picks]
    sfreq = raw.info["sfreq"]
    X = raw.get_data(picks=picks) * 1e6  # MNE volts -> EEGLAB microvolts
    del raw

    print(f"{set_file.stem}: {len(labels)} channels")
    X = highpass(X, sfreq)
    X -= X.mean(axis=1, keepdims=True)
    P = rank_projection(X)
    Xr = P @ X
    del X

    model = AMICA(n_models=1, n_mix=3, device="cuda")
    t0 = time.perf_counter()
    model.fit(Xr, max_iter=MAX_ITER)
    elapsed = time.perf_counter() - t0
    print(f"  AMICA fit took {elapsed / 60:.1f} min ({elapsed:.1f} s)")

    A = model.get_mixing_matrix()
    # Confirm get_mixing_matrix is sensor space, i.e. sources == pinv(A) @ data.
    ref, chk = model.transform(Xr[:, :5000]), np.linalg.pinv(A) @ Xr[:, :5000]
    r = np.array([np.corrcoef(a, b)[0, 1] for a, b in zip(ref, chk)])
    if np.abs(r).min() < 0.999:
        print(f"WARNING {set_file.stem}: mixing-matrix convention check failed ({np.abs(r).min():.3f})")

    A = A[:, model.variance_order()]  # EEGLAB order: IC1 = highest variance
    out_dir.mkdir(parents=True, exist_ok=True)
    model.write_amica_output(str(out_dir / "amicaout"))  # NB: in the rank-reduced space
    sio.savemat(
        out_dir / "amica_eeglab.mat",
        {
            "icaweights": np.linalg.pinv(A),
            "icasphere": P,
            "icawinv": P.T @ A,
            "chanlabels": np.array(labels, dtype=object),
            "setfile": str(set_file),
            "srate": sfreq,
            "final_ll": model.final_ll_,
        },
    )
    with open(out_dir / "amica_runtime.json", "w") as f:
        json.dump({"elapsed_seconds": elapsed}, f, indent=2)
    return model.final_ll_


if __name__ == "__main__":
    bids = load_bids(BIDS_MAT)
    rawdata = Path(bids["pth"])  # BIDS.pth is <project_root>\rawdata; derivatives live inside it
    deriv_in = rawdata / "derivatives" / DERIV_IN_DIR
    deriv_out = rawdata / "derivatives" / DERIV_OUT_DIR

    for set_file in deriv_paths(bids, DESC, deriv_in):
        out_dir = deriv_out / set_file.stem
        if (out_dir / "amica_eeglab.mat").is_file():
            print(f"skip {set_file.stem} (already done)")
            continue
        print(f"{set_file.stem}: final LL = {run_amica(set_file, out_dir):.5f}")