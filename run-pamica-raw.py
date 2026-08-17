r"""
AMICA (pamica) directly on BIDS "filt" derivative recordings (BrainVision
.vhdr) -- the file run_gedai_bids.m itself reads before GEDAI's stage-specific
artefact rejection. Skips GEDAI entirely and instead does only the same
minimal prep run_gedai_bids.m does ahead of GEDAI: drop the non-EEG channel
range, remove bad channels (from the cached *_badchans.mat), average
re-reference -- then Chebyshev-II high-pass, slim the recording down to
wake/REM/N1 epochs plus every second N2/N3 epoch (using the matched sleep
scoring, to keep AMICA's input to a manageable size), rank projection, AMICA,
and unmixing matrices written in the form EEGLAB expects in
EEG.icaweights / EEG.icasphere / EEG.icawinv.

.venv\Scripts\Activate.ps1
python run-pamica-raw.py
"""

import csv
import json
import re
import time
import xml.etree.ElementTree as ET
from datetime import datetime
from pathlib import Path

import mne
import numpy as np
import scipy.io as sio
import torch
from pamica import AMICA
from scipy.signal import cheb2ord, cheby2, sosfiltfilt

BIDS_MAT      = Path(__file__).parent / "BidsFiles" / "BIDS_DROP.mat"
DERIV_IN_DIR  = "prep-ged"    # derivatives subfolder holding the desc-filt .vhdr files (run_gedai_bids.m's own input) and *_badchans.mat
DERIV_OUT_DIR = "pamica-raw"  # derivatives subfolder to write AMICA output under
DESC          = "filt"        # opts.inputdesc default in run_gedai_bids.m
SUBJECTS      = ["drop0001"]  # None = all subjects
SESSIONS      = ["t1"]  # None = all sessions
TASKS         = ["Sleep", "sleep"]  # None = all tasks

SCORING_SUBDIR   = Path("scoring") / "scores" / "Manual_Checked"  # opts.scoringpath default, under <BIDS root>/derivatives
EPOCH_LENGTH_SEC = 30                # opts.epochlength default in run_gedai_bids.m
NOTEEG_CHANNELS  = range(257, 301)   # opts.noteegchannels default (1-based, MATLAB-inclusive)

STAGE_WAKE, STAGE_REM, STAGE_N1, STAGE_N2, STAGE_N3 = 1, 0, -1, -2, -3  # qol/scoreloader.m's relabeled digits

MAX_ITER = 2000  # EEGLAB-AMICA's usual budget; pamica's fit default is lower
DO_REJECT = False  # Fortran-style outlier rejection; costs ~2x GPU memory (pamica clones the full array
                   # every iteration under do_reject, core.py:1710) -- OOMs on this recording size/GPU
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
    """Map each raw .vhdr recording in the BIDS struct onto its desc-filt derivative
    .vhdr file (run_gedai_bids.m's own input, ahead of GEDAI).

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
            paths.append(deriv_in.joinpath(*folders, f"{file_id}_desc-{desc}_eeg.vhdr"))

    missing = [p for p in paths if not p.is_file()]
    for p in missing:
        print(f"missing: {p}")
    return [p for p in paths if p.is_file()]


def parse_entities(file_id):
    """Recover the BIDS entities dict from a '<key>-<val>_<key>-<val>...' file_id, for
    scoring-file matching (mirrors bids.internal.parse_filename's .entities)."""
    return dict(token.split("-", 1) for token in file_id.split("_") if "-" in token)


def collect_scoring_files(scoring_dir):
    """Port of GEDAI/+gedai/collectScoringFiles.m: xml preferred, then json, then csv
    (excluding *events.csv)."""
    base = Path(scoring_dir)
    if not base.is_dir():
        return []
    for pattern in ("*.xml", "*.json", "*.csv"):
        files = sorted(base.rglob(pattern))
        if pattern == "*.csv":
            files = [f for f in files if not f.name.lower().endswith("events.csv")]
        if files:
            return files
    return []


def match_scoring_file(entities, candidates):
    """Port of GEDAI/+gedai/matchScoringFile.m: first candidate whose path contains
    every entity value (case-insensitive substring, after its 'key-' marker), ignoring
    the acq entity."""
    fields = [k for k in entities if "acq" not in k.lower()]
    for candidate in candidates:
        text = str(candidate)
        ok = True
        for key in fields:
            m = re.search(rf"(?<={re.escape(key)}-)[^_]+", text, re.IGNORECASE)
            if not (m and entities[key].lower() in m.group(0).lower()):
                ok = False
                break
        if ok:
            return candidate
    return None


_SCOREMAP = {5: 0, 0: 1, 1: -1, 2: -2, 3: -3}  # qol/scoreloader.m's scoremap (From -> To)


def _relabel(digits):
    """Apply qol/scoreloader.m's scoremap (raw stage code -> Wake=1/REM=0/N1=-1/N2=-2/N3=-3)."""
    out = digits.astype(float).copy()
    finite = ~np.isnan(out)
    out[finite] = [_SCOREMAP.get(int(v), v) for v in out[finite]]
    return out


def _read_csv_scoring(path):
    """First column of a Sleeptrip-style scoring .csv; non-numeric rows (header) are
    skipped, matching MATLAB readmatrix's automatic header detection."""
    digits = []
    with open(path, newline="") as f:
        for row in csv.reader(f):
            if not row:
                continue
            try:
                digits.append(float(row[0]))
            except ValueError:
                continue
    return np.array(digits)


def _read_json_scoring(path):
    """Scoringhero-style .json ({"digit": [...]} nested one level in a list); already in
    the target Wake=1/REM=0/N1=-1/N2=-2/N3=-3 encoding, no scoremap relabeling."""
    with open(path) as f:
        data = json.load(f)
    return np.array(data[0]["digit"], dtype=float)


def _read_xml_scoring(path):
    """Compumedics/Nox EventExport .xml: SLEEP-* events -> per-epoch hypnogram, port of
    scoreloader.m's extract_sleep_stages/stageCode."""
    stage_code = {"SLEEP-S0": 0, "SLEEP-S1": 1, "SLEEP-S2": 2, "SLEEP-S3": 3, "SLEEP-S4": 3, "SLEEP-REM": 5}
    fmt = "%Y-%m-%dT%H:%M:%S.%f"
    starts, codes = [], []
    for ev in ET.parse(path).getroot().iter("Event"):
        type_el = ev.find("Type")
        if type_el is None or type_el.text not in stage_code:
            continue
        starts.append(datetime.strptime(ev.find("StartTime").text, fmt))
        codes.append(stage_code[type_el.text])
    if not starts:
        raise ValueError(f"no SLEEP-* events found in {path}")
    order = np.argsort(starts)
    starts = [starts[i] for i in order]
    codes = np.array(codes)[order]
    epoch_sec = np.median([(b - a).total_seconds() for a, b in zip(starts[:-1], starts[1:])])
    epoch_idx = np.round([(t - starts[0]).total_seconds() / epoch_sec for t in starts]).astype(int)
    digits = np.full(epoch_idx.max() + 1, np.nan)
    digits[epoch_idx] = codes
    return digits


def scoreloader(path):
    """Port of qol/scoreloader.m: load a sleep-scoring file (.csv/.json/.xml) into
    per-30s-epoch stage digits (Wake=1, REM=0, N1=-1, N2=-2, N3=-3)."""
    path = Path(path)
    ext = path.suffix.lower()
    if ext == ".csv":
        return _relabel(_read_csv_scoring(path))
    if ext == ".json":
        return _read_json_scoring(path)
    if ext == ".xml":
        return _relabel(_read_xml_scoring(path))
    raise ValueError(f"unsupported scoring file type: {ext} ({path})")


def stage_per_epoch(digits, n_epochs):
    """Per-epoch stage digit for epochs 0..n_epochs-1 (NaN where the scoring is shorter
    than the recording), mirroring run_gedai_bids.m:159-160's length reconciliation."""
    out = np.full(n_epochs, np.nan)
    n = min(len(digits), n_epochs)
    out[:n] = digits[:n]
    return out


def select_slim_epochs(stage):
    """Keep every Wake/REM/N1 epoch, every second N2 epoch and every second N3 epoch
    (the 1st, 3rd, 5th, ... of each); drop everything else, including any unscored tail."""
    keep = np.isin(stage, (STAGE_WAKE, STAGE_REM, STAGE_N1))
    for s in (STAGE_N2, STAGE_N3):
        idx = np.flatnonzero(stage == s)
        keep[idx[::2]] = True
    return keep


def apply_epoch_mask(X, sfreq, keep_mask, epoch_len_sec=EPOCH_LENGTH_SEC):
    """Concatenate the kept epochs' time samples in original order; samples past the
    last full epoch are dropped."""
    ep_samples = int(round(epoch_len_sec * sfreq))
    segments = [X[:, i * ep_samples:(i + 1) * ep_samples] for i in np.flatnonzero(keep_mask)]
    kept = np.concatenate(segments, axis=1) if segments else X[:, :0]
    print(f"  epoch selection: kept {int(keep_mask.sum())}/{len(keep_mask)} epochs"
          f" ({kept.shape[1] / sfreq / 60:.1f} / {X.shape[1] / sfreq / 60:.1f} min)")
    return kept


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


def drop_noneeg_channels(labels, X):
    """Mirror run_gedai_bids.m:163's pop_select(EEG,'nochannel',intersect(1:nbchan,257:300))
    -- badchans.mat's removed_channels mask is indexed against what remains after this."""
    drop = {i - 1 for i in NOTEEG_CHANNELS if i - 1 < len(labels)}
    keep = [i for i in range(len(labels)) if i not in drop]
    print(f"  dropped {len(drop)} non-EEG channel(s)"
          f" (index range {NOTEEG_CHANNELS.start}-{NOTEEG_CHANNELS.stop - 1}, 1-based)")
    return [labels[i] for i in keep], X[keep]


def load_bad_channel_mask(badchans_path, n_channels):
    """Load run_gedai_bids.m's cached clean_channels() output
    (badchans.mat: removed_channels/corr/znoise)."""
    if not badchans_path.is_file():
        raise FileNotFoundError(f"badchans file not found: {badchans_path}")
    mat = sio.loadmat(badchans_path, simplify_cells=True)
    removed = np.atleast_1d(mat["removed_channels"]).astype(bool).ravel()
    if removed.size != n_channels:
        raise ValueError(
            f"{badchans_path.name}: removed_channels has {removed.size} entries, expected {n_channels}"
        )
    return removed


def average_reference(X):
    """EEGLAB-side average reference, exactly run_gedai_bids.m:218's
    EEG.data - sum(EEG.data,1)/(size(EEG.data,1)+1) -- note the +1 divisor, kept as-is."""
    return X - X.sum(axis=0, keepdims=True) / (X.shape[0] + 1)


def run_amica(data_file, out_dir, scoring_candidates):
    raw = mne.io.read_raw_brainvision(data_file, preload=True)
    labels = list(raw.ch_names)
    sfreq = raw.info["sfreq"]
    X = raw.get_data() * 1e6  # MNE volts -> EEGLAB microvolts
    del raw

    print(f"{data_file.stem}: {len(labels)} channels")
    labels, X = drop_noneeg_channels(labels, X)

    file_id = data_file.stem.split("_desc-")[0]
    badchans_path = data_file.with_name(f"{file_id}_badchans.mat")
    removed = load_bad_channel_mask(badchans_path, len(labels))
    print(f"  {int(removed.sum())} bad channel(s) removed (of {len(labels)}), per {badchans_path.name}")
    labels = [l for l, bad in zip(labels, removed) if not bad]
    X = X[~removed]

    X = average_reference(X)
    X = highpass(X, sfreq)

    scoring_file = match_scoring_file(parse_entities(file_id), scoring_candidates)
    if scoring_file is None:
        raise RuntimeError(f"no scoring file matched for {file_id} -- required to slim the recording")
    print(f"  scoring -> {scoring_file}")
    digits = scoreloader(scoring_file)
    n_ep = int(X.shape[1] // (EPOCH_LENGTH_SEC * sfreq))
    stage = stage_per_epoch(digits, n_ep)
    keep_mask = select_slim_epochs(stage)
    X = apply_epoch_mask(X, sfreq, keep_mask)

    X -= X.mean(axis=1, keepdims=True)
    P = rank_projection(X)  # kept in float64: cheap, and wants a clean rank cut
    Xr = (P @ X).astype(np.float32)
    del X

    model = AMICA(n_models=1, n_mix=3, device="cuda")
    t0 = time.perf_counter()
    model.fit(Xr, max_iter=MAX_ITER, block_size=8192, dtype=torch.float32, do_reject=DO_REJECT)
    elapsed = time.perf_counter() - t0
    print(f"  AMICA fit took {elapsed / 60:.1f} min ({elapsed:.1f} s)")

    A = model.get_mixing_matrix()
    # Confirm get_mixing_matrix is sensor space, i.e. sources == pinv(A) @ data.
    ref, chk = model.transform(Xr[:, :5000]), np.linalg.pinv(A) @ Xr[:, :5000]
    r = np.array([np.corrcoef(a, b)[0, 1] for a, b in zip(ref, chk)])
    if np.abs(r).min() < 0.999:
        print(f"WARNING {data_file.stem}: mixing-matrix convention check failed ({np.abs(r).min():.3f})")

    A = A[:, model.variance_order()]  # EEGLAB order: IC1 = highest variance
    out_dir.mkdir(parents=True, exist_ok=True)
    model.write_amica_output(str(out_dir / "amicaout"))  # NB: in the rank-reduced, epoch-slimmed space
    sio.savemat(
        out_dir / "amica_eeglab.mat",
        {
            "icaweights": np.linalg.pinv(A),
            "icasphere": P,
            "icawinv": P.T @ A,
            "chanlabels": np.array(labels, dtype=object),
            "setfile": str(data_file),
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
    scoring_dir = rawdata / "derivatives" / SCORING_SUBDIR
    scoring_candidates = collect_scoring_files(scoring_dir)
    if not scoring_candidates:
        print(f"warning: no scoring files found under {scoring_dir}")

    for data_file in deriv_paths(bids, DESC, deriv_in):
        out_dir = deriv_out / data_file.stem
        if (out_dir / "amica_eeglab.mat").is_file():
            print(f"skip {data_file.stem} (already done)")
            continue
        print(f"{data_file.stem}: final LL = {run_amica(data_file, out_dir, scoring_candidates):.5f}")
