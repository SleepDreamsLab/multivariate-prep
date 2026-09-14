# CLAUDE.md
This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo is
Artefact subtraction for high-density sleep EEG: usually 256 channels, sometimes 128, usually at 250 Hz. The methods are
*multivariate*: they exploit spatial oversampling instead of treating electrodes independently (ICA, GED, GEDAI).
Data are organised as BIDS and outputs go to BIDS derivatives. All projects are in active development (DM, DROP, ERCP,
KCs, PM); **DROP is the current focus**.

Typical tasks: bug fixes, performance optimisation, reducing RAM use, new features.

## Start here
1. `Pipelines/PrepPipelineDROP.m` is the reference pipeline. It defines the step order and the parameters for each step.
   The other `Pipelines/PrepPipeline*.m` scripts follow the same pattern.
2. Open the `bidsfun_*` stage that the task touches. Each stage walks every recording in the BIDS layout:
   - `GEDAI/bidsfun_detect_badchans.m`: flat and bad channels (`desc` hp)
   - `GEDAI/bidsfun_hp_zap_cleanline.m`: high-pass, bad-channel removal, zapline-plus, cleanline (`desc` hpzc)
   - `GEDAI/bidsfun_gedai.m`: average reference plus stage-specific GEDAI (`desc` hpzcged; helpers in `+gedai`, `+run`)
   - `ICA/run-pamica.py`: AMICA and ICLabel in Python, **GPU machine only** (`desc` pamica, writes `*_ica.mat`)
   - `ICA/bidsfun_iclabel.m`: MATLAB ICLabel on that `*_ica.mat`. Use this *or* the Python ICLabel pass for a recording, never both.
   - `GEDAI/bidsfun_evalfigs.m` (plots in `+evalplots`), `GEDAI/bidsfun_check_outputs.m`: evaluation and QC
   - `Leadfield/`, `GED/`, `SleepOsci/`, `viewer/`: leadfields, GED, oscillation analyses, data viewer
3. Shared helpers: `qol/` (e.g. `bidswizard`, `smartcache`, `claimFile`, `sidecarjson`, `scoreloader`).
   `patches/` holds modified copies of EEGLAB, zapline, and cleanline functions that shadow the originals on the path.
4. `dependancies.m` puts everything on the path. External repos are **sibling checkouts** in `..` (cloned by
   `setup_dependencies.ps1`). `colormaps/` is vendored (slanCM); don't edit it.

## Stage function contract (required for new and modified stages)
- Signature `failures = bidsfun_<name>(BIDS, opts)` with a validated `arguments` block. Shared options come in via the
  pipeline's `defaults{:}` (`subjectfilter`, `sessionfilter`, `derivfolder`, `acqlabel`, `refresh`).
- Outputs are BIDS derivatives named `<fileID>_desc-<desc>_<suffix>.<ext>`, each with a JSON sidecar (`sidecarjson`).
- Respect `refresh`: skip a recording whose outputs already exist unless `refresh` is true.
- **Parallel machines.** Several machines run the same pipeline over the same share at once. Claim each recording
  with `claimFile` before doing any work, and release the claim on success, on error, and on Ctrl-C.
- Wrap each recording in try/catch, collect failures, and return them. One bad recording must never stop the batch.

## Environment
- MATLAB **R2025b or newer** (`C:\Program Files\MATLAB\R2025b\bin\matlab.exe`). Data live on the Windows share
  `\\vs03.herseninstituut.knaw.nl` (DROP: `VS03-SandC-1\data\nin\data-drop\rawdata`).
- **Batch machines have 64 GB RAM.** A 2 TB machine exists, but code must fit in 64 GB for 256-channel overnight
  recordings. Think about peak memory (chunking, `single`, clearing large intermediates, `ramsaver` options,
  `+utils/freemem`, pool size from `gedai.autoPoolSize`).
- Python: always `uv` with `.venv`, set up via `setup-venv.ps1` (`-Cuda` on GPU machines). `pamica` is an editable
  install from `../pAMICA`.
- `dependancies.m` is shared by all machines. Don't add machine-specific paths. Editing it is fine when needed.

## Testing
- **Test only on `sub-drop0001`, `ses-t1`.** Agents may not access other subjects' files. Run the relevant stage with
  `'subjectfilter', {'sub-drop0001'}, 'sessionfilter', {'ses-t1'}`. Write test outputs to a separate `derivfolder`
  (e.g. `prep-test-claude`) so production derivatives aren't overwritten, unless told otherwise.
- Where `sub-drop0001`'s `ses-t1` files actually sit (all under `\\vs03.herseninstituut.knaw.nl\VS03-SandC-1\data\nin\data-drop`,
  matching the default/pipeline paths in `bidsfun_gedai_benchmark.m` and `Pipelines/PrepPipelineDROP.m`):
  - EEG (raw, BrainVision): `rawdata\sub-drop0001\ses-t1\eeg\sub-drop0001_ses-t1_task-sleep_run-01_eeg.vhdr` (+ `.eeg`/`.vmrk`)
  - Scoring (`derivatives\scores\final`, the `scoringpath` the pipeline scripts pass):
    `derivatives\scores\final\sub-drop0001_ses-t1_task-sleep_run-01_eeg.csv`
  - Leadfield (Brainstorm head model, on the separate `VS03-SandD-4` share, the pipeline's `leadfieldpath`):
    `\\vs03.herseninstituut.knaw.nl\VS03-SandD-4\PM\Data_Analysis\Brainstorm_db\DROP_Leadfields2\data\sub-drop0001\ses-t1\headmodel_surf_openmeeg.mat`
  - SFP montage (resolved by `gedai.matchSfpFile`'s DROP branch, `sourcedata\gps`, not under `rawdata`):
    `sourcedata\gps\sub-drop0001\solved\sub-drop0001_ses-t1_task-sleep_run-01_acq-domesolved_eeg.sfp`
- Keep figures invisible in verification runs. They pop up over the user's work, even under `-batch`. Start every run
  with `set(0, 'DefaultFigureVisible', 'off')`; `exportgraphics` and `saveas` still work. Run long jobs in the background.
- Syntax checks: `checkcode` for MATLAB, `.venv\Scripts\python.exe -m py_compile ICA\run-pamica.py` for Python.

## Coding style
- Follow the [MATLAB Coding Guidelines](https://github.com/mathworks/MATLAB-Coding-Guidelines/blob/main/MATLAB-Coding-Guidelines.md):
  lowerCamelCase variables and functions (verb phrases), 4-space indent, 120-character lines, `arguments` blocks,
  preallocation, an `otherwise` in every `switch`, no `eval`/`assignin`/globals, and no float equality tests.
- Apply the guidelines to code you write or touch. Don't mass-reformat existing files, and don't rename existing
  name-value options that pipeline scripts pass.
- Match the existing help-block style: a detailed header that documents every option, its default, and the
  input and output paths.

## Git
- Build a new worktree for your work. Commit when a change looks final. **Never push.** Ask before creating a new branch.
- Editing sibling repos (`../GEDAI-master`, `../pAMICA`, `../eeg-oscillations`, etc.) is allowed when needed. Put fixes
  to third-party EEGLAB-ecosystem functions in `patches/` rather than the upstream checkout.

## Known pitfalls
- **RAM:** out-of-memory crashes happen on 64 GB machines. Treat memory as a first-class constraint in every change.
- **Plotting:** figures have come out as black images, and MATLAB has crashed during plotting (cause unknown). Close
  figures after saving, and be careful with large and parallel or background figure rendering.
- **GEDAI branch:** `../GEDAI-master` has two remotes: `origin` is the fork (SvennoNito) and `upstream` is the
  original (neurotuning). The pipeline runs on the fork's newest branch. As of 2026-09-14 that is
  `experiment/whole-night-physiology`, which builds on `sleep-fast`; `dependancies.m` and `setup_dependencies.ps1`
  still name `sleep-fast`. Branch names change, so don't rely on a prefix. First run
  `git -C ../GEDAI-master fetch --all`, then
  `git -C ../GEDAI-master for-each-ref --sort=-committerdate refs/remotes/origin`, and use the newest `origin` branch.
  Confirm with the user if it's unclear. `upstream` moves independently; don't merge from it unless asked.
- **ICLabel (MATLAB)** needs the sibling `ICLabel` cloned with submodules (matconvnet) and the sibling `firfilt`.
