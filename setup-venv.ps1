<#
.SYNOPSIS
    Sets up the Python environment for multivariate-prep (run-pamica.py and friends).

.DESCRIPTION
    - Installs uv (the Python package/venv manager) if it isn't already on PATH.
    - Clones the pAMICA repo as a sibling of this repo if it isn't already there
      (pyproject.toml pulls it in as an editable dependency via a relative path,
      so it must live at ..\pAMICA relative to this repo).
    - Runs `uv sync` to create .venv and install mne/numpy/scipy/pamica exactly
      as pinned in uv.lock.
    - Installs the CPU build of torch by default. Pass -Cuda on a machine that
      has an NVIDIA GPU to install the CUDA 13.0 build instead (pyproject.toml's
      pytorch-cu130 index) so pamica's AMICA fit actually runs on the GPU.

.PARAMETER Cuda
    Install the CUDA build of torch (`uv sync --extra cuda`) instead of the
    CPU build. Only useful on a machine with an NVIDIA GPU and a driver new
    enough for CUDA 13 (check `nvidia-smi`'s "CUDA Version" field).

.NOTES
    Run from anywhere; paths are resolved relative to this script's location, not cwd.
#>

param(
    [switch]$Cuda
)

$ErrorActionPreference = "Stop"

$repoRoot   = Split-Path -Parent $MyInvocation.MyCommand.Path
$pamicaRoot = Join-Path (Split-Path -Parent $repoRoot) "pAMICA"
$pamicaUrl  = "https://github.com/sccn/pAMICA.git"

# 1. uv itself
if (-not (Get-Command uv -ErrorAction SilentlyContinue)) {
    Write-Host "uv not found, installing..."
    powershell -ExecutionPolicy ByPass -Command "irm https://astral.sh/uv/install.ps1 | iex"
    $env:Path = "$env:USERPROFILE\.local\bin;$env:Path"
    if (-not (Get-Command uv -ErrorAction SilentlyContinue)) {
        throw "uv install finished but 'uv' still isn't on PATH. Open a new shell and re-run this script."
    }
}
Write-Host "uv: $(uv --version)"

# 2. sibling pAMICA checkout (editable dependency, see pyproject.toml [tool.uv.sources])
if (-not (Test-Path (Join-Path $pamicaRoot "pyproject.toml"))) {
    Write-Host "pAMICA not found at $pamicaRoot, cloning..."
    git clone $pamicaUrl $pamicaRoot
} else {
    Write-Host "pAMICA found at $pamicaRoot"
}

# 3. sync this repo's venv (.venv) from pyproject.toml / uv.lock
if ($Cuda) {
    if (-not (Get-Command nvidia-smi -ErrorAction SilentlyContinue)) {
        Write-Warning "nvidia-smi not found -- no NVIDIA driver detected on this machine. Continuing with -Cuda anyway, but torch.cuda.is_available() will likely come back False."
    }
    Write-Host "Syncing .venv for multivariate-prep (CUDA build of torch)..."
} else {
    Write-Host "Syncing .venv for multivariate-prep (CPU build of torch)..."
}
Push-Location $repoRoot
try {
    if ($Cuda) {
        uv sync --extra cuda
    } else {
        uv sync --extra cpu
    }
} finally {
    Pop-Location
}

Write-Host ""
Write-Host "Done. Activate with:"
Write-Host "  $repoRoot\.venv\Scripts\Activate.ps1"
Write-Host "or run scripts directly with:"
Write-Host "  uv run --project `"$repoRoot`" run-pamica.py"
