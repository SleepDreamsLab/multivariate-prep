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
    - Installs the CUDA build of torch by default (pyproject.toml's
      pytorch-cu126 index) so pamica's AMICA fit runs on the GPU. Pass -Cpu on
      a machine without an NVIDIA GPU to install the CPU build instead.

.PARAMETER Cpu
    Install the CPU build of torch (`uv sync --extra cpu`) instead of the
    CUDA build. Use this on a machine without an NVIDIA GPU, or whose driver
    is too old for CUDA 12.6 (check `nvidia-smi`'s "CUDA Version" field).

.NOTES
    Run from anywhere; paths are resolved relative to this script's location, not cwd.
#>

param(
    [switch]$Cpu
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
if ($Cpu) {
    Write-Host "Syncing .venv for multivariate-prep (CPU build of torch)..."
} else {
    if (-not (Get-Command nvidia-smi -ErrorAction SilentlyContinue)) {
        Write-Warning "nvidia-smi not found -- no NVIDIA driver detected on this machine. Continuing with the CUDA build anyway, but torch.cuda.is_available() will likely come back False. Pass -Cpu to install the CPU build instead."
    }
    Write-Host "Syncing .venv for multivariate-prep (CUDA build of torch)..."
}
Push-Location $repoRoot
try {
    if ($Cpu) {
        uv sync --extra cpu
    } else {
        uv sync --extra cuda
    }
} finally {
    Pop-Location
}

Write-Host ""
Write-Host "Done. Activate with:"
Write-Host "  $repoRoot\.venv\Scripts\Activate.ps1"
Write-Host "or run scripts directly with:"
Write-Host "  uv run --project `"$repoRoot`" run-pamica.py"
