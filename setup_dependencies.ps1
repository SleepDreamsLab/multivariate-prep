<#
.SYNOPSIS
    Clones the external GitHub dependencies listed in dependancies.m into the
    folder one level up from this repo, and registers each one with GitHub
    Desktop (via its `github` CLI).

.DESCRIPTION
    Repos already present at the target path are left untouched (no pull/
    overwrite) - they are just registered with GitHub Desktop. Only missing
    repos are freshly cloned. Folders that exist but are not git repos are
    skipped with a warning rather than being clobbered.
#>

$ErrorActionPreference = 'Stop'

$repoRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$parent   = Split-Path -Parent $repoRoot

# Mirrors the external dependencies listed in dependancies.m
$deps = @(
    [pscustomobject]@{ Name = 'GEDAI-master';    Url = 'https://github.com/SvennoNito/GEDAI-master.git';      Branch = 'sleep-fast' }
    [pscustomobject]@{ Name = 'bids-matlab';      Url = 'https://github.com/bids-standard/bids-matlab.git';    Branch = $null }
    [pscustomobject]@{ Name = 'cleanline';        Url = 'https://github.com/sccn/cleanline.git';               Branch = $null }
    [pscustomobject]@{ Name = 'zapline-plus';     Url = 'https://github.com/MariusKlug/zapline-plus.git';      Branch = $null }
    [pscustomobject]@{ Name = 'bva-io';           Url = 'https://github.com/sccn/bva-io.git';                  Branch = $null }
    [pscustomobject]@{ Name = 'pAMICA';           Url = 'https://github.com/sccn/pAMICA.git';                  Branch = $null }
    # ICLabel carries matconvnet as a submodule and run_ICL calls into it (vl_setupnn),
    # so a plain clone gives you a plugin that fails at the inference step. Recurse = $true
    # both clones with --recurse-submodules and back-fills submodules in an existing clone.
    [pscustomobject]@{ Name = 'ICLabel';          Url = 'https://github.com/sccn/ICLabel.git';                 Branch = $null; Recurse = $true }
    [pscustomobject]@{ Name = 'firfilt';          Url = 'https://github.com/sccn/firfilt.git';                 Branch = $null }
    [pscustomobject]@{ Name = 'brainstorm3';      Url = 'https://github.com/brainstorm-tools/brainstorm3.git'; Branch = $null }
    [pscustomobject]@{ Name = 'eeglab';           Url = 'https://github.com/sccn/eeglab.git';                  Branch = $null }
    [pscustomobject]@{ Name = 'eeg-oscillations'; Url = 'https://github.com/SvennoNito/eeg-oscillations.git';  Branch = $null }
)

# PATH is only refreshed for new processes, so a terminal opened before Git
# was installed won't see it yet even though `git --version` works elsewhere.
# Fall back to the registry (Git for Windows always records its InstallPath
# there) and then to every drive letter, since Git isn't always on C:.
$git = (Get-Command git -ErrorAction SilentlyContinue).Source

if (-not $git) {
    $regPaths = @(
        'HKLM:\SOFTWARE\GitForWindows'
        'HKLM:\SOFTWARE\WOW6432Node\GitForWindows'
        'HKCU:\SOFTWARE\GitForWindows'
    )
    foreach ($regPath in $regPaths) {
        $installPath = (Get-ItemProperty -Path $regPath -Name InstallPath -ErrorAction SilentlyContinue).InstallPath
        if ($installPath) {
            $candidate = Join-Path $installPath 'cmd\git.exe'
            if (Test-Path $candidate) { $git = $candidate; break }
        }
    }
}

if (-not $git) {
    $drives = (Get-PSDrive -PSProvider FileSystem).Root
    $suffixes = @('Program Files\Git\cmd\git.exe', 'Program Files (x86)\Git\cmd\git.exe')
    $candidates = foreach ($drive in $drives) { foreach ($suffix in $suffixes) { Join-Path $drive $suffix } }
    $candidates += "$env:LOCALAPPDATA\Programs\Git\cmd\git.exe"
    $git = $candidates | Where-Object { Test-Path $_ } | Select-Object -First 1
}

if (-not $git) {
    throw "git.exe not found on PATH, in the registry, or on any drive's 'Program Files\Git'. Close and reopen your terminal after installing Git (PATH is only re-read on new sessions), or install Git for Windows from https://git-scm.com/download/win."
}

$githubCli = Get-Command github -ErrorAction SilentlyContinue
if (-not $githubCli) {
    Write-Warning "GitHub Desktop CLI ('github' command) not found on PATH. Repos will still be cloned, but won't be auto-added to GitHub Desktop - add them manually via File > Add Local Repository."
}

foreach ($dep in $deps) {
    $target = Join-Path $parent $dep.Name
    Write-Host "== $($dep.Name) ==" -ForegroundColor Cyan

    if (-not (Test-Path $target)) {
        Write-Host "  Cloning into $target"
        $cloneArgs = @('clone')
        if ($dep.Branch)  { $cloneArgs += @('--branch', $dep.Branch) }
        if ($dep.Recurse) { $cloneArgs += '--recurse-submodules' }
        & $git @cloneArgs $dep.Url $target
    }
    elseif (Test-Path (Join-Path $target '.git')) {
        Write-Host "  Already present as a git repo - leaving contents untouched."
        # Exception to "untouched": submodules that are still empty. Filling them in only
        # adds files the clone was always meant to have, and without them the plugin is
        # broken in a way that surfaces late, deep inside a run.
        if ($dep.Recurse) {
            Write-Host "  Initialising submodules"
            & $git -C $target submodule update --init --recursive
        }
    }
    else {
        Write-Warning "  '$target' exists but is not a git repository. Skipping - resolve manually if you want it tracked (e.g. rename it aside and re-run)."
        continue
    }

    if ($githubCli) {
        Write-Host "  Registering with GitHub Desktop"
        github open $target
    }
}

Write-Host "`nDone." -ForegroundColor Green
