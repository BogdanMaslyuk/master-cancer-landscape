param(
    [switch]$Recreate
)

$ErrorActionPreference = "Stop"
$Root = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$Venv = Join-Path $Root ".venv"
$Python = Join-Path $Venv "Scripts\python.exe"
$Constraints = Join-Path $Root "constraints\python-3.13.txt"
$Explorer = Join-Path $Root "apps\explorer"

Set-Location $Root

function Run([string]$Label, [scriptblock]$Command) {
    Write-Host "`n=== $Label ===" -ForegroundColor Cyan
    & $Command
    if ($LASTEXITCODE -ne 0) {
        throw "$Label failed with exit code $LASTEXITCODE"
    }
}

if ($Recreate -and (Test-Path $Venv)) {
    Write-Host "Removing existing .venv ..." -ForegroundColor Yellow
    Remove-Item $Venv -Recurse -Force
}

if (-not (Test-Path $Python)) {
    $Uv = Get-Command uv.exe -ErrorAction SilentlyContinue
    if ($Uv) {
        # uv creates minimal environments without pip unless --seed is requested.
        # Architecture v1 deliberately seeds pip because the rest of the bootstrap
        # and verification workflow uses `python -m pip` consistently on Windows
        # and in CI.
        Run "Create Python 3.13 environment with uv" {
            & uv.exe venv --python 3.13 --seed $Venv
        }
    } else {
        $Launcher = Get-Command py.exe -ErrorAction SilentlyContinue
        if ($Launcher) {
            $Detected = & py.exe -3.13 -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')" 2>$null
            if ($LASTEXITCODE -eq 0 -and $Detected.Trim() -eq "3.13") {
                Run "Create Python 3.13 environment" { & py.exe -3.13 -m venv $Venv }
            } else {
                throw "Python 3.13 was not found by py.exe and uv.exe is unavailable. Install Python 3.13 or uv, then re-run bootstrap."
            }
        } else {
            $SystemPython = Get-Command python.exe -ErrorAction SilentlyContinue
            if (-not $SystemPython) {
                throw "Python 3.13 was not found. Install Python 3.13 or uv."
            }
            $Version = & python.exe -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')"
            if ($Version.Trim() -ne "3.13") {
                throw "python.exe is Python $Version; MCL Architecture v1 requires Python 3.13 for the locked environment."
            }
            Run "Create Python 3.13 environment" { & python.exe -m venv $Venv }
        }
    }
}

$VenvVersion = & $Python -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')"
if ($VenvVersion.Trim() -ne "3.13") {
    throw ".venv uses Python $VenvVersion. Re-run .\scripts\bootstrap.ps1 -Recreate with Python 3.13 available."
}

if (-not (Test-Path $Constraints)) {
    throw "Dependency constraints are missing: $Constraints"
}

Run "Pin pip" { & $Python -m pip install "pip==26.2.1" }
Run "Install locked Python dependencies" {
    & $Python -m pip install -c $Constraints -e ".[dev,analysis]" -e ".\apps\api[dev]"
}
Run "Check Python dependency graph" { & $Python -m pip check }
Run "Verify Python constraints" { & $Python .\scripts\verify_python_constraints.py }

Push-Location $Explorer
try {
    Run "Install locked frontend dependencies" { & npm.cmd ci }
} finally {
    Pop-Location
}

Write-Host "`nMCL bootstrap completed successfully." -ForegroundColor Green
Write-Host "Python: $Python"
Write-Host "Next: .\scripts\verify.ps1"
