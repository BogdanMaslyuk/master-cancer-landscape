param(
    [switch]$ListOnly,
    [switch]$SkipDownload
)

$ErrorActionPreference = "Stop"
$Root = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$Python = Join-Path $Root ".venv\Scripts\python.exe"

if (-not (Test-Path $Python)) {
    throw "Python environment not found: $Python"
}

Set-Location $Root

function Run-Step([string]$Label, [string[]]$Arguments) {
    Write-Host "`n=== $Label ===" -ForegroundColor Cyan
    & $Python @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "$Label failed with exit code $LASTEXITCODE"
    }
}

if ($ListOnly) {
    Run-Step "PRISM 24Q2 download plan" @(".\scripts\fetch_prism_24q2.py", "--list-only")
    Write-Host "`nNo files were downloaded. Re-run without -ListOnly to build the layer." -ForegroundColor Yellow
    exit 0
}

if (-not $SkipDownload) {
    Run-Step "Download PRISM 24Q2 minimal source set" @(".\scripts\fetch_prism_24q2.py")
} else {
    Write-Host "`n=== Download PRISM 24Q2 ===" -ForegroundColor Cyan
    Write-Host "Skipped by -SkipDownload; using files already present under data\raw\pharmacology\prism_24q2."
}

Run-Step "Normalize PRISM 24Q2" @(".\scripts\ingest_prism_24q2.py")
Run-Step "Build MCL Pharmacology Layer v1" @(".\scripts\build_pharmacology_layer.py")

Write-Host "`nMCL Pharmacology Layer v1 build completed." -ForegroundColor Green
Write-Host "Open a model card and inspect the new 'Фармакология модели' section after restarting backend/frontend."
