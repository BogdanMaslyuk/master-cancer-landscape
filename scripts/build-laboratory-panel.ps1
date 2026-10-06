param()

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

Run-Step "Map department cell-line collection to DepMap" @(".\scripts\build_laboratory_panel.py")
Run-Step "Intersect Candidate v2 with laboratory panel" @(".\scripts\build_laboratory_candidates.py")
Run-Step "Build mechanism-aware experimental panels" @(".\scripts\build_laboratory_mechanism_panels_v1_1.py")

Write-Host "`nMCL Laboratory Panel v1.1 build completed." -ForegroundColor Green
Write-Host "Review data\processed\laboratory_panel.tsv for identity mappings before using laboratory recommendations."
Write-Host "A549 and PC3 use curated DepMap identity decisions; U251 remains unresolved until the physical source is confirmed."
Write-Host "Restart backend/frontend and open /laboratory."
