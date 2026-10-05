param(
    [switch]$SkipAudit,
    [switch]$SkipBenchmark
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

if (-not $SkipAudit) {
    Run-Step "Audit pharmacology layer" @(".\scripts\audit_pharmacology_layer.py")
}

Run-Step "Build transparent candidate hypotheses" @(".\scripts\build_candidate_hypotheses.py")

if (-not $SkipBenchmark) {
    Run-Step "Known mechanism benchmark" @(".\scripts\build_known_mechanism_benchmark.py")
}

Write-Host "`nMCL Candidate Prioritization v1 build completed." -ForegroundColor Green
Write-Host "Restart backend/frontend and open /hypotheses."
Write-Host "Read outputs\qc\known_mechanism_benchmark.tsv before using novel candidates as laboratory priorities."
