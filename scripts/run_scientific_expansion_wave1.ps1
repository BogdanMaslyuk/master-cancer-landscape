$ErrorActionPreference = "Stop"

$root = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$python = Join-Path $root ".venv\Scripts\python.exe"

if (-not (Test-Path $python)) {
    throw "Python environment not found: $python"
}

Write-Host "=== Scientific Expansion Wave 1: cohort materialization ===" -ForegroundColor Cyan
& $python (Join-Path $PSScriptRoot "materialize_scientific_expansion_wave1.py")
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

Write-Host "`n=== Scientific Expansion Wave 1: genome-wide CRISPR analysis ===" -ForegroundColor Cyan
& $python (Join-Path $PSScriptRoot "analyze_scientific_expansion_wave1.py")
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

Write-Host "`nScientific Expansion Wave 1 completed successfully." -ForegroundColor Green
