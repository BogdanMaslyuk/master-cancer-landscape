param(
    [switch]$ForceStructureRefresh
)

$ErrorActionPreference = "Stop"
$Root = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$Python = Join-Path $Root ".venv\Scripts\python.exe"

if (-not (Test-Path $Python)) {
    throw "Python environment not found: $Python"
}

Set-Location $Root

Write-Host "=== Curated oncology reference compounds ===" -ForegroundColor Cyan
& $Python ".\scripts\build_oncology_reference_compounds.py"
if ($LASTEXITCODE -ne 0) {
    throw "Oncology reference catalog build failed with exit code $LASTEXITCODE"
}

Write-Host "`n=== Resolve oncology reference structures ===" -ForegroundColor Cyan
$ResolveArgs = @(".\scripts\resolve_oncology_reference_structures.py")
if ($ForceStructureRefresh) { $ResolveArgs += "--force" }
& $Python @ResolveArgs
if ($LASTEXITCODE -ne 0) {
    throw "Oncology reference structure resolution failed with exit code $LASTEXITCODE"
}

Write-Host "`nOncology reference build completed." -ForegroundColor Green
Write-Host "Review data\runtime\pharmacology\oncology_reference_target_summary.tsv and oncology_reference_structures.tsv."
