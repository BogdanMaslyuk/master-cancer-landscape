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

Write-Host "=== Validate approved own-compound registry: PYZ-001...PYZ-065 ===" -ForegroundColor Cyan
& $Python ".\scripts\validate_pyz65_registry.py"
if ($LASTEXITCODE -ne 0) {
    throw "Approved PYZ registry validation failed with exit code $LASTEXITCODE"
}

Write-Host "`n=== Build strict oncology reference universe ===" -ForegroundColor Cyan
$ReferenceArgs = @(".\scripts\build-oncology-reference-compounds.ps1")
if ($ForceStructureRefresh) { $ReferenceArgs += "-ForceStructureRefresh" }
& powershell.exe -NoProfile -ExecutionPolicy Bypass -File @ReferenceArgs
if ($LASTEXITCODE -ne 0) {
    throw "Oncology reference build failed with exit code $LASTEXITCODE"
}

Write-Host "`n=== MCL Own Compound Structural Similarity v1.3 ===" -ForegroundColor Cyan
& $Python ".\scripts\build_own_compound_similarity_v1_3.py"
if ($LASTEXITCODE -ne 0) {
    throw "Own compound strict oncology-reference similarity build failed with exit code $LASTEXITCODE"
}

Write-Host "`nOwn compound similarity build completed for PYZ-001...PYZ-065." -ForegroundColor Green
Write-Host "Primary result: data\runtime\own_compounds\oncology_reference_target_similarity_summary.parquet"
Write-Host "Reference catalog: data\runtime\pharmacology\oncology_reference_compounds.tsv"
Write-Host "Reference structures: data\runtime\pharmacology\oncology_reference_structures.tsv"
Write-Host "The earlier broad target-linked comparison is preserved in its old output files but is no longer the primary result."
