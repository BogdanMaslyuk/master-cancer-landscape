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

Write-Host "=== Resolve structures for compounds linked to priority MCL targets ===" -ForegroundColor Cyan
$EnrichmentArgs = @(".\scripts\enrich_priority_compound_structures_pubchem.py")
if ($ForceStructureRefresh) {
    $EnrichmentArgs += "--force"
}
& $Python @EnrichmentArgs
if ($LASTEXITCODE -ne 0) {
    throw "Known compound structure enrichment failed with exit code $LASTEXITCODE"
}

Write-Host "`n=== MCL Own Compound Structural Similarity v1.1 ===" -ForegroundColor Cyan
& $Python ".\scripts\build_own_compound_similarity_v1_1.py"
if ($LASTEXITCODE -ne 0) {
    throw "Own compound structural similarity build failed with exit code $LASTEXITCODE"
}

Write-Host "`nOwn compound similarity build completed." -ForegroundColor Green
Write-Host "Review data\runtime\pharmacology\compound_structure_registry.tsv, data\runtime\own_compounds\target_similarity_summary.tsv and QC outputs."
