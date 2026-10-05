param()

$ErrorActionPreference = "Stop"
$Root = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$Python = Join-Path $Root ".venv\Scripts\python.exe"

if (-not (Test-Path $Python)) {
    throw "Python environment not found: $Python"
}

Set-Location $Root

Write-Host "=== MCL Own Compound Structural Similarity v1 ===" -ForegroundColor Cyan
& $Python ".\scripts\build_own_compound_similarity.py"
if ($LASTEXITCODE -ne 0) {
    throw "Own compound structural similarity build failed with exit code $LASTEXITCODE"
}

Write-Host "`nOwn compound similarity build completed." -ForegroundColor Green
Write-Host "Review data\runtime\own_compounds\target_similarity_summary.tsv and outputs\qc\own_compound_similarity_qc.tsv."
