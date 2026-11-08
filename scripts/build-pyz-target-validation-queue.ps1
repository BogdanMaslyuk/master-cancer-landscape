$ErrorActionPreference = "Stop"
$Root = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$Python = Join-Path $Root ".venv\Scripts\python.exe"

if (-not (Test-Path $Python)) {
    throw "Python environment not found: $Python"
}

Set-Location $Root

$Required = @(
    "data\runtime\pyz_target_evidence\pyz_target_validation_shortlist.parquet",
    "config\pyz_target_validation_structures.tsv"
)
$Missing = @($Required | Where-Object { -not (Test-Path $_) })
if ($Missing.Count -gt 0) {
    Write-Host "Missing required inputs:" -ForegroundColor Red
    $Missing | ForEach-Object { Write-Host "  $_" -ForegroundColor Red }
    throw "Run .\scripts\build-pyz-target-evidence-matrix.ps1 first."
}

Write-Host "=== MCL PYZ Target Validation Queue v1 ===" -ForegroundColor Cyan
& $Python ".\scripts\build_pyz_target_validation_queue.py"
if ($LASTEXITCODE -ne 0) {
    throw "PYZ Target Validation Queue v1 failed with exit code $LASTEXITCODE"
}

Write-Host "`nPYZ Target Validation Queue build completed." -ForegroundColor Green
Write-Host "Queue: data\runtime\pyz_target_evidence\pyz_target_validation_queue.tsv"
Write-Host "Docking structure plan: data\runtime\pyz_target_evidence\pyz_target_validation_structures.tsv"
Write-Host "QC: outputs\qc\pyz_target_validation_queue_qc.tsv"
