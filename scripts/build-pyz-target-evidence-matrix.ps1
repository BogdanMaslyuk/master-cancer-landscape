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

$Required = @(
    "data\runtime\own_compounds\target_ligand_space_similarity_summary.parquet",
    "data\runtime\own_compounds\target_ligand_space_similarity_top_hits.parquet",
    "data\runtime\cellular_evidence\target_ligand_cellular_summary.parquet",
    "data\runtime\cellular_evidence\pyz_cellular_evidence_summary.parquet"
)

$Missing = @($Required | Where-Object { -not (Test-Path $_) })
if ($Missing.Count -gt 0) {
    Write-Host "Missing required inputs:" -ForegroundColor Red
    $Missing | ForEach-Object { Write-Host "  $_" -ForegroundColor Red }
    throw "Run Target Ligand Space and Cellular Evidence Layer v1.1 first."
}

Run-Step "MCL PYZ Target Evidence Matrix v1" @(
    ".\scripts\build_pyz_target_evidence_matrix.py"
)

Write-Host "`nPYZ Target Evidence Matrix build completed." -ForegroundColor Green
Write-Host "Full matrix: data\runtime\pyz_target_evidence\pyz_target_evidence_matrix.tsv"
Write-Host "Exploratory validation shortlist: data\runtime\pyz_target_evidence\pyz_target_validation_shortlist.tsv"
Write-Host "QC: outputs\qc\pyz_target_evidence_matrix_qc.tsv"
