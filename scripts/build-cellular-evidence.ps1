param(
    [switch]$RebuildTargetLigandSpace,
    [switch]$ForceTargetLigandRefresh
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

Run-Step "Validate approved own-compound registry: PYZ-001...PYZ-065" @(
    ".\scripts\validate_pyz65_registry.py"
)

if ($RebuildTargetLigandSpace) {
    Write-Host "`n=== Rebuild Target Ligand Space v1 ===" -ForegroundColor Cyan
    if ($ForceTargetLigandRefresh) {
        & powershell.exe -NoProfile -ExecutionPolicy Bypass -File ".\scripts\build-target-ligand-space.ps1" -ForceRefresh
    } else {
        & powershell.exe -NoProfile -ExecutionPolicy Bypass -File ".\scripts\build-target-ligand-space.ps1"
    }
    if ($LASTEXITCODE -ne 0) {
        throw "Target Ligand Space rebuild failed with exit code $LASTEXITCODE"
    }
}

$Required = @(
    "data\runtime\target_ligand_space\ligand_catalog.parquet",
    "data\runtime\own_compounds\target_ligand_space_similarity_top_hits.parquet",
    "data\runtime\pharmacology\compounds.parquet",
    "data\runtime\pharmacology\responses.parquet",
    "data\processed\depmap_model_target_context.parquet",
    "data\processed\depmap_crispr_model_atlas.parquet"
)

$Missing = @($Required | Where-Object { -not (Test-Path $_) })
if ($Missing.Count -gt 0) {
    Write-Host "Missing required inputs:" -ForegroundColor Red
    $Missing | ForEach-Object { Write-Host "  $_" -ForegroundColor Red }
    throw "Build Target Ligand Space, PRISM pharmacology and Candidate v2 before Cellular Evidence Layer v1."
}

Run-Step "MCL Cellular Evidence Layer v1: direct ligand -> PRISM cells -> DepMap CRISPR" @(
    ".\scripts\build_cellular_evidence_layer.py"
)

Write-Host "`nCellular Evidence Layer build completed." -ForegroundColor Green
Write-Host "Primary ligand summary: data\runtime\cellular_evidence\target_ligand_cellular_summary.tsv"
Write-Host "Model-level evidence: data\runtime\cellular_evidence\target_ligand_cellular_models.parquet"
Write-Host "PYZ summary: data\runtime\cellular_evidence\pyz_cellular_evidence_summary.tsv"
Write-Host "PYZ annotated neighbors: data\runtime\cellular_evidence\pyz_cellular_evidence_top_hits.tsv"
