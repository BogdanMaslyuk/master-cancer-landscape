param(
    [double]$CutoffNm = 10000,
    [int]$TopPerTarget = 20,
    [switch]$ForceRefresh
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

Write-Host "`n=== MCL Target Ligand Space v1: ChEMBL + BindingDB ===" -ForegroundColor Cyan
$LigandArgs = @(
    ".\scripts\build_target_ligand_space.py",
    "--cutoff-nm", $CutoffNm
)
if ($ForceRefresh) {
    $LigandArgs += "--force-refresh"
}
& $Python @LigandArgs
if ($LASTEXITCODE -ne 0) {
    throw "Target Ligand Space build failed with exit code $LASTEXITCODE"
}

Write-Host "`n=== PYZ-001...PYZ-065 against experimental target ligand space ===" -ForegroundColor Cyan
& $Python ".\scripts\build_pyz_target_ligand_similarity.py" --top-per-target $TopPerTarget
if ($LASTEXITCODE -ne 0) {
    throw "PYZ Target Ligand Space similarity failed with exit code $LASTEXITCODE"
}

Write-Host "`nTarget Ligand Space build completed." -ForegroundColor Green
Write-Host "Review data\runtime\target_ligand_space\target_summary.tsv"
Write-Host "Primary PYZ result: data\runtime\own_compounds\target_ligand_space_similarity_summary.tsv"
Write-Host "Nearest experimental ligands: data\runtime\own_compounds\target_ligand_space_similarity_top_hits.tsv"
