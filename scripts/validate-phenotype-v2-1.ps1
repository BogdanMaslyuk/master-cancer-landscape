# Isolated Scientific Core v2.1 P0/P1 validation. Does not modify source data.
param([string]$SourceRoot = (Join-Path $env:USERPROFILE "Desktop\master-cancer-landscape"))
$ErrorActionPreference = "Stop"
$Review = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$Source = (Resolve-Path $SourceRoot).Path
$Python = Join-Path $Source ".venv\Scripts\python.exe"
if (-not (Test-Path $Python)) { throw "Python environment missing: $Python" }

$Previous = (Get-Location).Path
$OldPath = [Environment]::GetEnvironmentVariable("PYTHONPATH", "Process")
$OldData = [Environment]::GetEnvironmentVariable("MCL_SCIENTIFIC_DATA_ROOT", "Process")
$OldRoot = [Environment]::GetEnvironmentVariable("MCL_ROOT", "Process")
try {
    Set-Location $Review
    $env:PYTHONPATH = Join-Path $Review "apps\api"
    Remove-Item Env:\MCL_SCIENTIFIC_DATA_ROOT -ErrorAction SilentlyContinue
    Remove-Item Env:\MCL_ROOT -ErrorAction SilentlyContinue
    Write-Host ""
    Write-Host "=== Active source and phenotype data safety ===" -ForegroundColor Cyan
    @'
from pathlib import Path
import sys
from mcl_api import phenotype
from mcl_api.state import phenotype_store

review = Path(sys.argv[1]).resolve()
source = Path(sys.argv[2]).resolve()
expected = (review / "apps" / "api" / "mcl_api").resolve()
actual = Path(phenotype.__file__).resolve().parent
runtime = phenotype_store.runtime.resolve()
if actual != expected:
    raise SystemExit(f"WRONG SOURCE: {actual} != {expected}")
if not runtime.is_relative_to(review) or runtime.is_relative_to(source):
    raise SystemExit("UNSAFE phenotype runtime outside review checkout")
assert phenotype_store.compounds()["total"] == 65
models = {r["cell_line_id"]: r for r in phenotype_store.models()["items"]}
assert {"HEK293", "HEK293T", "HEK293FT", "MRC-5", "BJ5ta"} <= set(models)
assert models["MRC-5"]["assay_ready"] is False
print("Code:", actual)
print("Phenotype runtime:", runtime)
print("65 PYZ and separate HEK/MRC models: PASS")
'@ | & $Python - $Review $Source
    if ($LASTEXITCODE -ne 0) { throw "Source/runtime safety check failed" }

    Write-Host ""
    Write-Host "=== Phenotype and scientific regressions ===" -ForegroundColor Cyan
    & $Python -m pytest .\apps\api\tests\test_phenotype.py .\apps\api\tests\test_architecture_boundaries.py .\apps\api\tests\test_openapi_contract.py .\apps\api\tests\test_scientific_reference_context.py .\apps\api\tests\test_uploaded_molecules.py -q
    if ($LASTEXITCODE -ne 0) { throw "Tests failed" }
    Write-Host ""
    Write-Host "=== Generated API types ===" -ForegroundColor Cyan
    & $Python .\scripts\generate_frontend_api_types.py --check
    if ($LASTEXITCODE -ne 0) { throw "OpenAPI type contract failed" }
    Write-Host ""
    Write-Host "Scientific Core v2.1 P0/P1: PASS" -ForegroundColor Green
}
finally {
    if ($null -eq $OldPath) { Remove-Item Env:\PYTHONPATH -ErrorAction SilentlyContinue }
    else { $env:PYTHONPATH = $OldPath }
    if ($null -eq $OldData) { Remove-Item Env:\MCL_SCIENTIFIC_DATA_ROOT -ErrorAction SilentlyContinue }
    else { $env:MCL_SCIENTIFIC_DATA_ROOT = $OldData }
    if ($null -eq $OldRoot) { Remove-Item Env:\MCL_ROOT -ErrorAction SilentlyContinue }
    else { $env:MCL_ROOT = $OldRoot }
    Set-Location $Previous
}
