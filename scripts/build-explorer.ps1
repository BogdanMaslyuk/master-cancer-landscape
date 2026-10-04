param(
    [switch]$RefreshReference,
    [switch]$AllowPartialMultiomics
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

Run-Step "Build DepMap model profiles" @(".\scripts\build_depmap_model_profiles.py")

$MultiomicsArgs = @(".\scripts\build_depmap_multiomics.py")
if ($AllowPartialMultiomics) {
    $MultiomicsArgs += "--allow-partial"
}
Run-Step "Build DepMap multi-omics indexes" $MultiomicsArgs
Run-Step "Build fast DepMap runtime arrays" @(".\scripts\build_depmap_runtime_arrays.py")

$Reference = Join-Path $Root "data\processed\gene_explorer\gene_reference.parquet"
if ($RefreshReference -or -not (Test-Path $Reference)) {
    Run-Step "Build human gene reference snapshot" @(".\scripts\build_gene_reference_snapshot.py")
} else {
    Write-Host "`n=== Human gene reference snapshot ===" -ForegroundColor Cyan
    Write-Host "Existing reference snapshot retained. Use -RefreshReference to rebuild it."
}

Run-Step "Build Gene Explorer materialized indexes" @(".\scripts\build_gene_explorer_index.py")
Run-Step "Audit genome-wide expansion readiness" @(".\scripts\audit_genomewide_expansion.py")
Run-Step "Verify Explorer runtime indexes" @(".\scripts\verify_runtime_indexes.py")

Write-Host "`nExplorer build completed." -ForegroundColor Green
