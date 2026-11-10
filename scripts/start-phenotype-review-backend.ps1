# Run review code while existing DepMap/PRISM are read from the source checkout.
# PhenotypeStore writes only to the review checkout, independently of MCL_ROOT.
param(
    [string]$SourceRoot = (Join-Path $env:USERPROFILE "Desktop\master-cancer-landscape"),
    [int]$Port = 8001
)
$ErrorActionPreference = "Stop"
$Review = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$Source = (Resolve-Path $SourceRoot).Path
$Python = Join-Path $Source ".venv\Scripts\python.exe"
if (-not (Test-Path $Python)) { throw "Python not found: $Python" }
$env:PYTHONPATH = Join-Path $Review "apps\api"
$env:MCL_ROOT = $Source
$env:MCL_SCIENTIFIC_DATA_ROOT = $Source
Set-Location $Review
Write-Host "MCL Scientific Core v2.1 backend http://127.0.0.1:$Port" -ForegroundColor Cyan
Write-Host "Existing scientific data: $Source (read-only)"
Write-Host "Raw phenotype uploads: $Review\data\runtime\phenotype"
Write-Host "Keep this PowerShell window open while using the interface."
& $Python -m uvicorn mcl_api.main:app --app-dir (Join-Path $Review "apps\api") --host 127.0.0.1 --port $Port
exit $LASTEXITCODE
