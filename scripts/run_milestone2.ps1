$ErrorActionPreference = "Stop"

Write-Host "[1/3] Validate four Open Targets disease mappings"
mcl validate-ot-diseases

Write-Host "[2/3] Fetch Wave 1 direct Open Targets associations + tractability"
mcl fetch-opentargets

Write-Host "[3/3] QC"
mcl qc-opentargets
