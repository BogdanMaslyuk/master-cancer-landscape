# E1-P0-B3: read-only evidence register for the completed full BindingDB audit.
# Centralized MCL_Project; no additional BindingDB download, unzip, UI or API changes.
param(
    [string]$Root = (Join-Path $env:USERPROFILE 'MCL_Project'),
    [string]$SourceZip = ''
)
$ErrorActionPreference = 'Stop'

$Repo = Join-Path $Root 'repositories\master-cancer-landscape'
$Python = Join-Path $Repo '.venv\Scripts\python.exe'
$Audit = Join-Path $PSScriptRoot 'audit_pyz_e1_p0b3_evidence.py'
$Reports = Join-Path $Root 'reports'

foreach ($Path in @($Root, $Repo, $Python, $Audit, $Reports)) {
    if (-not (Test-Path -LiteralPath $Path)) { throw "Required MCL path not found: $Path" }
}

if (-not $SourceZip) {
    $Candidates = @(Get-ChildItem -LiteralPath $Reports -File -Filter 'MCL_E1_P0B2_RESULT_*.zip' |
        Sort-Object LastWriteTime -Descending)
    if ($Candidates.Count -eq 0) {
        throw 'MCL_E1_P0B2_RESULT_*.zip not found in MCL_Project\reports'
    }
    $SourceZip = $Candidates[0].FullName
}
if (-not (Test-Path -LiteralPath $SourceZip -PathType Leaf)) {
    throw "The E1-P0-B2 result ZIP does not exist: $SourceZip"
}
$SourceZip = (Resolve-Path -LiteralPath $SourceZip).Path
$SourceHashBefore = (Get-FileHash -LiteralPath $SourceZip -Algorithm SHA256).Hash
$Stamp = Get-Date -Format 'yyyyMMdd_HHmmss_fff'
$Output = Join-Path $Reports ("E1_P0B3_EVIDENCE_" + $Stamp)

Write-Host '=== MCL / E1-P0-B3 / Scientific evidence traceability ===' -ForegroundColor Cyan
Write-Host "Project: $Root"
Write-Host "Input: $SourceZip"
& $Python $Audit --source-zip $SourceZip --out-dir $Output
if ($LASTEXITCODE -ne 0) { throw 'B3 evidence audit failed. No scientific sources were changed.' }

$SourceHashAfter = (Get-FileHash -LiteralPath $SourceZip -Algorithm SHA256).Hash
if ($SourceHashAfter -ne $SourceHashBefore) {
    throw 'Source B2 result ZIP changed during analysis; halt and investigate.'
}
$ResultPath = Join-Path $Output 'evidence_review_report.json'
$Report = Get-Content -LiteralPath $ResultPath -Raw -Encoding UTF8 | ConvertFrom-Json
if ($Report.scientific_status -ne 'NO_PRIMARY_ASSAYS_APPROVED' -or $Report.statistics.primary_assays_approved -ne 0) {
    throw 'Scientific validation guard failed; result not archived.'
}
$Files = @(
    $ResultPath,
    (Join-Path $Output 'assay_evidence_register.tsv'),
    (Join-Path $Output 'publication_review_queue.tsv'),
    (Join-Path $Output 'manual_source_review_checklist.tsv')
)
foreach ($File in $Files) {
    if (-not (Test-Path -LiteralPath $File -PathType Leaf)) {
        throw "Audit output file is missing: $File"
    }
}
$ResultZip = Join-Path $Reports ("MCL_E1_P0B3_RESULT_" + $Stamp + '.zip')
Compress-Archive -LiteralPath $Files -DestinationPath $ResultZip -CompressionLevel Fastest
if (-not (Test-Path -LiteralPath $ResultZip -PathType Leaf)) {
    throw 'B3 result ZIP could not be created.'
}
Write-Host '=== B3 completed ===' -ForegroundColor Green
Write-Host "Records: $($Report.statistics.records)"
Write-Host "Papers: $($Report.statistics.papers)"
Write-Host "Source InChIKey differences: $($Report.statistics.source_inchikey_discrepancy)"
Write-Host 'Primary assays approved: 0'
Write-Host "Result: $ResultZip" -ForegroundColor Green
# Do not launch explorer.exe: shell access may be disabled on Windows machines.