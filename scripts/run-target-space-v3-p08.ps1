param(
    [string]$Root = "$env:USERPROFILE\MCL_Project",
    [string]$P0Zip = '',
    [string]$P06Zip = '',
    [string]$P07Zip = '',
    [string]$ChEMBLSqlite = '',
    [string]$BindingDBTsv = '',
    [string]$ChEMBLRelease = '',
    [string]$BindingDBRelease = '',
    [ValidateSet('focus6', 'sample40')][string]$Panel = 'focus6',
    [ValidateRange(1, 10000)][int]$MaxRowsPerTarget = 1000
)

$ErrorActionPreference = 'Stop'
$Root = [System.IO.Path]::GetFullPath($Root)
$Repo = Join-Path $Root 'repositories\master-cancer-landscape'
$Python = Join-Path $Repo '.venv\Scripts\python.exe'
$Auditor = Join-Path $PSScriptRoot 'audit_target_space_v3_p08.py'
$Reports = Join-Path $Root 'reports'
if (-not $P0Zip) { $P0Zip = Join-Path $Reports 'TLS_V3_P0_20261119T133253114Z_e9ccd63f.zip' }
if (-not $P06Zip) { $P06Zip = Join-Path $Reports 'TLS_V3_P06_20261119T133323433Z_e9ccd63f.zip' }
if (-not $P07Zip) { $P07Zip = Join-Path $Reports 'TLS_V3_P07_20261119T154010621Z_6b853654.zip' }
foreach ($path in @($Python, $Auditor, $P0Zip, $P06Zip, $P07Zip)) {
    if (-not (Test-Path -LiteralPath $path -PathType Leaf)) {
        throw "Required read-only input missing: $path"
    }
}
if (-not (Test-Path -LiteralPath $Reports -PathType Container)) {
    New-Item -Path $Reports -ItemType Directory -Force | Out-Null
}
$Commit = & git -C $PSScriptRoot rev-parse --short=8 HEAD 2>$null
if ($LASTEXITCODE -ne 0 -or -not $Commit) { $Commit = 'unknown' }
$Stamp = (Get-Date).ToUniversalTime().ToString('yyyyMMddTHHmmssfffZ')
$Output = Join-Path $Reports "TLS_V3_P08_${Stamp}_${Commit}"
$Zip = "$Output.zip"
$Log = Join-Path $Reports "TLS_V3_P08_DIAGNOSTIC_${Stamp}_${Commit}.log"
$PythonArgs = @(
    $Auditor, '--p0-zip', $P0Zip, '--p06-zip', $P06Zip,
    '--p07-zip', $P07Zip, '--out', $Output, '--panel', $Panel,
    '--max-rows-per-target', "$MaxRowsPerTarget"
)
if ($ChEMBLSqlite) {
    if (-not (Test-Path -LiteralPath $ChEMBLSqlite -PathType Leaf)) { throw "ChEMBL SQLite missing: $ChEMBLSqlite" }
    if (-not $ChEMBLRelease) { throw 'Provide -ChEMBLRelease (e.g., chembl_36)' }
    $PythonArgs += @('--chembl-sqlite', $ChEMBLSqlite, '--chembl-release', $ChEMBLRelease)
}
if ($BindingDBTsv) {
    if (-not (Test-Path -LiteralPath $BindingDBTsv -PathType Leaf)) { throw "BindingDB TSV/ZIP missing: $BindingDBTsv" }
    if (-not $BindingDBRelease) { throw 'Provide -BindingDBRelease (e.g., 202611)' }
    $PythonArgs += @('--bindingdb-tsv', $BindingDBTsv, '--bindingdb-release', $BindingDBRelease)
}
Write-Host '=== MCL TLS V3 P0.8: OFFLINE local experimental-source validation ===' -ForegroundColor Cyan
Write-Host "Panel: $Panel; Python: $Python"
Write-Host "ChEMBL SQLite: $(if ($ChEMBLSqlite) { $ChEMBLSqlite } else { 'NOT_PROVIDED' })"
Write-Host "BindingDB TSV: $(if ($BindingDBTsv) { $BindingDBTsv } else { 'NOT_PROVIDED' })"
Write-Host 'NO HTTP. P0/P0.6/P0.7, PYZ 65x13, scientific source files: READ-ONLY.'
$OldErrorAction = $ErrorActionPreference
$ExitCode = 1
try {
    $ErrorActionPreference = 'Continue'
    & $Python @PythonArgs 2>&1 | Tee-Object -FilePath $Log
    $ExitCode = $LASTEXITCODE
} finally {
    $ErrorActionPreference = $OldErrorAction
}
if ($ExitCode -ne 0) {
    Write-Host "Original Python diagnostic: $Log" -ForegroundColor Yellow
    Get-Content -LiteralPath $Log -Tail 90 | Out-Host
    throw "TLS V3 P0.8 failed (exit code $ExitCode)"
}
if (-not (Test-Path -LiteralPath (Join-Path $Output 'p08_summary.json') -PathType Leaf)) {
    throw 'P0.8 summary absent despite zero exit code'
}
Compress-Archive -Path (Join-Path $Output '*') -DestinationPath $Zip -CompressionLevel Optimal
if (-not (Test-Path -LiteralPath $Zip -PathType Leaf)) { throw 'P0.8 archive missing' }
Write-Host "Report ZIP: $Zip" -ForegroundColor Green
Write-Host 'scientific_acceptance=NOT_APPROVED; merge_permitted=false'