param(
    [string]$Root = "$env:USERPROFILE\MCL_Project",
    [string]$P0Zip = '',
    [string]$P06Zip = '',
    [switch]$EnableTargetedRetry,
    [ValidateRange(1,24)] [int]$MaxRequests = 12
)
$ErrorActionPreference = 'Stop'
$Root = [IO.Path]::GetFullPath($Root)
$Reports = Join-Path $Root 'reports'
$Repo = Join-Path $Root 'repositories\master-cancer-landscape'
$Python = Join-Path $Repo '.venv\Scripts\python.exe'
$Script = Join-Path $PSScriptRoot 'audit_target_space_v3_p07.py'
if (-not $P0Zip) { $P0Zip = Join-Path $Reports 'TLS_V3_P0_20261119T133253114Z_e9ccd63f.zip' }
if (-not $P06Zip) { $P06Zip = Join-Path $Reports 'TLS_V3_P06_20261119T133323433Z_e9ccd63f.zip' }
foreach ($File in @($Python, $Script, $P0Zip, $P06Zip)) {
    if (-not (Test-Path -LiteralPath $File -PathType Leaf)) { throw "File not found: $File" }
}
$Commit = & git -C $PSScriptRoot rev-parse --short=8 HEAD 2>$null
if ($LASTEXITCODE -ne 0 -or -not $Commit) { $Commit='unknown' }
$Stamp = (Get-Date).ToUniversalTime().ToString('yyyyMMddTHHmmssfffZ')
$Output = Join-Path $Reports "TLS_V3_P07_${Stamp}_${Commit}"
$Archive = "$Output.zip"
$Log = Join-Path $Reports "TLS_V3_P07_DIAGNOSTIC_${Stamp}_${Commit}.log"
if ((Test-Path -LiteralPath $Output) -or (Test-Path -LiteralPath $Archive)) {
    throw 'Refusing to overwrite old audit results'
}
$ArgsList = @($Script,'--p0-zip',$P0Zip,'--p06-zip',$P06Zip,'--out',$Output)
if ($EnableTargetedRetry) {
    $ArgsList += @('--enable-targeted-retry','--max-requests',"$MaxRequests")
    Write-Host 'P0.7: explicitly approved limited ChEMBL retries' -ForegroundColor Yellow
} else {
    Write-Host 'P0.7: OFFLINE - replay previous archives, no network access' -ForegroundColor Cyan
}
$PrevPreference = $ErrorActionPreference
$ExitCode = 1
try {
    $ErrorActionPreference = 'Continue'
    & $Python @ArgsList 2>&1 | Tee-Object -FilePath $Log
    $ExitCode = $LASTEXITCODE
} finally {
    $ErrorActionPreference = $PrevPreference
}
if ($ExitCode -ne 0) {
    Write-Host "Python error log: $Log" -ForegroundColor Red
    Get-Content -LiteralPath $Log -Tail 100 | Out-Host
    throw "P0.7 failed with exit code $ExitCode"
}
if (-not (Test-Path -LiteralPath (Join-Path $Output 'p07_summary.json'))) { throw 'Missing P0.7 summary' }
Compress-Archive -Path (Join-Path $Output '*') -DestinationPath $Archive -CompressionLevel Optimal
if (-not (Test-Path -LiteralPath $Archive -PathType Leaf)) { throw 'P0.7 ZIP missing' }
Write-Host "P0.7 ZIP: $Archive" -ForegroundColor Green
Write-Host 'scientific_acceptance=NOT_APPROVED; merge_permitted=false'