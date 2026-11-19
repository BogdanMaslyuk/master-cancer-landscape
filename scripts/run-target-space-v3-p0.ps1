param(
    [string]$Root = "$env:USERPROFILE\MCL_Project",
    [string]$ScienceRoot = '',
    [string]$ReferenceRoot = '',
    [string]$ReleaseSnapshot = '26Q1'
)

$ErrorActionPreference = 'Stop'
$Root = [IO.Path]::GetFullPath($Root)
if (-not $ScienceRoot) { $ScienceRoot = Join-Path $Root 'repositories\mcl-pyz-v2-review' }
if (-not $ReferenceRoot) { $ReferenceRoot = Join-Path $Root 'repositories\master-cancer-landscape' }
$ScienceRoot = [IO.Path]::GetFullPath($ScienceRoot)
$ReferenceRoot = [IO.Path]::GetFullPath($ReferenceRoot)
$Repo = Join-Path $Root 'repositories\master-cancer-landscape'
$Python = Join-Path $Repo '.venv\Scripts\python.exe'
$Script = Join-Path $PSScriptRoot 'audit_target_space_v3_p0.py'

foreach ($File in @($Python, $Script)) {
    if (-not (Test-Path -LiteralPath $File -PathType Leaf)) { throw "Required file not found: $File" }
}
foreach ($Dir in @($ScienceRoot, $ReferenceRoot)) {
    if (-not (Test-Path -LiteralPath $Dir -PathType Container)) { throw "Source folder not found: $Dir" }
}
$GitCommit = (& git -C $PSScriptRoot rev-parse --short=8 HEAD 2>$null)
if ($LASTEXITCODE -ne 0 -or -not $GitCommit) {
    $GitCommit = 'unknown'
}
$Stamp = (Get-Date).ToUniversalTime().ToString('yyyyMMddTHHmmssfffZ')
$Out = Join-Path $Root "reports\TLS_V3_P0_${Stamp}_${GitCommit}"
$Zip = "$Out.zip"
if ((Test-Path -LiteralPath $Out) -or (Test-Path -LiteralPath $Zip)) {
    throw "Refusing to overwrite existing report: $Out"
}
$Parent = Split-Path $Out -Parent
if (-not (Test-Path $Parent)) { New-Item -ItemType Directory -Force -Path $Parent | Out-Null }

Write-Host '=== MCL TLS V3 / P0: read-only target coverage inventory ===' -ForegroundColor Cyan
Write-Host "Scientific root: $ScienceRoot"
Write-Host "Reference root: $ReferenceRoot"
Write-Host 'Network disabled. No changes to scientific data, PYZ pages, or current 65x13 matrix.'
$PythonArgs = @(
    $Script, '--science-root', $ScienceRoot,
    '--reference-root', $ReferenceRoot,
    '--out', $Out, '--offline', '--release-snapshot', $ReleaseSnapshot
)
# Capture the *original Python exception* even in Windows PowerShell 5.1.
# Native stderr can become a PowerShell error record if ErrorActionPreference=Stop,
# so temporarily disable terminating on STDERR while preserving Python exit code.
$Log = Join-Path $Parent ("TLS_V3_P0_DIAGNOSTIC_${Stamp}_${GitCommit}.log")
$OriginalErrorActionPreference = $ErrorActionPreference
$PythonExitCode = 1
try {
    $ErrorActionPreference = 'Continue'
    & $Python @PythonArgs 2>&1 | Tee-Object -FilePath $Log
    $PythonExitCode = $LASTEXITCODE
} finally {
    $ErrorActionPreference = $OriginalErrorActionPreference
}
if ($PythonExitCode -ne 0) {
    Write-Host "Python error log: $Log" -ForegroundColor Yellow
    Write-Host '=== LAST 90 LINES OF ORIGINAL PYTHON OUTPUT ===' -ForegroundColor Red
    Get-Content -LiteralPath $Log -Tail 90 | Out-Host
    throw "Inventory failed with exit code $PythonExitCode; see the Python error log above"
}
Compress-Archive -Path (Join-Path $Out '*') -DestinationPath $Zip -CompressionLevel Optimal
if (-not (Test-Path -LiteralPath $Zip -PathType Leaf)) { throw 'Report archive not created' }
Write-Host "Report folder: $Out" -ForegroundColor Green
Write-Host "Report ZIP: $Zip" -ForegroundColor Green
Write-Host 'scientific_acceptance=NOT_APPROVED; merge_permitted=false'