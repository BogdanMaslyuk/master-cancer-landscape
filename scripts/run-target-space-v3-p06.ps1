param(
    [string]$Root = "$env:USERPROFILE\MCL_Project",
    [string]$P0Report = '',
    [switch]$EnableExternalPilot,
    [ValidateRange(1,50)] [int]$MaxTargets = 40
)

$ErrorActionPreference = 'Stop'
$Root = [IO.Path]::GetFullPath($Root)
$Repo = Join-Path $Root 'repositories\master-cancer-landscape'
$Python = Join-Path $Repo '.venv\Scripts\python.exe'
$Script = Join-Path $PSScriptRoot 'pilot_target_space_v3_p06.py'
$P0Launcher = Join-Path $PSScriptRoot 'run-target-space-v3-p0.ps1'
$Reports = Join-Path $Root 'reports'
foreach ($File in @($Python, $Script, $P0Launcher)) {
    if (-not (Test-Path -LiteralPath $File -PathType Leaf)) { throw "Required file missing: $File" }
}
if (-not (Test-Path -LiteralPath $Reports)) { New-Item -ItemType Directory -Force -Path $Reports | Out-Null }

if (-not $P0Report) {
    Write-Host '=== 1. P0: свежая локальная техническая инвентаризация без сети ===' -ForegroundColor Cyan
    $Before = @((Get-ChildItem $Reports -Directory -Filter 'TLS_V3_P0_*' | ForEach-Object { $_.FullName }))
    powershell.exe -NoProfile -ExecutionPolicy Bypass -File $P0Launcher -Root $Root
    if ($LASTEXITCODE -ne 0) { throw 'P0 inventory failed; inspect P0 diagnostic log' }
    $Created = @(Get-ChildItem $Reports -Directory -Filter 'TLS_V3_P0_*' | Where-Object { $Before -notcontains $_.FullName })
    if ($Created.Count -ne 1) { throw "Expected exactly one new P0 report folder; got $($Created.Count)" }
    $P0Report = $Created[0].FullName
} else {
    $P0Report = [IO.Path]::GetFullPath($P0Report)
    Write-Host "Using explicitly selected P0 report: $P0Report" -ForegroundColor Cyan
}

Write-Host '=== 2. P0.6: выборка и проверка внешних источников ===' -ForegroundColor Cyan
$Commit = & git -C $PSScriptRoot rev-parse --short=8 HEAD 2>$null
if ($LASTEXITCODE -ne 0 -or -not $Commit) { $Commit='unknown' }
$Stamp = (Get-Date).ToUniversalTime().ToString('yyyyMMddTHHmmssfffZ')
$Output = Join-Path $Reports "TLS_V3_P06_${Stamp}_${Commit}"
$Zip = "$Output.zip"
$Log = Join-Path $Reports "TLS_V3_P06_DIAGNOSTIC_${Stamp}_${Commit}.log"
$ArgsList = @($Script, '--p0-report', $P0Report, '--out', $Output, '--max-targets', "$MaxTargets")
if ($EnableExternalPilot) {
    $ArgsList += '--enable-external-pilot'
    Write-Host 'EXPLICIT NETWORK OPT-IN: official UniProt, ChEMBL, BindingDB; max 50 targets / 320 requests / 50 MB' -ForegroundColor Yellow
} else {
    $ArgsList += '--offline'
    Write-Host 'OFFLINE: no network requests; only reproducible panel planning' -ForegroundColor Yellow
}
$OldPreference = $ErrorActionPreference
try {
    $ErrorActionPreference = 'Continue'
    & $Python @ArgsList 2>&1 | Tee-Object -FilePath $Log
    $ExitCode = $LASTEXITCODE
} finally {
    $ErrorActionPreference = $OldPreference
}
if ($ExitCode -ne 0) {
    Write-Host "Python error log: $Log" -ForegroundColor Yellow
    Get-Content -LiteralPath $Log -Tail 100 | Out-Host
    throw "TLS V3 P0.6 failed: $ExitCode"
}
if (-not (Test-Path -LiteralPath (Join-Path $Output 'p06_summary.json'))) { throw 'Missing pilot summary' }
Compress-Archive -Path (Join-Path $Output '*') -DestinationPath $Zip -CompressionLevel Optimal
if (-not (Test-Path -LiteralPath $Zip -PathType Leaf)) { throw 'P0.6 ZIP not created' }
Write-Host "P0 REPORT: $P0Report" -ForegroundColor Green
Write-Host "P0.6 ZIP: $Zip" -ForegroundColor Green
Write-Host 'scientific_acceptance=NOT_APPROVED; merge_permitted=false'