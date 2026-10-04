param(
    [string]$ApiUrl = "http://127.0.0.1:8000",
    [switch]$ForceBuild
)

$ErrorActionPreference = "Stop"
$Root = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$Explorer = Join-Path $Root "apps\explorer"

Set-Location $Explorer
$env:NEXT_PUBLIC_MCL_API_URL = $ApiUrl

if (-not (Test-Path (Join-Path $Explorer "node_modules"))) {
    throw "Frontend dependencies are missing. Run 'npm.cmd install' once in apps\explorer."
}

$BuildId = Join-Path $Explorer ".next\BUILD_ID"
if ($ForceBuild -or -not (Test-Path $BuildId)) {
    Write-Host "No production build found (or -ForceBuild requested). Building frontend..." -ForegroundColor Cyan
    if (Test-Path (Join-Path $Explorer ".next")) {
        Remove-Item -Recurse -Force (Join-Path $Explorer ".next")
    }
    & npm.cmd run build
    if ($LASTEXITCODE -ne 0) {
        throw "Frontend build failed. Do not run 'next start' until the build is fixed."
    }
}

Write-Host "MCL frontend" -ForegroundColor Cyan
Write-Host "UI:  http://127.0.0.1:3000"
Write-Host "API: $ApiUrl"
Write-Host "Keep this PowerShell window open while using MCL Explorer." -ForegroundColor Yellow

& npm.cmd run start
exit $LASTEXITCODE
