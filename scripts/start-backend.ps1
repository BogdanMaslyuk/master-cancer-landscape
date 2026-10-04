param(
    [string]$HostAddress = "127.0.0.1",
    [int]$Port = 8000
)

$ErrorActionPreference = "Stop"
$Root = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$Python = Join-Path $Root ".venv\Scripts\python.exe"

if (-not (Test-Path $Python)) {
    throw "Python environment not found: $Python. Create/sync the root .venv first."
}

Set-Location $Root
Write-Host "MCL backend" -ForegroundColor Cyan
Write-Host "Root: $Root"
Write-Host "API:  http://${HostAddress}:$Port"
Write-Host "Entrypoint: mcl_api.main:app"
Write-Host "Keep this PowerShell window open while using MCL Explorer." -ForegroundColor Yellow

& $Python -m uvicorn mcl_api.main:app --app-dir (Join-Path $Root "apps\api") --host $HostAddress --port $Port
exit $LASTEXITCODE
