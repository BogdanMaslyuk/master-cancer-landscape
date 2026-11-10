param([int]$BackendPort = 8001, [int]$FrontendPort = 3000)
$ErrorActionPreference = "Stop"
$Review = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$Explorer = Join-Path $Review "apps\explorer"
if (-not (Test-Path (Join-Path $Explorer "node_modules"))) {
    throw "Node dependencies missing. In review checkout, run: cd apps\explorer; npm.cmd ci"
}
$env:NEXT_PUBLIC_MCL_API_URL = "http://127.0.0.1:$BackendPort"
Set-Location $Explorer
Write-Host "MCL v2.1 phenotype UI: http://127.0.0.1:$FrontendPort/phenotype" -ForegroundColor Cyan
Write-Host "API: $env:NEXT_PUBLIC_MCL_API_URL"
& npm.cmd run dev -- -p $FrontendPort
exit $LASTEXITCODE
