param(
    [switch]$SkipFrontendBuild
)

$ErrorActionPreference = "Stop"
$Root = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$Python = Join-Path $Root ".venv\Scripts\python.exe"

if (-not (Test-Path $Python)) {
    throw "Python environment not found: $Python"
}

Set-Location $Root

function Run([string]$Label, [scriptblock]$Command) {
    Write-Host "`n=== $Label ===" -ForegroundColor Cyan
    & $Command
    if ($LASTEXITCODE -ne 0) {
        throw "$Label failed with exit code $LASTEXITCODE"
    }
}

Run "Architecture v1 repository contract" { & $Python .\scripts\verify_architecture_v1.py }
Run "Python dependency constraints" { & $Python .\scripts\verify_python_constraints.py }
Run "Python dependency graph" { & $Python -m pip check }
if (Test-Path (Join-Path $Root "tests")) {
    Run "Scientific core tests" { & $Python -m pytest .\tests -q }
}
Run "API tests" { & $Python -m pytest .\apps\api\tests -q }
Run "Generated API type contract" { & $Python .\scripts\generate_frontend_api_types.py --check }
Run "Runtime index contract" { & $Python .\scripts\verify_runtime_indexes.py }
Run "Explorer API smoke checks" { & $Python .\scripts\smoke_explorer_api.py }

$Explorer = Join-Path $Root "apps\explorer"
Push-Location $Explorer
try {
    if (-not (Test-Path (Join-Path $Explorer "node_modules"))) {
        throw "Frontend dependencies are missing. Run .\scripts\bootstrap.ps1 from the repository root."
    }
    Run "Frontend typecheck" { & npm.cmd run typecheck }
    if (-not $SkipFrontendBuild) {
        $env:NEXT_PUBLIC_MCL_API_URL = "http://127.0.0.1:8000"
        Run "Frontend production build" { & npm.cmd run build }
    }
} finally {
    Pop-Location
}

Write-Host "`nMCL verification completed successfully." -ForegroundColor Green
