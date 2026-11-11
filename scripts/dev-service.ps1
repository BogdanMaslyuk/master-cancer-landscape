param(
    [Parameter(Mandatory = $true)]
    [ValidateSet('api', 'web')]
    [string]$Service
)

$ErrorActionPreference = 'Stop'
$Root = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))
$env:MCL_ROOT = $Root
$env:PYTHONPATH = (Join-Path $Root 'apps\api')

if ($Service -eq 'api') {
    $Host.UI.RawUI.WindowTitle = 'MCL PYZ 2.0 | PS1 API :8100'
    Write-Host '=== PS1 | MCL PYZ 2.0 | API ===' -ForegroundColor Cyan

    $PythonCandidates = @(
        (Join-Path $Root '.venv\Scripts\python.exe'),
        (Join-Path (Split-Path $Root -Parent) 'master-cancer-landscape\.venv\Scripts\python.exe'),
        (Join-Path (Split-Path $Root -Parent) 'mcl-vnext-review\.venv\Scripts\python.exe')
    )
    $Python = $null
    foreach ($Candidate in $PythonCandidates) {
        if (Test-Path -LiteralPath $Candidate -PathType Leaf) {
            $Python = $Candidate
            break
        }
    }
    if (-not $Python) {
        throw 'Python environment not found. Install/restore the venv and retry.'
    }

    Set-Location $Root
    Write-Host "Code/data: $Root" -ForegroundColor Gray
    Write-Host "Python: $Python" -ForegroundColor Gray
    Write-Host 'API: http://127.0.0.1:8100' -ForegroundColor Green
    Write-Host 'This window belongs to MCL_CONTROL. To restart, use the control menu.' -ForegroundColor Yellow

    & $Python -m uvicorn mcl_api.main:app --app-dir (Join-Path $Root 'apps\api') --host 127.0.0.1 --port 8100
    $Code = $LASTEXITCODE
    Write-Host "API process finished, exit code $Code. Control menu can restart it." -ForegroundColor Yellow
}
else {
    $Host.UI.RawUI.WindowTitle = 'MCL PYZ 2.0 | PS2 WEB :3100'
    Write-Host '=== PS2 | MCL PYZ 2.0 | SITE ===' -ForegroundColor Cyan

    $Frontend = Join-Path $Root 'apps\explorer'
    if (-not (Test-Path -LiteralPath (Join-Path $Frontend 'package.json') -PathType Leaf)) {
        throw "Frontend files missing: $Frontend"
    }
    Set-Location $Frontend
    $env:NEXT_PUBLIC_MCL_API_URL = 'http://127.0.0.1:8100'
    Write-Host "Code: $Frontend" -ForegroundColor Gray
    Write-Host 'Site: http://127.0.0.1:3100/pyz/PYZ-035' -ForegroundColor Green
    Write-Host 'API: http://127.0.0.1:8100' -ForegroundColor Gray
    Write-Host 'This window belongs to MCL_CONTROL. To restart, use the control menu.' -ForegroundColor Yellow

    if (-not (Test-Path -LiteralPath (Join-Path $Frontend 'node_modules') -PathType Container)) {
        Write-Host 'node_modules missing: installing dependencies once (npm ci) ...' -ForegroundColor Yellow
        & npm.cmd ci
        if ($LASTEXITCODE -ne 0) { throw 'npm ci failed; see output above.' }
    }
    & npm.cmd run dev -- -p 3100
    $Code = $LASTEXITCODE
    Write-Host "Website process finished, exit code $Code. Control menu can restart it." -ForegroundColor Yellow
}
