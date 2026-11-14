param(
    [string]$MclRoot = "$env:USERPROFILE\Desktop\mcl-pyz-v2-review",
    [switch]$Offline,
    [switch]$SkipPerformance
)
$ErrorActionPreference = 'Stop'
$Host.UI.RawUI.WindowTitle = 'PS3 | АУДИТ PYZ Э0 | ТОЛЬКО ЧТЕНИЕ'
$Script = Join-Path $PSScriptRoot 'audit_pyz_e0.py'
$Repo = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$Candidates = @(
    (Join-Path $Repo '.venv\Scripts\python.exe'),
    (Join-Path $MclRoot '.venv\Scripts\python.exe'),
    (Join-Path $env:USERPROFILE 'Desktop\master-cancer-landscape\.venv\Scripts\python.exe')
)
$Python = $Candidates | Where-Object { Test-Path -LiteralPath $_ -PathType Leaf } | Select-Object -First 1
if (-not $Python) { throw 'Python окружение не найдено. Создайте .venv или укажите доступную установленную среду.' }
if (-not (Test-Path (Join-Path $MclRoot 'config\own_compounds.tsv'))) { throw "Нет реестра PYZ в $MclRoot" }
Write-Host '=== Э0 | АУДИТ НАУЧНЫХ ДАННЫХ / ТОЛЬКО ЧТЕНИЕ ===' -ForegroundColor Cyan
Write-Host "Аудируемый MCL: $MclRoot"
Write-Host "Python: $Python"
$Arguments = @($Script, '--root', $MclRoot)
if ($Offline) { $Arguments += '--offline' }
if ($SkipPerformance) { $Arguments += '--skip-performance' }
& $Python @Arguments
$ExitCode = $LASTEXITCODE
Write-Host "Код завершения аудитора: $ExitCode" -ForegroundColor Cyan
Write-Host 'Отчёты: outputs\audits\pyz_e0 внутри аудируемого MCL.'
if ($ExitCode -eq 2) { Write-Host 'Найдены блокирующие проблемы: см. defects.md. Э0 не принят.' -ForegroundColor Yellow }
if ($ExitCode -ne 0 -and $ExitCode -ne 2) { throw "Аудитор завершился с ошибкой $ExitCode" }