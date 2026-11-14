param(
    [string]$MclRoot = "$env:USERPROFILE\Desktop\mcl-pyz-v2-review",
    [switch]$SkipBrowser,
    [switch]$SkipColdApi,
    [switch]$NoZip
)
$ErrorActionPreference = 'Stop'
$AuditRoot = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))
$MclRoot = [System.IO.Path]::GetFullPath($MclRoot)
$Reports = Join-Path $AuditRoot 'outputs\audits\pyz_e0'
$PythonFile = Join-Path $PSScriptRoot 'close_pyz_e0.py'
if (-not (Test-Path $PythonFile)) { throw 'Не найден close_pyz_e0.py: обновите ветку Э0' }
if (-not (Test-Path (Join-Path $MclRoot 'config\own_compounds.tsv'))) { throw 'Не найден реестр молекул в рабочем дереве MCL' }
if (-not (Test-Path $Reports)) { throw 'Нет исходных отчётов Э0. Сначала запустите основной аудитор.' }
$Baseline = Get-ChildItem $Reports -Directory | Where-Object { Test-Path (Join-Path $_.FullName 'E0_ACCEPTANCE.md') } |
    Sort-Object Name -Descending | Select-Object -First 1
if (-not $Baseline) { throw 'Не найдена папка отчёта с E0_ACCEPTANCE.md' }
$Candidates = @(
    (Join-Path $AuditRoot '.venv\Scripts\python.exe'),
    (Join-Path $MclRoot '.venv\Scripts\python.exe'),
    (Join-Path $env:USERPROFILE 'Desktop\master-cancer-landscape\.venv\Scripts\python.exe')
)
$Python = $Candidates | Where-Object { Test-Path -LiteralPath $_ -PathType Leaf } | Select-Object -First 1
if (-not $Python) { throw 'Python .venv не найден' }
Write-Host '=== MCL PYZ / Э0 — ЗАКРЫТИЕ АВТОМАТИЧЕСКИХ ПРОВЕРОК ===' -ForegroundColor Cyan
Write-Host "Исходники: $MclRoot"
Write-Host "Отчёт: $($Baseline.FullName)"
Write-Host "Python: $Python"
$Tools = Join-Path $Reports '_ui_tools'
$HasBrowser = -not $SkipBrowser
if ($HasBrowser -and -not (Test-Path (Join-Path $Tools 'node_modules\playwright-core'))) {
    $Npm = Get-Command npm.cmd -ErrorAction SilentlyContinue
    if (-not $Npm) {
        Write-Warning 'npm.cmd не найден: браузерные снимки не будут выполнены.'
        $HasBrowser = $false
    } else {
        Write-Host 'Устанавливаем playwright-core ТОЛЬКО в отдельную папку отчётов Э0...' -ForegroundColor Cyan
        & npm.cmd install --prefix $Tools --no-save --ignore-scripts --no-audit --no-fund playwright-core
        if ($LASTEXITCODE -ne 0) {
            Write-Warning 'Установка playwright-core не прошла; файловые/API проверки продолжатся.'
            $HasBrowser = $false
        }
    }
}
$Arg = @($PythonFile, '--root', $MclRoot, '--baseline', $Baseline.FullName, '--ui-tools', $Tools, '--overwrite')
if (-not $HasBrowser) { $Arg += '--skip-ui' }
if ($SkipColdApi) { $Arg += '--skip-cold-api' }
& $Python @Arg
$Exit = $LASTEXITCODE
if ($Exit -eq 3) { throw 'Скрипт Э0 остановился из-за ошибочного пути/состояния. Проверьте консоль.' }
$Closeout = Join-Path $Baseline.FullName 'closeout'
$Acceptance = Join-Path $Closeout 'E0_ACCEPTANCE_CLOSEOUT.md'
if (Test-Path $Acceptance) {
    Write-Host "`n=== РЕЗУЛЬТАТ ===" -ForegroundColor Cyan
    Get-Content $Acceptance -Encoding UTF8 | Select-Object -First 24
}
if (-not $NoZip) {
    $Zip = Join-Path $env:USERPROFILE 'Desktop\MCL_PYZ_E0_CLOSEOUT.zip'
    if (Test-Path $Zip) { Remove-Item -LiteralPath $Zip -Force }
    Compress-Archive -Path (Join-Path $Baseline.FullName '*') -DestinationPath $Zip -CompressionLevel Optimal
    Write-Host "`nАрхив Э0: $Zip" -ForegroundColor Green
    explorer.exe /select,"$Zip"
}
if ($Exit -eq 2) { Write-Warning 'Есть расхождения научных файлов / API. Не принимать Э0 до разбора.' }
Write-Host "Проверки завершены. Э0 не принимается автоматически. Код: $Exit" -ForegroundColor Cyan