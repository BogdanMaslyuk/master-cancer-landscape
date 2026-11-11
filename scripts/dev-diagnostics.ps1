param()

$ErrorActionPreference = 'Stop'
$Root = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..'))
$Host.UI.RawUI.WindowTitle = 'PS3 | ДИАГНОСТИКА MCL | только проверка'

Write-Host ''
Write-Host '======================================================' -ForegroundColor Yellow
Write-Host '  PS3 — ДИАГНОСТИКА MCL (НЕ СЕРВЕР)' -ForegroundColor Yellow
Write-Host '  Проверка PS1, PS2 и научных файлов без их изменения' -ForegroundColor Yellow
Write-Host '======================================================' -ForegroundColor Yellow
Write-Host "Проект: $Root" -ForegroundColor Gray

$Missing = @()
$Required = @(
    'config\own_compounds.tsv',
    'data\runtime\pyz_target_evidence\pyz_target_evidence_matrix.parquet',
    'data\runtime\own_compounds\target_ligand_space_similarity_top_hits.parquet',
    'data\runtime\explorer\gene_catalog.parquet',
    'data\runtime\explorer\gene_context_metrics.parquet',
    'data\runtime\explorer\gene_annotations.parquet',
    'data\runtime\explorer\gene_dependency_summary.parquet',
    'data\runtime\explorer\manifest.json'
)

Write-Host ''
Write-Host '[1] НАУЧНЫЕ ФАЙЛЫ' -ForegroundColor Cyan
foreach ($Relative in $Required) {
    $Path = Join-Path $Root $Relative
    $Found = Test-Path -LiteralPath $Path -PathType Leaf
    if (-not $Found -and $Relative -like '*.parquet' -and $Relative -notlike '*\explorer\*') {
        $Fallback = $Path -replace '\.parquet$', '.tsv'
        $Found = Test-Path -LiteralPath $Fallback -PathType Leaf
    }
    if ($Found) { Write-Host "OK: $Relative" -ForegroundColor Green }
    else {
        Write-Host "НЕТ: $Relative" -ForegroundColor Yellow
        $Missing += $Relative
    }
}

function Test-LocalPort([int]$Port, [string]$Title) {
    $Listener = @(Get-NetTCPConnection -LocalPort $Port -State Listen -ErrorAction SilentlyContinue)
    if ($Listener.Count -gt 0) {
        Write-Host "РАБОТАЕТ: $Title (порт $Port)" -ForegroundColor Green
        return $true
    }
    Write-Host "НЕ ЗАПУЩЕН: $Title (порт $Port)" -ForegroundColor Red
    return $false
}

Write-Host ''
Write-Host '[2] СЕРВЕРЫ' -ForegroundColor Cyan
$ApiOn = Test-LocalPort 8100 'PS1 БЭКЕНД / API'
$WebOn = Test-LocalPort 3100 'PS2 ФРОНТЕНД / САЙТ'

if ($ApiOn) {
    Write-Host ''
    Write-Host '[3] НАУЧНЫЙ API' -ForegroundColor Cyan
    try {
        $Health = Invoke-RestMethod 'http://127.0.0.1:8100/health' -TimeoutSec 10
        Write-Host "API OK; корневая папка: $($Health.mcl_root)" -ForegroundColor Green
        $ExpectedRoot = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
        $ActualRoot = [System.IO.Path]::GetFullPath([string]$Health.mcl_root).TrimEnd('\')
        if (-not [string]::Equals($ExpectedRoot, $ActualRoot, [System.StringComparison]::OrdinalIgnoreCase)) {
            Write-Host 'ВНИМАНИЕ: API читает научные файлы ДРУГОЙ папки!' -ForegroundColor Red
        }
    }
    catch { Write-Host "Ошибка /health: $($_.Exception.Message)" -ForegroundColor Red }

    try {
        $Summary = Invoke-RestMethod 'http://127.0.0.1:8100/api/pyz/summary' -TimeoutSec 20
        Write-Host "PYZ: молекул $($Summary.compounds_n); мишеней $($Summary.targets_n); пар $($Summary.matrix_rows_n)" -ForegroundColor Green
        $Card = Invoke-RestMethod 'http://127.0.0.1:8100/api/pyz/PYZ-035' -TimeoutSec 25
        Write-Host "Карточка PYZ-035: $($Card.profile_v2.schema_version)" -ForegroundColor Green
    }
    catch { Write-Host "Ошибка PYZ API: $($_.Exception.Message)" -ForegroundColor Red }
}
if ($WebOn) {
    Write-Host ''
    Write-Host '[4] САЙТ' -ForegroundColor Cyan
    try {
        $Response = Invoke-WebRequest 'http://127.0.0.1:3100/pyz/PYZ-035' -UseBasicParsing -TimeoutSec 25
        Write-Host "Сайт отвечает: HTTP $($Response.StatusCode)" -ForegroundColor Green
    }
    catch { Write-Host "Ошибка страницы PYZ-035: $($_.Exception.Message)" -ForegroundColor Red }
}
Write-Host ''
Write-Host '------------------------ ИТОГ ------------------------' -ForegroundColor Cyan
if ($Missing.Count) {
    Write-Host "Отсутствующих научных файлов: $($Missing.Count). Перезапуск не восстановит их." -ForegroundColor Yellow
} else {
    Write-Host 'Основные научные файлы на месте.' -ForegroundColor Green
}
if ($ApiOn -and $WebOn) {
    Write-Host 'PS1 и PS2 слушают свои порты. Проверьте результаты API выше.' -ForegroundColor Green
} else {
    Write-Host 'Проверяйте запуск/ошибки в соответствующих PS1 или PS2.' -ForegroundColor Yellow
}
Write-Host 'Диагностика завершена. Данные не изменялись.' -ForegroundColor Gray
