# E0 evidence assessment: one MCL_Project workspace, no Git merge, no science modifications.
param(
    [string]$Root = (Join-Path $env:USERPROFILE 'MCL_Project'),
    [string]$InputPath = '',
    [string]$ScienceRoot = ''
)
$ErrorActionPreference = 'Stop'
$Repo = Join-Path $Root 'repositories\master-cancer-landscape'
$Python = Join-Path $Repo '.venv\Scripts\python.exe'
$Audit = Join-Path $PSScriptRoot 'audit_pyz_e0_acceptance_gate.py'
$Reports = Join-Path $Root 'reports'
foreach ($p in @($Root, $Repo, $Python, $Audit, $Reports)) {
    if (-not (Test-Path -LiteralPath $p)) { throw "Не найден обязательный путь MCL: $p" }
}

if (-not $InputPath) {
    $PossibleZips = @()
    $Locations = @($Reports, (Join-Path $env:USERPROFILE 'Desktop'))
    foreach ($loc in $Locations) {
        if (Test-Path -LiteralPath $loc) {
            $PossibleZips += @(Get-ChildItem -LiteralPath $loc -File -Filter 'MCL_PYZ_E0_CLOSEOUT*.zip' -ErrorAction SilentlyContinue)
        }
    }
    $PossibleZips = @($PossibleZips | Sort-Object FullName -Unique)
    if ($PossibleZips.Count -eq 1) {
        $InputPath = $PossibleZips[0].FullName
    } elseif ($PossibleZips.Count -gt 1) {
        $Names = ($PossibleZips | ForEach-Object { $_.FullName }) -join [Environment]::NewLine
        throw "Найдено несколько архивов Э0. Укажите -InputPath (ни один не выбран): $Names"
    } else {
        $AuditDirs = @(
            (Join-Path $Root 'repositories\mcl-pyz-e0-audit\outputs\audits\pyz_e0'),
            (Join-Path $env:USERPROFILE 'Desktop\mcl-pyz-e0-audit\outputs\audits\pyz_e0')
        )
        $Candidates = @()
        foreach ($dir in $AuditDirs) {
            if (Test-Path -LiteralPath $dir) {
                $Candidates += @(Get-ChildItem -LiteralPath $dir -Directory -ErrorAction SilentlyContinue |
                    Where-Object { Test-Path -LiteralPath (Join-Path $_.FullName 'E0_ACCEPTANCE.md') })
            }
        }
        $Candidates = @($Candidates | Sort-Object LastWriteTime -Descending)
        if ($Candidates.Count -eq 0) {
            throw 'Не найден исходный отчёт Э0. Укажите -InputPath к MCL_PYZ_E0_CLOSEOUT.zip либо к папке отчёта Э0.'
        }
        $InputPath = $Candidates[0].FullName
        Write-Host "Архив CLOSEOUT не обнаружен; выбран последний отчёт Э0: $InputPath" -ForegroundColor Yellow
    }
}
if (-not (Test-Path -LiteralPath $InputPath)) { throw "Не найден источник Э0: $InputPath" }
$Stamp = Get-Date -Format 'yyyyMMdd_HHmmss_fff'
$Out = Join-Path $Reports ('E0_ACCEPTANCE_GATE_' + $Stamp)
Write-Host '=== MCL / Э0: СВОДНАЯ ПРОВЕРКА ПРИЁМКИ ===' -ForegroundColor Cyan
Write-Host "Данные Э0: $InputPath"
if (-not $ScienceRoot) {
    foreach ($candidate in @(
        (Join-Path $Root 'repositories\mcl-pyz-v2-review'),
        (Join-Path $env:USERPROFILE 'Desktop\mcl-pyz-v2-review')
    )) {
        if (Test-Path -LiteralPath (Join-Path $candidate 'config\own_compounds.tsv')) {
            $ScienceRoot = $candidate
            break
        }
    }
}
$Arguments = @($Audit, '--e0', $InputPath, '--out-dir', $Out)
if ($ScienceRoot) {
    Write-Host "Повторная SHA-проверка исходников: $ScienceRoot"
    $Arguments += @('--science-root', $ScienceRoot)
} else {
    Write-Warning 'Исходная научная папка не найдена; перепроверены только SHA из отчётов Э0.'
}
& $Python @Arguments
if ($LASTEXITCODE -ne 0) { throw 'Ошибка проверки отчёта Э0 (исходные данные не менялись)' }
$ReportPath = Join-Path $Out 'gate_report.json'
$Report = Get-Content -LiteralPath $ReportPath -Raw -Encoding UTF8 | ConvertFrom-Json
if ($Report.merge_permitted -ne $false -or $Report.criteria_count -ne 9) {
    throw 'Нарушен консервативный контракт аудита Э0'
}
$Files = @(
    $ReportPath,
    (Join-Path $Out 'acceptance_matrix.tsv'),
    (Join-Path $Out 'findings.tsv'),
    (Join-Path $Out 'source_hash_reconciliation.tsv'),
    (Join-Path $Out 'E0_ACCEPTANCE_GATE_RU.md')
)
$ZipPath = Join-Path $Reports ('MCL_E0_ACCEPTANCE_GATE_' + $Stamp + '.zip')
Compress-Archive -LiteralPath $Files -DestinationPath $ZipPath -CompressionLevel Fastest
if (-not (Test-Path -LiteralPath $ZipPath -PathType Leaf)) { throw 'Не сформирован ZIP отчёта Э0' }
Write-Host "Критериев AC: $($Report.criteria_count)" -ForegroundColor Green
Write-Host "Технический итог: $($Report.final_status)" -ForegroundColor Yellow
Write-Host 'Научная приёмка: НЕ УТВЕРЖДЕНА; слияние PR запрещено'
Write-Host "Отчёт создан: $ZipPath" -ForegroundColor Green