# Э1-P0-B2 — единый рабочий каталог MCL_Project

Весь код остаётся в существующем GitHub-репозитории MCL. Не требуется распаковывать отдельный архив с программой на рабочем столе. PR №7 черновой; Э0 и Э1-P0-A не приняты.

Структура:
- MCL_Project/repositories/master-cancer-landscape — основной репозиторий, Python .venv
- MCL_Project/repositories/mcl-pyz-e1-p0b-review — отдельное рабочее дерево кода аудита (внутри того же проекта)
- MCL_Project/data/raw/bindingdb — полная исходная BindingDB
- MCL_Project/reports — прежний neighbors_audit.tsv и все новые результаты

## Одна команда PowerShell (блок)

~~~powershell
$ErrorActionPreference = "Stop"
$Root = Join-Path $env:USERPROFILE "MCL_Project"
$Repo = Join-Path $Root "repositories\master-cancer-landscape"
$Review = Join-Path $Root "repositories\mcl-pyz-e1-p0b-review"
$Branch = "feature/pyz-e1-p0b-source-inventory"

if (-not (Test-Path -LiteralPath $Repo)) { throw "Не найден репозиторий MCL в MCL_Project" }
git -C $Repo fetch origin "refs/heads/$($Branch):refs/remotes/origin/$($Branch)"
if ($LASTEXITCODE -ne 0) { throw "Ошибка Git fetch" }
if (-not (Test-Path -LiteralPath $Review)) {
    git -C $Repo worktree add --detach $Review "refs/remotes/origin/$Branch"
} else {
    $Changes = @(git -C $Review status --porcelain)
    if ($LASTEXITCODE -ne 0 -or $Changes.Count -gt 0) {
        throw "Проверяемая копия содержит изменения. Обновление остановлено без удаления файлов."
    }
    git -C $Review switch --detach "refs/remotes/origin/$Branch"
}
if ($LASTEXITCODE -ne 0) { throw "Не удалось обновить проверяемую копию" }
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "$Review\scripts\run_pyz_e1_p0b_full_bindingdb.ps1" -Root $Root
if ($LASTEXITCODE -ne 0) { throw "Аудит BindingDB завершился с ошибкой" }
~~~

Сценарий сам ищет BindingDB_All*.zip или BindingDB_All*.tsv в Загрузках, MCL_Project/documents и MCL_Project/data/raw/bindingdb. Если найдены несколько выпусков, просит указать конкретный файл с параметром -BindingDBPath. Архив из Загрузок копируется в MCL_Project/data/raw/bindingdb с проверкой SHA-256; исходный файл не удаляется. При одинаковом имени и разных хэшах копирование останавливается без перезаписи.

Старый neighbors_audit.tsv находится автоматически в MCL_Project/reports. Сканер требует ровно 5 аналогов для PYZ-035/CDK4 и 5 для PYZ-055/BCL2L1. ZIP базы читается потоково; распаковывать его не требуется. Итог записывается в MCL_Project/reports/E1_P0B2_FULL_<timestamp> и MCL_E1_P0B2_RESULT_<timestamp>.zip.

Исходные научные Parquet/TSV, рабочие серверы и сайт не изменяются. Каждый эксперимент относится к известному референсному лиганду, а не к PYZ; совпадения структур требуют проверки стереохимии, комплексов и первичных публикаций. Научный статус остаётся NO_PRIMARY_ASSAYS_APPROVED.

Тесты: tests/test_pyz_e1_p0b_full_bindingdb.py. Проверка на реальной полной BindingDB выполняется только на компьютере пользователя, поскольку исходная загрузка не предоставлена в этом чате.
