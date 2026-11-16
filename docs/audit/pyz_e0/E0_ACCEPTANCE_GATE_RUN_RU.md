# МCL / Э0 — запуск сводного отчёта без новых папок вне MCL_Project

## Назначение

Сценарий только читает старый отчёт Э0 и при наличии научные файлы для SHA-проверки, не перезапускает сайт/API и не меняет Parquet/TSV. **Он не утверждает Э0.** Источники на рабочем столе также могут быть прочитаны без их перемещения.

## Команда PowerShell

~~~powershell
$ErrorActionPreference = 'Stop'
$Root = "$env:USERPROFILE\MCL_Project"
$Repo = "$Root\repositories\master-cancer-landscape"
$Review = "$Root\repositories\mcl-pyz-e0-acceptance-review"
$Branch = 'feature/pyz-e0-acceptance-gate-v1'
if (-not (Test-Path $Repo)) { throw 'Главный репозиторий MCL не обнаружен' }
git -C $Repo fetch origin "refs/heads/${Branch}:refs/remotes/origin/${Branch}"
if ($LASTEXITCODE -ne 0) { throw 'Не удалось скачать ветку проверки Э0' }
if (-not (Test-Path $Review)) {
    git -C $Repo worktree add --detach $Review "refs/remotes/origin/$Branch"
} else {
    $Changes = @(git -C $Review status --porcelain)
    if ($LASTEXITCODE -ne 0 -or $Changes.Count -gt 0) {
        throw 'Аудиторская копия содержит изменения. Обновление остановлено без удаления файлов.'
    }
    git -C $Review switch --detach "refs/remotes/origin/$Branch"
}
if ($LASTEXITCODE -ne 0) { throw 'Ошибка подготовки отдельной аудиторской копии' }
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "$Review\scripts\run-pyz-e0-acceptance-gate.ps1" -Root $Root
if ($LASTEXITCODE -ne 0) { throw 'Не удалось создать отчёт Э0' }
~~~

Скрипт ищет `MCL_PYZ_E0_CLOSEOUT*.zip` сначала в `MCL_Project/reports` и на старом рабочем столе. Если zip не найден, проверяет существующие каталоги отчётов Э0 в отдельных `mcl-pyz-e0-audit` внутри проекта или старого рабочего стола. При нескольких ZIP или отсутствии отчёта завершает работу безопасно, не выбирая версию наугад.

При необходимости можно передать существующий путь: `-InputPath 'C:\путь\MCL_PYZ_E0_CLOSEOUT.zip'`. Если исходники были перенесены, можно указать `-ScienceRoot 'C:\путь\mcl-pyz-v2-review'`. Скачивать BindingDB заново или распаковывать программу не требуется.

## Ожидаемый результат

`MCL_Project/reports/MCL_E0_ACCEPTANCE_GATE_<timestamp>.zip` и соответствующая папка содержат `gate_report.json`, `acceptance_matrix.tsv`, `findings.tsv`, `source_hash_reconciliation.tsv`, `E0_ACCEPTANCE_GATE_RU.md`.

Значения `BLOCKED`, `INCOMPLETE_EVIDENCE` и `PENDING_OWNER_ACCEPTANCE` не означают неудачу формирования отчёта. Код 0 означает, что протокол создан. **Научная приёмка всегда остаётся ручной.**

## Дополнение: полный пакет свидетельств для ручной приёмки (Э0 v1.1)

Начиная с исправления после первичного анализа от 16.11.2026, ZIP содержит не только пять итоговых файлов, но и:
- `review_evidence/closeout/ui/screenshots/*.png` — девять исходных снимков;
- `review_evidence/closeout/ui/browser_results.json` и `ui_diagnostics.tsv` — постраничные причины предупреждений, HTTP, заголовок, прокрутка, клавиатура и JS;
- `review_evidence/input_inventory.tsv`, `review_evidence/performance.tsv`, `review_evidence/closeout/cold_api.json` и `review_evidence/closeout/cache_diagnostics.json` — исходные данные для AC-02/06;
- `review_evidence_manifest.json` — список вложенных материалов, размеры и SHA-256 для независимой проверки целостности.

Для пересборки **не требуется повторно запускать аудитор Э0 или браузер**: используется существующая локальная папка `mcl-pyz-e0-audit/outputs/audits/pyz_e0/...`. Повторяется только консолидатор в отдельной новой папке `MCL_Project/reports`. Исходные файлы неизменны. Предыдущий архив на 4 КБ сохраняется.
