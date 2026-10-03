$ErrorActionPreference = "Stop"

Write-Host "[MCL] Milestone 3 — DepMap Public 26Q1"
Write-Host "1/3 Проверка входных файлов"
mcl validate-depmap-inputs

Write-Host "2/3 Расчёт функциональных зависимостей"
mcl analyze-depmap

Write-Host "3/3 Контроль качества"
mcl qc-depmap

Write-Host "Готово. Сначала просмотрите data/processed/depmap_context_audit.tsv, затем depmap_evidence.tsv."
