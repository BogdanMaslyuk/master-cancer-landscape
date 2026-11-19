# MCL TLS V3 / P0.7 — краткая инструкция

1. Обновить репозиторий в **отдельную** рабочую копию `mcl-tls-v3-p07-review` из ветки `feature/tls-v3-p07-source-reliability`.
2. Запустить `scripts/run-target-space-v3-p07.ps1` в Windows PowerShell 5.1 с `-Root "$env:USERPROFILE\MCL_Project"`. **По умолчанию 0 сетевых запросов.**
3. Получить `MCL_Project\reports\TLS_V3_P07_...zip` и изучить `p07_QC.tsv`, `p07_chembl_observations.tsv`, `p07_bindingdb_protein_names.tsv`, `p07_retry_plan.tsv`.
4. Только после отдельного решения оператора допускается **новый** запуск с `-EnableTargetedRetry -MaxRequests 12` для строго ограниченных запросов ChEMBL. Никогда не нужен повтор P0 / P0.6.

Используются именно существующие исходные файлы:

- `TLS_V3_P0_20261119T133253114Z_e9ccd63f.zip`
- `TLS_V3_P06_20261119T133323433Z_e9ccd63f.zip`

Если имена иные, можно передать `-P0Zip '...' -P06Zip '...'` и соответствующая контрольная сумма все равно проверяется.

**Никаких изменений 65×13, PYZ API/сайта, научных Parquet, PR merge.**