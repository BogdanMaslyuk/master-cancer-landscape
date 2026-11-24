# MCL TLS V3 / P0.9 — автономный запуск в Windows PowerShell 5.1

P0.9 читает только валидированный ZIP P0.8. ChEMBL SQLite и BindingDB TSV/ZIP повторно не открываются.

```powershell
$Root = "$env:USERPROFILE\MCL_Project"
$Repo = "$Root\repositories\mcl-codex"
$P08 = "$Root\reports\TLS_V3_P08_COMBINED_VALIDATED_20261123T170023017175Z_e1661977.zip"
$P08Sha256 = '4bafe542b3599db87833c3c98bd8421bdbefe51a8bc82d2dd699dda37e6176f8'
$Python = "$Root\repositories\master-cancer-landscape\.venv\Scripts\python.exe"

powershell.exe -NoProfile -ExecutionPolicy Bypass `
  -File "$Repo\scripts\run-known-ligand-catalog-v3-p09.ps1" `
  -Root $Root `
  -Python $Python `
  -P08Zip $P08 `
  -P08ExpectedSHA256 $P08Sha256
if ($LASTEXITCODE -ne 0) { throw 'P0.9 failed' }
```

Результат создаётся в новой папке `MCL_Project\reports\TLS_V3_P09_<UTC>_<commit>\` и в одноимённом ZIP. Существующий отчёт не перезаписывается.

Перед научным использованием проверить `p09_input_verification.json`, `p09_QC.tsv`, группы дублей, очередь публикаций и ограничения в `P09_REPORT_RU.md`. Успешный запуск подтверждает воспроизводимость преобразования, но не научную корректность первичных экспериментов.

`scientific_acceptance=NOT_APPROVED`

`merge_permitted=false`
