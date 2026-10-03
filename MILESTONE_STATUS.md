# Master Cancer Landscape — состояние реализации

## Этап 0 — каркас проекта
**ЗАВЕРШЕНО**

## Этап 1 — нормализация мишеней
**ЗАВЕРШЕНО И ПРОВЕРЕНО**

- 31/31 мишеней;
- HGNC / Ensembl / NCBI Gene / UniProt;
- контроль качества без ошибок.

## Этап 2A — сопоставление заболеваний Open Targets
**ЗАВЕРШЕНО И ПРОВЕРЕНО**

- 4/4 опухолевых контекста;
- актуальные MONDO-идентификаторы релиза Open Targets 26.09;
- 0 ошибок, 0 предупреждений.

## Этап 2B — доказательные данные Open Targets
**ЗАВЕРШЕНО И ПРОВЕРЕНО**

- 34/34 пары «контекст × мишень»;
- association score, типы и источники доказательств, аннотации пригодности;
- исходные ответы JSON, происхождение данных и контрольные суммы;
- контроль качества без ошибок и предупреждений.

## Этап 3 — DepMap: функциональная зависимость
**ЗАВЕРШЕНО И ПРОВЕРЕНО НА DEPMAP 26Q1**

- 34/34 пары «опухолевый контекст × мишень» рассчитаны;
- 4 опухолевых контекста сформированы и проаудированы;
- Chronos Gene Effect и вероятность зависимости;
- критерий Манна—Уитни, дельта Клиффа, поправка Бенджамини—Хохберга;
- QC: 0 ошибок; предупреждения отражают отсутствие корректного comparator для GBM и малую FLT3-mutant AML группу.

Команды основного этапа:

```powershell
mcl validate-depmap-inputs
mcl analyze-depmap
mcl qc-depmap
```

## Этап 3.1 — чувствительность KRAS-контроля
**КОД РЕАЛИЗОВАН; ГОТОВ К ЗАПУСКУ НА ТЕХ ЖЕ ФАЙЛАХ 26Q1**

Цель — не заменять первичный анализ, а проверить устойчивость выводов для двух KRAS-контекстов:

- LUAD KRAS p.G12C;
- pancreatic adenocarcinoma KRAS p.G12D.

Для каждой из 17 соответствующих пар «контекст × мишень» выполняются два дополнительных сравнения:

1. точный KRAS-вариант vs KRAS WT-прокси;
2. точный KRAS-вариант vs другие KRAS driver/hotspot варианты.

KRAS WT-прокси определяется только среди моделей с WES/WGS и без записи о KRAS-варианте. Неоднозначные KRAS-варианты не считаются WT и исключаются.

Команды:

```powershell
mcl analyze-depmap-kras-sensitivity
mcl qc-depmap-kras-sensitivity
```

Ожидаемый основной результат: `data/processed/depmap_kras_sensitivity.tsv` (34 строки при наличии обоих типов comparator).

## Следующий научный узел

После интерпретации M3.1 — переход к структурной пригодности мишени (RCSB PDB) и известной химии (ChEMBL/BindingDB).

## v0.3.6 — DepMap sync layer

Added a reproducible synchronization gate before DepMap validation/analysis:

- `mcl sync-depmap` checks the pinned release cache;
- official metadata discovery uses `https://depmap.org/portal/api/no-captcha/download/files`;
- newer releases are reported but never silently replace the pinned release;
- five required files receive fast header checks, size and SHA256 recording;
- synchronization inventory is written to `data/processed/depmap_sync_inventory.tsv`;
- manifest is written to `outputs/reports/depmap_sync_manifest.json`;
- `--offline` validates the local cache without network access;
- optional `DEPMAP_BEARER_TOKEN` support is ready for a future official URL-bearing authenticated DepMap download catalog;
- CAPTCHA-protected downloads are never scraped or bypassed.

Execution order for DepMap is now:

`sync-depmap -> validate-depmap-inputs -> analyze-depmap -> qc-depmap -> analyze-depmap-kras-sensitivity -> qc-depmap-kras-sensitivity`.

## Milestone 3.2 — Genome-wide Dependency Explorer

**IMPLEMENTED in v0.4.0.**

The module expands context-specific DepMap analysis from the curated Wave 1 target list to every gene present in the pinned CRISPR matrices, while retaining explicit cohort rules, effect sizes, uncertainty, FDR, broad-dependency screening and provenance. It also produces an offline interactive HTML explorer.
