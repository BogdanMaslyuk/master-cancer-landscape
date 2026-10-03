# Master Cancer Landscape Computational Pipeline v0.4.0

Reproducible Python pipeline for the **Master Cancer Landscape** project.

## Current implementation status

This repository implements **Milestone 0 (project skeleton)** and the code for
**Milestone 1 (HGNC-first target identity normalization)**.

The pipeline intentionally separates:

`official source -> immutable raw snapshot -> normalized table -> QC/provenance -> Excel export`

The Excel workbook is a human-readable product, not the raw-data store.

## Scientific guardrails

- Canonical human target key: **HGNC approved symbol**; HGNC ID is the stable technical ID.
- Never merge targets by fuzzy text similarity alone.
- Ambiguous aliases are marked `manual_review_required`.
- Missing values are never silently converted to `0`.
- Existing DrugDev CMP / PPB3 / SwissTargetPrediction results are not used for target selection.
- Manual scientific interpretation fields are not owned by the target-normalization module.

## Quick start

```bash
python -m venv .venv
# Windows: .venv\\Scripts\\activate
# Linux/macOS: source .venv/bin/activate
pip install -e ".[dev,analysis]"

mcl status
mcl fetch-hgnc
mcl normalize-targets
mcl qc-targets
```

`mcl fetch-hgnc` downloads the official current HGNC approved and withdrawn TSV snapshots
into a date-stamped immutable directory and records SHA256 metadata.

`mcl normalize-targets` reads `data/input/targets_seed.tsv`, resolves each input using the
HGNC snapshot, writes normalized identifiers to `data/processed/`, writes field-level
provenance, and emits QC/manifest outputs.

## Reproducible/offline mode

To normalize against a pinned snapshot rather than downloading the current files:

```bash
mcl normalize-targets \
  --hgnc-tsv data/raw/hgnc/2026-09-26/hgnc_complete_set.txt \
  --withdrawn-tsv data/raw/hgnc/2026-09-26/withdrawn.txt
```

This is the preferred mode for publication/diploma reproducibility.

## Main outputs of Milestone 1

- `data/processed/target_identifiers.parquet`
- `data/processed/target_identifiers.tsv`
- `data/processed/provenance.parquet`
- `data/processed/provenance.tsv`
- `outputs/qc/target_normalization_qc.json`
- `outputs/qc/target_normalization_qc.tsv`
- `outputs/reports/run_manifest.json`

## Important

Milestone 1 normalizes **identity**, not biological priority. It does not assign target scores,
clinical value, tumour selectivity, or medicinal-chemistry priority.

## Milestone 2 — Open Targets Wave 1

Milestone 2 uses Open Targets only as an evidence aggregation layer. It does **not** convert association scores into an overall target priority and does not claim that disease-level Open Targets associations are molecular-subtype-specific.

### 1. Validate disease ontology mappings first

```powershell
mcl validate-ot-diseases
```

Expected: 4 disease mappings checked. `ERROR` must be 0 before continuing. A disease-name `WARNING` means the canonical ontology name changed and requires review; it must not be silently auto-corrected.

### 2. Fetch direct associations and tractability

```powershell
mcl fetch-opentargets
```

The command:
- requires completed Milestone 1 `target_identifiers.tsv`;
- queries direct disease-target associations (`enableIndirect: false`);
- paginates the full disease association table before filtering the 34 Wave 1 pairs;
- stores all raw API responses under `data/raw/opentargets/<release>/`;
- preserves overall score, datatype scores, datasource scores, and target tractability;
- records absent associations as `association_found=False` and **not** as score `0`;
- appends Open Targets provenance to the canonical provenance table.

### 3. QC

```powershell
mcl qc-opentargets
```

`ERROR: 0` is required before interpretation. `INFO association_absent` is not an error: absence in Open Targets must not be interpreted as experimental evidence of no biological relationship.

### Wave 1 ontology rule

The Open Targets disease entity is a disease-level context only. Molecular contexts (`KRAS G12C`, `KRAS G12D`, `IDH-wildtype`, `FLT3-mutant`) remain separate project annotations and must later be tested with subtype-aware evidence such as DepMap and primary literature.


## Open Targets disease identifier drift

Disease identifiers can change between Open Targets releases. A configured EFO ID that was valid historically may return `null` in a newer release. `mcl validate-ot-diseases` therefore:

1. checks the configured ID;
2. if absent, resolves the curated disease name with Open Targets `mapIds`;
3. accepts only a single case-insensitive exact-name disease hit;
4. round-trips the resolved ID through `disease(efoId:)`;
5. preserves configured and resolved IDs plus raw responses in provenance.

Warnings for identifier drift are expected and are not treated as mapping errors when the exact disease name is resolved uniquely.

## Этап 3 — DepMap: функциональная зависимость

Следующий вычислительный этап работает только с **зафиксированным официальным релизом DepMap Public 26Q1**, скачанным локально. Интерфейс портала не очищается и не скрейпится.

Нужны пять файлов одного релиза в каталоге `data/raw/depmap/26Q1/`:

- `Model.csv` или `Models.csv`;
- `CRISPRGeneEffect.csv`;
- `CRISPRGeneDependency.csv`;
- `OmicsSomaticMutations.csv`;
- `OmicsProfiles.csv`.

`OmicsProfiles.csv` нужен принципиально: отсутствие строки мутации не считается автоматически отсутствием мутации, если для модели не подтверждено WES/WGS-секвенирование.

Порядок запуска:

```powershell
mcl validate-depmap-inputs
mcl analyze-depmap
mcl qc-depmap
```

Перед расчётом обязательно требуется `ERROR: 0` на этапе `validate-depmap-inputs`.

Основные результаты:

- `data/processed/depmap_context_audit.tsv` — прозрачный состав основной и контрольной групп с причиной включения/исключения каждой модели;
- `data/processed/depmap_evidence.tsv` — количественные показатели для всех 34 пар первой волны;
- `outputs/qc/depmap_qc.tsv` — контроль качества;
- `outputs/reports/run_manifest_depmap.json` — релиз, хэши, конфигурация и воспроизводимость запуска.

Статистическая схема первой версии:

- основной показатель — Chronos Gene Effect;
- бинарная зависимость — `CRISPRGeneDependency > 0.5`;
- контроль — та же опухоль без определяющего молекулярного признака;
- двусторонний критерий Манна—Уитни;
- величина эффекта — дельта Клиффа;
- поправка Бенджамини—Хохберга по всем вычислимым сравнениям Wave 1;
- предупреждение о малой выборке при `n < 5` хотя бы в одной группе.

Подробная русская методология этапа: `docs/MILESTONE3_DEPMAP_RU.md`.

## DepMap 26Q1: важное замечание по глиобластоме

В релизе 26Q1 модели ЦНС были переклассифицированы по обновлённому OncoTree. Для CANCER-011 используется точная текущая аннотация `Glioblastoma, IDH-Wildtype` / код `GB` / линия `CNS/Brain`. IDH1 R132 и IDH2 R172 проверяются как контроль согласованности. Если внутри этого точного подтипа нет IDH-мутантных моделей, программа не создаёт искусственную контрольную группу и не заявляет мутационно-специфическую статистику.



## DepMap M3.1 — KRAS comparator sensitivity

M3.1 is a secondary sensitivity analysis and does **not** replace the primary M3 same-disease comparator. For LUAD KRAS p.G12C and pancreatic adenocarcinoma KRAS p.G12D, it separates the non-target disease cohort into:

- sequencing-backed `kras_wildtype_proxy`: default WES/WGS profile and no KRAS variant record;
- `other_kras_driver_hotspot`: non-target KRAS variant with a recognized DepMap driver/hotspot annotation;
- ambiguous non-target KRAS variants are excluded rather than silently called WT.

Run:

```powershell
mcl analyze-depmap-kras-sensitivity
mcl qc-depmap-kras-sensitivity
```

Outputs:

- `data/processed/depmap_kras_sensitivity.tsv`
- `data/processed/depmap_kras_sensitivity_audit.tsv`
- `outputs/qc/depmap_kras_sensitivity_qc.tsv`
- `outputs/qc/depmap_kras_sensitivity_meta.json`
- `outputs/reports/run_manifest_depmap_kras_sensitivity.json`

## DepMap synchronization (`mcl sync-depmap`)

The project now has a conservative synchronization layer for the pinned DepMap release:

```powershell
mcl sync-depmap
```

It:
- queries the official no-CAPTCHA DepMap metadata catalog when online;
- keeps the configured release pinned and only reports if a newer release exists;
- checks the five required local files;
- performs fast header/schema sanity checks;
- records file size and SHA256;
- writes `data/processed/depmap_sync_inventory.tsv` and `outputs/reports/depmap_sync_manifest.json`;
- never scrapes or bypasses CAPTCHA-protected file downloads.

If raw-file download URLs become available through official bearer-token authentication, set `DEPMAP_BEARER_TOKEN` and the same command can use those authorized URLs. Until then, only missing files must be downloaded once through the official DepMap portal and placed in the pinned cache directory.

Offline cache verification is available with:

```powershell
mcl sync-depmap --offline
```

## Этап 3.2 — Genome-wide Dependency Explorer

M3.2 расширяет DepMap-анализ от заранее выбранных мишеней до **всех генов CRISPR-матрицы релиза**. Это отдельный discovery-слой и он не заменяет M3/M3.1.

Для одного опухолевого контекста программа:

1. воспроизводимо формирует те же молекулярные когорты, что M3/M3.1;
2. одним проходом по `CRISPRGeneEffect.csv` сохраняет только нужные клеточные модели, но анализирует все гены;
3. одним проходом по `CRISPRGeneDependency.csv` одновременно рассчитывает зависимости выбранных моделей и broad-dependency fraction по всей панели;
4. для каждого гена сохраняет отдельные наблюдаемые оси: медианный Gene Effect, долю зависимых моделей, Δ между группами, Cliff's delta, p/q и broad-dependency flag;
5. применяет Benjamini–Hochberg внутри одного genome-wide сравнения;
6. не создаёт общего «target score»;
7. строит автономный офлайн HTML Explorer.

### Запуск

Основной контроль M3:

```powershell
mcl analyze-depmap-genome-wide --cancer-id CANCER-001 --comparison primary
```

Для точных KRAS-контекстов доступны дополнительные сравнения M3.1:

```powershell
mcl analyze-depmap-genome-wide --cancer-id CANCER-001 --comparison kras-wt
mcl analyze-depmap-genome-wide --cancer-id CANCER-001 --comparison other-kras
mcl analyze-depmap-genome-wide --cancer-id CANCER-004 --comparison kras-wt
mcl analyze-depmap-genome-wide --cancer-id CANCER-004 --comparison other-kras
```

QC:

```powershell
mcl qc-depmap-genome-wide --cancer-id CANCER-001 --comparison kras-wt
```

### Основные выходы

Для ключа `<Cancer_ID>__<comparison>`:

- `data/processed/depmap_genomewide/<key>_genes.tsv` — человекочитаемая genome-wide таблица;
- `data/processed/depmap_genomewide/<key>_genes.parquet` — канонический аналитический формат;
- `data/processed/depmap_genomewide/<key>_cohort.tsv` — аудит состава групп;
- `data/processed/depmap_genomewide/<key>_gene_effect_matrix.parquet` — Gene Effect по выбранным моделям и всем генам;
- `data/processed/depmap_genomewide/<key>_dependency_matrix.parquet` — Probability of Dependency по выбранным моделям и всем генам;
- `outputs/reports/depmap_explorer_<key>.html` — автономная интерактивная визуализация;
- `outputs/reports/depmap_genomewide_<key>_meta.json` — релиз, хэши и параметры анализа;
- `outputs/qc/depmap_genomewide_<key>_qc.tsv` — контроль качества.

### Что показывает Explorer

- scatter «медианный Gene Effect в контроле vs контексте» по всем генам;
- тепловую карту наиболее контекстно-селективных зависимостей по отдельным клеточным линиям;
- карточку выбранного гена;
- поиск и фильтрацию по всем генам;
- FDR и broad-dependency как отдельные признаки.

**Интерпретационное ограничение:** CRISPR-нокаут — это генетическая зависимость, а не доказательство того, что белок уже является пригодной, безопасной и фармакологически достижимой мишенью малой молекулы.

## Milestone 3.3 — pathway enrichment

After the required genome-wide DepMap comparisons have been generated:

```powershell
mcl analyze-pathways
mcl qc-pathways
```

M3.3 reconstructs per-comparison top candidates from the canonical genome-wide tables, derives recurrent/core candidate sets, builds eligibility-matched custom statistical backgrounds, and queries g:Profiler for GO Biological Process, Reactome, KEGG, and CORUM enrichment. Source-version metadata and raw API responses are retained for provenance. Before enrichment, Entrez Gene IDs are explicitly normalized with g:Convert to one canonical Ensembl gene per MCL/HGNC symbol; g:GOSt receives only ENSG query/background identifiers. The resulting `data/processed/pathways/identifier_resolution.tsv` preserves the Entrez→ENSG decision and any unresolved ambiguity for audit.
