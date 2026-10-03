# Changelog

## 0.3.5 — DepMap M3.1 KRAS comparator sensitivity

- Added a secondary, non-destructive sensitivity analysis for the two exact KRAS contexts: LUAD KRAS p.G12C and pancreatic adenocarcinoma KRAS p.G12D.
- The original M3 comparator is preserved as the primary analysis. M3.1 separately compares the exact KRAS variant against (1) a sequencing-backed KRAS WT proxy and (2) other KRAS driver/hotspot variants.
- KRAS WT proxy requires a default WES/WGS profile and no KRAS variant record; non-target KRAS variants lacking a recognized DepMap driver/hotspot annotation are excluded as ambiguous rather than called WT.
- Added Mann–Whitney U, Cliff's delta, dependency fractions, BH q-values, sample-size warnings, per-model audit, QC and run metadata for M3.1.
- Added CLI commands `analyze-depmap-kras-sensitivity` and `qc-depmap-kras-sensitivity`.
- Primary M3 evidence and QC files are not overwritten.
- Test suite: 30 passed; synthetic CLI end-to-end run completed with QC ERROR=0.

## v0.3.4
- Исправлена семантика low_sample_size: отсутствие контрольной группы больше не маркирует большую исследуемую группу как малую выборку.
- Отсутствие comparator остаётся отдельным предупреждением comparator_nonempty.

# Changelog

## 0.3.3
- Fixed `mcl analyze-depmap`: CLI now calls `depmap_qc` with the correct `expected_pairs_n` argument.
- Added regression coverage for the CLI/QC interface mismatch that previously raised `TypeError` after successful DepMap input validation.

# v0.3.2 — DepMap 26Q1 CNS ontology correction

- Updated CANCER-011 to the current DepMap 26Q1 CNS/Brain annotation: OncoTree code `GB`, subtype `Glioblastoma, IDH-Wildtype`, lineage `CNS/Brain`.
- Retained IDH1 R132 / IDH2 R172 mutation checking as a concordance safeguard rather than pretending that an IDH-mutant glioblastoma is a current same-disease comparator.
- Validation now emits explicit pre-specified low-sample warnings for context or comparator cohorts below n=5.
- No broadening of CANCER-011 to all CNS/Brain or all diffuse glioma models.
- Test suite: 26 passed.

# v0.3.1 — DepMap input diagnostics hotfix

- Recognize official `IsDefaultEntryForModel` in `OmicsProfiles.csv`.
- Print cancer-context IDs and observed counts during `validate-depmap-inputs`.
- Add per-context diagnostic breakdown: disease metadata matches, sequencing coverage, alteration statuses, context/comparator/excluded counts.
- No scientific context rules were broadened automatically.

# Changelog

## 0.3.0 — DepMap Wave 1 functional-dependency module

- Added pinned local DepMap Public 26Q1 workflow; portal scraping is not used.
- Added required input inventory for Model/Models, CRISPRGeneEffect, CRISPRGeneDependency, OmicsSomaticMutations and OmicsProfiles.
- Added WES/WGS coverage gate so missing mutation rows are not silently treated as wild type.
- Added explicit context rules for KRAS G12C LUAD, KRAS G12D pancreatic adenocarcinoma, IDH-wildtype glioblastoma and FLT3-mutant AML.
- Added auditable per-model context/comparator/exclusion table with reasons and relevant variants.
- Added Chronos and dependency-probability summaries for all 34 Wave 1 Cancer×Target pairs.
- Added same-disease comparator statistics: median difference, Mann–Whitney U, Cliff's delta, p-value and global Wave-1 Benjamini–Hochberg q-value.
- Added pre-specified low-sample warning at n < 5 in either comparison group.
- Added broad dependency fraction as a screening warning only; it is not an official common-essential or safety label.
- Added DepMap provenance, input SHA256 hashes, manifest and QC.
- Added `validate-depmap-inputs`, `analyze-depmap` and `qc-depmap` commands.
- Test suite: 24 passed.

## 0.2.2 — Wave-1 Open Targets ontology mapping correction

- Replaced stale/over-broad Wave-1 disease mappings with current MONDO disease-level entities.
- CANCER-001: `MONDO_0005061` — lung adenocarcinoma (API-confirmed in Open Targets 26.09).
- CANCER-004: `MONDO_0005184` — pancreatic ductal adenocarcinoma; KRAS G12D remains separate molecular context.
- CANCER-011: `MONDO_0018177` — glioblastoma; IDH-wildtype remains separate molecular context.
- CANCER-013: `MONDO_0018874` — acute myeloid leukemia (API-confirmed in Open Targets 26.09).
- Legacy EFO IDs retained in config for provenance/documentation.
- PDAC mapping upgraded from broader pancreatic-carcinoma proxy to exact disease-level mapping.

## 0.2.1 — Open Targets disease-ID drift hotfix

- Added Open Targets API/data release metadata query.
- Added runtime disease-name resolution via `mapIds` when a configured disease ID returns `null`.
- Resolver accepts only one case-insensitive exact-name disease hit; ambiguous/fuzzy hits require review.
- Preserves configured historical disease ID and release-resolved disease ID in provenance.
- `validate-ot-diseases` now reports configured ID => resolved ID and the live API/data release.
- `fetch-opentargets` now uses release-resolved disease IDs instead of stale configured IDs.
- Replaced batch `targets(ensemblIds: ...)` tractability dependency with cached single-target `target(ensemblId: ...)` calls.
- Added identifier-drift and ambiguity regression tests.
- Test suite: 17 passed.


## 0.1.0 — 2026-09-26
- Created reproducible project skeleton.
- Added HGNC current-snapshot downloader with hashes and immutable raw-data policy.
- Added local HGNC TSV resolver for approved symbols, previous symbols, aliases and withdrawn/merged symbols.
- Added typed target/provenance/QC models.
- Added target-normalization pipeline, QC, manifest and CLI.
- Added regression tests including `COX2 -> PTGS2` and ambiguous-alias behavior.

## 0.2.0 — Milestone 2A + 2B
- Added explicit Open Targets disease mappings for the four Wave 1 cancer contexts.
- Added `cancer_target_pairs_wave1.tsv` with the 34 curated Wave 1 pairs.
- Added `validate-ot-diseases`, `fetch-opentargets`, and `qc-opentargets` CLI commands.
- Added Open Targets GraphQL source adapter with retries and immutable raw JSON snapshots.
- Added direct disease-target association extraction (`enableIndirect: false`).
- Added datatype and datasource score preservation, target tractability, provenance, and QC.
- Molecular subtype context is deliberately kept separate from Open Targets disease ontology mappings.
- Added 7 Open Targets tests including an offline end-to-end parser/raw/QC roundtrip; total tests: 15.

## 0.3.6 — DepMap synchronization layer
- Added `mcl sync-depmap`.
- Uses the official DepMap no-CAPTCHA metadata endpoint for release/file discovery.
- Keeps the configured release pinned; newer releases are reported but never substituted silently.
- Checks the five required local release files, basic CSV schema, size and SHA256.
- Writes a machine-readable synchronization manifest and TSV inventory.
- Added optional future `DEPMAP_BEARER_TOKEN` support for an official URL-bearing authenticated catalog.
- Does not scrape, solve, or bypass DepMap CAPTCHA protections.

## 0.4.0 — Genome-wide Dependency Explorer

- added M3.2 genome-wide DepMap analysis for all CRISPR genes in one molecular context;
- added memory-conscious row-chunk scanning of very wide DepMap matrices;
- added all-model broad dependency fraction without retaining the full dependency matrix;
- added primary, KRAS-WT-proxy and other-KRAS comparator modes;
- added per-gene Mann–Whitney, Cliff's delta and within-comparison Benjamini–Hochberg FDR;
- added compact genome-wide QC;
- added self-contained offline HTML explorer with scatter, heatmap, gene card and searchable table;
- no opaque aggregate target score was introduced.
