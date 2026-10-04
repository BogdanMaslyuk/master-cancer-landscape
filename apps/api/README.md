# MCL Explorer API

Read-only FastAPI layer over Master Cancer Landscape scientific outputs and materialized Explorer indexes.

The API must not recompute genome-wide DepMap statistics, pathway enrichment, or full ontology projections during interactive requests.

## Preferred local run

From the repository root:

```powershell
.\scripts\start-backend.ps1
```

Keep that PowerShell window open while using Explorer.

Open `http://127.0.0.1:8000/docs` for interactive API documentation.

Set `MCL_ROOT` only if the API is launched outside the repository layout.

## Current runtime entrypoint

The current compatibility entrypoint is:

```text
mcl_api.main_fast:app
```

It uses the same API routes as `main.py` but replaces the Gene Explorer store with a materialized-index runtime so `/api/genes/search`, `/api/genes/facets` and `/api/gene-matrix` do not rebuild ontology mappings on request.

After the clean-build verification gate is green, this behavior should be folded into the standard `mcl_api.main:app` entrypoint and `main_fast.py` removed.

## Explorer runtime indexes

Build them with:

```powershell
.\scripts\build-explorer.ps1
```

Validate only the Gene Explorer runtime contract with:

```powershell
.\.venv\Scripts\python.exe .\scripts\verify_runtime_indexes.py
```

Current contract:

```text
mcl-gene-explorer-runtime-v1
schema 1.0
```

Required materialized files:

- `data/processed/gene_explorer/gene_catalog.parquet`
- `data/processed/gene_explorer/gene_context_metrics.parquet`
- `data/processed/gene_explorer/gene_annotations.parquet`
- `data/processed/gene_explorer/manifest.json`

## Cell-model molecular profile index

The Cancer Atlas distinguishes three evidence levels:

1. patient tumour evidence;
2. molecular context used to stratify models;
3. molecular profile of an individual experimental model.

The raw DepMap mutation file is intentionally not read on every API request. `scripts/build_depmap_model_profiles.py` creates compact local processed model metadata and mutation indexes. `scripts/build-explorer.ps1` runs that builder as part of the standard Explorer build.

Until the index is present, Explorer deliberately shows only the variants that were used to assign a model to an MCL molecular context. It does not pretend those variants are the complete genetics of the cell line.

Patient cohort genomics is a separate evidence layer. Statistics across DepMap cell lines must not be interpreted as mutation frequencies in patients.

## Full verification

From the repository root:

```powershell
.\scripts\verify.ps1
```

This runs API tests, runtime-index validation, API smoke checks, TypeScript checks and the production frontend build.
