# Master Cancer Landscape — Architecture v1

Status: Architecture v1 complete on `mcl-explorer-v0.1` as of 2026-10-04.

The Architecture v1 baseline is intended to be tagged `mcl-architecture-v1` after the final local verification on the target workstation.

## 1. Core rule

MCL separates scientific computation from interactive serving:

```text
official sources
    -> immutable/raw snapshots
    -> normalized scientific tables
    -> QC + provenance
    -> canonical processed scientific outputs
    -> runtime serving snapshot
    -> read-only API
    -> MCL Explorer UI
```

Interactive HTTP requests must not rebuild ontologies, recompute enrichment, re-run genome-wide statistics, or parse very wide scientific matrices merely to render a page.

## 2. Repository layers

### Scientific core — `src/mcl/`

Owns reproducible analysis and domain logic:

- source acquisition and pinned releases;
- normalization and identifier resolution;
- cohort construction;
- DepMap statistics;
- pathway analysis;
- QC and provenance;
- scientific exports.

This is the source of scientific truth.

### Configuration — `config/`

Owns explicit project assumptions and pinned settings:

- cancer contexts;
- thresholds;
- source versions;
- pathway comparisons;
- MCL functional taxonomy;
- field ownership.

Configuration changes are scientific changes and must be reviewable in Git.

### Data

```text
data/input/       manually curated project inputs
data/raw/         immutable source snapshots; local/ignored when large
data/interim/     rebuildable intermediate calculations
data/processed/   canonical scientific outputs and build-stage derived tables
data/runtime/     rebuildable serving artifacts optimized for Explorer latency
```

Architecture v1 defines a hard semantic boundary:

- `processed/` answers **what did the scientific pipeline produce?**
- `runtime/` answers **how can Explorer serve that result quickly?**

Runtime artifacts are derived, disposable and ignored by Git. Deleting `data/runtime/` must never delete scientific truth; it only requires rebuilding Explorer serving artifacts.

## 3. API architecture

The canonical API entrypoint is:

```text
mcl_api.main:app
```

The architectural invariant for major Explorer domains is:

```text
router -> service -> repository -> store/runtime data
```

Routers own HTTP concerns only. Services own use-case orchestration and request-level caching. Repositories own data-access composition. Stores own concrete file/data semantics and scientific guardrails.

Architecture regression tests prevent routers from importing stores directly, owning `lru_cache`, or importing heavy scientific dataframe/statistics libraries for request-time computation.

## 4. Build-time vs run-time boundary

### Build time

Examples:

- DepMap model profile indexing;
- multi-omics extraction;
- human-gene reference snapshot;
- MCL Functional Domain projection;
- Gene Explorer catalog generation;
- Gene x comparison metrics;
- annotation provenance materialization;
- transposition of wide model x gene Parquet matrices into fast gene x model runtime arrays.

Canonical command:

```powershell
.\scripts\build-explorer.ps1
```

### Run time

Examples:

- gene search;
- facets;
- Gene x Cancer matrix;
- gene card;
- pathway browsing;
- model browsing;
- explicitly requested deep analyses.

Run-time endpoints read prebuilt serving artifacts. If a required runtime artifact is missing or invalid, the API fails fast with an actionable rebuild message instead of silently falling back to expensive scientific matrices.

## 5. Explorer runtime contract

Canonical runtime root:

```text
data/runtime/explorer/
```

Current contract:

```text
index_contract = mcl-gene-explorer-runtime-v1
schema_version = 1.0
```

Required serving snapshot:

```text
data/runtime/explorer/
├── gene_catalog.parquet
├── gene_context_metrics.parquet
├── gene_annotations.parquet
├── manifest.json
├── gene_reference.parquet
├── gene_reference_terms.parquet
└── model_layers/
    ├── gene_effect.npy
    ├── gene_effect.json
    ├── expression.npy
    ├── expression.json
    ├── copy_number.npy
    └── copy_number.json
```

The `.npy` matrices are `gene x model` float32 arrays designed for memory-mapped access. They replace request-time access to extremely wide scientific Parquet matrices. Scientific source matrices remain under `data/processed/`.

Validate with:

```powershell
.\.venv\Scripts\python.exe .\scripts\verify_runtime_indexes.py
```

## 6. Processed -> runtime publication

Gene Explorer is built in two stages:

```text
processed scientific/build-stage outputs
        ↓
build_gene_explorer_index.py
        ↓
data/processed/gene_explorer/
        ↓
publish serving snapshot
        ↓
data/runtime/explorer/
```

Model-level fast arrays are generated directly from processed DepMap multi-omics matrices into `data/runtime/explorer/model_layers/`.

This preserves provenance and rebuildability while preventing UI optimizations from becoming scientific source data.

## 7. Typed API contract

Stable Explorer GET endpoints publish named Pydantic response models.

The contract chain is:

```text
Pydantic
   -> FastAPI OpenAPI
   -> generated TypeScript
   -> Next.js
```

Generated frontend types live at:

```text
apps/explorer/lib/generated/api-types.ts
```

Regenerate/check with:

```powershell
.\.venv\Scripts\python.exe .\scripts\generate_frontend_api_types.py
.\.venv\Scripts\python.exe .\scripts\generate_frontend_api_types.py --check
```

## 8. Scientific guardrails

These rules are architectural invariants:

- HGNC-approved symbols/IDs are canonical gene identity anchors.
- Missing values are not silently converted to zero.
- CRISPR knockout dependency is not equivalent to pharmacological inhibition.
- RNA expression is not equivalent to dependency or protein activity.
- Relative copy number is not automatically a clinical amplification/deletion call.
- Cell-line observations are not patient prevalence.
- Descriptive model summaries are separated from pre-specified target-vs-comparator inference.
- Broad dependency and low sample size are visible flags, not hidden penalties.
- No opaque global target score.
- Functional categories preserve provenance to formal source terms.
- Heavy exploratory mutation screens are not part of the default gene-page critical path.

## 9. Performance rule

Anything that can be deterministically derived once after data refresh should be considered for build-time materialization rather than repeated request-time computation.

Current examples:

- materialized gene catalog;
- materialized Gene x Cancer metrics;
- materialized annotation mappings;
- memory-mapped model-level multi-omics arrays.

Performance optimizations must not alter scientific meaning.

## 10. Dependency reproducibility

Architecture v1 uses two dependency locks:

- Python 3.13: `constraints/python-3.13.txt`;
- frontend: `apps/explorer/package-lock.json`.

Windows bootstrap is standardized through:

```powershell
.\scripts\bootstrap.ps1
```

For a deliberate clean rebuild of the Python environment:

```powershell
.\scripts\bootstrap.ps1 -Recreate
```

The Windows bootstrap supports `uv` and creates a seeded Python 3.13 environment with `pip` before installing the constrained dependency graph.

## 11. Local operation

Use three stable roles when debugging manually:

```text
PowerShell #1  backend
PowerShell #2  frontend
PowerShell #3  diagnostics / Git / builds
```

Normal startup:

```powershell
# PowerShell #1
.\scripts\start-backend.ps1

# PowerShell #2
.\scripts\start-frontend.ps1
```

`start-backend.ps1` launches the single canonical API entrypoint `mcl_api.main:app`.

## 12. Verification gate

Before accepting changes to the Architecture v1 baseline, run:

```powershell
.\scripts\verify.ps1
```

The gate includes:

- Architecture v1 repository contract;
- Python dependency constraints and `pip check`;
- scientific-core tests;
- API and architecture-boundary tests;
- generated OpenAPI -> TypeScript contract check;
- runtime-index contract verification;
- Explorer API smoke checks including gene detail latency;
- frontend TypeScript check;
- frontend production build.

A baseline is not considered stable until the full gate passes locally and CI is green.

## 13. Architecture v1 completion criteria

All Architecture v1 criteria are satisfied:

1. canonical router -> service -> repository boundaries for major API domains;
2. runtime serving files are isolated under `data/runtime/explorer/` from the API perspective;
3. stable API responses use explicit Pydantic schemas;
4. frontend types are mechanically generated from OpenAPI;
5. Node and Python dependency resolution is reproducible;
6. Windows clean bootstrap is documented and verified;
7. runtime indexes are rebuildable from processed scientific outputs;
8. tests prevent heavy scientific rebuilds and wide-matrix fallbacks during ordinary HTTP requests;
9. local verification and GitHub CI are green at the accepted baseline;
10. the baseline is ready for the `mcl-architecture-v1` tag before major scientific expansion resumes.

## 14. Long-term architecture

```text
Cancer context
  -> cell models
  -> genetic dependencies
  -> statistical robustness
  -> biological modules
  -> gene/target evidence
  -> druggability / safety / structure
  -> chemical tractability
  -> molecule x target x context hypotheses
  -> validation planning
```

MCL remains an evidence-navigation and hypothesis-generation system, not a black-box ranking engine.
