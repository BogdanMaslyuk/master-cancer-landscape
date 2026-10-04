# Master Cancer Landscape — Architecture v1

Status: Architecture v1 migration in progress on `mcl-explorer-v0.1`.

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

Request flow:

```text
FastAPI router
    -> service
    -> repository
    -> scientific/runtime stores
    -> materialized data
```

Current Gene Explorer path follows this pattern. Other domains are migrated incrementally behind the existing API contract.

### Routers

Own HTTP concerns only:

- path/query parameters;
- HTTP status mapping;
- response model binding;
- no scientific orchestration.

### Services

Own use-case orchestration and request-level caching:

- combine repository results;
- define lightweight interactive workflows;
- keep heavy analyses out of the page critical path.

### Repositories

Own access to stores/materialized data:

- runtime Gene Explorer indexes;
- processed scientific evidence;
- no HTTP knowledge.

### Stores

Own concrete file/data semantics and scientific guardrails.

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

Run-time endpoints read prebuilt serving artifacts. If a required runtime artifact is missing, the API fails fast with an actionable rebuild message instead of recomputing scientific layers inside the HTTP request.

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
├── gene_reference.parquet              # when reference snapshot is available
├── gene_reference_terms.parquet        # when reference terms are available
└── model_layers/
    ├── gene_effect.npy
    ├── gene_effect.json
    ├── expression.npy
    ├── expression.json
    ├── copy_number.npy
    └── copy_number.json
```

The `.npy` matrices are `gene x model` float32 arrays designed for memory-mapped access. They exist because opening an arbitrary gene from a ~20k-column scientific Parquet matrix caused tens-of-seconds cold latency on Windows. Scientific source matrices remain under `data/processed/`.

Validate with:

```powershell
.\.venv\Scripts\python.exe .\scripts\verify_runtime_indexes.py
```

The manifest is part of the runtime API contract, not decorative metadata.

## 6. Processed -> runtime publication

Gene Explorer is built in two stages:

```text
processed scientific/build-stage outputs
        ↓
build_gene_explorer_index.py
        ↓
data/processed/gene_explorer/     build-stage materialization/reference snapshot
        ↓
publish serving snapshot
        ↓
data/runtime/explorer/            canonical API runtime
```

Model-level fast arrays are generated directly from processed DepMap multi-omics matrices into `data/runtime/explorer/model_layers/`.

This preserves provenance and rebuildability while preventing UI optimizations from becoming scientific source data.

## 7. Scientific guardrails

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

## 8. Performance rule

Anything that can be deterministically derived once after data refresh should be considered for build-time materialization rather than repeated request-time computation.

Current examples:

- materialized gene catalog;
- materialized Gene x Cancer metrics;
- materialized annotation mappings;
- memory-mapped model-level multi-omics arrays.

Performance optimizations must not alter scientific meaning.

## 9. Local operation

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

## 10. Verification gate

Before accepting an Architecture v1 migration step, run:

```powershell
.\scripts\verify.ps1
```

The gate includes:

- scientific-core tests;
- API tests;
- runtime-index contract verification;
- API smoke checks including gene detail latency;
- frontend TypeScript check;
- frontend production build.

A migration step is not baseline-stable until this gate passes locally.

## 11. Dependency reproducibility

Frontend dependencies are locked with `apps/explorer/package-lock.json` and CI uses:

```text
npm ci
```

Python dependency locking remains an Architecture v1 task. The declared Python version range is currently 3.12–3.13.

## 12. Architecture v1 completion criteria

Architecture v1 is complete when all of the following hold:

1. canonical router -> service -> repository boundaries for major API domains;
2. runtime serving files live only under `data/runtime/explorer/` from the API perspective;
3. API responses use explicit Pydantic schemas for stable public contracts;
4. frontend types are generated or mechanically synchronized from OpenAPI;
5. Node and Python dependency resolution is reproducible;
6. a clean checkout can rebuild runtime artifacts through documented commands;
7. tests prevent heavy scientific rebuilds from occurring during ordinary HTTP requests;
8. local verification and CI are green;
9. the baseline receives an Architecture v1 tag before major scientific expansion resumes.

## 13. Long-term architecture

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
