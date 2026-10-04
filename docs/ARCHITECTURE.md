# Master Cancer Landscape — Architecture

Status: clean-baseline architecture for `mcl-explorer-v0.1`.

## 1. Core rule

MCL separates scientific computation from interactive serving:

```text
official sources
    -> immutable/raw snapshots
    -> normalized scientific tables
    -> QC + provenance
    -> materialized Explorer indexes
    -> read-only API
    -> MCL Explorer UI
```

Interactive HTTP requests must not rebuild ontologies, recompute enrichment, or re-run genome-wide statistics.

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
data/processed/   canonical scientific outputs
```

Explorer materialized indexes currently live under:

```text
data/processed/gene_explorer/
```

They are rebuildable runtime artifacts and are ignored by Git. A future migration may move them to `data/runtime/explorer/`; until then the contract, not the folder name, defines their role.

### API — `apps/api/`

Read-only presentation layer over processed outputs and materialized indexes.

The API may:

- read Parquet/TSV/JSON/YAML;
- filter and paginate;
- join already-materialized lightweight tables;
- perform explicitly designed small model-level exploratory calculations.

The API must not:

- download external resources;
- rebuild GO/Reactome/KEGG/CORUM projections during requests;
- run full genome-wide analyses during requests;
- mutate scientific source tables.

### Explorer — `apps/explorer/`

Next.js scientific decision-support interface.

The UI presents separate evidence axes rather than a hidden target score. Missing data remain missing.

## 3. Build-time vs run-time boundary

### Build time

Examples:

- DepMap model profile indexing;
- multi-omics extraction;
- human-gene reference snapshot;
- MCL Functional Domain projection;
- Gene Explorer catalog generation;
- Gene x comparison metrics;
- annotation provenance materialization.

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
- model browsing.

Run-time endpoints must prefer materialized indexes and fail clearly when required indexes are missing or incompatible.

## 4. Gene Explorer runtime contract

Current contract:

```text
index_contract = mcl-gene-explorer-runtime-v1
schema_version = 1.0
```

Required files:

```text
gene_catalog.parquet
gene_context_metrics.parquet
gene_annotations.parquet
manifest.json
```

Validate with:

```powershell
.\.venv\Scripts\python.exe .\scripts\verify_runtime_indexes.py
```

The manifest is part of the runtime API contract, not decorative metadata.

## 5. Scientific guardrails

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

## 6. Local operation

Use three stable roles when debugging manually:

```text
PowerShell #1  backend
PowerShell #2  frontend
PowerShell #3  diagnostics
```

Normal startup should use scripts instead of memorizing long commands:

```powershell
# PowerShell #1
.\scripts\start-backend.ps1

# PowerShell #2
.\scripts\start-frontend.ps1
```

`start-frontend.ps1` automatically builds the production frontend if `.next/BUILD_ID` is absent.

## 7. Verification gate

Before adding a new major biological layer, run:

```powershell
.\scripts\verify.ps1
```

The gate includes:

- scientific-core tests when present;
- API tests;
- runtime-index schema verification;
- API smoke checks;
- frontend TypeScript check;
- frontend production build.

A feature is not considered baseline-stable until this gate passes locally.

## 8. Current compatibility layer and cleanup rule

The fast materialized-index runtime currently enters through `mcl_api.main_fast:app`. This is a compatibility layer introduced after interactive Gene Explorer requests were found to rebuild expensive annotation mappings.

It must remain until the local clean-build gate is green. After that verification, the fast runtime will be folded into the standard `mcl_api.main:app` entrypoint and `main_fast.py` removed.

Likewise, duplicate/legacy Gene Explorer modules must not be deleted until their import/call sites and tests have been audited.

## 9. Planned structural cleanup

After a green clean-build gate:

1. fold the materialized-index runtime into `mcl_api.main:app`;
2. remove `main_fast.py`;
3. remove duplicate Gene Explorer implementations after call-site audit;
4. split API into repositories, services, schemas, and routers;
5. replace broad frontend `Record<string, any>` contracts with explicit/generated API types;
6. split large gene pages into focused components;
7. pin Python/Node dependency locks;
8. make CI a required clean-build gate.

Do not perform destructive consolidation before local verification.

## 10. Long-term architecture

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

MCL should remain an evidence-navigation and hypothesis-generation system, not a black-box ranking engine.
