# Master Cancer Landscape — Architecture v1 acceptance audit

Date: 2026-10-04
Branch: `mcl-explorer-v0.1`
Target tag: `mcl-architecture-v1`

## Result

Architecture v1 is accepted as the technical baseline for the next scientific expansion stage of MCL, subject to the final tagged commit passing the same verification gate.

## Accepted invariants

- Scientific computation is separated from interactive serving.
- `data/processed/` remains scientific/build-stage output; `data/runtime/` is disposable serving state.
- Explorer runtime root is `data/runtime/explorer/`.
- Major API domains follow `router -> service -> repository -> store/runtime data`.
- HTTP routers do not directly own scientific dataframe/statistics computation or request caches.
- Missing runtime artifacts fail fast instead of silently rebuilding scientific layers.
- Wide DepMap model-by-gene Parquet matrices are not used as an ordinary arbitrary-gene runtime fallback.
- Stable Explorer GET endpoints expose named Pydantic response models.
- Frontend API types are generated mechanically from FastAPI OpenAPI.
- Python 3.13 dependencies are constrained by `constraints/python-3.13.txt`.
- Frontend dependencies are locked by `apps/explorer/package-lock.json` and installed with `npm ci`.
- Windows clean bootstrap is supported by `scripts/bootstrap.ps1`, including `uv`-managed Python 3.13 with seeded `pip`.

## Verification evidence

The accepted local gate includes:

```text
Architecture v1 repository contract
Python dependency constraints
pip check
scientific core tests
API tests
OpenAPI -> TypeScript contract check
runtime index contract
Explorer API smoke checks
frontend typecheck
Next.js production build
```

During Architecture v1 validation, arbitrary-gene cold loading was reduced from tens of seconds to sub-second/low-single-second API latency by introducing memory-mapped `gene x model` runtime arrays. This optimization changes serving mechanics only, not scientific semantics.

The Windows clean bootstrap path was tested from a recreated `.venv`. Two Windows-specific bootstrap defects were found and fixed during acceptance testing: Python Launcher discovery and `uv venv` environments without seeded `pip`.

GitHub CI on the pre-finalization architecture checkpoint was green for both Python and frontend jobs. The final tagged commit must also pass CI before the tag is treated as the frozen baseline.

## Scientific guardrails preserved

- CRISPR knockout dependency is not pharmacological inhibition.
- RNA expression is not dependency or protein activity.
- Relative copy number is not a clinical amplification/deletion call.
- Cell-line observations are not patient prevalence.
- Missing values are not zero.
- Correlation is not causation.
- Mutation association is exploratory and does not establish mutation-caused dependency.
- Broad dependency and low sample size remain explicit flags.
- No opaque global target score is introduced by Architecture v1.

## Baseline freeze rule

After the final full local verification and green CI:

```powershell
git tag -a mcl-architecture-v1 -m "MCL Architecture v1 baseline"
git push origin mcl-architecture-v1
```

Do not merge `mcl-explorer-v0.1` into `main` automatically. Tagging freezes the architecture checkpoint on the development branch; scientific expansion may continue from that baseline.

## Next stage after freeze

Architecture work becomes maintenance rather than the main project focus. The next major workstream is scientific expansion:

1. cancer-context expansion;
2. Gene x Cancer evidence matrix;
3. co-dependency and mechanism layers;
4. patient-level evidence;
5. target safety/druggability;
6. integration with molecule-level discovery workflows.
