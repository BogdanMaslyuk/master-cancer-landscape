# MCL / PYZ E1-P0-A — per-drug model observations

Scope: code-only, additive changes based on immutable `bc5a835b1395168ba960aae2ec77250f11dbca0c`. Draft implementation. E0 is not formally accepted.

## Scientific meaning

- `cell_context.items[].models[].role` is a **legacy aggregate** from Candidate v2, NOT a universal property of a cell line and NOT a drug-response outcome of PYZ. Kept unchanged for backward compatibility.
- New `legacy_aggregated_role` explicitly marks the same historical value.
- New `reference_drug_observations[]` retains **each Candidate v2 row** linked to a known compound in that model: `compound_id`, `preferred_name`, `hypothesis_id`, `candidate_priority_status`, `reference_model_role`, `prism_lfc`, `dependency_probability`, `gene_effect`, `source_layer`, `source_scope` and `measurement_owner`.
- `measurement_owner=known_reference_compound_not_pyz`. Existing `scientific_summary.own_cellular_activity` and `own_target_binding` remain `not_measured`.
- `source_scope=candidate_v2_top_five_per_role`: truncated research candidate list. Absence of an observation here cannot be read as a negative PRISM result. It does not establish primary assay linkage; no `assay_id` is fabricated.
- CRISPR dependency/gene-effect values remain a distinct genetic context and must not be averaged or merged with PRISM reference-drug responses.
- Missing/nonfinite `prism_lfc` is `null`, never zero.

## Files changed

- `apps/api/mcl_api/own_compounds.py` — explicit Candidate v2 provenance for already-present raw observations.
- `apps/api/mcl_api/pyz_profile_v2.py` — transfer complete observation array to the additive PYZ profile without role aggregation or Top-N re-filtering.
- `apps/api/mcl_api/schemas/own_compounds.py` — typed observation contract, legacy role.
- `apps/explorer/lib/pyz-profile.ts` — matching TypeScript contract.
- `apps/explorer/app/pyz/[id]/page.tsx` — per-compound expandable records under each cell model, independent CRISPR columns.
- `apps/explorer/app/pyz/[id]/profile-v2.module.css` — scoped presentation.
- `apps/api/tests/test_pyz_e1_p0a_observations.py` — synthetic opposite-drug-response and null-value tests.

## Acceptance procedure

1. `python -m pytest apps/api/tests/test_pyz_profile_v2.py apps/api/tests/test_pyz_e1_p0a_observations.py -q`.
2. `python -m pytest apps/api/tests tests -q` on a checkout with required dependencies.
3. `cd apps/explorer && npm ci && npm run typecheck && npm run build`.
4. On local **real Parquet/TSV**, compare one BCL2L1 cellular model with at least two distinct known compounds and opposing `reference_model_role` values. Check profile JSON retains both records, each drug and any `prism_lfc`, with no PYZ experimental-activity claim.
5. Review screenshots and keyboard behaviour on desktop and 390 px. Do not count this automated test as a complete UI acceptance.
6. Rehash all scientific sources against E0 `source_sha_after.tsv`. No scientific source file should differ.
7. Preserve the E0 draft PR. Do not merge E1 or approve E0 automatically.

## Explicit exclusions

Primary BindingDB/ChEMBL assay-id reconciliation, measured Ki/Kd/IC50, independent Ames/DILI, cache invalidation, mobile-wide layout fixes and performance tuning are **outside** this P0-A change. They remain separate tasks E1-P0-B, E1-P0-C, CACHE-01, UI-01 and PERF-01.
