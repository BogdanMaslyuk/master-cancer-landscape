# Target Ligand Space v2 — implementation status (P0.1/P0.2)

Branch: feature/expanded-target-ligand-space-v2, forked from feature/vnext-cell-recommendation-hardening (f6a905ee).

## Status and scientific boundary

The first deliverable is an auditable registry of uniquely mapped human proteins independent of the 13 legacy targets. This does **not** yet constitute an expanded validated ligand-binding search.

P0.3 (ChEMBL/BindingDB direct assay curation), P1.1 (indexed target search), P1.2 (new-target DepMap cell-line selection), and P1.3 (UI) remain unfinished.

The original v1 Target Ligand Space, Candidate v2, 65 PYZ records and phenotype measurements are not overwritten.

## P0.1 — baseline audit

scripts/audit_target_space_v2_baseline.py reads existing data, SHA-256 hashes pinned local input files, and counts distinct core targets and PYZ IDs. It writes only to outputs/reports/target_ligand_space_v2/.

--data-root points to the original scientific dataset. --strict fails if a required input or expected core/PYZ count is absent. Missing scientific data stays explicitly missing. Fixed regression SMILES are recorded, but this audit does not run full molecule analysis.

## P0.2 — independent human protein registry

scripts/build_expanded_target_registry.py reads the previously checked UniProt mapping and optionally an official frozen UniProtKB reviewed-human TSV (with exact release supplied by the operator). Neither Candidate v2 nor priority_for_in_vitro is used as the selection gate.

Eligibility:
- valid UniProt accession; reviewed human record and unique primary gene;
- reject duplicate conflicting gene-to-accession mappings;
- reject non-human, unreviewed, isoform, ambiguous identifiers;
- store all rejected source rows with reason_code.

A catalogue entry does **not** mean known direct ligands, cancer relevance or DepMap coverage. These three fields stay not_evaluated. No claim of proteome-wide completeness is made.

Output directory: data/processed/target_ligand_space_v2/
- protein_target_universe.tsv
- protein_target_universe_qc.tsv
- protein_target_universe_manifest.json

Every input is SHA-256 hashed in the manifest. No raw data are modified.

## Windows commands

Execute in a dedicated checkout of the new branch, with the old scientific data available separately:

~~~powershell
$Review = "$env:USERPROFILE\Desktop\mcl-vnext-review"
$Source = "$env:USERPROFILE\Desktop\master-cancer-landscape"
$Python = "$Source\.venv\Scripts\python.exe"

Set-Location $Review
$env:PYTHONPATH = "$Review\apps\api;$Review\src"
& $Python ".\scripts\audit_target_space_v2_baseline.py" --data-root $Source
& $Python -m pytest ".\tests\test_expanded_target_registry.py" -q

# Partial local seed without an external UniProt TSV:
& $Python ".\scripts\build_expanded_target_registry.py" --v1-registry "$Source\data\processed\target_registry\protein_targets.parquet"

# Pinned OFFICIAL reviewed-human TSV (operator must supply actual release):
# & $Python ".\scripts\build_expanded_target_registry.py" --v1-registry "$Source\data\processed\target_registry\protein_targets.parquet" --uniprot-tsv "C:\data\uniprot_reviewed_human.tsv" --uniprot-release "ACTUAL_RELEASE"
~~~

The official UniProt export is NOT bundled. Do not invent its release number, and do not call a MyGene mapping file an official human Swiss-Prot snapshot.

## P0.3 publication gate

Next steps before activating v2 search:
1. Download official pinned ChEMBL/BindingDB releases; record hashes/licenses.
2. Separate true independent assays from mirrored entries across providers.
3. Preserve Ki/Kd/IC50 separately with raw units and censored relation signs.
4. Reject cell-based measurements from the direct-binding evidence class.
5. Preserve duplicate structures when backed by independent measurements.
6. Precompute RDKit fingerprints under a fixed algorithm/version and benchmark memory and latency.
7. Run core-13 regression and quality controls before switching the API.

The existing /molecule-analysis route and phenotype-data routes are unchanged in this phase.
