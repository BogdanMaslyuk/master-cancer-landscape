from __future__ import annotations

from pathlib import Path

import pandas as pd

from mcl_api.pharmacology import MCLPharmacologyStore
from mcl_api.pharmacology_catalog import MCLPharmacologyCatalogStore


def _runtime(root: Path) -> None:
    runtime = root / "data" / "runtime" / "pharmacology"
    processed = root / "data" / "processed"
    runtime.mkdir(parents=True, exist_ok=True)
    processed.mkdir(parents=True, exist_ok=True)

    pd.DataFrame([
        {
            "compound_id": "CMP-1",
            "preferred_name": "Example inhibitor",
            "canonical_smiles": "CCO",
            "inchikey": "TESTKEY",
            "chembl_id": "CHEMBL1",
            "broad_id": "BRD-1",
            "models_n": 2,
            "observations_n": 2,
            "target_evidence_n": 1,
            "targets_n": 1,
            "sources_json": '["TEST"]',
            "endpoints_json": '["LFC"]',
            "target_genes_json": '["EGFR"]',
        }
    ]).to_parquet(runtime / "compound_catalog.parquet", index=False)

    pd.DataFrame([
        {
            "target_id": "EGFR",
            "target_gene": "EGFR",
            "compounds_n": 1,
            "target_evidence_n": 1,
            "actions_json": '["inhibitor"]',
            "evidence_types_json": '["curated_MOA"]',
            "models_n": 2,
            "compound_model_pairs_n": 2,
            "crispr_models_n": 2,
            "dependent_models_n": 1,
            "strong_dependency_models_n": 1,
        }
    ]).to_parquet(runtime / "target_catalog.parquet", index=False)

    pd.DataFrame([
        {
            "target_gene": "EGFR",
            "compound_id": "CMP-1",
            "evidence_rows_n": 1,
            "actions_json": '["inhibitor"]',
            "evidence_types_json": '["curated_MOA"]',
            "models_n": 2,
            "crispr_models_n": 2,
            "dependent_models_n": 1,
            "strong_dependency_models_n": 1,
            "preferred_name": "Example inhibitor",
            "canonical_smiles": "CCO",
            "chembl_id": "CHEMBL1",
        }
    ]).to_parquet(runtime / "target_compound_catalog.parquet", index=False)

    pd.DataFrame([
        {"compound_id": "CMP-1", "preferred_name": "Example inhibitor", "canonical_smiles": "CCO"}
    ]).to_parquet(runtime / "compounds.parquet", index=False)
    pd.DataFrame([
        {"observation_id": "OBS-1", "model_id": "ACH-1", "compound_id": "CMP-1", "source": "TEST", "endpoint": "LFC", "value": -1.2, "unit": "log2"},
        {"observation_id": "OBS-2", "model_id": "ACH-2", "compound_id": "CMP-1", "source": "TEST", "endpoint": "LFC", "value": -0.1, "unit": "log2"},
    ]).to_parquet(runtime / "responses.parquet", index=False)
    pd.DataFrame([
        {"evidence_id": "TGT-1", "compound_id": "CMP-1", "target_gene": "EGFR", "action": "inhibitor", "evidence_type": "curated_MOA", "source": "TEST"}
    ]).to_parquet(runtime / "target_evidence.parquet", index=False)
    pd.DataFrame([
        {"model_id": "ACH-1", "compound_id": "CMP-1", "target_gene": "EGFR", "gene_effect": -1.1, "crispr_support_level": "strong_dependency"},
        {"model_id": "ACH-2", "compound_id": "CMP-1", "target_gene": "EGFR", "gene_effect": -0.1, "crispr_support_level": "weak_or_none"},
    ]).to_parquet(runtime / "model_compound_target_links.parquet", index=False)
    pd.DataFrame([
        {"model_id": "ACH-1", "cell_line_name": "MODEL1", "mcl_cancer_name": "Cancer A", "mcl_organ_ru": "Орган A"},
        {"model_id": "ACH-2", "cell_line_name": "MODEL2", "mcl_cancer_name": "Cancer B", "mcl_organ_ru": "Орган B"},
    ]).to_parquet(processed / "depmap_crispr_model_atlas.parquet", index=False)


def test_compound_catalog_search_links_target(tmp_path: Path):
    _runtime(tmp_path)
    base = MCLPharmacologyStore(tmp_path)
    store = MCLPharmacologyCatalogStore(tmp_path, base)
    payload = store.compound_search(search="Example", target_gene="EGFR")

    assert payload["available"] is True
    assert payload["total"] == 1
    assert payload["items"][0]["compound_id"] == "CMP-1"
    assert payload["items"][0]["target_genes"] == ["EGFR"]


def test_target_catalog_keeps_gene_link_and_dependency_fraction(tmp_path: Path):
    _runtime(tmp_path)
    base = MCLPharmacologyStore(tmp_path)
    store = MCLPharmacologyCatalogStore(tmp_path, base)
    payload = store.target_search(search="EGFR")

    assert payload["total"] == 1
    assert payload["items"][0]["target_gene"] == "EGFR"
    assert payload["items"][0]["dependency_fraction"] == 0.5


def test_target_detail_connects_compound_and_model_without_claiming_isoform(tmp_path: Path):
    _runtime(tmp_path)
    base = MCLPharmacologyStore(tmp_path)
    store = MCLPharmacologyCatalogStore(tmp_path, base)
    payload = store.target_detail("egfr")

    assert payload["coding_gene"] == "EGFR"
    assert payload["resolution"] == "gene_mapped_protein_target"
    assert payload["compounds"][0]["compound_id"] == "CMP-1"
    assert payload["model_examples"][0]["model_id"] == "ACH-1"
    assert "изоформа" in payload["resolution_note_ru"]
