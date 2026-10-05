from __future__ import annotations

from pathlib import Path

import pandas as pd

from mcl_api.pharmacology import MCLPharmacologyStore
from mcl_api.pharmacology_catalog import MCLPharmacologyCatalogStore


def _runtime(root: Path) -> None:
    runtime = root / "data" / "runtime" / "pharmacology"
    processed = root / "data" / "processed"
    target_registry = processed / "target_registry"
    runtime.mkdir(parents=True, exist_ok=True)
    processed.mkdir(parents=True, exist_ok=True)
    target_registry.mkdir(parents=True, exist_ok=True)

    pd.DataFrame([
        {
            "target_gene": "EGFR",
            "source_resolution": "gene_mapped",
            "protein_mapping_status": "unique_swissprot",
            "protein_preferred_name": "Epidermal growth factor receptor",
            "protein_name_raw": "Epidermal growth factor receptor (EC 2.7.10.1)",
            "uniprot_primary_accession": "P00533",
            "uniprot_accessions_json": '["P00533"]',
            "uniprot_entry_name": "EGFR_HUMAN",
            "protein_families": "Protein kinase superfamily",
            "protein_length": 1210,
            "mapping_source": "TEST",
        }
    ]).to_parquet(target_registry / "protein_targets.parquet", index=False)

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


def test_target_catalog_separates_protein_name_from_gene(tmp_path: Path):
    _runtime(tmp_path)
    base = MCLPharmacologyStore(tmp_path)
    store = MCLPharmacologyCatalogStore(tmp_path, base)
    payload = store.target_search(search="epidermal growth factor receptor")

    assert payload["total"] == 1
    item = payload["items"][0]
    assert item["target_gene"] == "EGFR"
    assert item["protein_preferred_name"] == "Epidermal growth factor receptor"
    assert item["uniprot_primary_accession"] == "P00533"
    assert item["source_resolution"] == "gene_mapped"
    assert item["dependency_fraction"] == 0.5


def test_target_detail_connects_protein_gene_compound_and_model_without_claiming_isoform(tmp_path: Path):
    _runtime(tmp_path)
    base = MCLPharmacologyStore(tmp_path)
    store = MCLPharmacologyCatalogStore(tmp_path, base)
    payload = store.target_detail("egfr")

    assert payload["coding_gene"] == "EGFR"
    assert payload["resolution"] == "gene_mapped"
    assert payload["identity"]["protein_preferred_name"] == "Epidermal growth factor receptor"
    assert payload["identity"]["uniprot_primary_accession"] == "P00533"
    assert payload["compounds"][0]["compound_id"] == "CMP-1"
    assert payload["model_examples"][0]["model_id"] == "ACH-1"
    assert "изоформ" in payload["resolution_note_ru"]


def test_compound_target_summary_uses_protein_name_but_keeps_target_gene(tmp_path: Path):
    _runtime(tmp_path)
    base = MCLPharmacologyStore(tmp_path)
    store = MCLPharmacologyCatalogStore(tmp_path, base)
    payload = store.compound_detail("CMP-1")

    target = payload["targets"][0]
    assert target["target_gene"] == "EGFR"
    assert target["protein_preferred_name"] == "Epidermal growth factor receptor"
    assert target["uniprot_primary_accession"] == "P00533"
