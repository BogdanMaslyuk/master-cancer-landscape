from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from mcl_api.pharmacology import MCLPharmacologyStore


def _write_runtime(root: Path) -> None:
    runtime = root / "data" / "runtime" / "pharmacology"
    runtime.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        [
            {
                "compound_id": "CMP-1",
                "preferred_name": "Example inhibitor",
                "canonical_smiles": "CCO",
                "chembl_id": "CHEMBL1",
            }
        ]
    ).to_parquet(runtime / "compounds.parquet", index=False)
    pd.DataFrame(
        [
            {
                "observation_id": "OBS-1",
                "model_id": "ACH-000001",
                "compound_id": "CMP-1",
                "source": "TEST",
                "assay_type": "viability",
                "endpoint": "IC50",
                "value": 0.25,
                "unit": "uM",
            },
            {
                "observation_id": "OBS-2",
                "model_id": "ACH-000001",
                "compound_id": "CMP-1",
                "source": "TEST",
                "assay_type": "viability",
                "endpoint": "AUC",
                "value": 0.41,
                "unit": "fraction",
            },
        ]
    ).to_parquet(runtime / "responses.parquet", index=False)
    pd.DataFrame(
        [
            {
                "evidence_id": "TGT-1",
                "compound_id": "CMP-1",
                "target_gene": "EGFR",
                "action": "inhibitor",
                "evidence_type": "curated_MOA",
                "source": "TEST-TARGET",
                "confidence": "curated",
            }
        ]
    ).to_parquet(runtime / "target_evidence.parquet", index=False)
    pd.DataFrame(
        [
            {
                "model_id": "ACH-000001",
                "compound_id": "CMP-1",
                "target_gene": "EGFR",
                "action": "inhibitor",
                "evidence_type": "curated_MOA",
                "confidence": "curated",
                "source": "TEST-TARGET",
                "gene_effect": -1.2,
                "crispr_support_level": "strong_dependency",
            }
        ]
    ).to_parquet(runtime / "model_compound_target_links.parquet", index=False)
    (runtime / "manifest.json").write_text(
        json.dumps(
            {
                "contract": "mcl-pharmacology-v1",
                "schema_version": "1.0",
                "status": "available",
                "compounds_n": 1,
                "responses_n": 2,
                "models_with_response_n": 1,
                "sources": ["TEST"],
            }
        ),
        encoding="utf-8",
    )


def test_model_pharmacology_preserves_incompatible_endpoints(tmp_path: Path):
    _write_runtime(tmp_path)
    store = MCLPharmacologyStore(tmp_path)
    payload = store.model("ACH-000001", limit=100)

    assert payload["available"] is True
    assert payload["observations_n"] == 2
    assert payload["compounds_n"] == 1
    assert {row["endpoint"] for row in payload["items"]} == {"IC50", "AUC"}
    assert {row["value"] for row in payload["items"]} == {0.25, 0.41}


def test_target_annotation_and_crispr_are_separate_evidence(tmp_path: Path):
    _write_runtime(tmp_path)
    store = MCLPharmacologyStore(tmp_path)
    payload = store.model("ACH-000001", limit=100)

    target = payload["items"][0]["target_hypotheses"][0]
    assert target["target_gene"] == "EGFR"
    assert target["evidence_type"] == "curated_MOA"
    assert target["gene_effect"] == -1.2
    assert target["crispr_support_level"] == "strong_dependency"
    assert "не доказывает" in payload["interpretation"]["crispr"]


def test_missing_runtime_is_explicit_not_negative_result(tmp_path: Path):
    payload = MCLPharmacologyStore(tmp_path).model("ACH-000001")
    assert payload["available"] is False
    assert payload["status"] == "not_built"
    assert payload["observations_n"] == 0
    assert "build_pharmacology_layer.py" in payload["build_command"]
