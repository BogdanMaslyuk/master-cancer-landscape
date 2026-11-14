"""Deterministic E0 auditor tests; fixtures are never confused with live PYZ data."""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd

MOD_PATH = Path(__file__).resolve().parents[1] / "scripts" / "audit_pyz_e0.py"
spec = importlib.util.spec_from_file_location("pyz_e0_audit", MOD_PATH)
audit = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(audit)


def registry():
    return pd.DataFrame({"own_compound_id": audit.EXPECTED_IDS, "standardized_smiles": ["CCO"] * 65})


def matrix():
    return pd.DataFrame([
        {"own_compound_id": cid, "target_gene": gene, "best_tanimoto_morgan_r2_2048": 0.4}
        for cid in audit.EXPECTED_IDS for gene in (f"GENE{i}" for i in range(13))
    ])


def test_full_65x13_is_observed_not_assumed():
    issues = []
    result, coverage = audit.inspect_data({"registry": registry(), "matrix": matrix()}, issues)
    assert len(coverage) == 845
    assert result["coverage"]["present"] == 845
    assert result["coverage"]["complete_65_by_13"] is True
    assert result["matrix"]["unique_pairs"] == 845
    assert issues == []


def test_missing_pair_and_duplicate_not_silently_zero():
    original = matrix()
    dropped = original.iloc[1:].copy()
    duplicate = pd.concat([dropped, dropped.iloc[[0]]], ignore_index=True)
    issues = []
    result, coverage = audit.inspect_data({"registry": registry(), "matrix": duplicate}, issues)
    assert result["coverage"]["complete_65_by_13"] is False
    assert result["coverage"]["missing"] == 1
    assert result["coverage"]["duplicate"] == 1
    assert any(i["code"] == "DUPLICATE_PAIRS" for i in issues)
    assert next(r for r in coverage if r["state"] == "missing")["tanimoto"] is None


def test_invalid_tanimoto_including_nan_and_out_of_range():
    m = matrix()
    for i, val in enumerate([float("nan"), -0.1, 1.1, None]):
        m.loc[i, "best_tanimoto_morgan_r2_2048"] = val
    issues = []
    result, coverage = audit.inspect_data({"registry": registry(), "matrix": m}, issues)
    assert result["matrix"]["invalid_similarity_rows"] == 4
    assert sum(x["state"] == "invalid" for x in coverage) == 4
    assert any(x["code"] == "INVALID_TANIMOTO" for x in issues)


def test_bad_registry_and_unknown_pyz_id():
    r = registry()
    r.loc[64, "own_compound_id"] = "PYZ-001"
    m = matrix()
    m.loc[0, "own_compound_id"] = "PYZ-999"
    issues = []
    result, _ = audit.inspect_data({"registry": r, "matrix": m}, issues)
    assert result["registry"]["duplicate_ids_rows"] == 2
    assert "PYZ-065" in result["registry"]["missing_ids"]
    assert "PYZ-999" in result["matrix"]["unknown_pyz_ids"]
    assert any(i["code"] == "REGISTRY_QC" for i in issues)


def test_missing_columns_does_not_become_success():
    issues = []
    result, cov = audit.inspect_data({"registry": registry(), "matrix": pd.DataFrame({"unknown": [1]})}, issues)
    assert cov == []
    assert result["matrix_status"] == "FAIL"
    assert any(i["code"] == "MISSING_COLUMN" for i in issues)


def test_conflict_of_two_reference_drugs_detected_without_averaging():
    hypotheses = pd.DataFrame({"hypothesis_id": ["H1", "H2"], "target_gene": ["CDK4", "CDK4"]})
    models = pd.DataFrame({"hypothesis_id": ["H1", "H2"],
                           "model_id": ["ACH-EXAMPLE", "ACH-EXAMPLE"],
                           "preferred_name": ["DRUG-A", "DRUG-B"]})
    issues = []
    result, _ = audit.inspect_data({"registry": registry(), "matrix": matrix(),
                                    "candidate": hypotheses, "candidate_models": models}, issues)
    assert result["candidate"]["models_with_multiple_reference_drugs_n"] == 1
    assert result["candidate"]["unknown_hypothesis_ids_n"] == 0


def test_reference_measurements_remain_uninterpreted_not_fake_pyz():
    response = {"profile_v2": {"compound_id": "PYZ-035", "schema_version": audit.PROFILE_VERSION,
                "identity": {"own_compound_id": "PYZ-035", "standardized_smiles": "CCO"},
                "scientific_summary": {"own_target_binding": "not_measured", "own_cellular_activity": "not_measured"},
                "target_hypotheses": {"items": [{"reference_endpoint": "Ki", "reference_value_nm": 10},
                                                {"reference_endpoint": "IC50", "reference_value_nm": 20}]}}}
    sample = audit.safe_api_sample("PYZ-035", response)
    assert sample["target_items_n"] == 2
    assert sample["own_target_binding"] == "not_measured"
    assert "reference_value_nm" not in sample


def test_percentile_is_deterministic():
    assert audit.percentile([1, 2, 3, 4, 5], .5) == 3
    assert audit.percentile([1, 2, 3, 4, 5], .95) == 4.8
    assert audit.percentile([], .95) is None


def test_only_audit_output_modified(tmp_path):
    root = tmp_path / "mcl"
    (root / "config").mkdir(parents=True)
    (root / "apps" / "api").mkdir(parents=True)
    p = root / "config" / "own_compounds.tsv"
    registry().to_csv(p, sep="\t", index=False)
    before = audit.sha256(p)
    from argparse import Namespace
    report_root = tmp_path / "reports_code"
    out, code = audit.audit(Namespace(root=str(root), report_root=str(report_root),
                                           offline=True, api_base="http://127.0.0.1:8100",
                                           ui_base="http://127.0.0.1:3100", skip_performance=True, timeout=1))
    assert code == 2  # Matrix deliberately absent: never report success.
    assert out.is_dir()
    assert out.is_relative_to(report_root)
    assert not (root / "outputs").exists()
    assert (out / "input_inventory.tsv").exists()
    assert (out / "E0_ACCEPTANCE.md").exists()
    assert "НЕ ПРИНЯТО" in (out / "E0_ACCEPTANCE.md").read_text(encoding="utf-8")
    assert audit.sha256(p) == before