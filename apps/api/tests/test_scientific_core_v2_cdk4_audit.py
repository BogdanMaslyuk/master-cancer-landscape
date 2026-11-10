from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd


SCRIPT = Path(__file__).resolve().parents[3] / "scripts" / "audit_scientific_core_v2_cdk4.py"
spec = importlib.util.spec_from_file_location("mcl_cdk4_audit", SCRIPT)
assert spec is not None and spec.loader is not None
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


def _setup():
    lab = pd.DataFrame([
        {"lab_name": "HCT116", "model_id": "M1"},
        {"lab_name": "MCF-7", "model_id": "M2"},
        {"lab_name": "T-47D", "model_id": None},
    ])
    atlas = pd.DataFrame([
        {"model_id": "M1", "cell_line_name": "HCT 116", "mcl_cancer_id": "C1", "mcl_cancer_name": "Colorectal"},
        {"model_id": "M2", "cell_line_name": "MCF-7", "mcl_cancer_id": "C2", "mcl_cancer_name": "Breast"},
    ])
    return lab, atlas


def test_audit_preserves_unmatched_and_unverified_lab_lines():
    lab, atlas = _setup()
    lines, issues = audit.identify_lines(lab, atlas)
    assert len(lines) == 6
    hct = lines.set_index("lab_name").loc["HCT 116"]
    assert hct["model_id"] == "M1"
    assert not hct["physical_stock_verified"]
    assert lines.set_index("lab_name").loc["T-47D", "laboratory_mapping"] == "not_mapped"
    assert lines.set_index("lab_name").loc["HuH-7", "laboratory_mapping"] == "not_in_lab_registry"
    assert any(x["kind"] == "laboratory_mapping" for x in issues)


def test_audit_keeps_prism_a_and_b_distinct_and_does_not_invent_pyz_response():
    lab, atlas = _setup()
    lines, _ = audit.identify_lines(lab, atlas)
    context = pd.DataFrame([
        {"model_id": "M1", "target_gene": "CDK4", "dependency_probability": 0.9,
         "gene_effect": -0.8, "depmap_release": "26Q1"},
        {"model_id": "M2", "target_gene": "CDK4", "dependency_probability": 0.95,
         "gene_effect": -0.65, "depmap_release": "26Q1"},
    ])
    responses = pd.DataFrame([
        {"observation_id": "O-A", "model_id": "M1", "compound_id": "A",
         "source": "PRISM", "endpoint": "LFC", "value": -2.0, "source_assay_id": "ASSAY-A"},
        {"observation_id": "O-B", "model_id": "M1", "compound_id": "B",
         "source": "PRISM", "endpoint": "LFC", "value": -0.1, "source_assay_id": "ASSAY-B"},
        {"observation_id": "O-PYZ", "model_id": "M1", "compound_id": "PYZ-035",
         "source": "OTHER", "endpoint": "LFC", "value": -3.0, "source_assay_id": "OTHER"},
    ])
    links = pd.DataFrame([
        {"model_id": "M1", "compound_id": "A", "target_gene": "CDK4"},
        {"model_id": "M1", "compound_id": "B", "target_gene": "CDK4"},
        {"model_id": "M1", "compound_id": "A", "target_gene": "CDK4"},
    ])
    compounds = pd.DataFrame([
        {"compound_id": "A", "preferred_name": "Known A"},
        {"compound_id": "B", "preferred_name": "Known B"},
    ])
    target_evidence = pd.DataFrame([
        {"target_gene": "CDK4", "compound_id": "A", "evidence_id": "T-A"},
        {"target_gene": "CDK4", "compound_id": "B", "evidence_id": "T-B"},
    ])
    hypotheses = pd.DataFrame([
        {"hypothesis_id": "H-A", "compound_id": "A", "target_gene": "CDK4",
         "mcl_cancer_id": "C1", "priority_status": "priority_for_in_vitro"},
        {"hypothesis_id": "H-B", "compound_id": "B", "target_gene": "CDK4",
         "mcl_cancer_id": "C1", "priority_status": "priority_for_in_vitro"},
    ])
    candidate_models = pd.DataFrame([
        {"hypothesis_id": "H-A", "model_id": "M1", "model_role": "positive"},
        {"hypothesis_id": "H-B", "model_id": "M1", "model_role": "discordant_dependency_without_activity"},
    ])
    out = audit.analyze_rows(
        lines, context, responses, links, compounds, target_evidence,
        hypotheses, candidate_models
    )
    prism = out["prism_measurements"]
    assert set(prism["compound_id"]) == {"A", "B"}
    assert len(prism) == 2
    assert set(prism["observation_id"]) == {"O-A", "O-B"}
    assert "PYZ-035" not in set(prism["compound_id"])
    hct = out["models"].set_index("lab_name").loc["HCT 116"]
    assert hct["genetic_status"] == "dependent"
    assert hct["dependency_probability"] == 0.9
    assert any(out["issues"]["kind"].eq("multi_role_across_reference_compounds"))


def test_audit_never_turns_missing_crispr_into_negative():
    lab, atlas = _setup()
    lines, _ = audit.identify_lines(lab, atlas)
    context = pd.DataFrame(columns=["model_id", "target_gene", "dependency_probability", "gene_effect"])
    responses = pd.DataFrame(columns=["model_id", "compound_id", "source", "endpoint"])
    links = pd.DataFrame(columns=["model_id", "compound_id", "target_gene"])
    compounds = pd.DataFrame(columns=["compound_id", "preferred_name"])
    evidence = pd.DataFrame(columns=["compound_id", "target_gene"])
    hypotheses = pd.DataFrame(columns=["hypothesis_id", "compound_id", "target_gene",
                                        "mcl_cancer_id", "priority_status"])
    candidates = pd.DataFrame(columns=["hypothesis_id", "model_id", "model_role"])
    out = audit.analyze_rows(lines, context, responses, links, compounds, evidence, hypotheses, candidates)
    assert set(out["models"]["genetic_status"]) == {"unknown"}
    assert out["prism_measurements"].empty
    assert not out["issues"].empty


def test_trace_role_filter_is_compound_specific_and_exposes_excluded_positive_rows():
    candidate = pd.DataFrame([
        {"model_id": "M1", "compound_id": "A", "model_role": "positive",
         "priority_status": "supported_hypothesis"},
        {"model_id": "M1", "compound_id": "B",
         "model_role": "discordant_dependency_without_activity",
         "priority_status": "priority_for_in_vitro"},
        {"model_id": "M2", "compound_id": "C", "model_role": "positive",
         "priority_status": "priority_for_in_vitro"},
    ])
    prism = pd.DataFrame([
        {"model_id": "M1", "compound_id": "A", "preferred_name": "Known A",
         "observation_id": "OA1"},
        {"model_id": "M1", "compound_id": "A", "preferred_name": "Known A",
         "observation_id": "OA2"},
        {"model_id": "M1", "compound_id": "B", "preferred_name": "Known B",
         "observation_id": "OB1"},
        {"model_id": "M2", "compound_id": "C", "preferred_name": "Known C",
         "observation_id": "OC1"},
    ])
    trace, impact = audit.trace_role_filter(candidate, prism)
    assert len(trace) == 3
    assert set(trace["compound_id"]) == {"A", "B", "C"}
    assert trace.set_index("compound_id").loc["A", "normalized_prism_records_n"] == 2
    assert not trace.set_index("compound_id").loc["A", "included_by_legacy_priority_filter"]
    assert trace.set_index("compound_id").loc["B", "included_by_legacy_priority_filter"]
    m1 = impact.set_index("model_id").loc["M1"]
    assert m1["candidate_positive_rows"] == 1
    assert m1["legacy_included_positive_rows"] == 0
    assert m1["legacy_excluded_positive_rows"] == 1


def test_trace_role_filter_does_not_equate_missing_prism_with_negative():
    candidate = pd.DataFrame([{
        "model_id": "M1", "compound_id": "A", "model_role": "positive",
        "priority_status": "priority_for_in_vitro",
    }])
    prism = pd.DataFrame(columns=["model_id", "compound_id", "observation_id"])
    trace, impact = audit.trace_role_filter(candidate, prism)
    assert int(trace.iloc[0]["normalized_prism_records_n"]) == 0
    assert int(impact.iloc[0]["legacy_included_positive_rows"]) == 1
    assert "negative" not in set(trace["model_role"])
