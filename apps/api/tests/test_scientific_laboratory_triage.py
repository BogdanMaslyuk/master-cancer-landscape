from __future__ import annotations

from mcl_api.scientific_laboratory_triage import lab_genetic_triage


def _model(model: str, cancer: str | None, p: float | None,
           *, registered: bool = True, obs: int = 0):
    return {
        "target_gene": "CDK4",
        "model_id": model,
        "registered_in_laboratory": registered,
        "cell_line_name": model,
        "cancer_id": cancer,
        "cancer_name": cancer,
        "dependency_probability": p,
        "gene_effect": -1.0 if p is not None and p > 0.5 else -0.1,
        "depmap_release": "26Q1",
        "reference_drug_observations_n": obs,
        "reference_drug_evidence_status": "normalized_prism_reference_observations",
        "laboratory_lines": [model] if registered else [],
    }


def test_exploratory_triage_pairs_only_same_cancer_measured_genetic_models():
    context = {"available": True, "target_gene": "CDK4", "items": [
        _model("MCF-7", "BREAST", 1.0, obs=17),
        _model("T-47D", "BREAST", 0.98, obs=0),
        _model("HCT116", "COLON", 0.80, obs=5),
        _model("COLON-WT", "COLON", 0.10),
        _model("JURKAT", "LYMPHOID", 0.2),
        _model("MISSING", "LYMPHOID", None),
        _model("UNRELATED", "COLON", 0.0, registered=False),
        _model("BJ5ta", None, 0.03),  # non-tumor is excluded by caller
    ]}
    selected = {x["model_id"] for x in context["items"] if x["model_id"] != "BJ5ta"}
    result = lab_genetic_triage(context, eligible_tumor_model_ids=selected)
    by_id = {x["model_id"]: x for x in result["items"]}
    cohorts = {x["cancer_id"]: x for x in result["cohorts"]}
    assert result["laboratory_registered_model_ids_n"] == 7
    assert result["laboratory_human_tumor_model_ids_n"] == 6
    assert result["lab_roles_not_classified_as_human_tumor_n"] == 1
    assert result["paired_cancer_cohorts_n"] == 1
    assert cohorts["BREAST"]["same_cancer_genetic_comparison_possible"] is False
    assert cohorts["COLON"]["same_cancer_genetic_comparison_possible"] is True
    assert cohorts["LYMPHOID"]["same_cancer_genetic_comparison_possible"] is False
    assert by_id["MCF-7"]["genetic_triage_group"] == "comparator_missing"
    assert by_id["T-47D"]["genetic_role"] == "genetically_dependent_hypothesis_model"
    assert by_id["HCT116"]["genetic_triage_group"] == "genetic_pair_available"
    assert by_id["COLON-WT"]["genetic_role"] == "same_cancer_genetic_comparator_candidate"
    assert by_id["JURKAT"]["genetic_role"] == "genetically_nondependent_context_only"
    assert by_id["MISSING"]["genetic_role"] == "missing_dependency"
    assert all(x["laboratory_ready"] is False for x in result["items"])
    assert all(x["compound_specific_activity"] == "not_measured" for x in result["items"])
    assert "BJ5ta" not in by_id
    assert "UNRELATED" not in by_id


def test_prism_known_drug_observations_do_not_change_genetic_triage_rank():
    a = _model("A", "BREAST", 0.9, obs=0)
    b = _model("B", "BREAST", 0.9, obs=400)
    c = _model("C", "BREAST", 0.1, obs=800)
    input_a = {"available": True, "target_gene": "CDK4", "items": [a, b, c]}
    first = lab_genetic_triage(input_a, eligible_tumor_model_ids={"A", "B", "C"})
    a["reference_drug_observations_n"] = 999
    b["reference_drug_observations_n"] = 0
    c["reference_drug_observations_n"] = 0
    second = lab_genetic_triage(input_a, eligible_tumor_model_ids={"A", "B", "C"})
    assert [x["model_id"] for x in first["items"]] == [x["model_id"] for x in second["items"]]
    assert [x["genetic_role"] for x in first["items"]] == [x["genetic_role"] for x in second["items"]]


def test_missing_gene_or_cancer_does_not_generate_matched_control():
    result = lab_genetic_triage({
        "available": True, "target_gene": "CDK4", "items": [
            _model("DEPENDENT", None, 0.95),
            _model("NEGATIVE", None, 0.05),
        ]}, eligible_tumor_model_ids={"DEPENDENT", "NEGATIVE"})
    assert result["paired_cancer_cohorts_n"] == 0
    assert result["cohorts"] == []
    assert not any(x["same_cancer_comparison_possible"] for x in result["items"])
    assert not any(x["genetic_role"] == "same_cancer_genetic_comparator_candidate" for x in result["items"])


def test_unavailable_gene_context_abstains():
    result = lab_genetic_triage({"available": False, "target_gene": "CDK4"},
                                eligible_tumor_model_ids={"M1"})
    assert result["available"] is False
    assert result["items"] == []
