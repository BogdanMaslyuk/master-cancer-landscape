"""P0.3 duplicate rules never mistake shared DOI for independent corroboration."""
from __future__ import annotations

from mcl.analysis.expanded_ligand_dedup import classify_duplicates, summarize_duplicates


def row(source, record, *, origin="", doi="", assay="A1", value="100",
        smiles="CCO", endpoint="Ki", target="CDK4"):
    return {
        "source": source,
        "source_record_id": record,
        "original_chembl_activity_id": origin,
        "doi": doi, "assay_id": assay,
        "numeric_value_nm": value, "relation": "=",
        "original_smiles": smiles, "target_gene": target, "endpoint": endpoint,
    }


def test_explicit_cross_source_mirror_only_counted_once_provisionally():
    records = [
        row("ChEMBL", "123"),
        row("BindingDB", "567", origin="123"),
    ]
    parsed = classify_duplicates(records)
    assert [r["duplicate_class"] for r in parsed] == ["confirmed_source_mirror"] * 2
    assert [r["count_as_independent"] for r in parsed] == ["primary_candidate", "no"]
    assert parsed[0]["duplicate_cluster_id"] == parsed[1]["duplicate_cluster_id"]
    assert summarize_duplicates(parsed)["confirmed_no_double_count_n"] == 1


def test_cross_source_mirror_conflict_stays_unresolved():
    parsed = classify_duplicates([
        row("ChEMBL", "123", value="100"),
        row("BindingDB", "567", origin="123", value="101"),
    ])
    assert all(r["duplicate_class"] == "conflicting_source_mirror" for r in parsed)
    assert all(r["count_as_independent"] == "pending_curation" for r in parsed)


def test_same_doi_different_assay_values_remain_distinct():
    parsed = classify_duplicates([
        row("ChEMBL", "1", doi="10.1234/paper", value="100", assay="A1"),
        row("ChEMBL", "2", doi="10.1234/paper", value="200", assay="A2"),
    ])
    assert all(r["duplicate_class"] == "unmatched_independent_candidate" for r in parsed)


def test_probable_publication_copy_needs_review_not_counted_as_confirmed():
    parsed = classify_duplicates([
        row("ChEMBL", "1", doi="10.1234/paper"),
        row("BindingDB", "2", doi="10.1234/paper"),
    ])
    assert all(r["duplicate_class"] == "probable_publication_duplicate_requires_review" for r in parsed)
    assert all(r["count_as_independent"] == "pending_curation" for r in parsed)


def test_same_source_record_duplicates_are_identified():
    parsed = classify_duplicates([row("ChEMBL", "1"), row("ChEMBL", "1")])
    assert [r["count_as_independent"] for r in parsed] == ["primary_candidate", "no"]


def test_one_molecule_multiple_targets_or_endpoints_not_collapsed():
    parsed = classify_duplicates([
        row("ChEMBL", "1", target="CDK4"),
        row("ChEMBL", "2", target="EGFR"),
        row("ChEMBL", "3", endpoint="IC50"),
    ])
    assert all(r["duplicate_class"] == "unmatched_independent_candidate" for r in parsed)
