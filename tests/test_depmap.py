from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import yaml

from mcl.analysis.depmap import (
    analyze_depmap_wave1,
    analyze_depmap_kras_sensitivity,
    build_context_membership,
    build_kras_sensitivity_membership,
    cliffs_delta,
)
from mcl.qc.depmap import depmap_qc, depmap_kras_sensitivity_qc
from mcl.sources.depmap import load_gene_matrix, load_sequenced_model_ids, parse_gene_label


def test_parse_depmap_gene_label():
    assert parse_gene_label("KRAS (3845)") == ("KRAS", "3845")
    assert parse_gene_label("PTPN11") == ("PTPN11", None)


def test_load_gene_matrix_maps_symbol_and_rejects_missing(tmp_path):
    path = tmp_path / "CRISPRGeneEffect.csv"
    path.write_text("ModelID,KRAS (3845),PTPN11 (5781)\nACH-1,-1.0,-0.2\n", encoding="utf-8")
    df = load_gene_matrix(path, ["KRAS", "PTPN11"])
    assert list(df.columns) == ["KRAS", "PTPN11"]
    assert df.loc["ACH-1", "KRAS"] == -1.0


def test_omics_profiles_only_counts_wes_wgs_default(tmp_path):
    path = tmp_path / "OmicsProfiles.csv"
    pd.DataFrame(
        [
            {"ModelID": "ACH-1", "Datatype": "wes", "is_default_entry": "True"},
            {"ModelID": "ACH-2", "Datatype": "rna", "is_default_entry": "True"},
            {"ModelID": "ACH-3", "Datatype": "wgs", "is_default_entry": "False"},
        ]
    ).to_csv(path, index=False)
    ids, schema = load_sequenced_model_ids(path)
    assert ids == {"ACH-1"}
    assert schema["default_entry_column"] == "is_default_entry"


def test_omics_profiles_recognizes_official_is_default_entry_for_model(tmp_path):
    path = tmp_path / "OmicsProfiles.csv"
    pd.DataFrame(
        [
            {"ModelID": "ACH-1", "Datatype": "wes", "IsDefaultEntryForModel": "Yes"},
            {"ModelID": "ACH-2", "Datatype": "wgs", "IsDefaultEntryForModel": "No"},
        ]
    ).to_csv(path, index=False)
    ids, schema = load_sequenced_model_ids(path)
    assert ids == {"ACH-1"}
    assert schema["default_entry_column"] == "IsDefaultEntryForModel"


def test_kras_context_assignment_requires_sequencing_and_exact_variant():
    models = pd.DataFrame(
        [
            {"ModelID": "ACH-1", "CellLineName": "A", "OncotreeCode": "PAAD"},
            {"ModelID": "ACH-2", "CellLineName": "B", "OncotreeCode": "PAAD"},
            {"ModelID": "ACH-3", "CellLineName": "C", "OncotreeCode": "PAAD"},
        ]
    )
    muts = pd.DataFrame(
        [
            {"ModelID": "ACH-1", "GeneSymbol": "KRAS", "ProteinChange": "p.G12D"},
            {"ModelID": "ACH-2", "GeneSymbol": "KRAS", "ProteinChange": "p.G12C"},
        ]
    )
    contexts = {
        "CANCER-004": {
            "depmap": {
                "metadata": {"oncotree_codes": ["PAAD"]},
                "mutation": {
                    "mode": "exact_protein_change_present",
                    "gene": "KRAS",
                    "protein_change_regexes": [r"^p\.G12D$"],
                    "other_variants_ambiguous": False,
                },
                "context_assignment": "defining_present",
                "comparator_assignment": "defining_absent",
            }
        }
    }
    memberships, audit = build_context_membership(
        models, muts, {"ACH-1", "ACH-2"}, contexts, {"driver_flags": [], "text_columns": []}
    )
    assert memberships["CANCER-004"]["context"] == ["ACH-1"]
    assert memberships["CANCER-004"]["comparator"] == ["ACH-2"]
    excluded = [x for x in audit if x.model_id == "ACH-3"][0]
    assert excluded.assigned_group == "excluded"
    assert excluded.alteration_status == "unsequenced"


def test_idh_wildtype_excludes_noncanonical_idh_variant_as_ambiguous():
    models = pd.DataFrame(
        [
            {"ModelID": "ACH-WT", "OncotreeCode": "GBM"},
            {"ModelID": "ACH-R132", "OncotreeCode": "GBM"},
            {"ModelID": "ACH-OTHER", "OncotreeCode": "GBM"},
        ]
    )
    muts = pd.DataFrame(
        [
            {"ModelID": "ACH-R132", "GeneSymbol": "IDH1", "ProteinChange": "p.R132H"},
            {"ModelID": "ACH-OTHER", "GeneSymbol": "IDH1", "ProteinChange": "p.A100T"},
        ]
    )
    contexts = {
        "CANCER-011": {
            "depmap": {
                "metadata": {"oncotree_codes": ["GBM"]},
                "mutation": {
                    "mode": "absence_of_defining_variants",
                    "defining_variants": [
                        {"gene": "IDH1", "protein_change_regexes": [r"^p\.R132"]},
                        {"gene": "IDH2", "protein_change_regexes": [r"^p\.R172"]},
                    ],
                    "other_relevant_variants_ambiguous": True,
                },
                "context_assignment": "defining_absent",
                "comparator_assignment": "defining_present",
            }
        }
    }
    memberships, _ = build_context_membership(
        models,
        muts,
        {"ACH-WT", "ACH-R132", "ACH-OTHER"},
        contexts,
        {"driver_flags": [], "text_columns": []},
    )
    assert memberships["CANCER-011"]["context"] == ["ACH-WT"]
    assert memberships["CANCER-011"]["comparator"] == ["ACH-R132"]
    assert memberships["CANCER-011"]["excluded"] == ["ACH-OTHER"]


def test_cliffs_delta_sign_matches_more_negative_context():
    # Context values are all lower (more dependent) than comparator values.
    assert cliffs_delta([-1.2, -1.0, -0.8], [-0.2, -0.1, 0.0]) == -1.0


def _write_minimal_depmap_project(root: Path) -> None:
    (root / "data/input").mkdir(parents=True)
    (root / "data/raw/depmap/TEST").mkdir(parents=True)
    (root / "config").mkdir(parents=True)
    (root / "data/input/cancer_target_pairs_wave1.tsv").write_text(
        "Evidence_ID\tCancer_ID\tTarget_ID\tCancer_name\tMolecular_context\tHGNC_symbol\n"
        "EVID-1\tCANCER-004\tTGT-KRAS\tPDAC KRAS G12D\tKRAS p.G12D\tKRAS\n",
        encoding="utf-8",
    )
    cfg = {
        "contexts": {
            "CANCER-004": {
                "name": "PDAC KRAS G12D",
                "depmap": {
                    "context_definition": "PDAC + KRAS G12D",
                    "comparator_definition": "PDAC without KRAS G12D",
                    "metadata": {"oncotree_codes": ["PAAD"]},
                    "mutation": {
                        "mode": "exact_protein_change_present",
                        "gene": "KRAS",
                        "protein_change_regexes": [r"^p\.G12D$"],
                        "other_variants_ambiguous": False,
                    },
                    "context_assignment": "defining_present",
                    "comparator_assignment": "defining_absent",
                },
            }
        }
    }
    (root / "config/cancer_contexts.yaml").write_text(
        yaml.safe_dump(cfg, sort_keys=False), encoding="utf-8"
    )
    (root / "config/thresholds.yaml").write_text(
        yaml.safe_dump(
            {
                "depmap": {
                    "minimum_context_n_warning": 2,
                    "minimum_comparator_n_warning": 2,
                    "minimum_n_for_statistical_test": 2,
                    "dependency_probability_threshold": 0.5,
                    "broad_dependency_fraction_warning": 0.8,
                }
            }
        ),
        encoding="utf-8",
    )
    raw = root / "data/raw/depmap/TEST"
    pd.DataFrame(
        [
            {"ModelID": "ACH-1", "CellLineName": "A", "OncotreeCode": "PAAD"},
            {"ModelID": "ACH-2", "CellLineName": "B", "OncotreeCode": "PAAD"},
            {"ModelID": "ACH-3", "CellLineName": "C", "OncotreeCode": "PAAD"},
            {"ModelID": "ACH-4", "CellLineName": "D", "OncotreeCode": "PAAD"},
        ]
    ).to_csv(raw / "Model.csv", index=False)
    pd.DataFrame(
        [
            {"ModelID": mid, "Datatype": "wes", "is_default_entry": True}
            for mid in ["ACH-1", "ACH-2", "ACH-3", "ACH-4"]
        ]
    ).to_csv(raw / "OmicsProfiles.csv", index=False)
    pd.DataFrame(
        [
            {"ModelID": "ACH-1", "HugoSymbol": "KRAS", "ProteinChange": "p.G12D", "Hotspot": "Y"},
            {"ModelID": "ACH-2", "HugoSymbol": "KRAS", "ProteinChange": "p.G12D", "Hotspot": "Y"},
            {"ModelID": "ACH-3", "HugoSymbol": "KRAS", "ProteinChange": "p.G12C", "Hotspot": "Y"},
        ]
    ).to_csv(raw / "OmicsSomaticMutations.csv", index=False)
    pd.DataFrame(
        {
            "ModelID": ["ACH-1", "ACH-2", "ACH-3", "ACH-4"],
            "KRAS (3845)": [-1.2, -1.0, -0.2, -0.1],
        }
    ).to_csv(raw / "CRISPRGeneEffect.csv", index=False)
    pd.DataFrame(
        {
            "ModelID": ["ACH-1", "ACH-2", "ACH-3", "ACH-4"],
            "KRAS (3845)": [0.95, 0.9, 0.2, 0.1],
        }
    ).to_csv(raw / "CRISPRGeneDependency.csv", index=False)


def test_depmap_end_to_end_statistics_and_audit(tmp_path):
    _write_minimal_depmap_project(tmp_path)
    rows, audit, provenance, inventory, meta = analyze_depmap_wave1(tmp_path, "TEST")
    assert len(rows) == 1
    row = rows[0]
    assert row.context_models_n == 2
    assert row.comparator_models_n == 2
    assert row.dependency_fraction == 1.0
    assert row.comparator_dependency_fraction == 0.0
    assert row.delta_gene_effect < 0
    assert row.cliffs_delta == -1.0
    assert row.p_value is not None
    assert row.q_value == row.p_value
    assert row.low_sample_size is False
    assert len(audit) == 4
    assert len(inventory) == 5
    assert provenance
    qc = depmap_qc(rows, audit, expected_pairs_n=1)
    assert not [x for x in qc if x.severity == "ERROR"]


def test_depmap_26q1_glioblastoma_idh_wildtype_metadata_is_recognized():
    models = pd.DataFrame(
        [
            {
                "ModelID": "ACH-GB-WT",
                "CellLineName": "GB1",
                "DepmapModelType": "GB",
                "OncotreeCode": "GB",
                "OncotreeLineage": "CNS/Brain",
                "OncotreePrimaryDisease": "Adult-Type Diffuse Glioma",
                "OncotreeSubtype": "Glioblastoma, IDH-Wildtype",
            },
            {
                "ModelID": "ACH-ASTR-IDH",
                "CellLineName": "ASTR1",
                "DepmapModelType": "ASTR",
                "OncotreeCode": "ASTR",
                "OncotreeLineage": "CNS/Brain",
                "OncotreePrimaryDisease": "Adult-Type Diffuse Glioma",
                "OncotreeSubtype": "Astrocytoma, IDH-Mutant",
            },
        ]
    )
    muts = pd.DataFrame(columns=["ModelID", "GeneSymbol", "ProteinChange"])
    contexts = {
        "CANCER-011": {
            "depmap": {
                "metadata": {
                    "oncotree_codes": ["GB"],
                    "depmap_model_types": ["GB"],
                    "oncotree_subtypes": ["Glioblastoma, IDH-Wildtype"],
                    "oncotree_lineages": ["CNS/Brain"],
                },
                "mutation": {
                    "mode": "absence_of_defining_variants",
                    "defining_variants": [
                        {"gene": "IDH1", "protein_change_regexes": [r"^p\.R132"]},
                        {"gene": "IDH2", "protein_change_regexes": [r"^p\.R172"]},
                    ],
                    "other_relevant_variants_ambiguous": True,
                },
                "context_assignment": "defining_absent",
                "comparator_assignment": "defining_present",
            }
        }
    }
    memberships, audit = build_context_membership(
        models,
        muts,
        {"ACH-GB-WT", "ACH-ASTR-IDH"},
        contexts,
        {"driver_flags": [], "text_columns": []},
    )
    assert memberships["CANCER-011"]["context"] == ["ACH-GB-WT"]
    assert memberships["CANCER-011"]["comparator"] == []
    # Astrocytoma belongs to the same parent disease but must not be silently treated
    # as the same exact glioblastoma context.
    assert not [x for x in audit if x.model_id == "ACH-ASTR-IDH"]


def test_analyze_depmap_cli_uses_current_qc_keyword():
    """Regression test for v0.3.2 CLI/QC keyword mismatch."""
    import ast
    import inspect
    import textwrap

    from mcl.cli import analyze_depmap_cmd

    tree = ast.parse(textwrap.dedent(inspect.getsource(analyze_depmap_cmd)))
    calls = [
        node
        for node in ast.walk(tree)
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Name)
        and node.func.id == "depmap_qc"
    ]
    assert len(calls) == 1
    keywords = {kw.arg for kw in calls[0].keywords}
    assert "expected_pairs_n" in keywords
    assert "expected_n" not in keywords


def test_missing_comparator_is_not_mislabeled_as_low_sample_when_context_is_large():
    """No comparator != small context; QC must keep these limitations distinct."""
    from mcl.models.depmap import DepMapEvidence
    from mcl.qc.depmap import depmap_qc, depmap_kras_sensitivity_qc

    row = DepMapEvidence(
        evidence_id="EVID-GBM",
        cancer_id="CANCER-011",
        target_id="TGT-EGFR",
        target_symbol="EGFR",
        depmap_release="26Q1",
        context_definition="Glioblastoma, IDH-Wildtype",
        comparator_definition="",
        context_models_n=66,
        comparator_models_n=0,
        context_gene_effect_n=66,
        comparator_gene_effect_n=0,
        context_dependency_probability_n=66,
        comparator_dependency_probability_n=0,
        dependent_cell_lines_n=10,
        dependency_fraction=10/66,
        median_gene_effect=-0.2,
        mean_gene_effect=-0.2,
        gene_effect_iqr=0.1,
        median_dependency_probability=0.2,
        min_gene_effect=-1.0,
        max_gene_effect=0.1,
        missing_context_gene_effect_n=0,
        comparator_dependent_cell_lines_n=0,
        comparator_dependency_fraction=None,
        comparator_median_gene_effect=None,
        comparator_median_dependency_probability=None,
        delta_gene_effect=None,
        cliffs_delta=None,
        statistical_test=None,
        p_value=None,
        q_value=None,
        low_sample_size=False,
        statistical_test_performed=False,
        all_models_gene_effect_n=1000,
        all_models_dependency_probability_n=1000,
        all_models_median_gene_effect=-0.1,
        broad_dependency_fraction=0.1,
        broad_dependency_warning=False,
        gene_effect_file="g.csv",
        dependency_probability_file="d.csv",
        model_file="m.csv",
        mutation_file="mut.csv",
        omics_profiles_file="o.csv",
        retrieved_at="2026-10-03T00:00:00+00:00",
        notes="",
    )
    qc = depmap_qc([row], [], expected_pairs_n=1)
    checks = [x.check for x in qc]
    assert "comparator_nonempty" in checks
    assert "low_sample_size" not in checks


def test_kras_sensitivity_membership_splits_target_wt_other_and_ambiguous():
    models = pd.DataFrame([
        {"ModelID": "ACH-T", "CellLineName": "target", "OncotreeCode": "PAAD"},
        {"ModelID": "ACH-W", "CellLineName": "wt", "OncotreeCode": "PAAD"},
        {"ModelID": "ACH-O", "CellLineName": "other", "OncotreeCode": "PAAD"},
        {"ModelID": "ACH-A", "CellLineName": "amb", "OncotreeCode": "PAAD"},
        {"ModelID": "ACH-U", "CellLineName": "unseq", "OncotreeCode": "PAAD"},
    ])
    muts = pd.DataFrame([
        {"ModelID": "ACH-T", "GeneSymbol": "KRAS", "ProteinChange": "p.G12D", "Hotspot": "Y"},
        {"ModelID": "ACH-O", "GeneSymbol": "KRAS", "ProteinChange": "p.G12V", "Hotspot": "Y"},
        {"ModelID": "ACH-A", "GeneSymbol": "KRAS", "ProteinChange": "p.M188V", "Hotspot": "False"},
    ])
    contexts = {
        "CANCER-004": {
            "depmap": {
                "metadata": {"oncotree_codes": ["PAAD"]},
                "mutation": {
                    "mode": "exact_protein_change_present",
                    "gene": "KRAS",
                    "protein_change_regexes": [r"^p\.G12D$"],
                },
            }
        }
    }
    memberships, audit = build_kras_sensitivity_membership(
        models,
        muts,
        {"ACH-T", "ACH-W", "ACH-O", "ACH-A"},
        contexts,
        {"driver_flags": ["Hotspot"]},
    )
    groups = memberships["CANCER-004"]
    assert groups["target_variant"] == ["ACH-T"]
    assert groups["kras_wildtype_proxy"] == ["ACH-W"]
    assert groups["other_kras_driver_hotspot"] == ["ACH-O"]
    assert groups["excluded"] == ["ACH-A", "ACH-U"]
    statuses = {x.model_id: x.kras_status for x in audit}
    assert statuses["ACH-A"] == "ambiguous_kras_variant"
    assert statuses["ACH-U"] == "unsequenced"


def _write_kras_sensitivity_project(root: Path) -> None:
    (root / "data/input").mkdir(parents=True)
    (root / "data/raw/depmap/TEST").mkdir(parents=True)
    (root / "config").mkdir(parents=True)
    (root / "data/input/cancer_target_pairs_wave1.tsv").write_text(
        "Evidence_ID\tCancer_ID\tTarget_ID\tCancer_name\tMolecular_context\tHGNC_symbol\n"
        "EVID-1\tCANCER-004\tTGT-KRAS\tPDAC KRAS G12D\tKRAS p.G12D\tKRAS\n",
        encoding="utf-8",
    )
    cfg = {
        "contexts": {
            "CANCER-004": {
                "name": "PDAC KRAS G12D",
                "depmap": {
                    "metadata": {"oncotree_codes": ["PAAD"]},
                    "mutation": {
                        "mode": "exact_protein_change_present",
                        "gene": "KRAS",
                        "protein_change_regexes": [r"^p\.G12D$"],
                        "other_variants_ambiguous": False,
                    },
                    "context_assignment": "defining_present",
                    "comparator_assignment": "defining_absent",
                },
            }
        }
    }
    (root / "config/cancer_contexts.yaml").write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf-8")
    (root / "config/thresholds.yaml").write_text(
        yaml.safe_dump({"depmap": {
            "minimum_context_n_warning": 2,
            "minimum_comparator_n_warning": 2,
            "minimum_n_for_statistical_test": 2,
            "dependency_probability_threshold": 0.5,
        }}),
        encoding="utf-8",
    )
    raw = root / "data/raw/depmap/TEST"
    pd.DataFrame([
        {"ModelID": "ACH-T1", "CellLineName": "T1", "OncotreeCode": "PAAD"},
        {"ModelID": "ACH-T2", "CellLineName": "T2", "OncotreeCode": "PAAD"},
        {"ModelID": "ACH-W1", "CellLineName": "W1", "OncotreeCode": "PAAD"},
        {"ModelID": "ACH-W2", "CellLineName": "W2", "OncotreeCode": "PAAD"},
        {"ModelID": "ACH-O1", "CellLineName": "O1", "OncotreeCode": "PAAD"},
        {"ModelID": "ACH-O2", "CellLineName": "O2", "OncotreeCode": "PAAD"},
    ]).to_csv(raw / "Model.csv", index=False)
    pd.DataFrame([
        {"ModelID": x, "Datatype": "wes", "IsDefaultEntryForModel": "Yes"}
        for x in ["ACH-T1", "ACH-T2", "ACH-W1", "ACH-W2", "ACH-O1", "ACH-O2"]
    ]).to_csv(raw / "OmicsProfiles.csv", index=False)
    pd.DataFrame([
        {"ModelID": "ACH-T1", "HugoSymbol": "KRAS", "ProteinChange": "p.G12D", "Hotspot": "Y"},
        {"ModelID": "ACH-T2", "HugoSymbol": "KRAS", "ProteinChange": "p.G12D", "Hotspot": "Y"},
        {"ModelID": "ACH-O1", "HugoSymbol": "KRAS", "ProteinChange": "p.G12V", "Hotspot": "Y"},
        {"ModelID": "ACH-O2", "HugoSymbol": "KRAS", "ProteinChange": "p.G12R", "HessDriver": "Y"},
    ]).to_csv(raw / "OmicsSomaticMutations.csv", index=False)
    pd.DataFrame({
        "ModelID": ["ACH-T1", "ACH-T2", "ACH-W1", "ACH-W2", "ACH-O1", "ACH-O2"],
        "KRAS (3845)": [-2.0, -1.8, -0.2, -0.1, -1.0, -0.9],
    }).to_csv(raw / "CRISPRGeneEffect.csv", index=False)
    pd.DataFrame({
        "ModelID": ["ACH-T1", "ACH-T2", "ACH-W1", "ACH-W2", "ACH-O1", "ACH-O2"],
        "KRAS (3845)": [0.99, 0.98, 0.1, 0.1, 0.9, 0.85],
    }).to_csv(raw / "CRISPRGeneDependency.csv", index=False)


def test_kras_sensitivity_end_to_end(tmp_path):
    _write_kras_sensitivity_project(tmp_path)
    rows, audit, meta = analyze_depmap_kras_sensitivity(tmp_path, "TEST")
    assert len(rows) == 2
    by_type = {x.comparison_type: x for x in rows}
    wt = by_type["vs_kras_wildtype_proxy"]
    other = by_type["vs_other_kras_driver_hotspot"]
    assert wt.context_gene_effect_n == 2 and wt.comparator_gene_effect_n == 2
    assert wt.delta_gene_effect < 0 and wt.cliffs_delta == -1.0
    assert wt.p_value is not None and wt.q_value is not None
    assert other.delta_gene_effect < 0 and other.cliffs_delta == -1.0
    assert meta["analysis"].startswith("M3.1")
    qc = depmap_kras_sensitivity_qc(rows, audit, expected_rows_n=2)
    assert not [x for x in qc if x.severity == "ERROR"]


def test_depmap_catalog_parser_and_latest_release():
    from mcl.sources.depmap_sync import parse_catalog_csv, latest_catalog_release, catalog_releases

    text = (
        "release,file_name,description\n"
        "26Q1,Model.csv,models\n"
        "26Q1,CRISPRGeneEffect.csv,gene effect\n"
        "26Q3,Model.csv,models\n"
    )
    rows = parse_catalog_csv(text)
    assert len(rows) == 3
    assert catalog_releases(rows) == ["26Q1", "26Q3"]
    assert latest_catalog_release(rows) == "26Q3"


def test_sync_depmap_offline_checks_cache_hashes_and_headers(tmp_path):
    from mcl.sources.depmap_sync import sync_depmap

    raw = tmp_path / "data/raw/depmap/TEST"
    raw.mkdir(parents=True)
    pd.DataFrame([{"ModelID": "ACH-1", "CellLineName": "A"}]).to_csv(raw / "Model.csv", index=False)
    pd.DataFrame({"ModelID": ["ACH-1"], "KRAS (3845)": [-1.0]}).to_csv(
        raw / "CRISPRGeneEffect.csv", index=False
    )
    pd.DataFrame({"ModelID": ["ACH-1"], "KRAS (3845)": [0.9]}).to_csv(
        raw / "CRISPRGeneDependency.csv", index=False
    )
    pd.DataFrame(
        [{"ModelID": "ACH-1", "HugoSymbol": "KRAS", "ProteinChange": "p.G12D"}]
    ).to_csv(raw / "OmicsSomaticMutations.csv", index=False)
    pd.DataFrame([{"ModelID": "ACH-1", "Datatype": "wes"}]).to_csv(
        raw / "OmicsProfiles.csv", index=False
    )

    result = sync_depmap(tmp_path, "TEST", check_catalog=False)
    assert result["ready"] is True
    assert len(result["files"]) == 5
    assert all(row["status"] == "ready" for row in result["files"])
    assert all(row["sha256"] for row in result["files"])


def test_sync_depmap_offline_reports_only_missing_files(tmp_path):
    from mcl.sources.depmap_sync import sync_depmap

    raw = tmp_path / "data/raw/depmap/TEST"
    raw.mkdir(parents=True)
    pd.DataFrame([{"ModelID": "ACH-1"}]).to_csv(raw / "Model.csv", index=False)
    result = sync_depmap(tmp_path, "TEST", check_catalog=False)
    statuses = {x["logical_name"]: x["status"] for x in result["files"]}
    assert statuses["models"] == "ready"
    assert statuses["gene_effect"] == "manual_required"
    assert result["ready"] is False


def test_sync_depmap_cli_offline_end_to_end(tmp_path):
    from typer.testing import CliRunner
    from mcl.cli import app

    (tmp_path / "config").mkdir(parents=True)
    (tmp_path / "config/source_versions.yaml").write_text(
        "depmap:\n  release: TEST\n", encoding="utf-8"
    )
    raw = tmp_path / "data/raw/depmap/TEST"
    raw.mkdir(parents=True)
    pd.DataFrame([{"ModelID": "ACH-1"}]).to_csv(raw / "Model.csv", index=False)
    pd.DataFrame({"ModelID": ["ACH-1"], "KRAS (3845)": [-1.0]}).to_csv(raw / "CRISPRGeneEffect.csv", index=False)
    pd.DataFrame({"ModelID": ["ACH-1"], "KRAS (3845)": [0.9]}).to_csv(raw / "CRISPRGeneDependency.csv", index=False)
    pd.DataFrame([{"ModelID": "ACH-1", "HugoSymbol": "KRAS"}]).to_csv(raw / "OmicsSomaticMutations.csv", index=False)
    pd.DataFrame([{"ModelID": "ACH-1", "Datatype": "wes"}]).to_csv(raw / "OmicsProfiles.csv", index=False)

    result = CliRunner().invoke(app, ["sync-depmap", "--root", str(tmp_path), "--offline"])
    assert result.exit_code == 0, result.output
    assert "Ready files: 5/5" in result.output
    assert (tmp_path / "outputs/reports/depmap_sync_manifest.json").exists()
    assert (tmp_path / "data/processed/depmap_sync_inventory.tsv").exists()


def test_genomewide_subset_loader_and_background_fraction(tmp_path):
    from mcl.sources.depmap import (
        load_genomewide_dependency_subset_with_background,
        load_genomewide_matrix_subset,
    )

    path = tmp_path / "matrix.csv"
    pd.DataFrame(
        {
            "ModelID": ["A", "B", "C"],
            "KRAS (3845)": [-1.0, -0.2, -0.1],
            "GENE2 (2)": [-0.5, -0.6, -0.7],
        }
    ).to_csv(path, index=False)
    subset = load_genomewide_matrix_subset(path, ["A", "C"])
    assert subset.shape == (2, 2)
    assert list(subset.columns) == ["KRAS", "GENE2"]

    dep_path = tmp_path / "dep.csv"
    pd.DataFrame(
        {
            "ModelID": ["A", "B", "C"],
            "KRAS (3845)": [0.9, 0.2, 0.1],
            "GENE2 (2)": [0.9, 0.8, 0.7],
        }
    ).to_csv(dep_path, index=False)
    selected, broad, n = load_genomewide_dependency_subset_with_background(
        dep_path, ["A", "C"], dependency_threshold=0.5
    )
    assert selected.shape == (2, 2)
    assert broad["KRAS"] == 1 / 3
    assert broad["GENE2"] == 1.0
    assert n["GENE2"] == 3


def test_genomewide_analysis_and_explorer_end_to_end(tmp_path):
    from mcl.analysis.depmap_genomewide import analyze_depmap_genomewide
    from mcl.qc.depmap_genomewide import depmap_genomewide_qc
    from mcl.visualization.depmap_explorer import build_depmap_explorer_html

    _write_minimal_depmap_project(tmp_path)
    raw = tmp_path / "data/raw/depmap/TEST"
    pd.DataFrame(
        {
            "ModelID": ["ACH-1", "ACH-2", "ACH-3", "ACH-4"],
            "KRAS (3845)": [-1.5, -1.3, -0.2, -0.1],
            "SOS1 (6654)": [-0.1, -0.2, -0.7, -0.6],
            "BROAD1 (999)": [-1.1, -1.0, -1.2, -1.1],
        }
    ).to_csv(raw / "CRISPRGeneEffect.csv", index=False)
    pd.DataFrame(
        {
            "ModelID": ["ACH-1", "ACH-2", "ACH-3", "ACH-4"],
            "KRAS (3845)": [0.99, 0.95, 0.1, 0.1],
            "SOS1 (6654)": [0.1, 0.2, 0.9, 0.8],
            "BROAD1 (999)": [0.99, 0.99, 0.99, 0.99],
        }
    ).to_csv(raw / "CRISPRGeneDependency.csv", index=False)

    results, cohort, ge, dep, meta = analyze_depmap_genomewide(
        tmp_path, "TEST", "CANCER-004", "primary"
    )
    assert len(results) == 3
    kras = results.set_index("gene_symbol").loc["KRAS"]
    assert kras["delta_gene_effect"] < 0
    assert kras["context_dependency_fraction"] == 1.0
    broad = results.set_index("gene_symbol").loc["BROAD1"]
    assert broad["broad_dependency_warning"]
    assert meta["context_models_n"] == 2
    assert meta["comparator_models_n"] == 2
    qc = depmap_genomewide_qc(results, meta)
    assert not [x for x in qc if x.severity == "ERROR"]

    out = tmp_path / "explorer.html"
    build_depmap_explorer_html(results, cohort, ge, meta, out, top_n_heatmap=3)
    text = out.read_text(encoding="utf-8")
    assert "Genome-wide Dependency Explorer" in text
    assert "KRAS" in text
    assert "BROAD1" in text


def test_genomewide_cli_end_to_end(tmp_path, monkeypatch):
    from typer.testing import CliRunner
    from mcl.cli import app

    _write_minimal_depmap_project(tmp_path)
    (tmp_path / "config/source_versions.yaml").write_text(
        yaml.safe_dump({"depmap": {"release": "TEST"}}), encoding="utf-8"
    )
    raw = tmp_path / "data/raw/depmap/TEST"
    pd.DataFrame(
        {
            "ModelID": ["ACH-1", "ACH-2", "ACH-3", "ACH-4"],
            "KRAS (3845)": [-1.5, -1.3, -0.2, -0.1],
            "SOS1 (6654)": [-0.1, -0.2, -0.7, -0.6],
        }
    ).to_csv(raw / "CRISPRGeneEffect.csv", index=False)
    pd.DataFrame(
        {
            "ModelID": ["ACH-1", "ACH-2", "ACH-3", "ACH-4"],
            "KRAS (3845)": [0.99, 0.95, 0.1, 0.1],
            "SOS1 (6654)": [0.1, 0.2, 0.9, 0.8],
        }
    ).to_csv(raw / "CRISPRGeneDependency.csv", index=False)

    # The production project declares pyarrow as a required dependency.  The
    # lightweight CI container used for this unit test may omit it, so stub only
    # the serialization call while still exercising CLI wiring and HTML creation.
    def _fake_to_parquet(self, path, index=False, **kwargs):
        Path(path).write_text("parquet-stub", encoding="utf-8")
    monkeypatch.setattr(pd.DataFrame, "to_parquet", _fake_to_parquet)

    runner = CliRunner()
    result = runner.invoke(
        app,
        [
            "analyze-depmap-genome-wide",
            "--cancer-id", "CANCER-004",
            "--comparison", "primary",
            "--root", str(tmp_path),
            "--release", "TEST",
            "--top-n-heatmap", "5",
        ],
    )
    assert result.exit_code == 0, result.output
    assert "Genes analyzed: 2" in result.output
    assert (tmp_path / "outputs/reports/depmap_explorer_CANCER-004__primary.html").exists()
    assert (tmp_path / "data/processed/depmap_genomewide/CANCER-004__primary_genes.parquet").exists()
