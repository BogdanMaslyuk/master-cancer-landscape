from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from mcl_api.matrix_gene_explorer import MatrixGeneExplorerStore
from mcl_api.store import MCLDataStore


def _root(tmp_path: Path) -> Path:
    (tmp_path / "config").mkdir(parents=True)
    (tmp_path / "data/processed/depmap_genomewide").mkdir(parents=True)
    (tmp_path / "data/processed/pathways_sensitivity").mkdir(parents=True)
    (tmp_path / "outputs/reports").mkdir(parents=True)
    (tmp_path / "outputs/qc").mkdir(parents=True)

    (tmp_path / "config/pathways.yaml").write_text(
        """organism: hsapiens
comparisons:
  - label: TEST_AHR_vs_control
    cancer_id: CANCER-X
    comparison: test
""",
        encoding="utf-8",
    )
    (tmp_path / "config/cancer_contexts.yaml").write_text(
        """contexts:
  CANCER-X:
    name: Test cancer context
    depmap:
      context_definition: Test-positive models
""",
        encoding="utf-8",
    )
    (tmp_path / "outputs/reports/depmap_genomewide_CANCER-X__test_meta.json").write_text(
        json.dumps({"context_models_n": 3, "comparator_models_n": 3, "genes_analyzed_n": 1, "depmap_release": "26Q1"}),
        encoding="utf-8",
    )
    (tmp_path / "outputs/qc/depmap_genomewide_CANCER-X__test_qc.json").write_text("[]", encoding="utf-8")
    pd.DataFrame([
        {
            "gene_symbol": "AHR", "context_median_gene_effect": -0.8,
            "comparator_median_gene_effect": -0.2, "delta_gene_effect": -0.6,
            "cliffs_delta": -0.7, "p_value": 0.01, "q_value": 0.02,
            "fdr_0_05": True, "context_models_n": 3, "comparator_models_n": 3,
            "broad_dependency_warning": False, "low_sample_size": False,
        }
    ]).to_csv(tmp_path / "data/processed/depmap_genomewide/CANCER-X__test_genes.tsv", sep="\t", index=False)

    pd.DataFrame([{"gene_symbol":"AHR","present_all_thresholds":False}]).to_csv(
        tmp_path / "data/processed/pathways_sensitivity/gene_stability.tsv", sep="\t", index=False
    )

    ids=[f"ACH-{i}" for i in range(1,7)]
    ge=[-1.0,-0.8,-0.6,-0.3,-0.2,-0.1]
    expression=[8.0,7.0,6.0,4.0,3.0,2.0]
    cn=[1.4,1.3,1.2,1.0,0.9,0.8]
    for layer,values in (("gene_effect",ge),("expression",expression),("copy_number",cn)):
        pd.DataFrame({"model_id":ids,"AHR":values}).to_parquet(tmp_path / f"data/processed/depmap_model_{layer}.parquet",index=False)
        (tmp_path / f"data/processed/depmap_model_{layer}_genes.json").write_text(json.dumps([{"gene":"AHR"}]),encoding="utf-8")

    pd.DataFrame({
        "model_id":ids,
        "cell_line_name":[f"MODEL{i}" for i in range(1,7)],
        "oncotree_subtype":["Test"]*6,
    }).to_parquet(tmp_path / "data/processed/depmap_model_metadata.parquet",index=False)
    pd.DataFrame({
        "model_id":ids,
        "cancer_id":["CANCER-X"]*6,
        "assigned_group":["context","context","context","comparator","comparator","comparator"],
    }).to_parquet(tmp_path / "data/processed/depmap_context_audit.parquet",index=False)

    mutation_rows=[]
    for model_id in ids[:3]:
        mutation_rows.append({"model_id":model_id,"gene":"KRAS","is_functional":True,"driver":True,"hotspot":True})
    for model_id in ids[3:]:
        mutation_rows.append({"model_id":model_id,"gene":"TP53","is_functional":True,"driver":True,"hotspot":False})
    pd.DataFrame(mutation_rows).to_parquet(tmp_path / "data/processed/depmap_model_mutations.parquet",index=False)
    return tmp_path


def test_descriptive_contexts_are_separate_from_inferential_comparison(tmp_path: Path) -> None:
    root=_root(tmp_path)
    explorer=MatrixGeneExplorerStore(root,MCLDataStore(root))
    payload=explorer.descriptive_contexts("AHR")
    assert payload["available"] is True
    assert payload["items"][0]["context"]["n"] == 3
    assert payload["items"][0]["context"]["median"] == -0.8
    assert payload["items"][0]["comparator"]["median"] == -0.2
    assert payload["items"][0]["interpretation_level"] == "descriptive_only"


def test_multiomics_correlations_return_points_and_spearman(tmp_path: Path) -> None:
    root=_root(tmp_path)
    explorer=MatrixGeneExplorerStore(root,MCLDataStore(root))
    payload=explorer.correlations("AHR")
    assert payload["available"] is True
    assert len(payload["points"]) == 6
    relationships={row["layer"]:row for row in payload["relationships"]}
    assert relationships["expression"]["n"] == 6
    assert relationships["expression"]["rho"] is not None
    assert relationships["expression"]["rho"] < 0


def test_mutation_association_finds_kras_link_to_gene_effect(tmp_path: Path) -> None:
    root=_root(tmp_path)
    explorer=MatrixGeneExplorerStore(root,MCLDataStore(root))
    payload=explorer.mutation_associations("AHR",min_mutated=3,min_wildtype=3,limit=10)
    assert payload["available"] is True
    rows={row["mutation_gene"]:row for row in payload["items"]}
    assert "KRAS" in rows
    assert rows["KRAS"]["mutated_n"] == 3
    assert rows["KRAS"]["wildtype_n"] == 3
    assert rows["KRAS"]["delta_median_gene_effect"] < 0


def test_identity_exposes_model_insights_bundle(tmp_path: Path) -> None:
    root=_root(tmp_path)
    explorer=MatrixGeneExplorerStore(root,MCLDataStore(root))
    payload=explorer.identity("AHR")
    assert "model_insights" in payload
    assert payload["model_insights"]["descriptive_contexts"]["available"] is True
