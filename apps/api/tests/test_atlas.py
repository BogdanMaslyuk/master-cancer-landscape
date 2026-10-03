from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from mcl_api.atlas import MCLAtlas
from mcl_api.store import MCLDataStore


def _write_tsv(path: Path, rows: list[dict]):
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def test_atlas_builds_disease_model_hierarchy(tmp_path: Path):
    (tmp_path / "config").mkdir(parents=True)
    (tmp_path / "config/cancer_contexts.yaml").write_text(
        """contexts:
  C1:
    name: Lung adenocarcinoma — KRAS G12C
    inclusion:
      lineage: lung
      disease: lung adenocarcinoma
      alteration: p.G12C
    opentargets:
      molecular_context: KRAS p.G12C
    depmap:
      context_definition: LUAD with KRAS G12C
      comparator_definition: LUAD without KRAS G12C
""",
        encoding="utf-8",
    )
    (tmp_path / "config/pathways.yaml").write_text(
        "comparisons:\n  - label: G12C vs WT\n    cancer_id: C1\n    comparison: kras-wt\n",
        encoding="utf-8",
    )
    (tmp_path / "outputs/reports").mkdir(parents=True)
    (tmp_path / "outputs/reports/depmap_genomewide_C1__kras-wt_meta.json").write_text(
        json.dumps({"context_models_n": 1, "comparator_models_n": 1, "genes_analyzed_n": 100, "depmap_release": "test"}),
        encoding="utf-8",
    )

    variant = json.dumps([{"GeneSymbol": "KRAS", "ProteinChange": "p.G12C", "Hotspot": "True"}])
    _write_tsv(
        tmp_path / "data/processed/depmap_context_audit.tsv",
        [
            {
                "cancer_id": "C1", "model_id": "ACH-1", "cell_line_name": "LINE-A",
                "depmap_model_type": "LUAD", "oncotree_lineage": "Lung",
                "oncotree_primary_disease": "Non-Small Cell Lung Cancer",
                "oncotree_subtype": "Lung Adenocarcinoma", "oncotree_code": "LUAD",
                "sequencing_available": True, "alteration_status": "defining_present",
                "assigned_group": "context", "assignment_reason": "present",
                "qualifying_variants_json": variant, "other_relevant_variants_json": "[]",
            },
            {
                "cancer_id": "C1", "model_id": "ACH-2", "cell_line_name": "LINE-B",
                "depmap_model_type": "LUAD", "oncotree_lineage": "Lung",
                "oncotree_primary_disease": "Non-Small Cell Lung Cancer",
                "oncotree_subtype": "Lung Adenocarcinoma", "oncotree_code": "LUAD",
                "sequencing_available": True, "alteration_status": "defining_absent",
                "assigned_group": "comparator", "assignment_reason": "absent",
                "qualifying_variants_json": "[]", "other_relevant_variants_json": "[]",
            },
        ],
    )

    store = MCLDataStore(tmp_path)
    atlas = MCLAtlas(tmp_path, store)

    overview = atlas.atlas()
    assert overview["organs_n"] == 1
    assert overview["contexts_n"] == 1
    assert overview["models_n"] == 2

    context = atlas.context("C1")
    assert context["analysis_available"] is True
    assert context["context_models_n"] == 1
    assert context["comparator_models_n"] == 1

    models = atlas.models(cancer_id="C1")
    assert models["total"] == 2
    assert models["items"][0]["variants"][0]["gene"] == "KRAS"
    assert models["items"][0]["variants"][0]["protein_change"] == "p.G12C"
