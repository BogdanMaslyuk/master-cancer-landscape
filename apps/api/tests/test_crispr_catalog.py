from __future__ import annotations

from pathlib import Path

import pandas as pd

from mcl_api.crispr_catalog import CRISPRModelCatalog


def _write_tsv(path: Path, rows: list[dict]):
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def test_crispr_catalog_builds_organ_cancer_model_hierarchy(tmp_path: Path):
    _write_tsv(
        tmp_path / "data/processed/depmap_crispr_model_atlas.tsv",
        [
            {
                "model_id": "ACH-1",
                "cell_line_name": "LINE-A",
                "oncotree_lineage": "Lung",
                "oncotree_primary_disease": "Non-Small Cell Lung Cancer",
                "oncotree_subtype": "Lung Adenocarcinoma",
                "oncotree_code": "LUAD",
                "mcl_system_ru": "Дыхательная система",
                "mcl_organ_id": "легкое",
                "mcl_organ_ru": "Лёгкое",
                "mcl_organ_icon": "lung",
                "mcl_cancer_id": "lung-nsclc",
                "mcl_cancer_name": "Non-Small Cell Lung Cancer",
                "mcl_subtype_name": "Lung Adenocarcinoma",
                "mcl_subtype_id": "lung-nsclc-luad",
                "classification_status": "classified",
                "classification_confidence": "high",
                "primary_or_metastasis": "Primary",
                "sample_collection_site": "Lung",
                "has_crispr": True,
                "depmap_release": "test",
            },
            {
                "model_id": "ACH-2",
                "cell_line_name": "LINE-B",
                "oncotree_lineage": "Lung",
                "oncotree_primary_disease": "Non-Small Cell Lung Cancer",
                "oncotree_subtype": "Lung Adenocarcinoma",
                "oncotree_code": "LUAD",
                "mcl_system_ru": "Дыхательная система",
                "mcl_organ_id": "легкое",
                "mcl_organ_ru": "Лёгкое",
                "mcl_organ_icon": "lung",
                "mcl_cancer_id": "lung-nsclc",
                "mcl_cancer_name": "Non-Small Cell Lung Cancer",
                "mcl_subtype_name": "Lung Adenocarcinoma",
                "mcl_subtype_id": "lung-nsclc-luad",
                "classification_status": "classified",
                "classification_confidence": "high",
                "primary_or_metastasis": "Metastasis",
                "sample_collection_site": "Liver",
                "has_crispr": True,
                "depmap_release": "test",
            },
        ],
    )

    catalog = CRISPRModelCatalog(tmp_path)
    atlas = catalog.atlas()

    assert atlas["canonical_index_available"] is True
    assert atlas["models_n"] == 2
    assert atlas["organs_n"] == 1
    assert atlas["cancers_n"] == 1
    assert atlas["organs"][0]["name_ru"] == "Лёгкое"
    assert atlas["organs"][0]["cancers"][0]["models_n"] == 2

    models = catalog.models(organ_id="легкое", cancer_id="lung-nsclc")
    assert models["total"] == 2
    assert len({row["model_id"] for row in models["items"]}) == 2

    metastatic = catalog.model("ACH-2")
    assert metastatic["mcl_organ_ru"] == "Лёгкое"
    assert metastatic["primary_or_metastasis"] == "Metastasis"
    assert metastatic["sample_collection_site"] == "Liver"


def test_crispr_catalog_falls_back_to_context_audit(tmp_path: Path):
    _write_tsv(
        tmp_path / "data/processed/depmap_context_audit.tsv",
        [
            {
                "cancer_id": "C1",
                "model_id": "ACH-1",
                "cell_line_name": "LINE-A",
                "oncotree_lineage": "Lung",
                "oncotree_primary_disease": "Non-Small Cell Lung Cancer",
                "oncotree_subtype": "Lung Adenocarcinoma",
                "assigned_group": "context",
                "sequencing_available": True,
            },
            {
                "cancer_id": "C2",
                "model_id": "ACH-1",
                "cell_line_name": "LINE-A",
                "oncotree_lineage": "Lung",
                "oncotree_primary_disease": "Non-Small Cell Lung Cancer",
                "oncotree_subtype": "Lung Adenocarcinoma",
                "assigned_group": "comparator",
                "sequencing_available": True,
            },
        ],
    )

    catalog = CRISPRModelCatalog(tmp_path)
    atlas = catalog.atlas()
    models = catalog.models()

    assert atlas["canonical_index_available"] is False
    assert atlas["status"] == "fallback_context_audit"
    assert atlas["models_n"] == 1
    assert models["total"] == 1
    assert models["items"][0]["curated_contexts_n"] == 2
