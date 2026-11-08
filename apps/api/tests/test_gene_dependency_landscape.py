from __future__ import annotations

import pandas as pd

from mcl_api.dependency_gene_explorer import DependencyAwareGeneExplorerStore


def test_relationship_strength_labels_are_interpretable():
    assert DependencyAwareGeneExplorerStore._strength_label(0.04) == "практически отсутствует"
    assert DependencyAwareGeneExplorerStore._strength_label(-0.18) == "очень слабая"
    assert DependencyAwareGeneExplorerStore._strength_label(-0.24) == "слабая"
    assert DependencyAwareGeneExplorerStore._strength_label(0.51) == "умеренная"
    assert DependencyAwareGeneExplorerStore._strength_label(-0.72) == "сильная"


def test_dependency_landscape_ranks_cancer_groups_by_dependency_fraction():
    store = DependencyAwareGeneExplorerStore.__new__(DependencyAwareGeneExplorerStore)
    catalog_row = pd.Series({
        "dependency_type": "selective",
        "dependency_type_ru": "Селективная зависимость",
        "specificity_score": 0.55,
        "specificity_label_ru": "очень высокая",
    })
    store._gene_or_raise = lambda symbol: ("GENE1", catalog_row)
    store._model_gene_layer = lambda layer, symbol: pd.DataFrame({
        "model_id": [f"M{i}" for i in range(12)],
        "gene_effect": [-1.2, -1.0, -0.8, -0.7, -0.6, -0.1, -0.2, -0.15, -0.05, -0.1, -0.2, -0.3],
    })
    store._crispr_atlas_frame = lambda: pd.DataFrame({
        "model_id": [f"M{i}" for i in range(12)],
        "mcl_cancer_id": ["A"] * 6 + ["B"] * 6,
        "mcl_cancer_name": ["Cancer A"] * 6 + ["Cancer B"] * 6,
        "mcl_organ_ru": ["Орган A"] * 6 + ["Орган B"] * 6,
        "mcl_system_ru": ["Система"] * 12,
    })

    payload = store.dependency_landscape("GENE1")

    assert payload["available"] is True
    assert payload["summary"]["models_n"] == 12
    assert payload["summary"]["dependent_models_n"] == 5
    assert payload["cancers"][0]["cancer_name"] == "Cancer A"
    assert payload["cancers"][0]["dependent_models_n"] == 5
    assert payload["cancers"][0]["eligible_for_label"] is True
