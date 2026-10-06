from __future__ import annotations

from pathlib import Path

import pandas as pd

from mcl_api.model_dependencies import ModelDependencyStore


def test_model_dependency_ranking_and_specificity(tmp_path: Path):
    processed = tmp_path / "data" / "processed"
    gene_dir = processed / "gene_explorer"
    gene_dir.mkdir(parents=True)
    (tmp_path / "config").mkdir(parents=True)

    model_ids = [f"ACH-{i:02d}" for i in range(1, 21)]
    matrix = pd.DataFrame({"model_id": model_ids})
    matrix["CORE"] = -1.0
    matrix["SELECTIVE"] = [ -1.2 ] + [0.0] * 19
    matrix["CANCER"] = [-0.9] * 5 + [0.0] * 15
    matrix["WEAK"] = -0.1
    matrix.to_parquet(processed / "depmap_model_gene_effect.parquet", index=False)

    atlas = pd.DataFrame(
        {
            "model_id": model_ids,
            "cell_line_name": [f"LINE-{i}" for i in range(1, 21)],
            "mcl_cancer_id": ["lung-cancer"] * 5 + ["other-cancer"] * 15,
            "mcl_cancer_name": ["Lung cancer"] * 5 + ["Other cancer"] * 15,
            "mcl_organ_ru": ["Лёгкое"] * 5 + ["Другие"] * 15,
            "depmap_release": ["test"] * 20,
        }
    )
    atlas.to_parquet(processed / "depmap_crispr_model_atlas.parquet", index=False)

    catalog = pd.DataFrame(
        {
            "gene_symbol": ["CORE", "SELECTIVE", "CANCER", "WEAK"],
            "gene_name": ["core gene", "selective gene", "cancer gene", "weak gene"],
            "mcl_domains_json": ['["cell_cycle"]', '["signaling"]', '["signaling"]', '[]'],
            "mcl_subdomains_json": ["[]", "[]", "[]", "[]"],
            "protein_classes_json": ["[]", "[]", "[]", "[]"],
        }
    )
    catalog.to_parquet(gene_dir / "gene_catalog.parquet", index=False)
    (tmp_path / "config" / "mcl_functional_domains.yaml").write_text(
        """version: test
domains:
  - id: cell_cycle
    label_ru: Клеточный цикл
  - id: signaling
    label_ru: Сигнальные пути
""",
        encoding="utf-8",
    )

    store = ModelDependencyStore(tmp_path)
    payload = store.model("ACH-01", limit=20)

    assert payload["available"] is True
    assert payload["genes_measured_n"] == 4
    by_gene = {row["gene"]: row for row in payload["items"]}
    assert by_gene["CORE"]["dependency_type"] == "broad_core"
    assert by_gene["SELECTIVE"]["dependency_type"] == "model_selective"
    assert by_gene["CANCER"]["dependency_type"] == "cancer_enriched"
    assert by_gene["WEAK"]["dependency_type"] == "weak_or_none"
    assert by_gene["SELECTIVE"]["rank"] == 1
    assert by_gene["CORE"]["domains"][0]["label_ru"] == "Клеточный цикл"


def test_model_dependency_filters_keep_original_rank(tmp_path: Path):
    processed = tmp_path / "data" / "processed"
    processed.mkdir(parents=True)
    pd.DataFrame(
        {"model_id": ["ACH-1", "ACH-2"], "A": [-1.0, -1.0], "B": [-0.8, 0.0], "C": [-0.1, -0.1]}
    ).to_parquet(processed / "depmap_model_gene_effect.parquet", index=False)
    pd.DataFrame(
        {
            "model_id": ["ACH-1", "ACH-2"],
            "mcl_cancer_id": ["x", "y"],
            "mcl_cancer_name": ["X", "Y"],
            "mcl_organ_ru": ["X", "Y"],
        }
    ).to_parquet(processed / "depmap_crispr_model_atlas.parquet", index=False)

    store = ModelDependencyStore(tmp_path)
    payload = store.model("ACH-1", search="B", limit=10)
    assert payload["filtered_total"] == 1
    assert payload["items"][0]["gene"] == "B"
    assert payload["items"][0]["rank"] == 2
