from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from mcl_api.atlas import MCLAtlas
from mcl_api.store import MCLDataStore


def _write_tsv(path: Path, rows: list[dict]):
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def _build_fixture(tmp_path: Path) -> MCLAtlas:
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

    return MCLAtlas(tmp_path, MCLDataStore(tmp_path))


def test_atlas_builds_disease_model_hierarchy(tmp_path: Path):
    atlas = _build_fixture(tmp_path)

    overview = atlas.atlas()
    assert overview["organs_n"] == 1
    assert overview["contexts_n"] == 1
    assert overview["models_n"] == 2
    overview_context = overview["organs"][0]["contexts"][0]
    assert overview_context["model_genetics"]["status"] == "deferred"

    context = atlas.context("C1")
    assert context["analysis_available"] is True
    assert context["context_models_n"] == 1
    assert context["comparator_models_n"] == 1
    assert context["patient_layer"]["status"] == "not_connected"
    assert context["model_genetics"]["available"] is False

    models = atlas.models(cancer_id="C1")
    assert models["total"] == 2
    assert models["items"][0]["variants"][0]["gene"] == "KRAS"
    assert models["items"][0]["variants"][0]["protein_change"] == "p.G12C"


def test_atlas_overview_does_not_load_full_mutation_profiles(tmp_path: Path, monkeypatch):
    atlas = _build_fixture(tmp_path)

    def fail_if_called(*args, **kwargs):
        raise AssertionError("Atlas overview must not load full mutation profiles")

    monkeypatch.setattr(atlas, "_context_model_genetics", fail_if_called)

    overview = atlas.atlas()
    context = overview["organs"][0]["contexts"][0]
    assert context["model_genetics"]["status"] == "deferred"


def test_model_uses_full_indexed_mutation_profile_when_available(tmp_path: Path):
    atlas = _build_fixture(tmp_path)
    _write_tsv(
        tmp_path / "data/processed/depmap_model_metadata.tsv",
        [
            {
                "model_id": "ACH-1",
                "patient_id": "PT-1",
                "primary_or_metastasis": "Primary",
                "sample_collection_site": "Lung",
                "sex": "Male",
                "age": 58,
            }
        ],
    )
    _write_tsv(
        tmp_path / "data/processed/depmap_model_mutations.tsv",
        [
            {
                "model_id": "ACH-1", "gene": "KRAS", "protein_change": "p.G12C",
                "dna_change": "c.34G>T", "variant_type": "SNP", "vep_impact": "MODERATE",
                "driver": True, "hotspot": True, "likely_lof": False, "is_functional": True,
                "priority_score": 10, "allele_fraction": 0.41,
            },
            {
                "model_id": "ACH-1", "gene": "STK11", "protein_change": "p.Q37*",
                "dna_change": "c.109C>T", "variant_type": "SNP", "vep_impact": "HIGH",
                "driver": False, "hotspot": False, "likely_lof": True, "is_functional": True,
                "priority_score": 6, "allele_fraction": 0.52,
            },
        ],
    )

    atlas._model_metadata.cache_clear()
    atlas._model_mutations.cache_clear()

    model = atlas.model("ACH-1")
    assert model["genetics"]["availability"] == "full"
    assert model["genetics"]["mutated_genes_n"] == 2
    assert {x["gene"] for x in model["genetics"]["priority_variants"]} == {"KRAS", "STK11"}
    assert model["metadata"]["primary_or_metastasis"] == "Primary"

    context = atlas.context("C1")
    assert context["model_genetics"]["available"] is True
    genes = {x["gene"] for x in context["model_genetics"]["top_genes"]}
    assert "KRAS" in genes
    assert "STK11" in genes
