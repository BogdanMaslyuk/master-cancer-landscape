"""P0 regression tests: no invented targets, no Candidate v2 filter, reproducible QC."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pandas as pd
import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "build_expanded_target_registry.py"
spec = importlib.util.spec_from_file_location("build_expanded_target_registry", SCRIPT)
module = importlib.util.module_from_spec(spec)
assert spec and spec.loader
spec.loader.exec_module(module)


def _v1():
    return pd.DataFrame([
        {
            "target_gene": "CDK4", "protein_mapping_status": "unique_swissprot",
            "uniprot_primary_accession": "P11802", "protein_preferred_name": "CDK4",
            "uniprot_release": "fixture-1",
        },
        {
            "target_gene": "AMBIG", "protein_mapping_status": "multiple_swissprot",
            "uniprot_primary_accession": "", "protein_preferred_name": "",
        },
    ])


def _uniprot():
    return pd.DataFrame([
        {"Entry": "P11802", "Reviewed": "reviewed", "Gene Names (primary)": "CDK4",
         "Organism (ID)": "9606", "Protein names": "Cyclin-dependent kinase 4"},
        {"Entry": "P00533", "Reviewed": "reviewed", "Gene Names (primary)": "EGFR",
         "Organism (ID)": "9606", "Protein names": "Epidermal growth factor receptor"},
        {"Entry": "P99999", "Reviewed": "unreviewed", "Gene Names (primary)": "NOTOK",
         "Organism (ID)": "9606", "Protein names": "Bad"},
        {"Entry": "P11111", "Reviewed": "reviewed", "Gene Names (primary)": "MOUSE",
         "Organism (ID)": "10090", "Protein names": "Bad organism"},
    ])


def test_expands_beyond_legacy_without_candidate_status():
    targets, qc = module.collect(_v1(), _uniprot())
    assert set(targets.target_gene) == {"CDK4", "EGFR"}
    assert targets.target_gene.is_unique
    assert targets.uniprot_accession.is_unique
    assert set(targets.ligand_coverage_status) == {"not_evaluated"}
    assert set(qc.reason_code) >= {"mapping_not_unique", "not_reviewed", "non_human"}


def test_conflicting_gene_to_protein_mapping_is_excluded():
    uniprot = _uniprot()
    extra = pd.DataFrame([{
        "Entry": "P00534", "Reviewed": "reviewed", "Gene Names (primary)": "EGFR",
        "Organism (ID)": "9606", "Protein names": "Conflicting protein"
    }])
    targets, qc = module.collect(_v1(), pd.concat([uniprot, extra], ignore_index=True))
    assert "EGFR" not in set(targets.target_gene)
    assert "conflicting_gene_protein_mapping" in set(qc.reason_code)


def test_noncanonical_isoform_not_promoted():
    extra = pd.DataFrame([{
        "Entry": "P00533-2", "Reviewed": "reviewed", "Gene Names (primary)": "EGFR",
        "Organism (ID)": "9606", "Protein names": "Isoform"
    }])
    targets, qc = module.collect(_v1(), extra)
    assert "EGFR" not in set(targets.target_gene)
    assert "invalid_identifier" in set(qc.reason_code)


def test_bad_uniprot_schema_fails_closed():
    with pytest.raises(ValueError, match="Organism"):
        module.collect(_v1(), _uniprot().drop(columns=["Organism (ID)"]))


def test_expanded_registry_writes_isolated_artifacts(tmp_path: Path):
    v1 = tmp_path / "source.tsv"
    raw = tmp_path / "uniprot.tsv"
    _v1().to_csv(v1, sep="\t", index=False)
    _uniprot().to_csv(raw, sep="\t", index=False)
    result = module.build(v1, raw, tmp_path / "v2", "snapshot-fixture")
    assert result["candidate_v2_filter_applied"] is False
    assert result["target_genes_n"] == 2
    assert result["universe_is_exhaustive"] is False
    assert result["ligand_or_depmap_coverage_confirmed"] is False
    manifest = json.loads((tmp_path / "v2" / "protein_target_universe_manifest.json").read_text())
    assert manifest["uniprot_release"] == "snapshot-fixture"
    assert len(manifest["source_files"]) == 2
    assert (tmp_path / "v2" / "protein_target_universe_qc.tsv").is_file()


def test_official_snapshot_requires_release(tmp_path: Path):
    v1 = tmp_path / "source.tsv"
    raw = tmp_path / "uniprot.tsv"
    _v1().to_csv(v1, sep="\t", index=False)
    _uniprot().to_csv(raw, sep="\t", index=False)
    with pytest.raises(ValueError, match="uniprot-release"):
        module.build(v1, raw, tmp_path / "out")
    assert not (tmp_path / "out").exists()


def test_uniprot_real_export_gene_column_lowercase_accepted():
    """UniProt live TSV may label its column 'Gene names (primary)'."""
    uniprot = _uniprot().rename(
        columns={"Gene Names (primary)": "Gene names (primary)"}
    )
    targets, qc = module.collect(_v1(), uniprot)
    assert {"EGFR", "CDK4"} == set(targets.target_gene)
    assert "ambiguous_primary_gene" not in set(qc.reason_code)
