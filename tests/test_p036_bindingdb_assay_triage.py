"""Synthetic P0.3.6 assay triage contract tests (not clinical or ligand evidence)."""
from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "triage_bindingdb_assays_p036.py"
SPEC = importlib.util.spec_from_file_location("triage_bindingdb_assays_p036", SCRIPT)
MOD = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(MOD)


def sample(tmp_path):
    folder = tmp_path / "links"
    folder.mkdir()
    queue = folder / "assay_review_queue.tsv"
    report = folder / "assay_link_report.json"
    names = [
        "scope", "target_gene", "uniprot_accession", "reactant_set_id",
        "monomer_id", "smiles", "article_doi", "Ki_nM_raw",
        "Kd_nM_raw", "IC50_nM_raw", "entryid_assayid",
        "assay_name", "assay_description", "link_status", "review_decision",
    ]
    rows = [
        ["core", "EGFR", "P00533", "100", "M1", "CCO", "10.1234/a", "10", "", "", "285_1",
         "binding assay", "Recombinant EGFR binding", MOD.LINKED, ""],
        ["core", "EGFR", "P00533", "101", "M2", "CCN", "10.1234/a", "<20", "", "", "285_1",
         "binding assay", "Recombinant EGFR binding", MOD.LINKED, ""],
        ["expanded", "PARP1", "P09874", "102", "M3", "CCC", "", "", "", "not_numeric", "300_2",
         "functional", "In vitro cell assay", MOD.LINKED, ""],
        ["core", "EGFR", "P00533", "103", "M4", "CCCC", "10.1234/b", "", "", "50", "",
         "", "", "", "NO_ASSAY_MAPPING", ""],
    ]
    with queue.open("w", encoding="utf-8", newline="") as file:
        writer = csv.writer(file, delimiter="\t")
        writer.writerow(names)
        writer.writerows(rows)
    report.write_text(json.dumps({
        "contract": "mcl-p036-bindingdb-assay-links-v1",
        "source_archives": {"articles": {"sha256": "x"*64}},
        "link_counts": {MOD.LINKED: 3, "NO_ASSAY_MAPPING": 1},
    }), encoding="utf-8")
    return queue, report


def test_assay_inventory_and_balanced_review_candidates(tmp_path):
    queue, evidence = sample(tmp_path)
    report, inventory, pairs, first_pass = MOD.prepare(queue, evidence, per_gene=1)
    assert report["input_rows"] == 4
    assert report["linked_rows"] == 3
    assert report["assay_ids_n"] == 2
    assert report["protein_assay_pairs_n"] == 2
    assert report["genes_with_linked_assays_n"] == 2
    assert report["first_pass_assay_protein_pairs_n"] == 2
    assert report["scientific_status"] == "UNREVIEWED_NOT_APPROVED"
    egfr = next(x for x in inventory if x["entryid_assayid"] == "285_1")
    parp = next(x for x in inventory if x["entryid_assayid"] == "300_2")
    assert egfr["raw_measurements"] == 2
    assert egfr["unique_monomer_ids"] == 2
    assert egfr["censored_values"] == 1
    assert parp["invalid_values"] == 1
    assert "INVALID_RAW_VALUE" in parp["flags"]
    assert "MISSING_ARTICLE_DOI" in parp["flags"]
    assert all(x["verification_status"] == "pending_manual_review" for x in inventory)
    assert all(x["scientific_class"] == "" and x["reviewer"] == "" for x in first_pass)


def test_report_queue_mismatch_blocks_output(tmp_path):
    queue, evidence = sample(tmp_path)
    data = json.loads(evidence.read_text(encoding="utf-8"))
    data["link_counts"][MOD.LINKED] = 4
    evidence.write_text(json.dumps(data), encoding="utf-8")
    with pytest.raises(ValueError, match="totals disagree"):
        MOD.prepare(queue, evidence)


def test_report_cli_preserves_input_and_never_overwrites(tmp_path, monkeypatch):
    queue, evidence = sample(tmp_path)
    output = tmp_path / "review"
    monkeypatch.setattr("sys.argv", [
        "triage_bindingdb_assays_p036.py",
        "--review-queue", str(queue),
        "--link-report", str(evidence),
        "--out-dir", str(output),
    ])
    MOD.main()
    assert (output / "assay_inventory.tsv").exists()
    assert (output / "protein_assay_pairs.tsv").exists()
    assert (output / "first_pass_scientific_review.tsv").exists()
    assert (output / "assay_triage_report.json").exists()
    with pytest.raises(FileExistsError, match="Immutable"):
        MOD.main()
    assert queue.is_file()
    assert evidence.is_file()


def test_no_output_in_upstream_link_dir(tmp_path, monkeypatch):
    queue, evidence = sample(tmp_path)
    monkeypatch.setattr("sys.argv", [
        "triage_bindingdb_assays_p036.py",
        "--review-queue", str(queue),
        "--link-report", str(evidence),
        "--out-dir", str(queue.parent / "unsafe"),
    ])
    with pytest.raises(ValueError, match="outside source folder"):
        MOD.main()
