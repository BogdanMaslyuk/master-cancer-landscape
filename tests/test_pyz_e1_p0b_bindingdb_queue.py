"""Synthetic regression tests for P0.3.6 structural crosswalk, not scientific approval."""
from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

from rdkit import Chem

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "audit_pyz_e1_p0b_bindingdb_queue.py"
SPEC = importlib.util.spec_from_file_location("audit_pyz_e1_p0b_bindingdb_queue", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(MODULE)


def table(path: Path, header: list[str], rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=header, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def fixtures(tmp_path: Path):
    a, b = MODULE.inchikey("CCO"), MODULE.inchikey("CCN")
    assert a == Chem.MolToInchiKey(Chem.MolFromSmiles("CCO"))
    assert b == Chem.MolToInchiKey(Chem.MolFromSmiles("CCN"))
    neighbors = [
        {"own_compound_id": "PYZ-035", "target_gene": "CDK4", "rank_in_target": "1",
         "ligand_id": "TLS-A", "ligand_inchikey": a, "tanimoto": "0.40"},
        {"own_compound_id": "PYZ-055", "target_gene": "BCL2L1", "rank_in_target": "1",
         "ligand_id": "TLS-B", "ligand_inchikey": b, "tanimoto": "0.41"},
    ]
    queue_rows = [
        {"scope": "core", "target_gene": "CDK4", "uniprot_accession": "P11802",
         "reactant_set_id": "100", "monomer_id": "M100", "smiles": "CCO",
         "article_doi": "10.1234/example-a", "Ki_nM_raw": "", "Kd_nM_raw": "<500",
         "IC50_nM_raw": "300", "entryid_assayid": "300_2",
         "assay_name": "Binding", "assay_description": "Protein binding method",
         "link_status": "LINKED_DESCRIPTION_PRESENT_UNREVIEWED",
         "review_decision": "", "reviewer": ""},
        {"scope": "core", "target_gene": "BCL2L1", "uniprot_accession": "Q07817",
         "reactant_set_id": "101", "monomer_id": "M101", "smiles": "CCN",
         "article_doi": "", "Ki_nM_raw": "2000", "Kd_nM_raw": "",
         "IC50_nM_raw": "", "entryid_assayid": "", "assay_name": "",
         "assay_description": "", "link_status": "NO_ASSAY_MAPPING",
         "review_decision": "", "reviewer": ""},
        # Wrong target with same SMILES must not contaminate CDK4.
        {"scope": "core", "target_gene": "EGFR", "uniprot_accession": "P00533",
         "reactant_set_id": "102", "monomer_id": "M102", "smiles": "CCO",
         "article_doi": "", "Ki_nM_raw": "1", "Kd_nM_raw": "",
         "IC50_nM_raw": "", "entryid_assayid": "300_2",
         "assay_name": "Binding", "assay_description": "Protein binding method",
         "link_status": "LINKED_DESCRIPTION_PRESENT_UNREVIEWED",
         "review_decision": "", "reviewer": ""},
        # Wrong ligand in correct target.
        {"scope": "core", "target_gene": "CDK4", "uniprot_accession": "P11802",
         "reactant_set_id": "103", "monomer_id": "M103", "smiles": "CCC",
         "article_doi": "", "Ki_nM_raw": "1", "Kd_nM_raw": "",
         "IC50_nM_raw": "", "entryid_assayid": "300_2",
         "assay_name": "Binding", "assay_description": "Protein binding method",
         "link_status": "LINKED_DESCRIPTION_PRESENT_UNREVIEWED",
         "review_decision": "", "reviewer": ""},
    ]
    neighbors_path = tmp_path / "neighbors.tsv"
    queue_path = tmp_path / "queue.tsv"
    link_path = tmp_path / "link.json"
    triage_path = tmp_path / "triage.json"
    table(neighbors_path, list(neighbors[0]), neighbors)
    table(queue_path, list(MODULE.REQUIRED_QUEUE), queue_rows)
    link_path.write_text(json.dumps({
        "contract": "mcl-p036-bindingdb-assay-links-v1",
        "source_archives": {"articles": {"sha256": "synthetic"}}
    }), encoding="utf-8")
    triage_path.write_text(json.dumps({
        "contract": "mcl-p036-bindingdb-assay-triage-v1",
        "upstream_report_sha256": MODULE.checksum(link_path),
        "assay_queue_sha256": MODULE.checksum(queue_path),
    }), encoding="utf-8")
    return neighbors_path, queue_path, link_path, triage_path, neighbors


def test_exact_target_and_inchikey_keep_independent_endpoint_rows(tmp_path):
    _, queue, _, _, neighbors = fixtures(tmp_path)
    summary, measurements, counts = MODULE.analyze(
        neighbors, queue, {"contract": "mcl-p036-bindingdb-assay-links-v1"}
    )
    assert len(summary) == 2
    assert len(measurements) == 3
    assert counts["neighbors_with_queue_match"] == 2
    assert counts["matched_assay_keys_unreviewed"] == 1
    assert counts["primary_assays_approved"] == 0
    a = [r for r in measurements if r["target_gene"] == "CDK4"]
    assert len(a) == 2
    assert {r["endpoint_type"]: r["raw_value_nm"] for r in a} == {"Kd": "<500", "IC50": "300"}
    assert {r["reactant_set_id"] for r in a} == {"100"}
    assert all(r["candidate_status"] == "EXACT_STRUCTURE_LINKED_ASSAY_UNREVIEWED" for r in a)
    assert all(r["measurement_owner"] == "known_bindingdb_ligand_not_pyz" for r in measurements)
    assert all(r["scientific_approval"] == "NOT_APPROVED_MANUAL_REVIEW_REQUIRED" for r in measurements)
    b = next(r for r in measurements if r["target_gene"] == "BCL2L1")
    assert b["candidate_status"] == "EXACT_STRUCTURE_SOURCE_LINK_INCOMPLETE"
    assert b["entryid_assayid"] == ""
    assert summary[0]["link_status"] == "STRUCTURE_FOUND_P036_PENDING_MANUAL_REVIEW"


def test_cli_respects_hashes_and_does_not_overwrite_existing_output(tmp_path):
    neighbors, queue, link, triage, _ = fixtures(tmp_path)
    output = tmp_path / "results" / "run1"
    command = [
        sys.executable, str(SCRIPT),
        "--neighbors", str(neighbors),
        "--review-queue", str(queue),
        "--link-report", str(link),
        "--triage-report", str(triage),
        "--out-dir", str(output),
    ]
    source_bytes = {p: p.read_bytes() for p in (neighbors, queue, link, triage)}
    completed = subprocess.run(command, capture_output=True, text=True)
    assert completed.returncode == 0, completed.stderr
    report = json.loads((output / "bindingdb_crosswalk_report.json").read_text(encoding="utf-8"))
    assert report["triage_link_integrity"] == "link_report_and_queue_sha256_match"
    assert report["stats"]["matched_endpoint_rows"] == 3
    assert report["scientific_status"] == "NO_PRIMARY_ASSAYS_APPROVED"
    assert (output / "bindingdb_neighbor_matches.tsv").is_file()
    assert (output / "bindingdb_neighbor_coverage.tsv").is_file()
    assert all(p.read_bytes() == data for p, data in source_bytes.items())
    repeat = subprocess.run(command, capture_output=True, text=True)
    assert repeat.returncode != 0
    assert "FileExistsError" in repeat.stderr


def test_wrong_queue_hash_blocks_crosswalk(tmp_path):
    neighbors, queue, link, triage, _ = fixtures(tmp_path)
    queue.write_text(queue.read_text(encoding="utf-8") + "\n", encoding="utf-8")
    out = tmp_path / "results" / "invalid"
    completed = subprocess.run([
        sys.executable, str(SCRIPT),
        "--neighbors", str(neighbors),
        "--review-queue", str(queue),
        "--link-report", str(link),
        "--triage-report", str(triage),
        "--out-dir", str(out),
    ], capture_output=True, text=True)
    assert completed.returncode != 0
    assert "different assay review queue SHA" in completed.stderr
    assert not out.exists()


def test_mismatched_neighbor_inchikey_not_selected(tmp_path):
    _, queue, _, _, neighbors = fixtures(tmp_path)
    neighbors[0]["ligand_inchikey"] = "NONMATCHING-KEY"
    summary, measurements, counts = MODULE.analyze(
        neighbors, queue, {"contract": "mcl-p036-bindingdb-assay-links-v1"}
    )
    assert counts["neighbors_with_queue_match"] == 1
    assert all(x["target_gene"] != "CDK4" for x in measurements)
    assert summary[0]["link_status"] == "NO_SAME_STRUCTURE_IN_P036_ARTICLES_PILOT"
