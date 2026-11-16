"""E1-P0-B: exact structural crosswalk to already-linked BindingDB P0.3.6 rows.

This is READ-ONLY and never approves assays. It joins (gene, standardized
InChIKey) from v1 top-hits to the existing BindingDB Articles assay review queue.
A structurally matched Articles row is NOT proof of an identical v1 REST record.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

FOCUS = {("PYZ-035", "CDK4"), ("PYZ-055", "BCL2L1")}
ASSAY_KEY = re.compile(r"^[0-9]+_[0-9]+$")
VALUE = re.compile(r"^\s*(<=|>=|<|>|=|≤|≥)?\s*(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?\s*$")
ENDPOINT_COLUMNS = (("Ki", "Ki_nM_raw"), ("Kd", "Kd_nM_raw"), ("IC50", "IC50_nM_raw"))
REQUIRED_NEIGHBORS = {"own_compound_id", "target_gene", "rank_in_target", "ligand_id", "ligand_inchikey"}
REQUIRED_QUEUE = {
    "scope", "target_gene", "uniprot_accession", "reactant_set_id", "monomer_id",
    "smiles", "article_doi", "Ki_nM_raw", "Kd_nM_raw", "IC50_nM_raw",
    "entryid_assayid", "assay_name", "assay_description", "link_status",
    "review_decision", "reviewer",
}
ROW_FIELDS = [
    "own_compound_id", "target_gene", "rank_in_target", "ligand_id",
    "tanimoto", "ligand_inchikey", "queue_row_number", "uniprot_accession",
    "reactant_set_id", "monomer_id", "article_doi", "endpoint_type",
    "raw_value_nm", "entryid_assayid", "assay_name", "assay_description",
    "link_status", "queue_review_decision", "candidate_status",
    "measurement_owner", "scientific_approval",
]
NEIGHBOR_FIELDS = [
    "own_compound_id", "target_gene", "rank_in_target", "ligand_id",
    "ligand_inchikey", "matching_queue_rows", "matching_assay_keys",
    "matching_sources_with_doi", "link_status",
]


def checksum(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for piece in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            h.update(piece)
    return h.hexdigest()


def inchikey(smiles: str) -> str:
    """Use the same chemical standardization as TLS v1 build_target_ligand_space."""
    if not smiles.strip():
        return ""
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=True)
        if mol is None:
            return ""
        cleaned = rdMolStandardize.Cleanup(mol)
        parent = rdMolStandardize.FragmentParent(cleaned)
        canonical = Chem.MolToSmiles(parent, canonical=True, isomericSmiles=True)
        repaired = Chem.MolFromSmiles(canonical, sanitize=True)
        if repaired is None:
            return ""
        repaired.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(repaired)
        Chem.GetSymmSSSR(repaired)
        return Chem.MolToInchiKey(repaired) or ""
    except Exception:
        return ""


def _read_neighbors(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if not REQUIRED_NEIGHBORS.issubset(reader.fieldnames or []):
            raise ValueError("neighbors missing required columns")
        rows = [r for r in reader if (r["own_compound_id"], r["target_gene"]) in FOCUS]
    # Multiple rows with same target+ligand can invalidate source linkage.
    keys = [(r["own_compound_id"], r["target_gene"], r["ligand_id"]) for r in rows]
    if len(keys) != len(set(keys)):
        raise ValueError("Duplicate neighbor identity")
    for r in rows:
        if not r["ligand_inchikey"]:
            raise ValueError("Missing neighbor InChIKey; cannot crosswalk")
    return rows


def analyze(neighbors: list[dict[str, str]], queue: Path, link_report: dict[str, Any]
           ) -> tuple[list[dict], list[dict], dict]:
    if link_report.get("contract") != "mcl-p036-bindingdb-assay-links-v1":
        raise ValueError("Unexpected P0.3.6 linkage report contract")
    source_by_key: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for neighbor in neighbors:
        source_by_key[(neighbor["target_gene"].upper(), neighbor["ligand_inchikey"])].append(neighbor)

    matched = []
    scanned = 0
    target_rows = 0
    invalid_smiles = 0
    computed_keys: dict[str, str] = {}
    candidate_link_status = Counter()
    with queue.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        missing = sorted(REQUIRED_QUEUE - set(reader.fieldnames or []))
        if missing:
            raise ValueError("assay review queue missing columns: " + str(missing))
        for number, row in enumerate(reader, start=2):
            scanned += 1
            if None in row:
                raise ValueError("Malformed P0.3.6 queue row at " + str(number))
            gene = (row.get("target_gene") or "").strip().upper()
            if gene not in {g for _, g in FOCUS}:
                continue
            target_rows += 1
            smiles = (row.get("smiles") or "").strip()
            if smiles not in computed_keys:
                computed_keys[smiles] = inchikey(smiles)
            chemical_id = computed_keys[smiles]
            if not chemical_id:
                invalid_smiles += 1
                continue
            for neighbor in source_by_key.get((gene, chemical_id), []):
                assay = (row["entryid_assayid"] or "").strip()
                linked = (
                    row["link_status"] == "LINKED_DESCRIPTION_PRESENT_UNREVIEWED"
                    and ASSAY_KEY.fullmatch(assay) is not None
                    and bool((row["assay_description"] or "").strip())
                )
                candidate_status = (
                    "EXACT_STRUCTURE_LINKED_ASSAY_UNREVIEWED" if linked
                    else "EXACT_STRUCTURE_SOURCE_LINK_INCOMPLETE"
                )
                candidate_link_status[candidate_status] += 1
                for endpoint, column in ENDPOINT_COLUMNS:
                    value = (row[column] or "").strip()
                    if not value:
                        continue
                    matched.append({
                        "own_compound_id": neighbor["own_compound_id"],
                        "target_gene": gene,
                        "rank_in_target": neighbor["rank_in_target"],
                        "ligand_id": neighbor["ligand_id"],
                        "tanimoto": neighbor.get("tanimoto", ""),
                        "ligand_inchikey": chemical_id,
                        "queue_row_number": number,
                        "uniprot_accession": (row["uniprot_accession"] or "").strip(),
                        "reactant_set_id": (row["reactant_set_id"] or "").strip(),
                        "monomer_id": (row["monomer_id"] or "").strip(),
                        "article_doi": (row["article_doi"] or "").strip(),
                        "endpoint_type": endpoint,
                        "raw_value_nm": value,
                        "entryid_assayid": assay if linked else "",
                        "assay_name": (row["assay_name"] or "").strip() if linked else "",
                        "assay_description": (row["assay_description"] or "").strip() if linked else "",
                        "link_status": row["link_status"],
                        "queue_review_decision": (row["review_decision"] or "").strip(),
                        "candidate_status": (
                            candidate_status if VALUE.fullmatch(value)
                            else "INVALID_RAW_VALUE_REQUIRES_REVIEW"
                        ),
                        "measurement_owner": "known_bindingdb_ligand_not_pyz",
                        "scientific_approval": "NOT_APPROVED_MANUAL_REVIEW_REQUIRED",
                    })
    n_rows = []
    match_by_ligand = defaultdict(list)
    for row in matched:
        match_by_ligand[(row["own_compound_id"], row["target_gene"], row["ligand_id"])].append(row)
    for original in neighbors:
        pair = (original["own_compound_id"], original["target_gene"], original["ligand_id"])
        found = match_by_ligand.get(pair, [])
        n_rows.append({
            "own_compound_id": original["own_compound_id"],
            "target_gene": original["target_gene"],
            "rank_in_target": original["rank_in_target"],
            "ligand_id": original["ligand_id"],
            "ligand_inchikey": original["ligand_inchikey"],
            "matching_queue_rows": len({m["queue_row_number"] for m in found}),
            "matching_assay_keys": len({m["entryid_assayid"] for m in found if m["entryid_assayid"]}),
            "matching_sources_with_doi": len({m["queue_row_number"] for m in found if m["article_doi"]}),
            "link_status": (
                "STRUCTURE_FOUND_P036_PENDING_MANUAL_REVIEW" if found
                else "NO_SAME_STRUCTURE_IN_P036_ARTICLES_PILOT"
            ),
        })
    stats = {
        "queue_rows_scanned": scanned,
        "queue_rows_in_focus_genes": target_rows,
        "queue_rows_invalid_smiles_in_focus_genes": invalid_smiles,
        "unique_neighbor_structures": len(neighbors),
        "neighbors_with_queue_match": sum(x["matching_queue_rows"] > 0 for x in n_rows),
        "matched_endpoint_rows": len(matched),
        "matched_reactant_set_ids": len({x["reactant_set_id"] for x in matched if x["reactant_set_id"]}),
        "matched_assay_keys_unreviewed": len({x["entryid_assayid"] for x in matched if x["entryid_assayid"]}),
        "matched_dois_unreviewed": len({x["article_doi"] for x in matched if x["article_doi"]}),
        "candidate_status_counts": dict(sorted(Counter(x["candidate_status"] for x in matched).items())),
        "source_link_status_counts_for_matched_records": dict(sorted(Counter(x["link_status"] for x in matched).items())),
        "primary_assays_approved": 0,
    }
    return n_rows, matched, stats


def write_table(path: Path, fields: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8-sig") as output:
        writer = csv.DictWriter(output, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--neighbors", required=True, type=Path)
    parser.add_argument("--review-queue", required=True, type=Path)
    parser.add_argument("--link-report", required=True, type=Path)
    parser.add_argument("--triage-report", type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args()
    files = [args.neighbors, args.review_queue, args.link_report]
    if args.triage_report:
        files.append(args.triage_report)
    for path in files:
        if not path.is_file():
            raise FileNotFoundError(path)
    out = args.out_dir.resolve()
    if out.exists():
        raise FileExistsError(out)
    for path in files:
        if out.is_relative_to(path.resolve().parent):
            raise ValueError("Output must not be written under source directories")

    source = json.loads(args.link_report.read_text(encoding="utf-8-sig"))
    triage_hash_status = "not_provided"
    if args.triage_report:
        triage = json.loads(args.triage_report.read_text(encoding="utf-8-sig"))
        if triage.get("contract") != "mcl-p036-bindingdb-assay-triage-v1":
            raise ValueError("Unrecognized triage contract")
        observed_link_hash = checksum(args.link_report)
        declared_link_hash = triage.get("upstream_report_sha256", "")
        if declared_link_hash.lower() != observed_link_hash.lower():
            raise ValueError("Triage report refers to different assay linkage report SHA")
        declared_queue_hash = triage.get("assay_queue_sha256", "")
        observed_queue_hash = checksum(args.review_queue)
        if declared_queue_hash.lower() != observed_queue_hash.lower():
            raise ValueError("Triage report refers to different assay review queue SHA")
        triage_hash_status = "link_report_and_queue_sha256_match"

    neighbors = _read_neighbors(args.neighbors)
    matched_neighbors, rows, stats = analyze(neighbors, args.review_queue, source)
    output = {
        "contract": "mcl-pyz-e1-p0b-bindingdb-articles-crosswalk-v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {str(path): checksum(path) for path in files},
        "triage_link_integrity": triage_hash_status,
        "source_archive_hashes_declared_in_upstream_report": source.get("source_archives", {}),
        "scientific_status": "NO_PRIMARY_ASSAYS_APPROVED",
        "identity_basis": "target_gene_and_v1_standardized_InChIKey",
        "identity_limitations": (
            "Matching a BindingDB Articles row to a v1 ligand structure does not "
            "prove this was the exact original v1 REST measurement. "
            "Do not equate the two experimental records without reactant-set provenance."
        ),
        "stats": stats,
        "limitations": [
            "This crosswalk is observational and does not approve binding assays, papers or direct target engagement.",
            "The original v1 ligand measurements/catalog remain unavailable and are not recreated.",
            "BindingDB Articles is a potentially incomplete pilot subset; a missing match is not proof of no experimental data.",
            "No assay is created for a v1 REST source row, and no assay identifier is transferred to PYZ.",
            "Source archive SHA values are inherited metadata; this run does not rehash all original zipped archives.",
            "Review every complex/construct, species, assay method, article and numeric qualifier before publication.",
        ],
    }
    out.mkdir(parents=True)
    write_table(out / "bindingdb_neighbor_matches.tsv", ROW_FIELDS, rows)
    write_table(out / "bindingdb_neighbor_coverage.tsv", NEIGHBOR_FIELDS, matched_neighbors)
    (out / "bindingdb_crosswalk_report.json").write_text(
        json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print("E1-P0-B BindingDB P0.3.6 structural crosswalk")
    print("Matched analogues: {} / {}".format(stats["neighbors_with_queue_match"], len(neighbors)))
    print("Candidate endpoint rows: {}".format(stats["matched_endpoint_rows"]))
    print("Unique mapped assay keys (unreviewed): {}".format(stats["matched_assay_keys_unreviewed"]))
    print("Primary-source approvals: 0")
    print("Output: " + str(out))


if __name__ == "__main__":
    main()
