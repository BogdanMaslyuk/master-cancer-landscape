#!/usr/bin/env python3
"""P0.3.6: inventory linked BindingDB assays before human scientific review.

One assay key can describe many measurement rows and target proteins.
This script creates an inspectable REVIEW QUEUE, not automatic A/B approvals.
Input is the immutable output of link_bindingdb_assays_p036.py.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
from pathlib import Path

csv.field_size_limit(1 << 28)
LINKED = "LINKED_DESCRIPTION_PRESENT_UNREVIEWED"
VALUE = re.compile(r"^\s*(?:<=|>=|<|>|=|≤|≥)?\s*(\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)\s*$")
KEY = re.compile(r"^\d+_\d+$")
ENDPOINTS = {"Ki": "Ki_nM_raw", "Kd": "Kd_nM_raw", "IC50": "IC50_nM_raw"}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for data in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(data)
    return h.hexdigest()


def inspect_value(raw: str) -> tuple[bool, bool]:
    raw = raw.strip()
    if not raw:
        return False, False
    match = VALUE.fullmatch(raw)
    if not match:
        return False, True
    try:
        valid = Decimal(match.group(1))
        if not valid.is_finite() or valid <= 0:
            return False, True
    except (InvalidOperation, ValueError):
        return False, True
    return raw.startswith(("<", ">", "≤", "≥")), False


def prepare(queue: Path, upstream: Path, *, per_gene: int = 1) -> tuple[dict, list[dict], list[dict], list[dict]]:
    if per_gene < 1 or per_gene > 3:
        raise ValueError("per_gene must be 1..3")
    evidence = json.loads(upstream.read_text(encoding="utf-8"))
    if evidence.get("contract") != "mcl-p036-bindingdb-assay-links-v1":
        raise ValueError("Unrecognized assay-link report contract")
    with queue.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        required = {
            "scope", "target_gene", "uniprot_accession", "reactant_set_id",
            "monomer_id", "smiles", "article_doi", *ENDPOINTS.values(),
            "entryid_assayid", "assay_name", "assay_description",
            "link_status", "review_decision",
        }
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"Missing required review queue columns: {sorted(required - set(reader.fieldnames or []))}")
        assays = {}
        pairs = {}
        status_count = Counter()
        rows_scanned = 0
        linked_rows = 0
        for row in reader:
            rows_scanned += 1
            if None in row:
                raise ValueError(f"Malformed row at {rows_scanned + 1}")
            status = row["link_status"].strip()
            status_count[status] += 1
            if status != LINKED:
                continue
            linked_rows += 1
            assay_id = row["entryid_assayid"].strip()
            if not KEY.fullmatch(assay_id) or not row["assay_description"].strip():
                raise ValueError(f"Invalid linked assay metadata at row {rows_scanned+1}")
            gene = row["target_gene"].strip()
            acc = row["uniprot_accession"].strip()
            if not gene or not acc or not row["reactant_set_id"].strip():
                raise ValueError(f"Missing protein/measurement identity at row {rows_scanned+1}")
            pair_key = (assay_id, gene, acc)
            data = assays.setdefault(assay_id, {
                "entryid_assayid": assay_id,
                "assay_names": set(), "descriptions": set(),
                "genes": set(), "accessions": set(),
                "monomers": set(), "dois": set(), "rsids": set(),
                "raw_measurements": 0, "Ki_rows": 0, "Kd_rows": 0,
                "IC50_rows": 0, "missing_doi_rows": 0,
                "censored_values": 0, "invalid_values": 0,
                "missing_smiles_rows": 0,
            })
            data["assay_names"].add(row["assay_name"].strip())
            data["descriptions"].add(row["assay_description"].strip())
            data["genes"].add(gene)
            data["accessions"].add(acc)
            data["monomers"].add(row["monomer_id"].strip())
            data["rsids"].add(row["reactant_set_id"].strip())
            data["raw_measurements"] += 1
            doi = row["article_doi"].strip().lower()
            if doi:
                data["dois"].add(doi)
            else:
                data["missing_doi_rows"] += 1
            if not row["smiles"].strip():
                data["missing_smiles_rows"] += 1
            pair = pairs.setdefault(pair_key, {
                "entryid_assayid": assay_id,
                "target_gene": gene,
                "uniprot_accession": acc,
                "scope": row["scope"].strip(),
                "raw_measurements": 0, "Ki_rows": 0, "Kd_rows": 0,
                "IC50_rows": 0, "monomers": set(), "dois": set(),
            })
            pair["raw_measurements"] += 1
            pair["monomers"].add(row["monomer_id"].strip())
            if doi:
                pair["dois"].add(doi)
            for endpoint, field in ENDPOINTS.items():
                raw = row[field].strip()
                if not raw:
                    continue
                data[endpoint + "_rows"] += 1
                pair[endpoint + "_rows"] += 1
                censored, invalid = inspect_value(raw)
                data["censored_values"] += int(censored)
                data["invalid_values"] += int(invalid)
    declared = {k: int(v) for k, v in evidence.get("link_counts", {}).items()}
    if dict(status_count) != declared:
        raise ValueError(f"Queue/report link-status totals disagree: observed={dict(status_count)}, declared={declared}")
    inventory = []
    for key in sorted(assays):
        item = assays[key]
        descriptions = sorted(item["descriptions"])
        names = sorted(item["assay_names"])
        flags = []
        if len(descriptions) != 1 or len(names) != 1:
            flags.append("CONFLICTING_METHOD_TEXT")
        if len(item["genes"]) > 1:
            flags.append("METHOD_LINKED_TO_MULTIPLE_GENES")
        if item["missing_doi_rows"]:
            flags.append("MISSING_ARTICLE_DOI")
        if item["invalid_values"]:
            flags.append("INVALID_RAW_VALUE")
        if item["censored_values"]:
            flags.append("CENSORED_VALUE")
        if item["missing_smiles_rows"]:
            flags.append("MISSING_SMILES")
        inventory.append({
            "entryid_assayid": key,
            "assay_name": names[0] if len(names) == 1 else "",
            "assay_description": descriptions[0] if len(descriptions) == 1 else "",
            "target_genes": ";".join(sorted(item["genes"])),
            "uniprot_accessions": ";".join(sorted(item["accessions"])),
            "target_genes_n": len(item["genes"]),
            "raw_measurements": item["raw_measurements"],
            "unique_reactant_set_ids": len(item["rsids"]),
            "unique_monomer_ids": len(item["monomers"] - {""}),
            "article_dois_n": len(item["dois"]),
            "article_dois": ";".join(sorted(item["dois"])),
            "Ki_rows": item["Ki_rows"], "Kd_rows": item["Kd_rows"],
            "IC50_rows": item["IC50_rows"],
            "censored_values": item["censored_values"],
            "invalid_values": item["invalid_values"],
            "missing_smiles_rows": item["missing_smiles_rows"],
            "flags": ";".join(flags),
            "scientific_class": "",
            "verification_status": "pending_manual_review",
            "reviewer": "",
            "evidence_url": "",
            "review_notes": "",
        })
    inventory_by_key = {row["entryid_assayid"]: row for row in inventory}
    by_gene = defaultdict(list)
    pair_records = []
    for (assay_id, gene, acc), row in sorted(pairs.items()):
        parent = inventory_by_key[assay_id]
        pair_row = {
            "scope": row["scope"], "target_gene": gene, "uniprot_accession": acc,
            "entryid_assayid": assay_id,
            "assay_name": parent["assay_name"],
            "assay_description_excerpt": parent["assay_description"][:350],
            "raw_measurements": row["raw_measurements"],
            "unique_monomer_ids": len(row["monomers"] - {""}),
            "article_dois": ";".join(sorted(row["dois"])),
            "Ki_rows": row["Ki_rows"], "Kd_rows": row["Kd_rows"],
            "IC50_rows": row["IC50_rows"],
            "flags": parent["flags"],
            "scientific_class": "",
            "verification_status": "pending_manual_review",
            "reviewer": "", "review_notes": "",
        }
        pair_records.append(pair_row)
        if "CONFLICTING_METHOD_TEXT" not in parent["flags"]:
            by_gene[gene].append(pair_row)
    first_pass = []
    for gene in sorted(by_gene):
        ranked = sorted(by_gene[gene],
                        key=lambda x: (-x["raw_measurements"], -x["unique_monomer_ids"], x["entryid_assayid"]))
        first_pass.extend(ranked[:per_gene])
    first_pass = sorted(first_pass, key=lambda x: (x["scope"] != "core", x["target_gene"],
                                                  -x["raw_measurements"], x["entryid_assayid"]))
    report = {
        "contract": "mcl-p036-bindingdb-assay-triage-v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "assay_queue_sha256": sha256(queue),
        "upstream_report_sha256": sha256(upstream),
        "upstream_source_archives": evidence.get("source_archives", {}),
        "input_rows": rows_scanned, "linked_rows": linked_rows,
        "assay_ids_n": len(inventory),
        "protein_assay_pairs_n": len(pair_records),
        "genes_with_linked_assays_n": len(by_gene),
        "first_pass_assay_protein_pairs_n": len(first_pass),
        "first_pass_limit_per_gene": per_gene,
        "linked_status_counts": dict(status_count),
        "scientific_status": "UNREVIEWED_NOT_APPROVED",
        "methodology": "First pass samples high-row-count assay per gene for inspection, NOT predicted experimental quality.",
        "limitations": [
            "Assay metadata is not a scientific assay audit or evidence of physical binding.",
            "Descriptions are unstructured, and numeric records/DOIs require primary-source verification.",
            "Multiple assay-protein pairs may share one assay key; this is not an independent experiment count.",
            "No A/B quality grade is assigned; all manual review and primary-paper verification remain pending.",
        ],
    }
    return report, inventory, pair_records, first_pass


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"No rows for {path.name}; incomplete source?")
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--review-queue", type=Path, required=True)
    parser.add_argument("--link-report", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--per-gene", type=int, default=1)
    args = parser.parse_args()
    out = args.out_dir.resolve()
    if out.exists():
        raise FileExistsError(f"Immutable report already exists: {out}")
    for source in (args.review_queue, args.link_report):
        if out.is_relative_to(source.resolve().parent):
            raise ValueError(f"Reports must be outside source folder: {source.parent}")
    report, assays, pairs, first = prepare(
        args.review_queue, args.link_report, per_gene=args.per_gene)
    out.mkdir(parents=True)
    write_tsv(out/"assay_inventory.tsv", assays)
    write_tsv(out/"protein_assay_pairs.tsv", pairs)
    write_tsv(out/"first_pass_scientific_review.tsv", first)
    (out/"assay_triage_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"Measurement rows: {report['input_rows']}")
    print(f"Rows linked to described assays: {report['linked_rows']}")
    print(f"Unique assay keys: {report['assay_ids_n']}")
    print(f"Protein-assay pairs: {report['protein_assay_pairs_n']}")
    print(f"Proteins with linked methods: {report['genes_with_linked_assays_n']}")
    print(f"First-pass protein-assay review items: {report['first_pass_assay_protein_pairs_n']}")
    print(f"Scientific status: {report['scientific_status']}")
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
