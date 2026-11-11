#!/usr/bin/env python3
"""P0.3.6: link published BindingDB measurements to exact assay descriptions.

Read ONLY the pinned BindingDB_BindingDB_Articles, rsid_eaids and Assays ZIPs.
Links are preliminary provenance, NOT proof of binding or A/B assay approval.
All source archives remain unmodified. Missing/ambiguous joins fail closed.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
import zipfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

csv.field_size_limit(1 << 28)

ARTICLES_FILE = "BindingDB_BindingDB_Articles.tsv"
MAP_FILE = "BindingDB_rsid_eaids.tsv"
ASSAYS_FILE = "BindingDB_Assays.tsv"
CHAIN_RE = re.compile(r"^UniProt \(SwissProt\) Primary ID of Target Chain(?: ([1-9]\d*))?$")
ASSAY_KEY_RE = re.compile(r"^(\d+)_(\d+)$")
ENDPOINTS = ("Ki (nM)", "Kd (nM)", "IC50 (nM)")
REQUIRED_ARTICLES = {
    "BindingDB Reactant_set_id", "BindingDB MonomerID", "Ligand SMILES",
    *ENDPOINTS,
}

def checksum(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for data in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            h.update(data)
    return h.hexdigest()


def tsv_rows(zip_path: Path, filename: str):
    """Yield positional rows from exactly one TSV member; never extract to disk."""
    with zipfile.ZipFile(zip_path) as archive:
        selected = [m for m in archive.infolist()
                    if Path(m.filename).name == filename and not m.is_dir()]
        if len(selected) != 1:
            raise ValueError(f"Expected exactly one {filename} in {zip_path.name}")
        with archive.open(selected[0]) as binary:
            with io.TextIOWrapper(binary, encoding="utf-8-sig", newline="") as text:
                reader = csv.reader(text, delimiter="\t", quoting=csv.QUOTE_NONE)
                yield from reader


def require_columns(header: list[str], required: set[str], filename: str) -> dict[str, int]:
    if len(set(header)) != len(header):
        raise ValueError(f"{filename} has duplicate column names")
    missing = required - set(header)
    if missing:
        raise ValueError(f"{filename} missing columns {sorted(missing)}")
    return {key: header.index(key) for key in required}


def extract_articles(path: Path, panel: dict[str, dict]) -> tuple[list[dict], dict]:
    rows = tsv_rows(path, ARTICLES_FILE)
    header = [label.strip() for label in next(rows)]
    index = require_columns(header, REQUIRED_ARTICLES, ARTICLES_FILE)
    chain_indexes = [i for i, name in enumerate(header) if CHAIN_RE.fullmatch(name)]
    if not chain_indexes:
        raise ValueError(f"Articles do not expose numbered/unnumbered SwissProt target chain: "
                         f"{[n for n in header if 'UniProt' in n][:20]}")
    count_index = next((i for i, name in enumerate(header)
                        if name.startswith("Number of Protein Chains in Target")), None)
    doi_index = header.index("Article DOI") if "Article DOI" in header else None
    out = []
    tally = Counter()
    for number, cells in enumerate(rows, start=2):
        tally["total_articles_rows"] += 1
        if len(cells) <= max(chain_indexes + list(index.values())):
            tally["short_rows"] += 1
            continue
        accession_values = [cells[i].strip() for i in chain_indexes if cells[i].strip()]
        accessions = set(accession_values)
        matches = sorted(accessions & panel.keys())
        if not matches:
            continue
        tally["matched_panel_rows"] += 1
        chain_count = (cells[count_index].strip()
                       if count_index is not None and count_index < len(cells) else "")
        if len(accessions) != 1 or len(accession_values) != 1 or chain_count not in ("", "1", "1.0"):
            tally["excluded_multichain_rows"] += 1
            continue
        if not any(cells[index[field]].strip() for field in ENDPOINTS):
            continue
        accession = matches[0]
        row = panel[accession]
        tally["eligible_raw_endpoint_rows"] += 1
        out.append({
            "scope": row["scope"],
            "target_gene": row["target_gene"],
            "uniprot_accession": accession,
            "reactant_set_id": cells[index["BindingDB Reactant_set_id"]].strip(),
            "monomer_id": cells[index["BindingDB MonomerID"]].strip(),
            "smiles": cells[index["Ligand SMILES"]].strip(),
            "article_doi": cells[doi_index].strip() if doi_index is not None
            and doi_index < len(cells) else "",
            "Ki_nM_raw": cells[index["Ki (nM)"]].strip(),
            "Kd_nM_raw": cells[index["Kd (nM)"]].strip(),
            "IC50_nM_raw": cells[index["IC50 (nM)"]].strip(),
            "article_line": number,
        })
    return out, dict(tally)


def load_reactant_assay_mapping(path: Path, ids: set[str]) -> tuple[dict[str, set[str]], dict]:
    rows = tsv_rows(path, MAP_FILE)
    cols = require_columns([x.strip() for x in next(rows)],
                           {"REACTANT_SET_ID", "ENTRYID_ASSAYID"}, MAP_FILE)
    mapping: dict[str, set[str]] = defaultdict(set)
    stats = Counter()
    for cells in rows:
        stats["mapping_rows_scanned"] += 1
        if len(cells) <= max(cols.values()):
            stats["malformed_mapping_rows"] += 1
            continue
        rsid = cells[cols["REACTANT_SET_ID"]].strip()
        if rsid not in ids:
            continue
        key = cells[cols["ENTRYID_ASSAYID"]].strip()
        if not ASSAY_KEY_RE.fullmatch(key):
            stats["invalid_assay_key_rows"] += 1
            mapping[rsid].add("INVALID:" + key)
        else:
            mapping[rsid].add(key)
            stats["matched_mapping_rows"] += 1
    return dict(mapping), dict(stats)


def load_assays(path: Path, keys: set[str]) -> tuple[dict[str, dict], set[str], dict]:
    rows = tsv_rows(path, ASSAYS_FILE)
    cols = require_columns([x.strip() for x in next(rows)],
                           {"ENTRYID", "ASSAYID", "ASSAY_NAME", "DESCRIPTION"}, ASSAYS_FILE)
    assays = {}
    conflicts = set()
    stats = Counter()
    for cells in rows:
        stats["assays_rows_scanned"] += 1
        if len(cells) <= max(cols.values()):
            stats["malformed_assay_rows"] += 1
            continue
        key = cells[cols["ENTRYID"]].strip() + "_" + cells[cols["ASSAYID"]].strip()
        if key not in keys:
            continue
        info = {
            "assay_name": cells[cols["ASSAY_NAME"]].strip(),
            "assay_description": cells[cols["DESCRIPTION"]].strip(),
        }
        if key in assays and assays[key] != info:
            conflicts.add(key)
        else:
            assays[key] = info
    return assays, conflicts, dict(stats)


def link(articles: Path, mapping: Path, assays: Path, panel_path: Path,
         expected: dict[str, str]) -> tuple[dict, list[dict], list[dict]]:
    inputs = {"articles": articles, "mapping": mapping, "assays": assays}
    file_hashes = {}
    for name, path in inputs.items():
        if not path.is_file():
            raise FileNotFoundError(path)
        actual = checksum(path)
        if name not in expected or actual.lower() != expected[name].lower():
            raise ValueError(f"Pinned {name} archive SHA-256 mismatch or absent")
        file_hashes[name] = {"path": str(path.resolve()), "sha256": actual}
    panel = pd.read_csv(panel_path, sep="\t", dtype=str, keep_default_na=False)
    required = {"scope", "target_gene", "uniprot_accession"}
    if not required.issubset(panel.columns) or len(panel) != 20:
        raise ValueError("Expected exact 20-protein pilot panel")
    if panel["uniprot_accession"].duplicated().any():
        raise ValueError("Duplicate panel UniProt accessions")
    panel_map = {r["uniprot_accession"]: r for r in panel.to_dict("records")}
    measurements, article_stats = extract_articles(articles, panel_map)
    ids = {m["reactant_set_id"] for m in measurements if m["reactant_set_id"]}
    mapping_by_id, map_stats = load_reactant_assay_mapping(mapping, ids)
    keys = {key for assoc in mapping_by_id.values() for key in assoc
            if ASSAY_KEY_RE.fullmatch(key)}
    assay_by_key, conflicts, assay_stats = load_assays(assays, keys)
    rsid_counts = Counter(row["reactant_set_id"] for row in measurements
                          if row["reactant_set_id"])
    status_counts = Counter()
    by_protein = defaultdict(Counter)
    for row in measurements:
        rsid = row["reactant_set_id"]
        possible = mapping_by_id.get(rsid, set())
        row.update({
            "entryid_assayid": "",
            "assay_name": "",
            "assay_description": "",
            "link_status": "",
            "review_decision": "",
            "reviewer": "",
            "review_notes": "",
        })
        if not rsid:
            status = "MISSING_REACTANT_SET_ID"
        elif rsid_counts[rsid] > 1:
            status = "DUPLICATE_REACTANT_SET_ID_IN_ARTICLES"
        elif not possible:
            status = "NO_ASSAY_MAPPING"
        elif any(x.startswith("INVALID:") for x in possible):
            status = "INVALID_ASSAY_MAPPING"
        elif len(possible) > 1:
            status = "MULTIPLE_ASSAY_MAPPINGS"
        else:
            key = next(iter(possible))
            row["entryid_assayid"] = key
            if key in conflicts:
                status = "CONFLICTING_ASSAY_METADATA"
            elif key not in assay_by_key:
                status = "NO_ASSAY_METADATA"
            else:
                metadata = assay_by_key[key]
                row["assay_name"] = metadata["assay_name"]
                row["assay_description"] = metadata["assay_description"]
                status = ("LINKED_DESCRIPTION_PRESENT_UNREVIEWED"
                          if row["assay_description"] else "ASSAY_DESCRIPTION_MISSING")
        row["link_status"] = status
        status_counts[status] += 1
        by_protein[row["target_gene"]]["raw_rows"] += 1
        by_protein[row["target_gene"]][status] += 1
    breakdown = []
    for gene in panel["target_gene"]:
        stats = by_protein[gene]
        breakdown.append({
            "target_gene": gene,
            "raw_measurements": stats["raw_rows"],
            "with_assay_description": stats["LINKED_DESCRIPTION_PRESENT_UNREVIEWED"],
            "missing_or_ambiguous_link": (
                stats["raw_rows"] - stats["LINKED_DESCRIPTION_PRESENT_UNREVIEWED"]
            ),
            "approved_measurements": 0,
        })
    report = {
        "contract": "mcl-p036-bindingdb-assay-links-v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "release_declared": "202611",
        "panel_sha256": checksum(panel_path),
        "source_archives": file_hashes,
        "source_scope": "BindingDB_BindingDB_Articles (not full BindingDB)",
        "article_counts": article_stats,
        "mapping_counts": map_stats,
        "assay_counts": assay_stats,
        "link_counts": dict(status_counts),
        "distinct_assay_ids_requested": len(keys),
        "distinct_assay_ids_with_metadata": len(assay_by_key),
        "conflicting_assay_keys": sorted(conflicts),
        "scientific_status": "NO_ASSAYS_APPROVED_PENDING_MANUAL_REVIEW",
        "limitations": [
            "An exact entry/assay join establishes metadata provenance only.",
            "No assessment of assay specificity, target isolation, units, numerical validity or independently duplicated results.",
            "Multichain rows are excluded from single-protein endpoint counts.",
            "Rows with incomplete or ambiguous joins remain pending, not filled with guessed metadata.",
            "Manual review fields are empty and cannot be treated as approval.",
        ],
    }
    return report, measurements, breakdown


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--articles-zip", type=Path, required=True)
    parser.add_argument("--mapping-zip", type=Path, required=True)
    parser.add_argument("--assays-zip", type=Path, required=True)
    parser.add_argument("--pilot-panel", type=Path, required=True)
    parser.add_argument("--articles-sha256", required=True)
    parser.add_argument("--mapping-sha256", required=True)
    parser.add_argument("--assays-sha256", required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    out = args.out_dir.resolve()
    if out.exists():
        raise FileExistsError(f"Immutable assay linkage output exists: {out}")
    for path in (args.articles_zip, args.mapping_zip, args.assays_zip, args.pilot_panel):
        if out.is_relative_to(path.resolve().parent):
            raise ValueError(f"Refusing to write inside source directory: {path.parent}")
    report, records, summary = link(
        args.articles_zip, args.mapping_zip, args.assays_zip, args.pilot_panel,
        {"articles": args.articles_sha256, "mapping": args.mapping_sha256,
         "assays": args.assays_sha256})
    out.mkdir(parents=True)
    fields = [
        "scope", "target_gene", "uniprot_accession", "reactant_set_id",
        "monomer_id", "smiles", "article_doi", "Ki_nM_raw", "Kd_nM_raw",
        "IC50_nM_raw", "article_line", "entryid_assayid", "assay_name",
        "assay_description", "link_status", "review_decision", "reviewer",
        "review_notes",
    ]
    write_tsv(out / "assay_review_queue.tsv", records, fields)
    write_tsv(out / "assay_link_coverage.tsv", summary,
              ["target_gene", "raw_measurements", "with_assay_description",
               "missing_or_ambiguous_link", "approved_measurements"])
    (out / "assay_link_report.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"Articles endpoints for pilot (raw): {len(records)}")
    print(f"Distinct requested assay keys: {len({r['entryid_assayid'] for r in records if r['entryid_assayid']})}")
    for name, amount in sorted(report["link_counts"].items()):
        print(f"  {name}: {amount}")
    print("Scientific status: NO ASSAYS APPROVED / MANUAL SCIENTIFIC REVIEW REQUIRED")
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
