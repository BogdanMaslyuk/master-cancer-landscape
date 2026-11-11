#!/usr/bin/env python3
"""P0.3.6 exploratory BindingDB Articles coverage from pinned ZIP files.

This is a strictly observational, read-only subset audit. It does NOT
approve A/B biochemical measurements, deduplicate papers, assert target
engagement, infer PYZ efficacy or replace the FULL BindingDB release.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import zipfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

csv.field_size_limit(1 << 28)

EXPECTED_ARCHIVES = {
    "articles": "BindingDB_BindingDB_Articles.tsv",
    "assays": "BindingDB_Assays.tsv",
    "mapping": "BindingDB_rsid_eaids.tsv",
}
NEEDED = (
    "BindingDB Reactant_set_id",
    "BindingDB MonomerID",
    "Ligand SMILES",
    "UniProt (SwissProt) Primary ID of Target Chain",
    "Ki (nM)",
    "Kd (nM)",
    "IC50 (nM)",
)
ENDPOINTS = ("Ki", "Kd", "IC50")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def archive_info(path: Path, expected_member: str, expected_hash: str | None = None) -> dict:
    if not path.is_file():
        raise FileNotFoundError(path)
    checksum = sha256(path)
    if expected_hash and checksum.lower() != expected_hash.lower():
        raise ValueError(f"ZIP checksum mismatch: {path.name}")
    with zipfile.ZipFile(path) as archive:
        matches = [x for x in archive.infolist()
                   if Path(x.filename).name == expected_member and not x.is_dir()]
        if len(matches) != 1:
            raise ValueError(f"Expected exactly one {expected_member} in {path}")
        member = matches[0]
        if member.file_size <= 0:
            raise ValueError(f"Empty TSV member: {expected_member}")
        with archive.open(member) as stream:
            with io.TextIOWrapper(stream, encoding="utf-8-sig", newline="") as text:
                reader = csv.reader(text, delimiter="\t", quoting=csv.QUOTE_NONE)
                header = next(reader, [])
                first = next(reader, [])
        if not header or len(set(header)) != len(header):
            # BindingDB bulk format may include repeated chain-annotation columns,
            # so keep exact positions instead of guessing a remapped schema.
            if not header:
                raise ValueError(f"No TSV header: {expected_member}")
        return {
            "path": str(path.resolve()),
            "zip_bytes": path.stat().st_size,
            "sha256": checksum,
            "inner_filename": member.filename,
            "inner_bytes": member.file_size,
            "column_count": len(header),
            "first_row_column_count": len(first),
            "columns": header,
            "duplicate_column_labels": sorted(k for k, n in Counter(header).items() if n > 1),
            "zip_tested_header": True,
        }


def _panel(path: Path) -> tuple[list[dict], dict[str, dict]]:
    panel = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    required = {"target_gene", "uniprot_accession", "scope"}
    if not required.issubset(panel.columns):
        raise ValueError(f"Panel missing {sorted(required - set(panel.columns))}")
    if len(panel) < 18 or len(panel) > 23:
        raise ValueError("Expected fixed 13 core + 5..10 expanded panel")
    if int((panel["scope"] == "core").sum()) != 13:
        raise ValueError("Panel must have exactly 13 core targets")
    if not 5 <= int((panel["scope"] == "expanded").sum()) <= 10:
        raise ValueError("Panel must have 5..10 expanded targets")
    if panel["uniprot_accession"].duplicated().any() or panel["target_gene"].duplicated().any():
        raise ValueError("Pilot contains ambiguous gene/UniProt mappings")
    records = panel.to_dict("records")
    return records, {row["uniprot_accession"]: row for row in records}


def screen_articles(zip_path: Path, panel: dict[str, dict]) -> dict:
    counts = defaultdict(Counter)
    ligands = defaultdict(set)
    reaction_sets = defaultdict(set)
    references = defaultdict(set)
    total_rows = 0
    pilot_rows = 0
    malformed = 0
    missing_ids = 0
    found_multichain = 0
    header_labels = set()
    with zipfile.ZipFile(zip_path) as archive:
        matches = [x for x in archive.infolist() if Path(x.filename).name == EXPECTED_ARCHIVES["articles"]]
        if len(matches) != 1:
            raise ValueError("Missing or duplicated curated-articles TSV in ZIP")
        with archive.open(matches[0]) as stream:
            with io.TextIOWrapper(stream, encoding="utf-8-sig", newline="") as text:
                reader = csv.reader(text, delimiter="\t", quoting=csv.QUOTE_NONE)
                columns = next(reader)
                positions = {name: columns.index(name) for name in NEEDED if name in columns}
                absent = set(NEEDED) - set(positions)
                if absent:
                    raise ValueError(f"BindingDB articles missing columns: {sorted(absent)}")
                chain_column = next((columns.index(name) for name in columns
                                     if name.startswith("Number of Protein Chains in Target")), None)
                doi_column = columns.index("Article DOI") if "Article DOI" in columns else None
                header_labels = set(columns)
                max_index = max(positions.values())
                for cells in reader:
                    total_rows += 1
                    if len(cells) <= max_index:
                        malformed += 1
                        continue
                    accession = cells[positions["UniProt (SwissProt) Primary ID of Target Chain"]].strip()
                    if accession not in panel:
                        continue
                    pilot_rows += 1
                    counter = counts[accession]
                    counter["source_rows"] += 1
                    if chain_column is not None and len(cells) > chain_column:
                        chain = cells[chain_column].strip()
                        if chain and chain not in {"1", "1.0"}:
                            found_multichain += 1
                            counter["multichain_rows"] += 1
                    rsid = cells[positions["BindingDB Reactant_set_id"]].strip()
                    monomer = cells[positions["BindingDB MonomerID"]].strip()
                    if not rsid or not monomer:
                        missing_ids += 1
                        counter["missing_reactant_or_monomer_id"] += 1
                    else:
                        reaction_sets[accession].add(rsid)
                        ligands[accession].add(monomer)
                    if doi_column is not None and len(cells) > doi_column:
                        doi = cells[doi_column].strip()
                        if doi:
                            references[accession].add(doi.lower())
                    measured = False
                    for endpoint in ENDPOINTS:
                        raw = cells[positions[endpoint + " (nM)"]].strip()
                        if raw:
                            measured = True
                            counter[endpoint + "_nonempty_raw"] += 1
                            if raw.startswith(("<", ">", "≤", "≥")):
                                counter[endpoint + "_censored_raw"] += 1
                    if measured:
                        counter["rows_with_Ki_Kd_or_IC50"] += 1
    return {
        "total_article_rows_scanned": total_rows,
        "pilot_articles_rows": pilot_rows,
        "unparseable_short_rows": malformed,
        "pilot_rows_missing_reactant_or_monomer_id": missing_ids,
        "pilot_rows_with_multichain_annotation": found_multichain,
        "coverage": {
            accession: {
                "source_rows": counts[accession]["source_rows"],
                "raw_Ki_rows": counts[accession]["Ki_nonempty_raw"],
                "raw_Kd_rows": counts[accession]["Kd_nonempty_raw"],
                "raw_IC50_rows": counts[accession]["IC50_nonempty_raw"],
                "rows_with_endpoint": counts[accession]["rows_with_Ki_Kd_or_IC50"],
                "censored_raw_Ki": counts[accession]["Ki_censored_raw"],
                "censored_raw_Kd": counts[accession]["Kd_censored_raw"],
                "censored_raw_IC50": counts[accession]["IC50_censored_raw"],
                "multichain_rows": counts[accession]["multichain_rows"],
                "distinct_BindingDB_monomer_ids": len(ligands[accession]),
                "distinct_reaction_set_ids": len(reaction_sets[accession]),
                "distinct_DOIs_in_articles_file": len(references[accession]),
            }
            for accession in panel
        },
        "source_has_doi_column": "Article DOI" in header_labels,
        "note": "These are raw rows from an ARTICLES-ONLY subset; no experiment has been approved.",
    }


def audit(articles: Path, assays: Path, mapping: Path, panel_file: Path,
          expected_hashes: dict[str, str] | None = None) -> tuple[dict, list[dict]]:
    expected_hashes = expected_hashes or {}
    panel_rows, panel = _panel(panel_file)
    sources = {}
    for name, path in (("articles", articles), ("assays", assays), ("mapping", mapping)):
        sources[name] = archive_info(path, EXPECTED_ARCHIVES[name], expected_hashes.get(name))
    observed = screen_articles(articles, panel)
    records = []
    for row in panel_rows:
        info = observed["coverage"][row["uniprot_accession"]]
        records.append({
            "scope": row["scope"],
            "target_gene": row["target_gene"],
            "uniprot_accession": row["uniprot_accession"],
            **info,
            "experimental_quality_status": "UNREVIEWED_SUBSET",
        })
    result = {
        "contract": "mcl-p036-bindingdb-articles-preliminary-coverage-v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "release_declared": "202611",
        "source_scope": "BindingDB_BindingDB_Articles only; NOT BindingDB_All",
        "pilot_panel_sha256": sha256(panel_file),
        "sources": sources,
        "counts": {k: v for k, v in observed.items() if k != "coverage"},
        "total_pilot_proteins": len(records),
        "proteins_with_any_articles_endpoint_rows": sum(r["rows_with_endpoint"] > 0 for r in records),
        "scientific_status": "EXPLORATORY_NOT_A_B_APPROVED",
        "limitations": [
            "Only BindingDB-curated articles subset, not the full BindingDB database.",
            "Raw source records do not establish independent evidence, assay identity or valid single-protein binding.",
            "Missing a ligand in this subset never means no ligand exists in ChEMBL or elsewhere.",
            "The ZIP filenames and their checksums prove local file identity, not official authenticity.",
            "The assay and entry-assay mapping archives are inventoried, not joined at this step.",
        ],
    }
    return result, records


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--articles-zip", type=Path, required=True)
    parser.add_argument("--assays-zip", type=Path, required=True)
    parser.add_argument("--mapping-zip", type=Path, required=True)
    parser.add_argument("--pilot-panel", type=Path, required=True)
    parser.add_argument("--expected-articles-sha256")
    parser.add_argument("--expected-assays-sha256")
    parser.add_argument("--expected-mapping-sha256")
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.out_dir.resolve().is_relative_to(args.articles_zip.resolve().parent):
        raise SystemExit("Audit reports must not be written inside the immutable source archive folder.")
    if args.out_dir.resolve().is_relative_to(args.pilot_panel.resolve().parent):
        raise SystemExit("Audit reports must not overwrite the fixed pilot panel directory.")
    expected = {k: v for k, v in {
        "articles": args.expected_articles_sha256,
        "assays": args.expected_assays_sha256,
        "mapping": args.expected_mapping_sha256,
    }.items() if v}
    report, records = audit(args.articles_zip, args.assays_zip, args.mapping_zip,
                            args.pilot_panel, expected_hashes=expected)
    if args.out_dir.exists():
        raise FileExistsError("Immutable exploratory report exists; choose a new report directory")
    args.out_dir.mkdir(parents=True)
    with (args.out_dir / "pilot_bindingdb_articles_coverage.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(records)
    (args.out_dir / "bindingdb_articles_schema_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(f"BindingDB Articles scanned: {report['counts']['total_article_rows_scanned']}")
    print(f"Pilot-matched article rows: {report['counts']['pilot_articles_rows']}")
    print(f"Proteins with raw Ki/Kd/IC50: {report['proteins_with_any_articles_endpoint_rows']}/{len(records)}")
    for row in records:
        print(f"{row['scope']:8} {row['target_gene']:7} "
              f"rows={row['rows_with_endpoint']:6} "
              f"Ki={row['raw_Ki_rows']:5} Kd={row['raw_Kd_rows']:5} "
              f"IC50={row['raw_IC50_rows']:6}")
    print(f"Report: {args.out_dir / 'bindingdb_articles_schema_report.json'}")
    print("Status: ARTICLES-ONLY / UNREVIEWED / NOT A/B VERIFIED")


if __name__ == "__main__":
    main()
