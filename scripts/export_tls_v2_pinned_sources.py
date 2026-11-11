#!/usr/bin/env python3
"""Export a joined ChEMBL SQLite snapshot + pin an official BindingDB TSV.

Unlike a live REST query, an explicitly versioned SQLite download makes the
activity/assay/target joins reproducible. This command does not fetch data.
It requires an operator-supplied official snapshot and exact source metadata.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sqlite3
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

SQL = """
SELECT
  CAST(act.activity_id AS TEXT) AS activity_id,
  a.chembl_id AS assay_chembl_id,
  td.chembl_id AS target_chembl_id,
  seq.accession AS uniprot_accession,
  CAST(td.tax_id AS TEXT) AS organism_id,
  td.target_type AS target_type,
  CAST(a.confidence_score AS TEXT) AS confidence_score,
  a.relationship_type AS relationship_type,
  a.assay_type AS assay_type,
  bao.label AS bao_label,
  act.standard_type AS standard_type,
  CAST(act.standard_value AS TEXT) AS standard_value,
  act.standard_units AS standard_units,
  act.standard_relation AS standard_relation,
  cs.canonical_smiles AS canonical_smiles,
  md.chembl_id AS molecule_chembl_id,
  d.chembl_id AS document_chembl_id,
  d.doi AS doi,
  CAST(d.pubmed_id AS TEXT) AS pubmed_id,
  act.data_validity_comment AS data_validity_comment,
  CAST(act.potential_duplicate AS TEXT) AS potential_duplicate,
  a.description AS assay_description
FROM activities act
JOIN assays a ON a.assay_id=act.assay_id
JOIN target_dictionary td ON td.tid=a.tid
JOIN (
 SELECT tid,MAX(component_id) AS component_id
 FROM target_components GROUP BY tid HAVING COUNT(*)=1
) tc ON tc.tid=td.tid
JOIN component_sequences seq ON seq.component_id=tc.component_id
JOIN temp.pinned_target_accessions pinned ON pinned.accession=seq.accession
JOIN molecule_dictionary md ON md.molregno=act.molregno
JOIN compound_structures cs ON cs.molregno=md.molregno
LEFT JOIN docs d ON d.doc_id=a.doc_id
LEFT JOIN bioassay_ontology bao ON bao.bao_id=a.bao_format
WHERE td.tax_id=9606
  AND td.target_type='SINGLE PROTEIN'
  AND a.variant_id IS NULL
  AND a.confidence_score=9
  AND a.relationship_type='D'
  AND a.assay_type IN ('B','F')
  AND act.standard_type IN ('Ki','Kd','IC50')
  AND act.standard_value IS NOT NULL
  AND cs.canonical_smiles IS NOT NULL
ORDER BY act.activity_id
"""

REQUIRED_SCHEMA = {
    "activities": {"activity_id","assay_id","molregno","standard_type",
                   "standard_value","standard_units","standard_relation",
                   "data_validity_comment","potential_duplicate"},
    "assays": {"assay_id","chembl_id","tid","doc_id","confidence_score",
               "relationship_type","assay_type","bao_format","variant_id","description"},
    "target_dictionary": {"tid","chembl_id","tax_id","target_type"},
    "target_components": {"tid","component_id"},
    "component_sequences": {"component_id","accession"},
    "molecule_dictionary": {"molregno","chembl_id"},
    "compound_structures": {"molregno","canonical_smiles"},
    "bioassay_ontology": {"bao_id","label"},
    "docs": {"doc_id","chembl_id","doi","pubmed_id"},
}


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024*1024), b""):
            h.update(block)
    return h.hexdigest()


def _check_tables(con: sqlite3.Connection) -> None:
    for name, columns in REQUIRED_SCHEMA.items():
        existing = {row[1] for row in con.execute(f"PRAGMA table_info({name})")}
        missing = columns-existing
        if missing:
            raise ValueError(f"ChEMBL SQLite incompatible: {name} lacks {sorted(missing)}")


def export_sqlite(
    database: Path, registry: Path, output_tsv: Path, *,
    release: str, expected_db_sha256: str,
    pilot_panel: Path | None = None,
) -> dict:
    if not release or not expected_db_sha256 or len(expected_db_sha256) != 64:
        raise ValueError("Pinned ChEMBL release and SQLite SHA-256 are mandatory")
    if not database.is_file():
        raise FileNotFoundError(database)
    if digest(database) != expected_db_sha256.lower():
        raise ValueError("Pinned ChEMBL database SHA-256 mismatch")
    mapping = pd.read_csv(registry, sep="\t", dtype=str, keep_default_na=False)
    required = {"uniprot_accession","mapping_status","organism_id"}
    if not required.issubset(mapping.columns):
        raise ValueError("Missing P0.2 registry columns")
    accessions = sorted(set(
        mapping.loc[
            mapping["mapping_status"].eq("unique_reviewed_human") &
            mapping["organism_id"].eq("9606"), "uniprot_accession"
        ]
    ))
    if pilot_panel is not None:
        panel = pd.read_csv(pilot_panel, sep="\t", dtype=str, keep_default_na=False)
        required_panel = {"uniprot_accession", "target_gene", "scope"}
        if not required_panel.issubset(panel.columns):
            raise ValueError("Pilot panel missing gene/accession/scope")
        if panel["uniprot_accession"].duplicated().any() or panel["target_gene"].duplicated().any():
            raise ValueError("Ambiguous duplicate pilot panel")
        mapping = dict(zip(mapping["uniprot_accession"], mapping["target_gene"]))
        if any(mapping.get(row["uniprot_accession"]) != row["target_gene"]
               for row in panel.to_dict("records")):
            raise ValueError("Pilot panel contains unmapped/mismatched UniProt accession")
        selected = set(panel["uniprot_accession"])
        if not selected.issubset(set(accessions)):
            raise ValueError("Pilot panel contains proteins outside reviewed registry")
        accessions = sorted(selected)
    if not accessions:
        raise ValueError("No verified human target accessions in P0.2 registry")
    if output_tsv.exists():
        raise FileExistsError(f"Output cannot overwrite previous extract: {output_tsv}")

    # Open the physical database in read-only mode. Temporary in-memory joins
    # are allowed, but the archived source bytes cannot be changed.
    con = sqlite3.connect(f"file:{database.resolve().as_posix()}?mode=ro", uri=True)
    temp = output_tsv.with_name(output_tsv.name + ".part")
    try:
        _check_tables(con)
        version_exists = con.execute("""
            SELECT COUNT(*) FROM sqlite_master
            WHERE type='table' AND name='version'
        """).fetchone()[0]
        version_snapshot = [
            list(map(str,row)) for row in con.execute("SELECT * FROM version LIMIT 10")
        ] if version_exists else []
        if not version_snapshot or not any(
            release.lower().replace("chembl_", "") in " ".join(row).lower()
            for row in version_snapshot
        ):
            raise ValueError(
                f"Declared release {release} not corroborated by SQLite VERSION rows"
            )
        con.execute("CREATE TEMP TABLE pinned_target_accessions (accession TEXT PRIMARY KEY)")
        con.executemany("INSERT INTO pinned_target_accessions VALUES (?)",
                        [(acc,) for acc in accessions])
        output_tsv.parent.mkdir(parents=True, exist_ok=True)
        count = 0
        with temp.open("w", encoding="utf-8", newline="") as file:
            cursor = con.execute(SQL)
            names = [entry[0] for entry in cursor.description]
            writer = csv.writer(file, delimiter="\t", lineterminator="\n")
            writer.writerow(names)
            for row in cursor:
                writer.writerow(["" if value is None else value for value in row])
                count += 1
        if count == 0:
            raise ValueError("No ChEMBL rows matched the reviewed protein/assay constraints")
        temp.replace(output_tsv)
        return {
            "rows": count, "output": str(output_tsv),
            "pilot_panel_sha256": digest(pilot_panel) if pilot_panel is not None else None,
            "output_sha256": digest(output_tsv),
            "database_sha256": expected_db_sha256.lower(),
            "registry_sha256": digest(registry),
            "sql_sha256": hashlib.sha256(SQL.encode("utf-8")).hexdigest(),
            "version_rows": version_snapshot,
        }
    finally:
        con.close()
        if temp.exists():
            temp.unlink()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chembl-sqlite", type=Path, required=True)
    parser.add_argument("--chembl-sqlite-sha256", required=True)
    parser.add_argument("--chembl-release", required=True)
    parser.add_argument("--chembl-url", required=True)
    parser.add_argument("--bindingdb-tsv", type=Path, required=True)
    parser.add_argument("--bindingdb-release", required=True)
    parser.add_argument("--bindingdb-url", required=True)
    parser.add_argument("--bindingdb-license", required=True, help="Verified terms for this exact source release")
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--pilot-panel", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    if not args.bindingdb_tsv.is_file():
        raise FileNotFoundError(args.bindingdb_tsv)
    if (not args.bindingdb_license.strip() or any(
        token in args.bindingdb_license.upper()
        for token in ("VERIFY", "UNKNOWN", "TODO", "REPLACE", "TBD", "ACTUAL_", "INSERT_")
    )):
        raise ValueError("BindingDB license/attribution must be verified for the selected source")
    if not args.chembl_url.startswith("https://") or not args.bindingdb_url.startswith("https://"):
        raise ValueError("Use official pinned HTTPS source URLs")
    output = args.output_dir/"chembl_joined.tsv"
    result = export_sqlite(
        args.chembl_sqlite,args.registry,output,
        release=args.chembl_release,expected_db_sha256=args.chembl_sqlite_sha256,
        pilot_panel=args.pilot_panel,
    )
    bdb_sha = digest(args.bindingdb_tsv)
    manifest = {
        "contract": "mcl-tls-v2-pinned-extracts-v1",
        "created_at_utc":datetime.now(timezone.utc).isoformat(),
        "source_extract_provenance": result,
        "sources": [
            {
                "source":"ChEMBL", "release":args.chembl_release,
                "path":output.name,"sha256":result["output_sha256"],
                "url":args.chembl_url,"license":"CC BY-SA 3.0",
                "format":"joined_activity_tsv",
            },
            {
                "source":"BindingDB", "release":args.bindingdb_release,
                "path":str(args.bindingdb_tsv.resolve()),
                "sha256":bdb_sha,"url":args.bindingdb_url,
                "license":args.bindingdb_license.strip(),
                "format":"bindingdb_bulk_tsv",
            },
        ],
    }
    target = args.output_dir/"source_manifest.json"
    if target.exists():
        raise FileExistsError(f"Cannot overwrite pinned manifest: {target}")
    target.write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding="utf-8")
    print(f"Pinned ChEMBL extracted rows: {result['rows']}")
    print(f"Source manifest: {target}")
    print("Next: ingest_expanded_target_ligands.py, then manual assay review, then build_expanded_target_index.py")


if __name__ == "__main__":
    main()
