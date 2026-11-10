#!/usr/bin/env python3
"""P0.3.5: stage a provenance-linked, assay-reviewed ligand catalogue and Morgan index.

This command builds a local STAGING release only; P0.3.6 must approve and
activate it. Never writes to v1, current MCL API, or phenotype outputs.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import sqlite3
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
from rdkit import rdBase

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from mcl.analysis.expanded_ligand_chemistry import (  # noqa: E402
    FP_BITS, FP_BYTES, FP_CHIRALITY, FP_RADIUS, characterize_smiles,
)

CHEM_FIELDS = (
    "exact_structure_id", "parent_connectivity_id", "canonical_smiles",
    "parent_smiles", "full_inchikey", "parent_inchikey", "chemistry_status",
    "chemistry_reasons", "rdkit_version",
)
EXTRA_FIELDS = (
    "curated_quality_tier", "curation_reason_codes", "review_method",
    "duplicate_status", "independence_status", "ligand_id",
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def reviewed_assays(path: Path) -> dict[tuple[str, str, str], dict]:
    if not path.is_file():
        raise FileNotFoundError(f"Assay-review source missing: {path}")
    required = {
        "source", "source_release", "assay_id", "method_class",
        "reviewer", "evidence_url", "verification_status",
    }
    approved = {}
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if not required.issubset(set(reader.fieldnames or [])):
            raise ValueError(f"Missing assay review columns: {sorted(required-set(reader.fieldnames or []))}")
        for line, row in enumerate(reader, 2):
            if None in row:
                raise ValueError(f"Malformed review TSV at line {line}")
            if row["verification_status"].strip().lower() != "approved":
                continue
            if row["method_class"] not in {"binding_direct", "biochemical_function"}:
                raise ValueError(f"Invalid method_class at line {line}")
            if not row["assay_id"].strip() or not row["reviewer"].strip() or not row["evidence_url"].startswith(("https://", "http://")):
                raise ValueError(f"Unverifiable assay approval at line {line}")
            key = tuple(row[k].strip() for k in ("source", "source_release", "assay_id"))
            if key in approved:
                raise ValueError(f"Duplicate reviewed assay: {key}")
            approved[key] = row
    return approved


def _setup(connection: sqlite3.Connection) -> None:
    connection.executescript("""
        PRAGMA journal_mode=DELETE;
        CREATE TABLE measurements (
            seq INTEGER PRIMARY KEY,
            source TEXT NOT NULL,
            source_id TEXT NOT NULL,
            endpoint TEXT NOT NULL,
            target_gene TEXT NOT NULL,
            exact_id TEXT NOT NULL,
            origin_activity_id TEXT,
            doi TEXT,
            relation TEXT NOT NULL,
            nm TEXT NOT NULL,
            tier TEXT NOT NULL,
            dup TEXT NOT NULL DEFAULT 'unmatched',
            independent TEXT NOT NULL DEFAULT 'pending_review',
            fingerprint BLOB,
            payload TEXT NOT NULL
        );
        CREATE INDEX idx_source ON measurements (source,source_id,endpoint);
        CREATE INDEX idx_mirror ON measurements (origin_activity_id,endpoint);
        CREATE INDEX idx_doi ON measurements (doi,target_gene,exact_id,endpoint,relation,nm);
        CREATE INDEX idx_catalog ON measurements (tier,target_gene,exact_id);
    """)


def _mark_group(con: sqlite3.Connection, indices: list[int], *,
                label: str, consistent: bool, canonical: int | None = None) -> None:
    for seq in indices:
        if not consistent:
            con.execute("UPDATE measurements SET tier='C', dup=?, independent='review' WHERE seq=?",
                        (label + "_conflict", seq))
        elif seq == canonical:
            con.execute("UPDATE measurements SET dup=?, independent='primary' WHERE seq=?",
                        (label, seq))
        else:
            con.execute("UPDATE measurements SET tier='C', dup=?, independent='mirror' WHERE seq=?",
                        (label, seq))


def _duplicate_control(con: sqlite3.Connection) -> None:
    # A provider/source activity ID identifies one endpoint; corrupted conflicting
    # duplicates quarantine ALL members, rather than silently keeping one.
    groups = con.execute("""
        SELECT source,source_id,endpoint
        FROM measurements
        GROUP BY source,source_id,endpoint HAVING COUNT(*)>1
    """).fetchall()
    for source, sid, endpoint in groups:
        rows = con.execute("""
            SELECT seq,target_gene,exact_id,relation,nm
            FROM measurements WHERE source=? AND source_id=? AND endpoint=?
            ORDER BY seq
        """, (source, sid, endpoint)).fetchall()
        signature = {(r[1],r[2],r[3],r[4]) for r in rows}
        _mark_group(con, [r[0] for r in rows],
                    label="same_source_record", consistent=len(signature)==1,
                    canonical=rows[0][0])

    # Explicit original ChEMBL Activity ID is the only automatic source
    # mirror linkage. DOIs alone do not prove the same experimental run.
    groups = con.execute("""
        SELECT origin_activity_id,endpoint
        FROM measurements WHERE origin_activity_id!=''
        GROUP BY origin_activity_id,endpoint HAVING COUNT(DISTINCT source)>1
    """).fetchall()
    for activity, endpoint in groups:
        rows = con.execute("""
            SELECT seq,source,target_gene,exact_id,relation,nm
            FROM measurements WHERE origin_activity_id=? AND endpoint=?
            ORDER BY CASE source WHEN 'ChEMBL' THEN 0 ELSE 1 END,seq
        """, (activity, endpoint)).fetchall()
        signature = {(r[2],r[3],r[4],r[5]) for r in rows}
        _mark_group(con, [r[0] for r in rows],
                    label="confirmed_source_mirror", consistent=len(signature)==1,
                    canonical=rows[0][0])

    # Assign provisional independent status BEFORE comparing publication-level
    # duplicates. Previously unclassified A/B rows would otherwise be skipped.
    con.execute("""
        UPDATE measurements SET independent='primary'
        WHERE tier IN ('A','B') AND independent='pending_review'
    """)

    # DOI+same reported value is only a possible duplicate. Mark for review,
    # even if both rows were provisionally eligible for A/B.
    groups = con.execute("""
        SELECT doi,target_gene,exact_id,endpoint,relation,nm
        FROM measurements
        WHERE doi!='' AND tier IN ('A','B') AND independent='primary'
        GROUP BY doi,target_gene,exact_id,endpoint,relation,nm HAVING COUNT(*)>1
    """).fetchall()
    for signature in groups:
        rows = con.execute("""
            SELECT seq FROM measurements
            WHERE doi=? AND target_gene=? AND exact_id=? AND endpoint=?
                  AND relation=? AND nm=? AND tier IN ('A','B')
                  AND independent='primary'
        """, signature).fetchall()
        for (seq,) in rows:
            con.execute("""
                UPDATE measurements SET tier='C',
                 dup='probable_publication_duplicate_requires_review',
                 independent='review' WHERE seq=?
            """, (seq,))
    con.execute("""
        UPDATE measurements SET independent='primary'
        WHERE tier IN ('A','B') AND independent='pending_review'
    """)


def _rows(con: sqlite3.Connection):
    for seq, tier, status, independent, payload in con.execute("""
        SELECT seq,tier,dup,independent,payload
        FROM measurements ORDER BY seq
    """):
        record = json.loads(payload)
        record["curated_quality_tier"] = tier
        if status != "unmatched":
            reasons = json.loads(record.get("curation_reason_codes") or "[]")
            record["curation_reason_codes"] = json.dumps(sorted(set(reasons + [status])))
        record["duplicate_status"] = status
        record["independence_status"] = independent
        yield record


def _write_parquet_records(path: Path, records, fields: list[str]) -> int:
    schema = pa.schema([(field, pa.string()) for field in fields])
    writer = pq.ParquetWriter(path, schema, compression="zstd")
    batch, count = [], 0
    try:
        for row in records:
            batch.append({field: str(row.get(field) or "") for field in fields})
            count += 1
            if len(batch) >= 4096:
                writer.write_table(pa.Table.from_pylist(batch, schema=schema))
                batch.clear()
        if batch:
            writer.write_table(pa.Table.from_pylist(batch, schema=schema))
    finally:
        writer.close()
    return count


def build(stage_dir: Path, registry: Path, review_file: Path, output_dir: Path) -> dict:
    manifest_path = stage_dir / "ingest_manifest.json"
    normalized_path = stage_dir / "measurements_normalized.parquet"
    if not manifest_path.is_file() or not normalized_path.is_file():
        raise FileNotFoundError("P0.3 normalized stage and ingest_manifest.json are required")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("contract") != "mcl-tls-v2-normalized-staging-v1" or manifest.get("status") != "staging_not_published":
        raise ValueError("Invalid or published P0.3 input")
    if manifest["normalized_sha256"] != sha256(normalized_path):
        raise ValueError("Normalized Parquet changed after its pinned ingest build")
    if not registry.is_file() or manifest["registry_sha256"] != sha256(registry):
        raise ValueError("P0.2 protein registry mismatch")
    reviews = reviewed_assays(review_file)
    if not reviews:
        raise ValueError("No approved assay methods; cannot build a VERIFIED ligand index")
    build_key = "|".join((manifest["normalized_sha256"], manifest["registry_sha256"],
                          sha256(review_file), "radius2-2048-no-chirality",
                          rdBase.rdkitVersion))
    build_id = hashlib.sha256(build_key.encode("utf-8")).hexdigest()[:20]
    stage = output_dir / "staging" / build_id
    if stage.exists():
        raise FileExistsError(f"Immutable index staging release already exists: {stage}")
    stage.mkdir(parents=True)
    con = None
    try:
        con = sqlite3.connect(str(stage / "working.sqlite3"))
        _setup(con)
        source_rows = pq.ParquetFile(normalized_path)
        fields = list(source_rows.schema_arrow.names)
        missing = set(("source","source_release","source_record_id","endpoint",
                       "target_gene","original_smiles","assay_id","quality_tier",
                       "relation","numeric_value_nm")) - set(fields)
        if missing:
            raise ValueError(f"Normalized input lacks columns: {sorted(missing)}")
        inserted = 0
        for batch in source_rows.iter_batches(batch_size=2048):
            for row in batch.to_pylist():
                chemical, fp = characterize_smiles(row.get("original_smiles") or "")
                prelim = row.get("quality_tier") or "D"
                key = (row.get("source") or "", row.get("source_release") or "",
                       row.get("assay_id") or "")
                review = reviews.get(key)
                final = "C" if prelim not in {"D"} else "D"
                reason = []
                if chemical["chemistry_status"] != "verified_pilot_scope":
                    reason.append(chemical["chemistry_reasons"] or "unverified_chemistry")
                if not review:
                    reason.append("assay_method_not_approved")
                if prelim == "candidate_A" and review and review["method_class"] == "binding_direct":
                    final = "A"
                elif prelim == "candidate_B" and review and review["method_class"] == "biochemical_function":
                    final = "B"
                else:
                    if prelim in {"candidate_A", "candidate_B"} and review:
                        reason.append("assay_endpoint_method_discordance")
                if chemical["chemistry_status"] != "verified_pilot_scope":
                    final = "C" if row.get("target_gene") else "D"
                if row.get("measurement_status") != "exact_eligible":
                    final = "C" if row.get("target_gene") else "D"
                    reason.append("non_exact_or_outside_strict_cutoff")
                if not row.get("target_gene") or not row.get("source_record_id"):
                    final = "D"
                    reason.append("missing_provenance_or_target")
                row.update(chemical)
                row["curation_reason_codes"] = json.dumps(sorted(set(reason)))
                row["review_method"] = review["method_class"] if review else ""
                # Parent group is for navigation only; exact group is for measurements.
                acc = (row.get("source") or "", row.get("source_record_id") or "",
                       row.get("endpoint") or "")
                origin = (row.get("original_chembl_activity_id") or "") if acc[0] != "ChEMBL" else acc[1]
                con.execute("""
                    INSERT INTO measurements
                    (source,source_id,endpoint,target_gene,exact_id,origin_activity_id,
                     doi,relation,nm,tier,fingerprint,payload)
                    VALUES (?,?,?,?,?,?,?,?,?,?,?,?)
                """, (acc[0],acc[1],acc[2],row.get("target_gene") or "",
                      chemical["exact_structure_id"],origin,
                      (row.get("doi") or "").lower().strip(),
                      row.get("relation") or "",row.get("numeric_value_nm") or "",
                      final,fp,json.dumps(row,ensure_ascii=False)))
                inserted += 1
            con.commit()
        if inserted != int(manifest["measurement_rows"]):
            raise ValueError("Normalized measurement row count differs from pinned manifest")
        _duplicate_control(con)
        con.commit()

        # No A/B index without actually reviewed methods. C/D retained for QC.
        eligible = int(con.execute(
            "SELECT count(*) FROM measurements WHERE tier IN ('A','B') AND independent='primary'"
        ).fetchone()[0])
        if not eligible:
            raise ValueError("Zero verified independent A/B measurements: no index published")

        extra = list(CHEM_FIELDS) + list(EXTRA_FIELDS)
        count_curated = _write_parquet_records(
            stage / "measurements_curated.parquet", _rows(con),
            list(dict.fromkeys(fields + extra))
        )
        duplicate_records = (
            {"source": source, "source_record_id": sid, "endpoint": end,
             "duplicate_status": dup, "independence_status": ind}
            for source, sid, end, dup, ind in con.execute("""
                SELECT source,source_id,endpoint,dup,independent
                FROM measurements WHERE dup!='unmatched' ORDER BY seq
            """)
        )
        duplicate_n = _write_parquet_records(
            stage / "duplicate_clusters.parquet", duplicate_records,
            ["source","source_record_id","endpoint","duplicate_status","independence_status"]
        )

        # Stream catalogue rows; never materialize a million ligand records in RAM.
        catalog_fields = [
            "ligand_id", "target_gene", "uniprot_accession", "exact_structure_id",
            "parent_connectivity_id", "canonical_smiles", "parent_smiles",
            "full_inchikey", "independent_measurements_n",
            "independent_sources_n", "endpoint_types_json", "binding_A_n",
            "biochemical_B_n",
        ]
        catalogue_schema = pa.schema([(k, pa.string()) for k in catalog_fields])
        links_schema = pa.schema([(k, pa.string()) for k in (
            "ligand_id","source","source_record_id","endpoint","tier"
        )])
        catalog_writer = pq.ParquetWriter(
            stage/"ligand_catalog.parquet", catalogue_schema, compression="zstd"
        )
        link_writer = pq.ParquetWriter(
            stage/"ligand_measurement_links.parquet", links_schema, compression="zstd"
        )
        fp_dir = stage/"fp_index"
        fp_dir.mkdir()
        catalog_buffer = []
        catalog_count = 0
        target_counts = Counter()
        con.execute("CREATE TABLE ligand_ids (ligand_id TEXT PRIMARY KEY)")
        try:
            with (fp_dir/"fingerprints.bin").open("wb") as data, \
                 (fp_dir/"entries.tsv").open("w",encoding="utf-8",newline="") as listing:
                table = csv.writer(listing,delimiter="\t")
                table.writerow(["ligand_id","target_gene","exact_structure_id",
                                "byte_offset","fingerprint_bytes"])
                groups = con.execute("""
                    SELECT target_gene,exact_id,COUNT(*)
                    FROM measurements
                    WHERE tier IN ('A','B') AND independent='primary'
                    GROUP BY target_gene,exact_id ORDER BY target_gene,exact_id
                """)
                for gene,exact,n in groups:
                    all_rows = con.execute("""
                        SELECT source,source_id,endpoint,tier,fingerprint,payload
                        FROM measurements
                        WHERE target_gene=? AND exact_id=? AND tier IN ('A','B')
                              AND independent='primary'
                        ORDER BY source,source_id,endpoint
                    """,(gene,exact)).fetchall()
                    first=json.loads(all_rows[0][5])
                    fp=all_rows[0][4]
                    if fp is None or len(fp)!=FP_BYTES:
                        raise ValueError(f"Missing fingerprint {gene}/{exact}")
                    if any(r[4]!=fp for r in all_rows):
                        raise ValueError(f"Inconsistent fingerprints {gene}/{exact}")
                    ligand_id="TLS2-"+hashlib.sha256(
                        (gene+"|"+exact).encode("utf-8")
                    ).hexdigest()[:18]
                    con.execute("INSERT INTO ligand_ids VALUES (?)",(ligand_id,))
                    offset=data.tell()
                    data.write(fp)
                    table.writerow([ligand_id,gene,exact,offset,FP_BYTES])
                    catalog_count+=1
                    target_counts[gene]+=1
                    catalog_buffer.append({
                        "ligand_id":ligand_id,"target_gene":gene,
                        "uniprot_accession":first.get("uniprot_accession") or "",
                        "exact_structure_id":exact,
                        "parent_connectivity_id":first["parent_connectivity_id"],
                        "canonical_smiles":first["canonical_smiles"],
                        "parent_smiles":first["parent_smiles"],
                        "full_inchikey":first["full_inchikey"],
                        "independent_measurements_n":str(n),
                        "independent_sources_n":str(len({r[0] for r in all_rows})),
                        "endpoint_types_json":json.dumps(sorted({r[2] for r in all_rows})),
                        "binding_A_n":str(sum(r[3]=="A" for r in all_rows)),
                        "biochemical_B_n":str(sum(r[3]=="B" for r in all_rows)),
                    })
                    if len(catalog_buffer)>=4096:
                        catalog_writer.write_table(pa.Table.from_pylist(
                            catalog_buffer,schema=catalogue_schema
                        ))
                        catalog_buffer.clear()
                    # Link table preserves every independently curated measurement.
                    for pos in range(0,len(all_rows),4096):
                        block=all_rows[pos:pos+4096]
                        links=[
                            {"ligand_id":ligand_id,"source":r[0],
                             "source_record_id":r[1],"endpoint":r[2],"tier":r[3]}
                            for r in block
                        ]
                        link_writer.write_table(pa.Table.from_pylist(
                            links,schema=links_schema
                        ))
                if catalog_buffer:
                    catalog_writer.write_table(pa.Table.from_pylist(
                        catalog_buffer,schema=catalogue_schema
                    ))
        finally:
            link_writer.close()
            catalog_writer.close()
        with (stage/"target_coverage.tsv").open("w",encoding="utf-8",newline="") as stream:
            writer=csv.writer(stream,delimiter="\t")
            writer.writerow(["target_gene","verified_ligands_n"])
            writer.writerows(sorted(target_counts.items()))
        if (fp_dir/"fingerprints.bin").stat().st_size!=catalog_count*FP_BYTES:
            raise AssertionError("Fingerprint index size/count mismatch")
        if catalog_count!=pq.read_metadata(stage/"ligand_catalog.parquet").num_rows:
            raise AssertionError("Catalogue row count mismatch")

        counts = dict(con.execute(
            "SELECT tier,count(*) FROM measurements GROUP BY tier"
        ).fetchall())
        dup_counts = dict(con.execute(
            "SELECT dup,count(*) FROM measurements GROUP BY dup"
        ).fetchall())
        output_files = [
            stage/"measurements_curated.parquet",
            stage/"duplicate_clusters.parquet",
            stage/"ligand_measurement_links.parquet",
            stage/"ligand_catalog.parquet",
            fp_dir/"fingerprints.bin", fp_dir/"entries.tsv",
            stage/"target_coverage.tsv",
        ]
        report = {
            "contract": "mcl-tls-v2-reviewed-ligand-index-staging-v1",
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "status": "staging_not_published",
            "build_id": build_id,
            "input_ingest_sha256": sha256(manifest_path),
            "input_normalized_sha256": manifest["normalized_sha256"],
            "input_registry_sha256": manifest["registry_sha256"],
            "assay_review_sha256": sha256(review_file),
            "approved_assays_n": len(reviews),
            "measurement_rows_n": count_curated,
            "verified_measurements_n": eligible,
            "unique_target_structure_pairs_n": catalog_count,
            "proteins_with_verified_ligands_n": len(target_counts),
            "duplicate_records_flagged_n": duplicate_n,
            "quality_tier_counts": counts,
            "duplicate_status_counts": dup_counts,
            "rdkit_version": rdBase.rdkitVersion,
            "fingerprint": {
                "algorithm":"Morgan", "radius":FP_RADIUS, "bits":FP_BITS,
                "use_chirality":FP_CHIRALITY, "bytes_per_entry":FP_BYTES,
            },
            "artifacts": [
                {"path":str(p.relative_to(stage)), "sha256":sha256(p),
                 "bytes":p.stat().st_size} for p in output_files
            ],
            "limitations": [
                "A/B records require explicit documented assay review, not only ChEMBL confidence score.",
                "Unresolved stereochemistry, multi-fragment structures and possible duplicate publications are quarantined.",
                "Fingerprints ignore chirality by legacy design; exact structure IDs do not.",
                "This is staging, not active MCL search; P0.3.6 scientific QC and rollback approval remain mandatory.",
            ],
        }
        con.close()
        con = None
        (stage/"working.sqlite3").unlink()
        (stage/"build_manifest.json").write_text(
            json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8"
        )
        return report
    except BaseException:
        if con is not None:
            con.close()
        shutil.rmtree(stage, ignore_errors=True)
        raise


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ingest-stage", type=Path, required=True)
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--reviewed-assays", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path,
                        default=ROOT/"data/runtime/target_ligand_space_v2")
    args = parser.parse_args()
    report = build(args.ingest_stage, args.registry, args.reviewed_assays,
                   args.output_dir)
    print(f"Staged TLS v2 index: {report['build_id']}")
    print(f"Verified ligand-target structures: {report['unique_target_structure_pairs_n']}")
    print(f"Status: {report['status']}; API NOT switched")


if __name__ == "__main__":
    main()
