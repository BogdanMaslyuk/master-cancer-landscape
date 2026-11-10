#!/usr/bin/env python3
"""P0.3.1/2: validate pinned ChEMBL and BindingDB extracts and stream measurements.

This script NEVER publishes a ligand catalog or fingerprint index. It expects
a verified, explicitly pre-joined ChEMBL export (not an unmodified SQL dump).
It does not download sources or guess release names.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from mcl.analysis.expanded_ligand_measurements import (  # noqa: E402
    PRIMARY_FIELDS, clean, normalize_bindingdb, normalize_chembl,
)

REQUIRED = {
    "ChEMBL": {
        "activity_id", "assay_chembl_id", "target_chembl_id",
        "uniprot_accession", "organism_id", "target_type",
        "confidence_score", "relationship_type", "assay_type", "bao_label",
        "standard_type", "standard_value", "standard_units",
        "standard_relation", "canonical_smiles", "molecule_chembl_id",
    },
    "BindingDB": {
        "BindingDB Reactant_set_id", "BindingDB MonomerID", "Ligand SMILES",
        "UniProt (SwissProt) Primary ID of Target Chain",
        "Target Source Organism According to Curator or DataSource",
        "Ki (nM)", "Kd (nM)", "IC50 (nM)",
    },
}


def sha256(path: Path) -> str:
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for part in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(part)
    return result.hexdigest()


def validate_manifest(manifest: dict, manifest_file: Path) -> list[dict]:
    if manifest.get("contract") != "mcl-tls-v2-pinned-extracts-v1":
        raise ValueError("Unrecognized source_manifest contract")
    sources = manifest.get("sources")
    if (not isinstance(sources, list) or len(sources) != len(REQUIRED)
            or not all(isinstance(s, dict) for s in sources)
            or sorted(s["source"] for s in sources) != sorted(REQUIRED)):
        raise ValueError("Exactly one ChEMBL and one BindingDB source are required")
    verified = []
    for entry in sources:
        source = entry["source"]
        if entry.get("format") != ("joined_activity_tsv" if source == "ChEMBL" else "bindingdb_bulk_tsv"):
            raise ValueError(f"Unsupported source extract format for {source}")
        if not all(clean(entry.get(k)) for k in ("release", "path", "sha256", "url", "license")):
            raise ValueError(f"Unpinned or unlicensed source: {source}")
        path = Path(entry["path"]).expanduser()
        if not path.is_absolute():
            path = manifest_file.parent / path
        path = path.resolve()
        if not path.is_file():
            raise FileNotFoundError(f"Source not available: {path}")
        expected = entry["sha256"].lower()
        if len(expected) != 64 or any(x not in "0123456789abcdef" for x in expected):
            raise ValueError(f"Invalid SHA-256 in manifest: {source}")
        actual = sha256(path)
        if actual != expected:
            raise ValueError(f"SHA-256 mismatch: {source}; no outputs published")
        verified.append({**entry, "resolved_path": path, "sha256": actual, "bytes": path.stat().st_size})
    return verified


def load_registry(path: Path) -> tuple[dict[str, str], str]:
    if not path.is_file():
        raise FileNotFoundError(f"Missing P0.2 registry: {path}")
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    required = {"target_gene", "uniprot_accession", "mapping_status", "organism_id"}
    if not required.issubset(frame.columns):
        raise ValueError(f"P0.2 registry columns missing: {sorted(required - set(frame.columns))}")
    valid = frame[
        frame["mapping_status"].eq("unique_reviewed_human")
        & frame["organism_id"].eq("9606")
    ]
    if valid.empty or valid["uniprot_accession"].duplicated().any() or valid["target_gene"].duplicated().any():
        raise ValueError("No unique reviewed-human P0.2 protein mapping")
    return dict(zip(valid["uniprot_accession"], valid["target_gene"])), sha256(path)


def stream_rows(source: dict):
    path: Path = source["resolved_path"]
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Missing TSV header: {path}")
        missing = REQUIRED[source["source"]] - set(reader.fieldnames)
        if missing:
            raise ValueError(f"{source['source']} TSV missing required columns: {sorted(missing)}")
        for row_number, row in enumerate(reader, start=2):
            if None in row:  # extra unexpected columns indicate malformed TSV
                raise ValueError(f"Malformed TSV extra columns in {path.name}:{row_number}")
            yield row_number, row


def build(manifest_path: Path, registry_path: Path, out_dir: Path, *, cutoff: str = "10000") -> dict:
    if float(cutoff) <= 0:
        raise ValueError("cutoff must be positive")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    verified = validate_manifest(manifest, manifest_path.resolve())
    registry, registry_sha = load_registry(registry_path)
    source_key = "|".join(
        f"{r['source']}:{r['release']}:{r['sha256']}" for r in sorted(verified, key=lambda a: a["source"])
    )
    build_id = hashlib.sha256(f"{source_key}|{registry_sha}|cutoff={cutoff}".encode()).hexdigest()[:20]
    staging = out_dir / "staging" / build_id
    if staging.exists():
        raise FileExistsError(f"Immutable staging run already exists: {staging}")
    staging.mkdir(parents=True)
    normalized = staging / "measurements_normalized.parquet"
    report_path = staging / "ingest_manifest.json"
    writer = None
    total = 0
    by_quality = {}
    by_source = {}
    try:
        schema = pa.schema([(name, pa.string()) for name in PRIMARY_FIELDS])
        writer = pq.ParquetWriter(normalized, schema, compression="zstd")
        batch = []
        for source in sorted(verified, key=lambda a: a["source"]):
            provider = source["source"]
            by_source[provider] = 0
            for line, row in stream_rows(source):
                items = []
                if provider == "ChEMBL":
                    items.append(normalize_chembl(
                        row, release=source["release"], sha=source["sha256"],
                        locator=f"{source['resolved_path'].name}:{line}",
                        registry=registry, cutoff_nm=cutoff,
                    ))
                else:
                    for endpoint in ("Ki", "Kd", "IC50"):
                        if clean(row.get(endpoint + " (nM)")):
                            items.append(normalize_bindingdb(
                                row, endpoint=endpoint, release=source["release"],
                                sha=source["sha256"],
                                locator=f"{source['resolved_path'].name}:{line}:{endpoint}",
                                registry=registry, cutoff_nm=cutoff,
                            ))
                for item in items:
                    batch.append(item)
                    total += 1
                    by_source[provider] += 1
                    label = item["quality_tier"]
                    by_quality[label] = by_quality.get(label, 0) + 1
                    if len(batch) >= 10000:
                        writer.write_table(pa.Table.from_pylist(batch, schema=schema))
                        batch.clear()
            if batch:
                writer.write_table(pa.Table.from_pylist(batch, schema=schema))
                batch.clear()
        writer.close()
        writer = None
        if not total:
            raise ValueError("No measurements found; refusing to publish empty staged output")
        report = {
            "contract": "mcl-tls-v2-normalized-staging-v1",
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "build_id": build_id,
            "status": "staging_not_published",
            "raw_sources": [{k: str(v) if isinstance(v, Path) else v for k, v in s.items()}
                            for s in verified],
            "registry_sha256": registry_sha,
            "cutoff_nm": cutoff,
            "measurement_rows": total,
            "source_rows_emitted": by_source,
            "candidate_quality_tiers": by_quality,
            "normalized_sha256": sha256(normalized),
            "limitations": [
                "Candidate quality tiers are provisional; structures and assay context require independent curation.",
                "No cross-source deduplication or Morgan search index has been built.",
                "No clinical/anti-cancer efficacy or activity of uploaded PYZ is inferred.",
            ],
        }
        report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
        return report
    except Exception:
        if writer is not None:
            writer.close()
        # A failed build leaves NO apparently valid reusable staging run.
        import shutil
        shutil.rmtree(staging, ignore_errors=True)
        raise


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-manifest", type=Path, required=True)
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, default=ROOT / "data/processed/target_ligand_space_v2")
    parser.add_argument("--cutoff-nm", default="10000")
    args = parser.parse_args()
    result = build(args.source_manifest, args.registry, args.out_dir, cutoff=args.cutoff_nm)
    print(f"P0.3 staged {result['measurement_rows']} rows (NOT PUBLISHED): {result['build_id']}")


if __name__ == "__main__":
    main()
