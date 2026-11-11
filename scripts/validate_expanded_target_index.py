#!/usr/bin/env python3
"""Independent P0.3.5 index consistency check. Read-only. Never activates MCL."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

import pyarrow.parquet as pq
from rdkit import Chem, DataStructs

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from mcl.analysis.expanded_ligand_chemistry import (  # noqa: E402
    FP_BITS,FP_BYTES,FP_RADIUS,FP_CHIRALITY,GENERATOR,characterize_smiles
)


def digest(path):
    h=hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda:stream.read(1024*1024),b""):
            h.update(block)
    return h.hexdigest()


def validate(stage: Path) -> dict:
    manifest_file=stage/"build_manifest.json"
    if not manifest_file.is_file():
        raise FileNotFoundError(manifest_file)
    manifest=json.loads(manifest_file.read_text(encoding="utf-8"))
    if manifest.get("contract")!="mcl-tls-v2-reviewed-ligand-index-staging-v1":
        raise ValueError("Unknown index manifest contract")
    if manifest.get("status")!="staging_not_published":
        raise ValueError("Index is not in pinned review stage")
    for artifact in manifest.get("artifacts",[]):
        file=stage/artifact["path"]
        if not file.is_file() or digest(file)!=artifact["sha256"]:
            raise ValueError(f"Index artifact missing or tampered: {file}")
    params=manifest["fingerprint"]
    if (params["radius"],params["bits"],params["use_chirality"],params["bytes_per_entry"]) != (
        FP_RADIUS,FP_BITS,FP_CHIRALITY,FP_BYTES
    ):
        raise ValueError("Fingerprint parameter mismatch")
    n=0
    path=stage/"fp_index"/"entries.tsv"
    data_path=stage/"fp_index"/"fingerprints.bin"
    total_bytes=data_path.stat().st_size
    catalog=pq.ParquetFile(stage/"ligand_catalog.parquet")
    with path.open(encoding="utf-8",newline="") as index, data_path.open("rb") as data:
        index_rows=csv.DictReader(index,delimiter="\t")
        def catalogue_rows():
            for batch in catalog.iter_batches(batch_size=2048):
                yield from batch.to_pylist()
        for row, entry in zip(catalogue_rows(), index_rows, strict=True):
            ligand_id=row["ligand_id"]
            offset=int(entry["byte_offset"])
            if row["target_gene"]!=entry["target_gene"] or ligand_id!=entry["ligand_id"]:
                raise ValueError(f"Catalogue/index ID mismatch: {ligand_id}")
            if row["exact_structure_id"]!=entry["exact_structure_id"]:
                raise ValueError(f"Exact chemical identity mismatch: {ligand_id}")
            if offset!=n*FP_BYTES or int(entry["fingerprint_bytes"])!=FP_BYTES:
                raise ValueError(f"Invalid fingerprint offset: {ligand_id}")
            if offset+FP_BYTES>total_bytes:
                raise ValueError(f"Out-of-range fingerprint: {ligand_id}")
            chem, expected=characterize_smiles(row["canonical_smiles"])
            if chem["exact_structure_id"]!=row["exact_structure_id"]:
                raise ValueError(f"Full chemical structure does not reproduce: {ligand_id}")
            if chem["parent_connectivity_id"]!=row["parent_connectivity_id"]:
                raise ValueError(f"Parent chemical structure mismatch: {ligand_id}")
            if chem["chemistry_status"]!="verified_pilot_scope" or expected is None:
                raise ValueError(f"Invalid catalogue chemistry: {ligand_id}")
            data.seek(offset)
            if data.read(FP_BYTES)!=expected:
                raise ValueError(f"Morgan fingerprint mismatch: {ligand_id}")
            n+=1
    if total_bytes!=n*FP_BYTES or n!=manifest["unique_target_structure_pairs_n"]:
        raise ValueError("Fingerprint index count or size mismatch")
    if pq.read_metadata(stage/"ligand_measurement_links.parquet").num_rows!=manifest["verified_measurements_n"]:
        raise ValueError("Measurement-to-ligand links count mismatch")
    return {
        "status":"verified_internal_consistency",
        "ligand_entries_recalculated":n,
        "experiment_rows_linked":manifest["verified_measurements_n"],
        "source_release_id":manifest["build_id"],
        "scientific_manual_review_complete":False,
        "note":"This validates file and chemical integrity, NOT biological efficacy or final publication.",
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage",type=Path,required=True)
    args=parser.parse_args()
    result=validate(args.stage)
    print(json.dumps(result,indent=2,ensure_ascii=False))


if __name__=="__main__":
    main()
