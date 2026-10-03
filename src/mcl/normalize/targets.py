from __future__ import annotations

import csv
from datetime import datetime, timezone
from pathlib import Path

from mcl.models import ProvenanceRecord, TargetIdentity, TargetSeed
from mcl.sources.hgnc import HGNCLocalResolver, PARSER_VERSION
from mcl.utils.hash import sha256_file


def _none_if_blank(value: str | None) -> str | None:
    value = (value or "").strip()
    return value or None


def load_target_seeds(path: str | Path) -> list[TargetSeed]:
    seeds: list[TargetSeed] = []
    with Path(path).open("r", encoding="utf-8-sig", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            row = dict(row)
            for key in ("existing_uniprot_id", "existing_ensembl_gene_id", "existing_ncbi_gene_id"):
                row[key] = _none_if_blank(row.get(key))
            seeds.append(TargetSeed.model_validate(row))
    return seeds


def normalize_targets(
    seed_path: str | Path,
    hgnc_tsv: str | Path,
    withdrawn_tsv: str | Path | None = None,
    source_release: str | None = None,
) -> tuple[list[TargetIdentity], list[ProvenanceRecord]]:
    resolver = HGNCLocalResolver(hgnc_tsv, withdrawn_tsv)
    seeds = load_target_seeds(seed_path)
    identities = [resolver.normalize(seed) for seed in seeds]

    retrieved_at = datetime.now(timezone.utc).isoformat()
    raw_hash = sha256_file(hgnc_tsv)
    provenance: list[ProvenanceRecord] = []
    fields = [
        "hgnc_id", "hgnc_symbol", "hgnc_approved_name", "hgnc_status",
        "hgnc_date_modified", "ensembl_gene_id", "ncbi_gene_id", "uniprot_id",
        "mapping_status", "mapping_basis",
    ]
    for identity in identities:
        record_id = identity.hgnc_id
        for field in fields:
            value = getattr(identity, field)
            if hasattr(value, "value"):
                value = value.value
            provenance.append(
                ProvenanceRecord(
                    entity_type="target",
                    entity_id=identity.target_id,
                    field_name=field,
                    value=None if value is None else str(value),
                    source_name="HGNC",
                    source_record_id=record_id,
                    source_release=source_release,
                    source_url_or_endpoint="https://hgnc.genenames.org/download/",
                    retrieved_at=retrieved_at,
                    evidence_nature="Observed database annotation",
                    parser_version=PARSER_VERSION,
                    raw_file=str(Path(hgnc_tsv)),
                    raw_record_hash=identity.raw_record_hash or raw_hash,
                )
            )
    return identities, provenance
