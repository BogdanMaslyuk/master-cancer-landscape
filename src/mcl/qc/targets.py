from __future__ import annotations

from collections import Counter

from mcl.models import MappingStatus, QCRecord, TargetIdentity, TargetSeed


def target_normalization_qc(
    seeds: list[TargetSeed], identities: list[TargetIdentity]
) -> list[QCRecord]:
    qc: list[QCRecord] = []

    target_ids = [x.target_id for x in identities]
    symbols = [x.hgnc_symbol for x in identities if x.hgnc_symbol]
    for value, count in Counter(target_ids).items():
        if count > 1:
            qc.append(QCRecord(severity="ERROR", check="unique_target_id", entity_id=value,
                               message=f"Duplicate Target_ID appears {count} times"))
    for value, count in Counter(symbols).items():
        if count > 1:
            qc.append(QCRecord(severity="ERROR", check="unique_hgnc_symbol", observed=value,
                               message=f"Canonical HGNC symbol appears {count} times"))

    seed_by_id = {s.target_id: s for s in seeds}
    for identity in identities:
        seed = seed_by_id[identity.target_id]
        if identity.mapping_status in {MappingStatus.NOT_FOUND, MappingStatus.MANUAL_REVIEW_REQUIRED}:
            qc.append(QCRecord(
                severity="ERROR",
                check="mapping_resolved",
                entity_id=identity.target_id,
                field="input_symbol",
                observed=identity.input_symbol,
                message=f"Mapping status: {identity.mapping_status.value}; {identity.notes or ''}".strip(),
            ))
            continue
        if identity.hgnc_status and identity.hgnc_status.lower() != "approved":
            qc.append(QCRecord(severity="ERROR", check="hgnc_status", entity_id=identity.target_id,
                               observed=identity.hgnc_status, expected="Approved",
                               message="Resolved HGNC record is not approved"))

        comparisons = [
            ("UniProt_ID", seed.existing_uniprot_id, identity.uniprot_id),
            ("Ensembl_gene_ID", seed.existing_ensembl_gene_id, identity.ensembl_gene_id),
            ("NCBI_Gene_ID", seed.existing_ncbi_gene_id, identity.ncbi_gene_id),
        ]
        for field, existing, normalized in comparisons:
            if existing and normalized and existing != normalized:
                qc.append(QCRecord(
                    severity="WARNING", check="existing_identifier_conflict",
                    entity_id=identity.target_id, field=field,
                    observed=existing, expected=normalized,
                    message="Existing workbook/seed identifier differs from current HGNC cross-reference; manual verification required.",
                ))

    if not qc:
        qc.append(QCRecord(severity="INFO", check="target_normalization", message="All checks passed"))
    return qc
