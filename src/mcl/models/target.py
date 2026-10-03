from __future__ import annotations

from enum import StrEnum
from pydantic import BaseModel, ConfigDict, Field


class MappingStatus(StrEnum):
    EXACT_APPROVED_SYMBOL = "exact_approved_symbol"
    PREVIOUS_SYMBOL = "previous_symbol"
    ALIAS_SYMBOL = "alias_symbol"
    WITHDRAWN_MERGED = "withdrawn_merged"
    MANUAL_REVIEW_REQUIRED = "manual_review_required"
    NOT_FOUND = "not_found"


class TargetSeed(BaseModel):
    model_config = ConfigDict(extra="forbid")

    target_id: str = Field(alias="Target_ID")
    input_symbol: str
    existing_uniprot_id: str | None = None
    existing_ensembl_gene_id: str | None = None
    existing_ncbi_gene_id: str | None = None


class TargetIdentity(BaseModel):
    model_config = ConfigDict(extra="forbid")

    target_id: str
    input_symbol: str
    mapping_status: MappingStatus
    mapping_basis: str
    hgnc_id: str | None = None
    hgnc_symbol: str | None = None
    hgnc_approved_name: str | None = None
    hgnc_status: str | None = None
    hgnc_date_modified: str | None = None
    hgnc_previous_symbols: list[str] = Field(default_factory=list)
    hgnc_aliases: list[str] = Field(default_factory=list)
    ensembl_gene_id: str | None = None
    ncbi_gene_id: str | None = None
    uniprot_id: str | None = None
    raw_record_hash: str | None = None
    source_snapshot: str | None = None
    notes: str | None = None
