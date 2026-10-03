from __future__ import annotations

from pydantic import BaseModel, ConfigDict


class ProvenanceRecord(BaseModel):
    model_config = ConfigDict(extra="forbid")

    entity_type: str
    entity_id: str
    field_name: str
    value: str | None
    source_name: str
    source_record_id: str | None
    source_release: str | None
    source_url_or_endpoint: str | None
    retrieved_at: str
    evidence_nature: str
    parser_version: str
    raw_file: str
    raw_record_hash: str
