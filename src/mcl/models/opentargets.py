from __future__ import annotations

from enum import StrEnum
from pydantic import BaseModel, ConfigDict, Field


class DiseaseMappingScope(StrEnum):
    EXACT_DISEASE = "exact_disease"
    BROADER_DISEASE_PROXY = "broader_disease_proxy"


class MappingConfidence(StrEnum):
    HIGH = "high"
    MEDIUM = "medium"
    LOW = "low"
    MANUAL_REVIEW_REQUIRED = "manual_review_required"


class OpenTargetsDiseaseMapping(BaseModel):
    model_config = ConfigDict(extra="forbid")

    cancer_id: str
    cancer_name: str
    molecular_context: str
    disease_id: str
    expected_disease_name: str
    mapping_scope: DiseaseMappingScope
    mapping_confidence: MappingConfidence
    molecular_context_encoded_in_ot: bool = False
    rationale: str
    mapping_source_url: str


class CancerTargetPairSeed(BaseModel):
    model_config = ConfigDict(extra="forbid")

    evidence_id: str
    cancer_id: str
    target_id: str
    cancer_name: str
    molecular_context: str
    hgnc_symbol: str


class OpenTargetsAssociation(BaseModel):
    model_config = ConfigDict(extra="forbid")

    evidence_id: str
    cancer_id: str
    target_id: str
    target_ensembl_id: str
    target_symbol: str
    disease_id: str
    disease_name: str
    expected_disease_name: str
    molecular_context: str
    mapping_scope: DiseaseMappingScope
    mapping_confidence: MappingConfidence
    molecular_context_encoded_in_ot: bool
    association_found: bool
    association_score: float | None = None
    direct_association: bool = True
    genetic_association_score: float | None = None
    somatic_mutation_score: float | None = None
    known_drug_score: float | None = None
    literature_score: float | None = None
    datatype_scores_json: str = "{}"
    datasource_scores_json: str = "{}"
    tractability_json: str = "[]"
    source_release: str
    api_endpoint: str
    retrieved_at: str
    raw_association_file: str
    raw_tractability_file: str | None = None
    raw_record_hash: str
    notes: str | None = None


class OpenTargetsDiseaseValidation(BaseModel):
    model_config = ConfigDict(extra="forbid")

    cancer_id: str
    configured_disease_id: str
    expected_disease_name: str
    returned_disease_id: str | None = None
    returned_disease_name: str | None = None
    id_match: bool = False
    name_match: bool = False
    mapping_scope: DiseaseMappingScope
    mapping_confidence: MappingConfidence
    status: str
    raw_file: str | None = None
    notes: str | None = None
