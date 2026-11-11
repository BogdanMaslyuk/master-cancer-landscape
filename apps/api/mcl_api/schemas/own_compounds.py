from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, ConfigDict, Field

from .common import ApiResponseModel


class PyzSummaryResponse(ApiResponseModel):
    pass


class PyzCatalogResponse(ApiResponseModel):
    pass


class PyzMatrixResponse(ApiResponseModel):
    pass


class ScientificClaim(BaseModel):
    rule_id: str
    title: str
    meaning: str
    source: str


class PyzBlock(BaseModel):
    status: Literal["available", "partial", "not_available", "pending_review"]
    data_available: bool


class PyzIdentity(PyzBlock):
    own_compound_id: str = Field(pattern=r"^PYZ-\d{3}$")
    standardized_smiles: str | None
    chemotype_ru: str | None
    ring_variant: str | None
    series_range: str | None
    source: str


class PyzProperties(PyzBlock):
    qed: float | None = Field(default=None, ge=0, le=1)
    source: str | None = None


class PyzAdmet(PyzBlock):
    integrated_risk: str | None = None
    note_ru: str | None = None
    source_report: str | None = None
    prediction_not_experiment: Literal[True]
    individual_ames: str | None = None
    individual_dili: str | None = None


class PyzNeighbor(BaseModel):
    ligand_id: str | None = None
    tanimoto: float | None = Field(default=None, ge=0, le=1)
    reference_endpoint: str | None = None
    reference_value_nm: float | None = None


class PyzTarget(BaseModel):
    target_gene: str
    uniprot_accession: str | None = None
    tanimoto: float | None = Field(default=None, ge=0, le=1)
    passes_navigation_cutoff: bool
    ligands_compared_n: int | None = None
    best_ligand_id: str | None = None
    reference_value_nm: float | None = None
    reference_endpoint: str | None = None
    reference_measurement_owner: Literal["known_ligand_not_pyz"]
    evidence_type: Literal["computational_prediction"]
    scientific_status: Literal["exploratory_not_directly_validated"]
    source: str
    neighbors: list[PyzNeighbor] = Field(default_factory=list)


class PyzTargetHypotheses(PyzBlock):
    scope: Literal["existing_13_core_targets"]
    cutoff_navigational: float
    items: list[PyzTarget]
    source: str


class PyzContextModel(BaseModel):
    model_id: str | None = None
    cell_line_name: str | None = None
    role: str | None = None
    dependency_probability: float | None = None
    gene_effect: float | None = None
    available_in_laboratory_registry: bool
    known_reference_compounds: list[str] = Field(default_factory=list)


class PyzContext(BaseModel):
    target_gene: str | None = None
    mcl_organ_ru: str | None = None
    mcl_cancer_name: str | None = None
    mcl_cancer_id: str | None = None
    tanimoto: float | None = Field(default=None, ge=0, le=1)
    models: list[PyzContextModel] = Field(default_factory=list)
    context_status: Literal["known_target_or_reference_context_not_pyz"]
    source: str


class PyzCellContext(PyzBlock):
    items: list[PyzContext] = Field(default_factory=list)


class PyzScientificSummary(BaseModel):
    established: list[ScientificClaim]
    hypotheses: list[ScientificClaim]
    gaps: list[ScientificClaim]
    next_step: ScientificClaim
    rule_ids: list[str]
    own_cellular_activity: Literal["not_measured"]
    own_target_binding: Literal["not_measured"]


class PyzProvenance(BaseModel):
    own_compounds_registry: str
    admet_report: str | None = None
    target_library: str
    candidate_context: str
    bindingdb_p036_status: Literal["NOT_INGESTED_UNREVIEWED"]
    docking_status: str


class PyzProfileV2(BaseModel):
    """Strictly validated additive profile. Scientific source rows are unchanged."""

    model_config = ConfigDict(extra="forbid")
    schema_version: Literal["mcl-pyz-profile-v2-p0"]
    compound_id: str = Field(pattern=r"^PYZ-\d{3}$")
    identity: PyzIdentity
    properties: PyzProperties
    admet: PyzAdmet
    target_hypotheses: PyzTargetHypotheses
    cell_context: PyzCellContext
    scientific_summary: PyzScientificSummary
    provenance: PyzProvenance


class PyzDetailResponse(ApiResponseModel):
    """The existing response is preserved, with a strictly typed nested view."""

    profile_v2: PyzProfileV2
