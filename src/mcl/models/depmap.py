from __future__ import annotations

from pydantic import BaseModel, ConfigDict


class DepMapEvidence(BaseModel):
    """Quantitative DepMap evidence for one Cancer_ID × Target_ID pair."""

    model_config = ConfigDict(extra="forbid")

    evidence_id: str
    cancer_id: str
    target_id: str
    target_symbol: str
    depmap_release: str

    context_definition: str
    comparator_definition: str
    context_models_n: int
    comparator_models_n: int
    context_gene_effect_n: int
    comparator_gene_effect_n: int
    context_dependency_probability_n: int
    comparator_dependency_probability_n: int

    dependent_cell_lines_n: int
    dependency_fraction: float | None
    median_gene_effect: float | None
    mean_gene_effect: float | None
    gene_effect_iqr: float | None
    median_dependency_probability: float | None
    min_gene_effect: float | None
    max_gene_effect: float | None
    missing_context_gene_effect_n: int

    comparator_dependent_cell_lines_n: int
    comparator_dependency_fraction: float | None
    comparator_median_gene_effect: float | None
    comparator_median_dependency_probability: float | None

    delta_gene_effect: float | None
    cliffs_delta: float | None
    statistical_test: str | None
    p_value: float | None
    q_value: float | None

    low_sample_size: bool
    statistical_test_performed: bool

    all_models_gene_effect_n: int
    all_models_dependency_probability_n: int
    all_models_median_gene_effect: float | None
    broad_dependency_fraction: float | None
    broad_dependency_warning: bool

    gene_effect_file: str
    dependency_probability_file: str
    model_file: str
    mutation_file: str
    omics_profiles_file: str
    retrieved_at: str
    notes: str


class DepMapContextAuditRow(BaseModel):
    """Auditable context/comparator assignment of one DepMap model."""

    model_config = ConfigDict(extra="forbid")

    cancer_id: str
    model_id: str
    cell_line_name: str | None = None
    depmap_model_type: str | None = None
    oncotree_lineage: str | None = None
    oncotree_primary_disease: str | None = None
    oncotree_subtype: str | None = None
    oncotree_code: str | None = None
    disease_metadata_match: bool
    sequencing_available: bool
    alteration_status: str
    assigned_group: str
    assignment_reason: str
    qualifying_variants_json: str
    other_relevant_variants_json: str


class DepMapKrasSensitivityEvidence(BaseModel):
    """Sensitivity analysis that separates KRAS WT proxy and other KRAS-driver comparators."""

    model_config = ConfigDict(extra="forbid")

    sensitivity_id: str
    evidence_id: str
    cancer_id: str
    target_id: str
    target_symbol: str
    depmap_release: str
    target_kras_variant: str
    comparison_type: str
    context_definition: str
    comparator_definition: str

    context_models_n: int
    comparator_models_n: int
    context_gene_effect_n: int
    comparator_gene_effect_n: int
    context_dependency_probability_n: int
    comparator_dependency_probability_n: int

    dependent_cell_lines_n: int
    dependency_fraction: float | None
    median_gene_effect: float | None
    median_dependency_probability: float | None
    comparator_dependent_cell_lines_n: int
    comparator_dependency_fraction: float | None
    comparator_median_gene_effect: float | None
    comparator_median_dependency_probability: float | None

    delta_gene_effect: float | None
    cliffs_delta: float | None
    statistical_test: str | None
    p_value: float | None
    q_value: float | None
    low_sample_size: bool
    statistical_test_performed: bool

    retrieved_at: str
    notes: str


class DepMapKrasSensitivityAuditRow(BaseModel):
    """Auditable KRAS status assignment for one same-disease DepMap model."""

    model_config = ConfigDict(extra="forbid")

    cancer_id: str
    target_kras_variant: str
    model_id: str
    cell_line_name: str | None = None
    depmap_model_type: str | None = None
    oncotree_lineage: str | None = None
    oncotree_primary_disease: str | None = None
    oncotree_subtype: str | None = None
    oncotree_code: str | None = None
    sequencing_available: bool
    kras_status: str
    assigned_group: str
    assignment_reason: str
    target_variants_json: str
    other_kras_variants_json: str
