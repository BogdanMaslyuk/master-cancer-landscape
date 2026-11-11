// AUTO-GENERATED FILE. DO NOT EDIT MANUALLY.
// Source of truth: FastAPI OpenAPI schema (`mcl_api.main:app`).
// Regenerate with: .\.venv\Scripts\python.exe .\scripts\generate_frontend_api_types.py

export type AtlasOrgan = {
  [key: string]: unknown;
  "id": string;
  "name_ru": string | null;
  "name_en": string | null;
  "icon": string | null;
  "contexts": Record<string, unknown>[];
  "contexts_n": number;
  "models_n": number;
  "analyses_n": number;
};

export type AtlasResponse = {
  [key: string]: unknown;
  "organs_n": number;
  "contexts_n": number;
  "models_n": number;
  "analyses_n": number;
  "organs": AtlasOrgan[];
};

export type CRISPRAtlasResponse = Record<string, unknown>;

export type CRISPRModelResponse = Record<string, unknown>;

export type CRISPRModelsResponse = Record<string, unknown>;

export type CancerContextResponse = {
  [key: string]: unknown;
  "id": string;
  "name": string;
  "organ_id": string;
  "organ_ru": string | null;
  "cancer_ru": string | null;
  "molecular_context": string | null;
  "models_n": number;
  "context_models_n": number;
  "comparator_models_n": number;
  "comparisons_n": number;
};

export type CohortResponse = Record<string, unknown>;

export type ComparisonGenesResponse = {
  [key: string]: unknown;
  "comparison_id": string;
  "page": number;
  "page_size": number;
  "total": number;
  "items": Record<string, unknown>[];
};

export type ComparisonListResponse = ComparisonSummary[];

export type ComparisonSummary = {
  [key: string]: unknown;
  "id": string;
  "label": string;
  "cancer_id": string;
  "comparison": string;
  "context_definition": string | null;
  "comparator_definition": string | null;
  "context_models_n": number | null;
  "comparator_models_n": number | null;
  "genes_analyzed_n": number | null;
  "depmap_release": string | null;
  "retrieved_at": string | null;
  "qc_status": string;
};

export type CompoundCatalogResponse = Record<string, unknown>;

export type CompoundPharmacologyResponse = Record<string, unknown>;

export type CoverageSummary = {
  [key: string]: unknown;
  "annotated_genes_n"?: number | null;
  "gene_universe_n"?: number | null;
  "coverage_fraction"?: number | null;
  "status"?: string | null;
  "reason"?: string | null;
};

export type GeneAnnotationsResponse = Record<string, unknown>;

export type GeneContextsResponse = Record<string, unknown>;

export type GeneDependencyLandscapeResponse = Record<string, unknown>;

export type GeneDetailResponse = {
  [key: string]: unknown;
  "identity": Record<string, unknown>;
  "model_insights"?: Record<string, unknown> | null;
  "comparisons"?: Record<string, unknown>[];
  "pathways"?: Record<string, unknown>[];
};

export type GeneFacetResponse = {
  [key: string]: unknown;
  "taxonomy_version"?: string | null;
  "status"?: string | null;
  "domains"?: Record<string, unknown>[];
  "protein_classes"?: Record<string, unknown>[];
  "compartments"?: Record<string, unknown>[];
  "hallmarks"?: Record<string, unknown>[];
  "sources"?: Record<string, unknown>[];
  "coverage"?: CoverageSummary | Record<string, unknown> | null;
  "reference_coverage"?: ReferenceCoverageSummary | Record<string, unknown> | null;
  "provenance_note"?: string | null;
};

export type GeneListResponse = Record<string, unknown>[];

export type GeneMatrixResponse = {
  [key: string]: unknown;
  "genes_total_after_filters": number;
  "genes_returned_n": number;
  "comparisons_n": number;
  "comparisons": Record<string, unknown>[];
  "rows": Record<string, unknown>[];
  "filters": Record<string, unknown>;
  "interpretation": Record<string, unknown>;
};

export type GeneModelsResponse = Record<string, unknown>;

export type GeneMutationAssociationsResponse = Record<string, unknown>;

export type GeneSearchResponse = {
  [key: string]: unknown;
  "page": number;
  "page_size": number;
  "total": number;
  "pages": number;
  "sort_by": string;
  "sort_order": string;
  "items": Record<string, unknown>[];
  "functional_coverage"?: CoverageSummary | null;
  "reference_coverage"?: ReferenceCoverageSummary | null;
};

export type GeneSuggestResponse = Record<string, unknown>[];

export type HypothesisCatalogResponse = Record<string, unknown>;

export type HypothesisDetailResponse = Record<string, unknown>;

export type HypothesisSummaryResponse = Record<string, unknown>;

export type LaboratoryCandidateDetailResponse = Record<string, unknown>;

export type LaboratoryCandidatesResponse = Record<string, unknown>;

export type LaboratoryLinesResponse = Record<string, unknown>;

export type LaboratoryMechanismPanelsResponse = Record<string, unknown>;

export type LaboratorySummaryResponse = Record<string, unknown>;

export type ModelDependenciesResponse = Record<string, unknown>;

export type ModelDetailResponse = Record<string, unknown>;

export type ModelListItem = {
  [key: string]: unknown;
  "cancer_id": string;
  "model_id": string;
  "cell_line_name"?: string | null;
  "cancer_ru"?: string | null;
  "organ_ru"?: string | null;
  "molecular_context"?: string | null;
  "assigned_group": string;
  "assigned_group_ru"?: string | null;
  "sequencing_available"?: boolean | null;
  "alteration_status"?: string | null;
  "assignment_reason"?: string | null;
  "variants"?: ModelVariant[];
};

export type ModelMultiomicsResponse = Record<string, unknown>;

export type ModelPharmacologyResponse = Record<string, unknown>;

export type ModelVariant = {
  [key: string]: unknown;
  "gene"?: string | null;
  "protein_change"?: string | null;
  "dna_change"?: string | null;
  "hotspot"?: unknown;
  "driver"?: unknown;
  "classification"?: string | null;
  "kind"?: string | null;
};

export type ModelsResponse = {
  [key: string]: unknown;
  "total": number;
  "items": ModelListItem[];
};

export type MoleculeAnalysisResponse = Record<string, unknown>;

export type MoleculeAnalyzeRequest = {
  "smiles": string;
  "name"?: string | null;
  "force_recompute"?: boolean;
};

export type MoleculeCellLinesResponse = Record<string, unknown>;

export type MoleculeDetailResponse = Record<string, unknown>;

export type MoleculeExperimentResponse = Record<string, unknown>;

export type MoleculeTargetsResponse = Record<string, unknown>;

export type MultiomicsResponse = Record<string, unknown>;

export type NetworkResponse = Record<string, unknown>;

export type OverviewResponse = {
  [key: string]: unknown;
  "genes_analyzed_n": number;
  "comparisons_n": number;
  "stable_recurrent_genes_n": number;
  "stable_pathways_n": number;
  "thresholds"?: number[];
  "per_threshold"?: Record<string, unknown>;
  "analysis_version"?: string | null;
  "generated_at"?: string | null;
  "data_release"?: string | null;
};

export type PathwayListResponse = Record<string, unknown>[];

export type PathwayStabilityResponse = Record<string, unknown>[];

export type PharmacologyConcordanceResponse = Record<string, unknown>;

export type PharmacologySummaryResponse = Record<string, unknown>;

export type PhenotypeCompoundsResponse = Record<string, unknown>;

export type PhenotypeExperimentCreateRequest = {
  "title": string;
  "protocol_id": string;
  "cell_line_ids": string[];
  "compound_ids": string[];
  "exposure_h": number;
  "assay_type"?: string;
};

export type PhenotypeExperimentResponse = Record<string, unknown>;

export type PhenotypeImportRequest = {
  "filename": string;
  "csv_text"?: string | null;
  "file_base64"?: string | null;
};

export type PhenotypeImportResponse = Record<string, unknown>;

export type PhenotypeModelsResponse = Record<string, unknown>;

export type PyzAdmet = {
  "status": "available" | "partial" | "not_available" | "pending_review";
  "data_available": boolean;
  "integrated_risk"?: string | null;
  "note_ru"?: string | null;
  "source_report"?: string | null;
  "prediction_not_experiment": boolean;
  "individual_ames"?: string | null;
  "individual_dili"?: string | null;
};

export type PyzCatalogResponse = Record<string, unknown>;

export type PyzCellContext = {
  "status": "available" | "partial" | "not_available" | "pending_review";
  "data_available": boolean;
  "items"?: PyzContext[];
};

export type PyzContext = {
  "target_gene"?: string | null;
  "mcl_organ_ru"?: string | null;
  "mcl_cancer_name"?: string | null;
  "mcl_cancer_id"?: string | null;
  "tanimoto"?: number | null;
  "models"?: PyzContextModel[];
  "context_status": string;
  "source": string;
};

export type PyzContextModel = {
  "model_id"?: string | null;
  "cell_line_name"?: string | null;
  "role"?: string | null;
  "dependency_probability"?: number | null;
  "gene_effect"?: number | null;
  "available_in_laboratory_registry": boolean;
  "known_reference_compounds"?: string[];
};

export type PyzDetailResponse = {
  [key: string]: unknown;
  "profile_v2": PyzProfileV2;
};

export type PyzIdentity = {
  "status": "available" | "partial" | "not_available" | "pending_review";
  "data_available": boolean;
  "own_compound_id": string;
  "standardized_smiles": string | null;
  "chemotype_ru": string | null;
  "ring_variant": string | null;
  "series_range": string | null;
  "source": string;
};

export type PyzMatrixResponse = Record<string, unknown>;

export type PyzNeighbor = {
  "ligand_id"?: string | null;
  "tanimoto"?: number | null;
  "reference_endpoint"?: string | null;
  "reference_value_nm"?: number | null;
};

export type PyzProfileV2 = {
  "schema_version": string;
  "compound_id": string;
  "identity": PyzIdentity;
  "properties": PyzProperties;
  "admet": PyzAdmet;
  "target_hypotheses": PyzTargetHypotheses;
  "cell_context": PyzCellContext;
  "scientific_summary": PyzScientificSummary;
  "provenance": PyzProvenance;
};

export type PyzProperties = {
  "status": "available" | "partial" | "not_available" | "pending_review";
  "data_available": boolean;
  "qed"?: number | null;
  "source"?: string | null;
};

export type PyzProvenance = {
  "own_compounds_registry": string;
  "admet_report"?: string | null;
  "target_library": string;
  "candidate_context": string;
  "bindingdb_p036_status": string;
  "docking_status": string;
};

export type PyzScientificSummary = {
  "established": ScientificClaim[];
  "hypotheses": ScientificClaim[];
  "gaps": ScientificClaim[];
  "next_step": ScientificClaim;
  "rule_ids": string[];
  "own_cellular_activity": string;
  "own_target_binding": string;
};

export type PyzSummaryResponse = Record<string, unknown>;

export type PyzTarget = {
  "target_gene": string;
  "uniprot_accession"?: string | null;
  "tanimoto"?: number | null;
  "passes_navigation_cutoff": boolean;
  "ligands_compared_n"?: number | null;
  "best_ligand_id"?: string | null;
  "reference_value_nm"?: number | null;
  "reference_endpoint"?: string | null;
  "reference_measurement_owner": string;
  "evidence_type": string;
  "scientific_status": string;
  "source": string;
  "neighbors"?: PyzNeighbor[];
};

export type PyzTargetHypotheses = {
  "status": "available" | "partial" | "not_available" | "pending_review";
  "data_available": boolean;
  "scope": string;
  "cutoff_navigational": number;
  "items": PyzTarget[];
  "source": string;
};

export type QCResponse = Record<string, unknown>;

export type ReferenceCoverageSummary = {
  [key: string]: unknown;
  "resolved_genes_n"?: number | null;
  "gene_universe_n"?: number | null;
  "coverage_fraction"?: number | null;
  "available"?: boolean | null;
};

export type ScientificClaim = {
  "rule_id": string;
  "title": string;
  "meaning": string;
  "source": string;
};

export type StableGenesResponse = Record<string, unknown>[];

export type TargetCatalogResponse = Record<string, unknown>;

export type TargetPharmacologyResponse = Record<string, unknown>;

