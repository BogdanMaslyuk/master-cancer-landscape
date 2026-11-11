/** Additive, strict presentation contract for GET /api/pyz/{id}. */
export type Availability = "available" | "partial" | "not_available" | "pending_review";

export interface EvidenceClaim {
  rule_id: string;
  title: string;
  meaning: string;
  source: string;
}
export interface DataBlock {
  status: Availability;
  data_available: boolean;
}
export interface PyzTargetNeighbor {
  ligand_id: string | null;
  tanimoto: number | null;
  reference_endpoint: string | null;
  reference_value_nm: number | null;
}
export interface PyzTarget {
  target_gene: string;
  uniprot_accession: string | null;
  tanimoto: number | null;
  passes_navigation_cutoff: boolean;
  ligands_compared_n: number | null;
  best_ligand_id: string | null;
  reference_value_nm: number | null;
  reference_endpoint: string | null;
  reference_measurement_owner: "known_ligand_not_pyz";
  evidence_type: "computational_prediction";
  scientific_status: "exploratory_not_directly_validated";
  source: string;
  neighbors: PyzTargetNeighbor[];
}
export interface PyzCellModel {
  model_id: string | null;
  cell_line_name: string | null;
  role: string | null;
  dependency_probability: number | null;
  gene_effect: number | null;
  available_in_laboratory_registry: boolean;
  known_reference_compounds: string[];
}
export interface PyzContext {
  target_gene: string | null;
  mcl_organ_ru: string | null;
  mcl_cancer_name: string | null;
  mcl_cancer_id: string | null;
  tanimoto: number | null;
  models: PyzCellModel[];
  context_status: "known_target_or_reference_context_not_pyz";
  source: string;
}
export interface PyzProfileV2 {
  schema_version: "mcl-pyz-profile-v2-p0";
  compound_id: string;
  identity: DataBlock & {
    own_compound_id: string;
    standardized_smiles: string | null;
    chemotype_ru: string | null;
    ring_variant: string | null;
    series_range: string | null;
    source: string;
  };
  properties: DataBlock & {
    qed: number | null;
    source: string | null;
  };
  admet: DataBlock & {
    integrated_risk: string | null;
    note_ru: string | null;
    source_report: string | null;
    prediction_not_experiment: true;
    individual_ames: string | null;
    individual_dili: string | null;
  };
  target_hypotheses: DataBlock & {
    scope: "existing_13_core_targets";
    cutoff_navigational: number;
    items: PyzTarget[];
    source: string;
  };
  cell_context: DataBlock & { items: PyzContext[] };
  scientific_summary: {
    established: EvidenceClaim[];
    hypotheses: EvidenceClaim[];
    gaps: EvidenceClaim[];
    next_step: EvidenceClaim;
    rule_ids: string[];
    own_cellular_activity: "not_measured";
    own_target_binding: "not_measured";
  };
  provenance: {
    own_compounds_registry: string;
    admet_report: string | null;
    target_library: string;
    candidate_context: string;
    bindingdb_p036_status: "NOT_INGESTED_UNREVIEWED";
    docking_status: string;
  };
}
export interface PyzDetailV2Response {
  available: boolean;
  profile_v2: PyzProfileV2;
}
