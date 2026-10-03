from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
PROCESSED = ROOT / "data" / "processed"
RAW_DEPMAP = ROOT / "data" / "raw" / "depmap"


def _truthy_series(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False).astype(bool)
    return (
        series.astype(str)
        .str.strip()
        .str.lower()
        .isin({"true", "1", "yes", "y", "t"})
    )


def _release_from_inventory() -> str | None:
    path = PROCESSED / "depmap_input_inventory.tsv"
    if not path.exists():
        return None
    frame = pd.read_csv(path, sep="\t", usecols=["depmap_release"], nrows=1)
    if frame.empty:
        return None
    value = str(frame.iloc[0]["depmap_release"]).strip()
    return value or None


def _resolve_release(explicit: str | None) -> str:
    if explicit:
        return explicit
    inventory_release = _release_from_inventory()
    if inventory_release:
        return inventory_release
    if RAW_DEPMAP.exists():
        candidates = sorted([p.name for p in RAW_DEPMAP.iterdir() if p.is_dir()], reverse=True)
        if candidates:
            return candidates[0]
    raise SystemExit("Cannot infer DepMap release. Pass --release, e.g. --release 26Q1")


def _load_audit_models() -> set[str]:
    parquet = PROCESSED / "depmap_context_audit.parquet"
    tsv = PROCESSED / "depmap_context_audit.tsv"
    if parquet.exists():
        audit = pd.read_parquet(parquet, columns=["model_id"])
    elif tsv.exists():
        audit = pd.read_csv(tsv, sep="\t", usecols=["model_id"], low_memory=False)
    else:
        raise SystemExit("Missing data/processed/depmap_context_audit.tsv/parquet. Build the DepMap context audit first.")
    return set(audit["model_id"].dropna().astype(str).str.strip())


def _build_metadata(model_path: Path, model_ids: set[str]) -> pd.DataFrame:
    header = pd.read_csv(model_path, nrows=0)
    wanted = [
        "ModelID", "PatientID", "CellLineName", "StrippedCellLineName", "DepmapModelType",
        "OncotreeLineage", "OncotreePrimaryDisease", "OncotreeSubtype", "OncotreeCode",
        "RRID", "Age", "AgeCategory", "Sex", "PatientRace", "PrimaryOrMetastasis",
        "SampleCollectionSite", "SourceType", "SourceDetail", "CatalogNumber", "ModelType",
        "TissueOrigin", "ModelDerivationMaterial", "ModelTreatment", "PatientTreatmentStatus",
        "PatientTreatmentType", "Stage", "StagingSystem", "PatientTumorGrade", "GrowthPattern",
        "CCLEName", "WTSIMasterCellID", "SangerModelID", "COSMICID", "ModelIDAlias",
    ]
    usecols = [c for c in wanted if c in header.columns]
    frame = pd.read_csv(model_path, usecols=usecols, low_memory=False)
    frame = frame[frame["ModelID"].astype(str).isin(model_ids)].copy()
    rename = {
        "ModelID": "model_id",
        "PatientID": "patient_id",
        "CellLineName": "cell_line_name",
        "StrippedCellLineName": "stripped_cell_line_name",
        "DepmapModelType": "depmap_model_type",
        "OncotreeLineage": "oncotree_lineage",
        "OncotreePrimaryDisease": "oncotree_primary_disease",
        "OncotreeSubtype": "oncotree_subtype",
        "OncotreeCode": "oncotree_code",
        "RRID": "rrid",
        "Age": "age",
        "AgeCategory": "age_category",
        "Sex": "sex",
        "PatientRace": "patient_race",
        "PrimaryOrMetastasis": "primary_or_metastasis",
        "SampleCollectionSite": "sample_collection_site",
        "SourceType": "source_type",
        "SourceDetail": "source_detail",
        "CatalogNumber": "catalog_number",
        "ModelType": "model_type",
        "TissueOrigin": "tissue_origin",
        "ModelDerivationMaterial": "model_derivation_material",
        "ModelTreatment": "model_treatment",
        "PatientTreatmentStatus": "patient_treatment_status",
        "PatientTreatmentType": "patient_treatment_type",
        "Stage": "stage",
        "StagingSystem": "staging_system",
        "PatientTumorGrade": "patient_tumor_grade",
        "GrowthPattern": "growth_pattern",
        "CCLEName": "ccle_name",
        "WTSIMasterCellID": "wtsi_master_cell_id",
        "SangerModelID": "sanger_model_id",
        "COSMICID": "cosmic_id",
        "ModelIDAlias": "model_id_alias",
    }
    return frame.rename(columns=rename).drop_duplicates("model_id")


def _build_mutations(mutation_path: Path, model_ids: set[str], chunksize: int) -> pd.DataFrame:
    header = pd.read_csv(mutation_path, nrows=0)
    wanted = [
        "SequencingID", "ModelID", "IsDefaultEntryForModel", "Chrom", "Pos", "Ref", "Alt",
        "AF", "DP", "VariantType", "DNAChange", "ProteinChange", "HugoSymbol",
        "EnsemblGeneID", "DbsnpRsID", "MolecularConsequence", "VepImpact", "VepClinSig",
        "OncogeneHighImpact", "TumorSuppressorHighImpact", "TranscriptLikelyLof", "CivicID",
        "CivicDescription", "LikelyLoF", "HessDriver", "HessSignature", "RevelScore",
        "AMClass", "AMPathogenicity", "Hotspot", "EntrezGeneID",
    ]
    usecols = [c for c in wanted if c in header.columns]
    chunks: list[pd.DataFrame] = []
    for chunk in pd.read_csv(mutation_path, usecols=usecols, chunksize=chunksize, low_memory=False):
        chunk = chunk[chunk["ModelID"].astype(str).isin(model_ids)].copy()
        if chunk.empty:
            continue
        if "IsDefaultEntryForModel" in chunk.columns:
            default_mask = _truthy_series(chunk["IsDefaultEntryForModel"])
            if default_mask.any():
                chunk = chunk[default_mask].copy()
        chunks.append(chunk)

    if not chunks:
        return pd.DataFrame(columns=["model_id", "gene"])

    frame = pd.concat(chunks, ignore_index=True)
    rename = {
        "SequencingID": "sequencing_id",
        "ModelID": "model_id",
        "Chrom": "chromosome",
        "Pos": "position",
        "Ref": "ref",
        "Alt": "alt",
        "AF": "allele_fraction",
        "DP": "depth",
        "VariantType": "variant_type",
        "DNAChange": "dna_change",
        "ProteinChange": "protein_change",
        "HugoSymbol": "gene",
        "EnsemblGeneID": "ensembl_gene_id",
        "DbsnpRsID": "dbsnp_rs_id",
        "MolecularConsequence": "molecular_consequence",
        "VepImpact": "vep_impact",
        "VepClinSig": "clin_sig",
        "OncogeneHighImpact": "oncogene_high_impact",
        "TumorSuppressorHighImpact": "tumor_suppressor_high_impact",
        "TranscriptLikelyLof": "transcript_likely_lof",
        "CivicID": "civic_id",
        "CivicDescription": "civic_description",
        "LikelyLoF": "likely_lof",
        "HessDriver": "driver",
        "HessSignature": "hess_signature",
        "RevelScore": "revel_score",
        "AMClass": "am_class",
        "AMPathogenicity": "am_pathogenicity",
        "Hotspot": "hotspot",
        "EntrezGeneID": "entrez_gene_id",
    }
    frame = frame.rename(columns=rename)
    frame = frame.drop(columns=["IsDefaultEntryForModel"], errors="ignore")

    for col in ["driver", "hotspot", "likely_lof", "oncogene_high_impact", "tumor_suppressor_high_impact", "transcript_likely_lof"]:
        if col in frame.columns:
            frame[col] = _truthy_series(frame[col])
        else:
            frame[col] = False

    impact = frame.get("vep_impact", pd.Series("", index=frame.index)).astype(str).str.upper()
    frame["is_functional"] = (
        impact.isin(["HIGH", "MODERATE"])
        | frame["driver"]
        | frame["hotspot"]
        | frame["likely_lof"]
        | frame["oncogene_high_impact"]
        | frame["tumor_suppressor_high_impact"]
    )
    frame["priority_score"] = (
        frame["driver"].astype(int) * 5
        + frame["hotspot"].astype(int) * 4
        + frame["likely_lof"].astype(int) * 3
        + impact.eq("HIGH").astype(int) * 3
        + impact.eq("MODERATE").astype(int)
        + frame["oncogene_high_impact"].astype(int) * 2
        + frame["tumor_suppressor_high_impact"].astype(int) * 2
    )

    dedup_cols = [c for c in ["model_id", "gene", "protein_change", "dna_change", "chromosome", "position", "ref", "alt"] if c in frame.columns]
    if dedup_cols:
        frame = frame.sort_values("priority_score", ascending=False).drop_duplicates(dedup_cols)
    return frame.reset_index(drop=True)


def _write_outputs(metadata: pd.DataFrame, mutations: pd.DataFrame, release: str) -> None:
    PROCESSED.mkdir(parents=True, exist_ok=True)
    metadata = metadata.copy()
    metadata["depmap_release"] = release
    mutations = mutations.copy()
    mutations["depmap_release"] = release

    metadata.to_parquet(PROCESSED / "depmap_model_metadata.parquet", index=False)
    metadata.to_csv(PROCESSED / "depmap_model_metadata.tsv", sep="\t", index=False)
    mutations.to_parquet(PROCESSED / "depmap_model_mutations.parquet", index=False)

    if mutations.empty:
        summary = pd.DataFrame(columns=["model_id", "mutations_n", "mutated_genes_n", "functional_variants_n", "driver_mutations_n", "hotspot_mutations_n", "likely_lof_n"])
    else:
        rows = []
        for model_id, group in mutations.groupby("model_id"):
            rows.append(
                {
                    "model_id": model_id,
                    "mutations_n": len(group),
                    "mutated_genes_n": group["gene"].dropna().astype(str).nunique(),
                    "functional_variants_n": int(group["is_functional"].sum()),
                    "driver_mutations_n": int(group["driver"].sum()),
                    "hotspot_mutations_n": int(group["hotspot"].sum()),
                    "likely_lof_n": int(group["likely_lof"].sum()),
                    "depmap_release": release,
                }
            )
        summary = pd.DataFrame(rows).sort_values("model_id")
    summary.to_parquet(PROCESSED / "depmap_model_profile_summary.parquet", index=False)
    summary.to_csv(PROCESSED / "depmap_model_profile_summary.tsv", sep="\t", index=False)

    print(f"DepMap release: {release}")
    print(f"Indexed models: {metadata['model_id'].nunique() if not metadata.empty else 0}")
    print(f"Mutation rows: {len(mutations)}")
    print(f"Models with mutation records: {mutations['model_id'].nunique() if not mutations.empty else 0}")
    print("Wrote data/processed/depmap_model_metadata.*")
    print("Wrote data/processed/depmap_model_mutations.parquet")
    print("Wrote data/processed/depmap_model_profile_summary.*")


def main() -> None:
    parser = argparse.ArgumentParser(description="Build compact per-cell-line molecular profiles for MCL Explorer.")
    parser.add_argument("--release", default=None, help="DepMap release directory, e.g. 26Q1. Inferred from inventory by default.")
    parser.add_argument("--chunksize", type=int, default=250_000, help="Rows per mutation CSV chunk.")
    args = parser.parse_args()

    release = _resolve_release(args.release)
    release_dir = RAW_DEPMAP / release
    model_path = release_dir / "Model.csv"
    mutation_path = release_dir / "OmicsSomaticMutations.csv"
    missing = [str(p) for p in [model_path, mutation_path] if not p.exists()]
    if missing:
        raise SystemExit("Missing DepMap source files:\n" + "\n".join(missing))

    model_ids = _load_audit_models()
    print(f"MCL audit models to index: {len(model_ids)}")
    metadata = _build_metadata(model_path, model_ids)
    mutations = _build_mutations(mutation_path, model_ids, max(10_000, args.chunksize))
    _write_outputs(metadata, mutations, release)


if __name__ == "__main__":
    main()
