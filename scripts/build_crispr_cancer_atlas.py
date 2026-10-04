from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
PROCESSED = ROOT / "data" / "processed"
RAW_DEPMAP = ROOT / "data" / "raw" / "depmap"
QC_DIR = ROOT / "outputs" / "qc"


# This is deliberately a compact display ontology, not a replacement for OncoTree.
# OncoTree remains the scientific classification source of truth; these labels only
# make the Atlas easier to navigate in Russian.
LINEAGE_DISPLAY: dict[str, dict[str, str]] = {
    "Lung": {"system_ru": "Дыхательная система", "organ_ru": "Лёгкое", "icon": "lung"},
    "Pancreas": {"system_ru": "Пищеварительная система", "organ_ru": "Поджелудочная железа", "icon": "pancreas"},
    "Breast": {"system_ru": "Молочная железа", "organ_ru": "Молочная железа", "icon": "tissue"},
    "Colorectal": {"system_ru": "Пищеварительная система", "organ_ru": "Толстая и прямая кишка", "icon": "tissue"},
    "Bowel": {"system_ru": "Пищеварительная система", "organ_ru": "Кишечник", "icon": "tissue"},
    "Stomach": {"system_ru": "Пищеварительная система", "organ_ru": "Желудок", "icon": "tissue"},
    "Esophagus": {"system_ru": "Пищеварительная система", "organ_ru": "Пищевод", "icon": "tissue"},
    "Liver": {"system_ru": "Пищеварительная система", "organ_ru": "Печень", "icon": "tissue"},
    "Biliary Tract": {"system_ru": "Пищеварительная система", "organ_ru": "Желчные пути", "icon": "tissue"},
    "Kidney": {"system_ru": "Мочевая система", "organ_ru": "Почка", "icon": "tissue"},
    "Bladder/Urinary Tract": {"system_ru": "Мочевая система", "organ_ru": "Мочевой пузырь и мочевые пути", "icon": "tissue"},
    "Prostate": {"system_ru": "Мужская репродуктивная система", "organ_ru": "Предстательная железа", "icon": "tissue"},
    "Ovary/Fallopian Tube": {"system_ru": "Женская репродуктивная система", "organ_ru": "Яичник и маточная труба", "icon": "tissue"},
    "Ovary": {"system_ru": "Женская репродуктивная система", "organ_ru": "Яичник", "icon": "tissue"},
    "Uterus": {"system_ru": "Женская репродуктивная система", "organ_ru": "Матка", "icon": "tissue"},
    "Endometrium": {"system_ru": "Женская репродуктивная система", "organ_ru": "Эндометрий", "icon": "tissue"},
    "Cervix": {"system_ru": "Женская репродуктивная система", "organ_ru": "Шейка матки", "icon": "tissue"},
    "Skin": {"system_ru": "Кожа", "organ_ru": "Кожа", "icon": "tissue"},
    "CNS/Brain": {"system_ru": "Нервная система", "organ_ru": "Головной мозг и ЦНС", "icon": "brain"},
    "Central Nervous System": {"system_ru": "Нервная система", "organ_ru": "Головной мозг и ЦНС", "icon": "brain"},
    "Peripheral Nervous System": {"system_ru": "Нервная система", "organ_ru": "Периферическая нервная система", "icon": "tissue"},
    "Myeloid": {"system_ru": "Кроветворная система", "organ_ru": "Миелоидные опухоли", "icon": "blood"},
    "Lymphoid": {"system_ru": "Лимфатическая и кроветворная система", "organ_ru": "Лимфоидные опухоли", "icon": "blood"},
    "Bone": {"system_ru": "Опорно-двигательная система", "organ_ru": "Кость", "icon": "tissue"},
    "Soft Tissue": {"system_ru": "Мягкие ткани", "organ_ru": "Мягкие ткани", "icon": "tissue"},
    "Head and Neck": {"system_ru": "Голова и шея", "organ_ru": "Опухоли головы и шеи", "icon": "tissue"},
    "Thyroid": {"system_ru": "Эндокринная система", "organ_ru": "Щитовидная железа", "icon": "tissue"},
    "Adrenal Gland": {"system_ru": "Эндокринная система", "organ_ru": "Надпочечник", "icon": "tissue"},
    "Pleura": {"system_ru": "Дыхательная система", "organ_ru": "Плевра", "icon": "lung"},
    "Eye": {"system_ru": "Орган зрения", "organ_ru": "Глаз", "icon": "tissue"},
}


def _slug(value: Any) -> str:
    text = str(value or "").strip().lower()
    slug = "-".join("".join(ch if ch.isalnum() else " " for ch in text).split())
    return slug or "unknown"


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
    inventory = _release_from_inventory()
    if inventory:
        return inventory
    if RAW_DEPMAP.exists():
        candidates = sorted((p.name for p in RAW_DEPMAP.iterdir() if p.is_dir()), reverse=True)
        if candidates:
            return candidates[0]
    raise SystemExit("Cannot infer DepMap release. Pass --release, e.g. --release 26Q1")


def _crispr_model_ids(path: Path) -> list[str]:
    # Only the first column is required. This avoids loading the ~18.5k-gene matrix.
    frame = pd.read_csv(path, usecols=[0], low_memory=False)
    if frame.empty:
        return []
    ids = frame.iloc[:, 0].dropna().astype(str).str.strip()
    ids = ids[ids != ""]
    return list(dict.fromkeys(ids.tolist()))


def _load_models(path: Path, ids: set[str]) -> pd.DataFrame:
    wanted = [
        "ModelID", "PatientID", "CellLineName", "StrippedCellLineName", "DepmapModelType",
        "OncotreeLineage", "OncotreePrimaryDisease", "OncotreeSubtype", "OncotreeCode",
        "PatientSubtypeFeatures", "RRID", "Age", "AgeCategory", "Sex", "PatientRace",
        "PrimaryOrMetastasis", "SampleCollectionSite", "SourceType", "SourceDetail",
        "CatalogNumber", "ModelType", "TissueOrigin", "ModelDerivationMaterial", "CCLEName",
        "WTSIMasterCellID", "SangerModelID", "COSMICID", "ModelIDAlias",
    ]
    header = pd.read_csv(path, nrows=0).columns.tolist()
    usecols = [c for c in wanted if c in header]
    frame = pd.read_csv(path, usecols=usecols, low_memory=False)
    frame["ModelID"] = frame["ModelID"].astype(str).str.strip()
    return frame[frame["ModelID"].isin(ids)].drop_duplicates("ModelID").copy()


def _display_for_lineage(lineage: Any) -> dict[str, str]:
    key = str(lineage or "").strip()
    if key in LINEAGE_DISPLAY:
        return LINEAGE_DISPLAY[key]
    return {
        "system_ru": "Прочие опухоли",
        "organ_ru": key or "Не классифицировано",
        "icon": "tissue",
    }


def build(release: str) -> pd.DataFrame:
    release_dir = RAW_DEPMAP / release
    model_path = release_dir / "Model.csv"
    crispr_path = release_dir / "CRISPRGeneEffect.csv"
    missing = [str(p) for p in (model_path, crispr_path) if not p.exists()]
    if missing:
        raise SystemExit("Missing DepMap source files:\n" + "\n".join(missing))

    ids = _crispr_model_ids(crispr_path)
    id_set = set(ids)
    models = _load_models(model_path, id_set)

    rename = {
        "ModelID": "model_id", "PatientID": "patient_id", "CellLineName": "cell_line_name",
        "StrippedCellLineName": "stripped_cell_line_name", "DepmapModelType": "depmap_model_type",
        "OncotreeLineage": "oncotree_lineage", "OncotreePrimaryDisease": "oncotree_primary_disease",
        "OncotreeSubtype": "oncotree_subtype", "OncotreeCode": "oncotree_code",
        "PatientSubtypeFeatures": "patient_subtype_features", "RRID": "rrid", "Age": "age",
        "AgeCategory": "age_category", "Sex": "sex", "PatientRace": "patient_race",
        "PrimaryOrMetastasis": "primary_or_metastasis", "SampleCollectionSite": "sample_collection_site",
        "SourceType": "source_type", "SourceDetail": "source_detail", "CatalogNumber": "catalog_number",
        "ModelType": "model_type", "TissueOrigin": "tissue_origin",
        "ModelDerivationMaterial": "model_derivation_material", "CCLEName": "ccle_name",
        "WTSIMasterCellID": "wtsi_master_cell_id", "SangerModelID": "sanger_model_id",
        "COSMICID": "cosmic_id", "ModelIDAlias": "model_id_alias",
    }
    models = models.rename(columns=rename)

    missing_ids = [model_id for model_id in ids if model_id not in set(models["model_id"].astype(str))]
    if missing_ids:
        filler = pd.DataFrame({"model_id": missing_ids})
        models = pd.concat([models, filler], ignore_index=True, sort=False)

    display = models.get("oncotree_lineage", pd.Series("", index=models.index)).map(_display_for_lineage)
    models["mcl_system_ru"] = display.map(lambda x: x["system_ru"])
    models["mcl_organ_ru"] = display.map(lambda x: x["organ_ru"])
    models["mcl_organ_icon"] = display.map(lambda x: x["icon"])
    models["mcl_system_id"] = models["mcl_system_ru"].map(_slug)
    models["mcl_organ_id"] = models["mcl_organ_ru"].map(_slug)

    primary = models.get("oncotree_primary_disease", pd.Series("", index=models.index)).fillna("").astype(str).str.strip()
    subtype = models.get("oncotree_subtype", pd.Series("", index=models.index)).fillna("").astype(str).str.strip()
    model_type = models.get("depmap_model_type", pd.Series("", index=models.index)).fillna("").astype(str).str.strip()
    models["mcl_cancer_name"] = primary.where(primary != "", subtype.where(subtype != "", model_type))
    models["mcl_cancer_id"] = (
        models.get("oncotree_lineage", pd.Series("", index=models.index)).fillna("").astype(str)
        + "|" + models["mcl_cancer_name"].fillna("").astype(str)
    ).map(_slug)
    models["mcl_subtype_name"] = subtype
    models["mcl_subtype_id"] = (
        models["mcl_cancer_id"].fillna("").astype(str) + "|" + subtype.fillna("").astype(str)
    ).map(_slug)

    lineage = models.get("oncotree_lineage", pd.Series("", index=models.index)).fillna("").astype(str).str.strip()
    cancer = models["mcl_cancer_name"].fillna("").astype(str).str.strip()
    models["classification_status"] = "classified"
    models.loc[(lineage == "") | (cancer == ""), "classification_status"] = "requires_review"
    models["classification_confidence"] = "high"
    models.loc[models["classification_status"] == "requires_review", "classification_confidence"] = "low"
    models["has_crispr"] = True
    models["depmap_release"] = release

    order = {model_id: i for i, model_id in enumerate(ids)}
    models["_order"] = models["model_id"].map(order)
    models = models.sort_values("_order").drop(columns="_order").reset_index(drop=True)
    return models


def _write(frame: pd.DataFrame, release: str) -> None:
    PROCESSED.mkdir(parents=True, exist_ok=True)
    QC_DIR.mkdir(parents=True, exist_ok=True)
    parquet = PROCESSED / "depmap_crispr_model_atlas.parquet"
    tsv = PROCESSED / "depmap_crispr_model_atlas.tsv"
    frame.to_parquet(parquet, index=False)
    frame.to_csv(tsv, sep="\t", index=False)

    qc = {
        "depmap_release": release,
        "crispr_models_n": int(frame["model_id"].nunique()),
        "organs_n": int(frame.loc[frame["classification_status"] == "classified", "mcl_organ_id"].nunique()),
        "cancers_n": int(frame.loc[frame["classification_status"] == "classified", "mcl_cancer_id"].nunique()),
        "subtypes_n": int(frame.loc[frame["mcl_subtype_name"].fillna("").astype(str).str.strip() != "", "mcl_subtype_id"].nunique()),
        "requires_review_n": int((frame["classification_status"] == "requires_review").sum()),
        "primary_models_n": int(frame.get("primary_or_metastasis", pd.Series("", index=frame.index)).fillna("").astype(str).str.lower().eq("primary").sum()),
        "metastatic_models_n": int(frame.get("primary_or_metastasis", pd.Series("", index=frame.index)).fillna("").astype(str).str.lower().eq("metastasis").sum()),
        "source": "DepMap Model.csv joined to model rows present in CRISPRGeneEffect.csv",
        "classification_source": "OncoTree fields from DepMap Model.csv",
        "note": "SampleCollectionSite is preserved as metadata and is never used as the primary tumor organ.",
    }
    (QC_DIR / "depmap_crispr_model_atlas_qc.json").write_text(
        json.dumps(qc, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    pd.DataFrame([qc]).to_csv(QC_DIR / "depmap_crispr_model_atlas_qc.tsv", sep="\t", index=False)

    print(f"DepMap release: {release}")
    print(f"CRISPR models: {qc['crispr_models_n']}")
    print(f"Atlas organs: {qc['organs_n']}")
    print(f"Atlas cancers: {qc['cancers_n']}")
    print(f"Requires review: {qc['requires_review_n']}")
    print(f"Wrote {parquet.relative_to(ROOT)}")
    print(f"Wrote {tsv.relative_to(ROOT)}")
    print("Wrote outputs/qc/depmap_crispr_model_atlas_qc.json/.tsv")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build the MCL Cancer Atlas from every DepMap model with a CRISPR Gene Effect profile."
    )
    parser.add_argument("--release", default=None, help="DepMap release, e.g. 26Q1. Inferred by default.")
    args = parser.parse_args()
    release = _resolve_release(args.release)
    frame = build(release)
    _write(frame, release)


if __name__ == "__main__":
    main()
