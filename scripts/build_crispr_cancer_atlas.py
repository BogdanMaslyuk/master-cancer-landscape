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


# Display-only navigation map. OncoTree fields remain the scientific source of truth.
LINEAGE_DISPLAY: dict[str, tuple[str, str, str]] = {
    "Lung": ("Дыхательная система", "Лёгкое", "lung"),
    "Pleura": ("Дыхательная система", "Плевра", "lung"),
    "Pancreas": ("Пищеварительная система", "Поджелудочная железа", "pancreas"),
    "Breast": ("Молочная железа", "Молочная железа", "tissue"),
    "Colorectal": ("Пищеварительная система", "Толстая и прямая кишка", "tissue"),
    "Bowel": ("Пищеварительная система", "Кишечник", "tissue"),
    "Stomach": ("Пищеварительная система", "Желудок", "tissue"),
    "Esophagus": ("Пищеварительная система", "Пищевод", "tissue"),
    "Liver": ("Пищеварительная система", "Печень", "tissue"),
    "Biliary Tract": ("Пищеварительная система", "Желчные пути", "tissue"),
    "Kidney": ("Мочевая система", "Почка", "tissue"),
    "Bladder/Urinary Tract": ("Мочевая система", "Мочевой пузырь и мочевые пути", "tissue"),
    "Prostate": ("Мужская репродуктивная система", "Предстательная железа", "tissue"),
    "Ovary/Fallopian Tube": ("Женская репродуктивная система", "Яичник и маточная труба", "tissue"),
    "Ovary": ("Женская репродуктивная система", "Яичник", "tissue"),
    "Uterus": ("Женская репродуктивная система", "Матка", "tissue"),
    "Endometrium": ("Женская репродуктивная система", "Эндометрий", "tissue"),
    "Cervix": ("Женская репродуктивная система", "Шейка матки", "tissue"),
    "Skin": ("Кожа", "Кожа", "tissue"),
    "CNS/Brain": ("Нервная система", "Головной мозг и ЦНС", "brain"),
    "Central Nervous System": ("Нервная система", "Головной мозг и ЦНС", "brain"),
    "Peripheral Nervous System": ("Нервная система", "Периферическая нервная система", "tissue"),
    "Myeloid": ("Кроветворная система", "Миелоидные опухоли", "blood"),
    "Lymphoid": ("Лимфатическая и кроветворная система", "Лимфоидные опухоли", "blood"),
    "Bone": ("Опорно-двигательная система", "Кость", "tissue"),
    "Soft Tissue": ("Мягкие ткани", "Мягкие ткани", "tissue"),
    "Head and Neck": ("Голова и шея", "Опухоли головы и шеи", "tissue"),
    "Thyroid": ("Эндокринная система", "Щитовидная железа", "tissue"),
    "Adrenal Gland": ("Эндокринная система", "Надпочечник", "tissue"),
    "Eye": ("Орган зрения", "Глаз", "tissue"),
}


def _text(value: Any) -> str:
    try:
        if value is None or pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    return str(value).strip()


def _slug(value: Any) -> str:
    text = _text(value).lower()
    return "-".join("".join(ch if ch.isalnum() else " " for ch in text).split()) or "unknown"


def _prepare_for_parquet(frame: pd.DataFrame) -> pd.DataFrame:
    """Stringify mixed object columns while preserving quantitative columns."""
    result = frame.copy()
    for col in result.columns:
        if result[col].dtype == object:
            result[col] = result[col].map(lambda x: None if _text(x) == "" else str(x))
    return result


def _resolve_release(explicit: str | None) -> str:
    if explicit:
        return explicit
    inventory = PROCESSED / "depmap_input_inventory.tsv"
    if inventory.exists():
        frame = pd.read_csv(inventory, sep="\t", usecols=["depmap_release"], nrows=1)
        if not frame.empty and _text(frame.iloc[0]["depmap_release"]):
            return _text(frame.iloc[0]["depmap_release"])
    if RAW_DEPMAP.exists():
        releases = sorted((p.name for p in RAW_DEPMAP.iterdir() if p.is_dir()), reverse=True)
        if releases:
            return releases[0]
    raise SystemExit("Cannot infer DepMap release. Pass --release, e.g. --release 26Q1")


def _crispr_model_ids(path: Path) -> list[str]:
    # Read only the row identifier column; never load the ~18.5k gene matrix here.
    frame = pd.read_csv(path, usecols=[0], low_memory=False)
    if frame.empty:
        return []
    ids = frame.iloc[:, 0].map(_text)
    return list(dict.fromkeys(x for x in ids.tolist() if x))


def _load_model_metadata(path: Path, model_ids: set[str]) -> pd.DataFrame:
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
    frame["ModelID"] = frame["ModelID"].map(_text)
    return frame[frame["ModelID"].isin(model_ids)].drop_duplicates("ModelID").copy()


def _display(lineage: Any) -> tuple[str, str, str]:
    key = _text(lineage)
    if key in LINEAGE_DISPLAY:
        return LINEAGE_DISPLAY[key]
    if not key:
        return ("Прочие опухоли", "Не классифицировано", "tissue")
    return ("Прочие опухоли", key, "tissue")


def build(release: str) -> pd.DataFrame:
    source = RAW_DEPMAP / release
    model_path = source / "Model.csv"
    crispr_path = source / "CRISPRGeneEffect.csv"
    missing = [str(p) for p in (model_path, crispr_path) if not p.exists()]
    if missing:
        raise SystemExit("Missing DepMap source files:\n" + "\n".join(missing))

    ordered_ids = _crispr_model_ids(crispr_path)
    models = _load_model_metadata(model_path, set(ordered_ids))
    models = models.rename(columns={
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
    })

    present = set(models["model_id"].map(_text)) if not models.empty else set()
    missing_ids = [x for x in ordered_ids if x not in present]
    if missing_ids:
        models = pd.concat([models, pd.DataFrame({"model_id": missing_ids})], ignore_index=True, sort=False)

    nav = models.get("oncotree_lineage", pd.Series("", index=models.index)).map(_display)
    models["mcl_system_ru"] = nav.map(lambda x: x[0])
    models["mcl_organ_ru"] = nav.map(lambda x: x[1])
    models["mcl_organ_icon"] = nav.map(lambda x: x[2])
    models["mcl_system_id"] = models["mcl_system_ru"].map(_slug)
    models["mcl_organ_id"] = models["mcl_organ_ru"].map(_slug)

    primary = models.get("oncotree_primary_disease", pd.Series("", index=models.index)).map(_text)
    subtype = models.get("oncotree_subtype", pd.Series("", index=models.index)).map(_text)
    model_type = models.get("depmap_model_type", pd.Series("", index=models.index)).map(_text)
    lineage = models.get("oncotree_lineage", pd.Series("", index=models.index)).map(_text)
    models["mcl_cancer_name"] = [p or s or m for p, s, m in zip(primary, subtype, model_type)]
    models["mcl_cancer_id"] = [ _slug(f"{l}|{c}") for l, c in zip(lineage, models["mcl_cancer_name"]) ]
    models["mcl_subtype_name"] = subtype
    models["mcl_subtype_id"] = [ _slug(f"{c}|{s}") for c, s in zip(models["mcl_cancer_id"], subtype) ]

    models["classification_status"] = ["classified" if l and c else "requires_review" for l, c in zip(lineage, models["mcl_cancer_name"])]
    models["classification_confidence"] = models["classification_status"].map({"classified": "high", "requires_review": "low"})
    models["has_crispr"] = True
    models["depmap_release"] = release

    order = {model_id: i for i, model_id in enumerate(ordered_ids)}
    models["_order"] = models["model_id"].map(order)
    return models.sort_values("_order").drop(columns="_order").reset_index(drop=True)


def _write(frame: pd.DataFrame, release: str) -> None:
    PROCESSED.mkdir(parents=True, exist_ok=True)
    QC_DIR.mkdir(parents=True, exist_ok=True)
    parquet = PROCESSED / "depmap_crispr_model_atlas.parquet"
    tsv = PROCESSED / "depmap_crispr_model_atlas.tsv"
    export = _prepare_for_parquet(frame)
    export.to_parquet(parquet, index=False)
    export.to_csv(tsv, sep="\t", index=False)

    subtype_mask = export["mcl_subtype_name"].fillna("").astype(str).str.strip() != ""
    origin = export.get("primary_or_metastasis", pd.Series("", index=export.index)).fillna("").astype(str).str.lower()
    qc = {
        "depmap_release": release,
        "crispr_models_n": int(export["model_id"].nunique()),
        "organs_n": int(export.loc[export["classification_status"] == "classified", "mcl_organ_id"].nunique()),
        "cancers_n": int(export.loc[export["classification_status"] == "classified", "mcl_cancer_id"].nunique()),
        "subtypes_n": int(export.loc[subtype_mask, "mcl_subtype_id"].nunique()),
        "requires_review_n": int((export["classification_status"] == "requires_review").sum()),
        "primary_models_n": int(origin.eq("primary").sum()),
        "metastatic_models_n": int(origin.eq("metastasis").sum()),
        "source": "DepMap Model.csv joined to rows present in CRISPRGeneEffect.csv",
        "classification_source": "OncoTree fields from DepMap Model.csv",
        "note": "SampleCollectionSite is metadata only and never defines the primary tumor organ.",
    }
    (QC_DIR / "depmap_crispr_model_atlas_qc.json").write_text(json.dumps(qc, ensure_ascii=False, indent=2), encoding="utf-8")
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
    parser = argparse.ArgumentParser(description="Build MCL Cancer Atlas from every DepMap model with a CRISPR Gene Effect profile.")
    parser.add_argument("--release", default=None, help="DepMap release, e.g. 26Q1. Inferred by default.")
    args = parser.parse_args()
    release = _resolve_release(args.release)
    _write(build(release), release)


if __name__ == "__main__":
    main()
