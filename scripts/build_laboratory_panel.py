from __future__ import annotations

import argparse
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config" / "laboratory_cell_lines.tsv"
IDENTITY_OVERRIDES = ROOT / "config" / "laboratory_identity_overrides.tsv"
PROCESSED = ROOT / "data" / "processed"
RAW_DEPMAP = ROOT / "data" / "raw" / "depmap"
ATLAS = PROCESSED / "depmap_crispr_model_atlas.parquet"
OUTPUT = PROCESSED / "laboratory_panel.parquet"
OUTPUT_TSV = PROCESSED / "laboratory_panel.tsv"
MANIFEST = PROCESSED / "laboratory_panel_manifest.json"
QC_JSON = ROOT / "outputs" / "qc" / "laboratory_panel_qc.json"
QC_TSV = ROOT / "outputs" / "qc" / "laboratory_panel_qc.tsv"


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _text(value: Any) -> str:
    try:
        if value is None or pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    return str(value).strip()


def _norm(value: Any) -> str:
    return "".join(ch for ch in _text(value).upper() if ch.isalnum())


def _bool(value: Any) -> bool:
    return _text(value).lower() in {"1", "true", "yes", "y"}


def _release(explicit: str | None) -> str:
    if explicit:
        return explicit
    multi = PROCESSED / "depmap_model_multiomics_manifest.json"
    if multi.exists():
        payload = json.loads(multi.read_text(encoding="utf-8"))
        value = _text(payload.get("depmap_release"))
        if value:
            return value
    atlas_qc = ROOT / "outputs" / "qc" / "depmap_crispr_model_atlas_qc.json"
    if atlas_qc.exists():
        payload = json.loads(atlas_qc.read_text(encoding="utf-8"))
        value = _text(payload.get("depmap_release"))
        if value:
            return value
    candidates = sorted((p.name for p in RAW_DEPMAP.iterdir() if p.is_dir()), reverse=True) if RAW_DEPMAP.exists() else []
    if candidates:
        return candidates[0]
    raise SystemExit("Cannot infer DepMap release. Pass --release 26Q1.")


def _candidate_tokens(row: pd.Series) -> set[str]:
    tokens: set[str] = set()
    for column in ("CellLineName", "StrippedCellLineName", "CCLEName", "ModelIDAlias"):
        value = _text(row.get(column))
        if not value:
            continue
        parts = [value]
        parts.extend(re.split(r"[|;,/_ ]+", value))
        for part in parts:
            token = _norm(part)
            if len(token) >= 2:
                tokens.add(token)
    return tokens


def _load_models(release: str) -> pd.DataFrame:
    path = RAW_DEPMAP / release / "Model.csv"
    if not path.exists():
        raise SystemExit(f"Missing {path.relative_to(ROOT)}")
    header = pd.read_csv(path, nrows=0).columns.tolist()
    wanted = [
        "ModelID", "CellLineName", "StrippedCellLineName", "CCLEName", "ModelIDAlias",
        "OncotreeLineage", "OncotreePrimaryDisease", "OncotreeSubtype", "OncotreeCode",
        "RRID", "CatalogNumber", "SourceType", "SourceDetail", "TissueOrigin",
    ]
    usecols = [c for c in wanted if c in header]
    frame = pd.read_csv(path, usecols=usecols, low_memory=False)
    frame["ModelID"] = frame["ModelID"].map(_text)
    frame["_tokens"] = frame.apply(_candidate_tokens, axis=1)
    return frame


def _load_identity_overrides() -> dict[str, dict[str, str]]:
    if not IDENTITY_OVERRIDES.exists():
        return {}
    frame = pd.read_csv(IDENTITY_OVERRIDES, sep="\t", dtype=str).fillna("")
    required = {"lab_name", "decision", "depmap_model_id", "status", "evidence_note_ru"}
    missing = required - set(frame.columns)
    if missing:
        raise SystemExit(
            "laboratory_identity_overrides.tsv is missing columns: " + ", ".join(sorted(missing))
        )
    return {
        _norm(row["lab_name"]): {str(k): _text(v) for k, v in row.items()}
        for _, row in frame.iterrows()
        if _norm(row["lab_name"])
    }


def _atlas_lookup() -> pd.DataFrame:
    if not ATLAS.exists():
        return pd.DataFrame(columns=["model_id"])
    frame = pd.read_parquet(ATLAS)
    keep = [c for c in (
        "model_id", "cell_line_name", "mcl_cancer_id", "mcl_cancer_name", "mcl_organ_ru",
        "mcl_system_ru", "oncotree_lineage", "oncotree_primary_disease", "oncotree_subtype",
        "oncotree_code", "depmap_release", "has_crispr",
    ) if c in frame.columns]
    return frame[keep].drop_duplicates("model_id")


def _payload_from_hit(hit: pd.Series, *, method: str) -> dict[str, Any]:
    payload: dict[str, Any] = {
        "match_status": "matched",
        "match_method": method,
        "match_candidates_n": 1,
        "model_id": _text(hit.get("ModelID")),
    }
    for source, target in (
        ("CellLineName", "depmap_cell_line_name"),
        ("StrippedCellLineName", "depmap_stripped_cell_line_name"),
        ("CCLEName", "ccle_name"),
        ("RRID", "rrid"),
        ("CatalogNumber", "catalog_number"),
        ("OncotreeLineage", "model_oncotree_lineage"),
        ("OncotreePrimaryDisease", "model_oncotree_primary_disease"),
        ("OncotreeSubtype", "model_oncotree_subtype"),
        ("OncotreeCode", "model_oncotree_code"),
    ):
        payload[target] = _text(hit.get(source)) or None
    return payload


def _match_one(
    row: pd.Series,
    models: pd.DataFrame,
    override: dict[str, str] | None,
) -> dict[str, Any]:
    species = _text(row.get("species")).lower()
    if species != "human":
        return {
            "match_status": "not_applicable_nonhuman",
            "match_method": None,
            "match_candidates_n": 0,
            "model_id": None,
            "identity_status": "not_applicable_nonhuman",
            "identity_note_ru": None,
            "preferred_candidate_model_id": None,
        }

    override = override or {}
    decision = _text(override.get("decision"))
    curated_id = _text(override.get("depmap_model_id"))
    identity_status = _text(override.get("status")) or "automatic_alias_match"
    identity_note = _text(override.get("evidence_note_ru")) or None

    if decision == "force_match" and curated_id:
        hits = models[models["ModelID"].astype(str) == curated_id]
        if len(hits) == 1:
            payload = _payload_from_hit(hits.iloc[0], method="curated_model_id_override")
            payload.update({
                "identity_status": identity_status,
                "identity_note_ru": identity_note,
                "preferred_candidate_model_id": curated_id,
            })
            return payload
        return {
            "match_status": "requires_review_override_missing",
            "match_method": "curated_model_id_override",
            "match_candidates_n": int(len(hits)),
            "model_id": None,
            "identity_status": "override_missing_in_pinned_release",
            "identity_note_ru": identity_note,
            "preferred_candidate_model_id": curated_id,
        }

    aliases = [_text(row.get("lab_name"))]
    aliases.extend(x.strip() for x in _text(row.get("aliases")).split("|") if x.strip())
    norms = {_norm(x) for x in aliases if _norm(x)}
    hit_mask = models["_tokens"].map(lambda tokens: bool(norms & tokens))
    hits = models.loc[hit_mask].copy()

    base_annotations = {
        "identity_status": identity_status,
        "identity_note_ru": identity_note,
        "preferred_candidate_model_id": curated_id or None,
    }

    if hits.empty:
        return {
            "match_status": "not_found",
            "match_method": None,
            "match_candidates_n": 0,
            "model_id": None,
            **base_annotations,
        }
    if hits["ModelID"].nunique() != 1:
        ids = sorted(set(hits["ModelID"].dropna().astype(str)))
        return {
            "match_status": "requires_review_ambiguous",
            "match_method": "exact_normalized_alias",
            "match_candidates_n": len(ids),
            "model_id": None,
            "candidate_model_ids": "|".join(ids[:20]),
            **base_annotations,
        }

    payload = _payload_from_hit(hits.iloc[0], method="exact_normalized_alias")
    payload.update(base_annotations)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description="Map the department's physical cell-line collection to pinned DepMap identities.")
    parser.add_argument("--release", default=None)
    args = parser.parse_args()

    if not CONFIG.exists():
        raise SystemExit(f"Missing {CONFIG.relative_to(ROOT)}")
    release = _release(args.release)
    source = pd.read_csv(CONFIG, sep="\t", dtype=str).fillna("")
    models = _load_models(release)
    atlas = _atlas_lookup()
    overrides = _load_identity_overrides()

    matched_rows: list[dict[str, Any]] = []
    for _, row in source.iterrows():
        base = row.to_dict()
        base["highlighted_on_source"] = _bool(base.get("highlighted_on_source"))
        base["laboratory_available"] = True
        override = overrides.get(_norm(row.get("lab_name")))
        base.update(_match_one(row, models, override))
        matched_rows.append(base)

    out = pd.DataFrame(matched_rows)
    if not atlas.empty:
        out = out.merge(atlas, on="model_id", how="left", suffixes=("", "_atlas"))
    out["has_crispr_atlas"] = out.get("mcl_cancer_id", pd.Series(index=out.index, dtype=object)).notna()
    out["laboratory_role"] = "other"
    out.loc[(out["species"].str.lower() == "human") & out["model_class"].str.startswith("tumor"), "laboratory_role"] = "human_tumor"
    out.loc[out["control_role"].eq("general_non_tumor_control"), "laboratory_role"] = "human_non_tumor_control"
    out.loc[out["model_class"].eq("technical_transformed"), "laboratory_role"] = "human_technical"
    out["depmap_release"] = release
    out["built_at"] = _now()

    PROCESSED.mkdir(parents=True, exist_ok=True)
    QC_JSON.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(OUTPUT, index=False, compression="zstd")
    out.to_csv(OUTPUT_TSV, sep="\t", index=False)

    human = out["species"].str.lower().eq("human")
    tumor = out["laboratory_role"].eq("human_tumor")
    ambiguous = out["match_status"].astype(str).str.startswith("requires_review")
    curated = out["match_method"].fillna("").astype(str).eq("curated_model_id_override") & out["match_status"].eq("matched")
    source_confirmation = out.get("identity_status", pd.Series("", index=out.index)).astype(str).eq("requires_source_confirmation")
    qc = {
        "contract": "mcl-laboratory-panel-v1.1",
        "built_at": _now(),
        "depmap_release": release,
        "laboratory_lines_n": int(len(out)),
        "human_lines_n": int(human.sum()),
        "human_tumor_lines_n": int(tumor.sum()),
        "human_non_tumor_controls_n": int(out["laboratory_role"].eq("human_non_tumor_control").sum()),
        "human_matched_depmap_n": int((human & out["match_status"].eq("matched")).sum()),
        "human_in_crispr_atlas_n": int((human & out["has_crispr_atlas"]).sum()),
        "curated_identity_overrides_n": int(curated.sum()),
        "requires_review_n": int(ambiguous.sum()),
        "requires_source_confirmation_n": int(source_confirmation.sum()),
        "not_found_human_n": int((human & out["match_status"].eq("not_found")).sum()),
        "non_tumor_control": "BJ5ta",
        "control_interpretation_ru": (
            "BJ5ta — единственный доступный человеческий неопухолевый контроль. Это hTERT-иммортализованные "
            "фибробласты, поэтому сравнение с ними оценивает доступную лабораторную селективность, но не доказывает "
            "безопасность для нормальной ткани соответствующего органа."
        ),
    }
    QC_JSON.write_text(json.dumps(qc, ensure_ascii=False, indent=2), encoding="utf-8")
    pd.DataFrame([qc]).to_csv(QC_TSV, sep="\t", index=False)
    manifest = {
        **qc,
        "source_registry": str(CONFIG.relative_to(ROOT)),
        "source_identity_overrides": str(IDENTITY_OVERRIDES.relative_to(ROOT)) if IDENTITY_OVERRIDES.exists() else None,
        "source_model_metadata": str((RAW_DEPMAP / release / "Model.csv").relative_to(ROOT)),
        "outputs": [str(OUTPUT.relative_to(ROOT)), str(OUTPUT_TSV.relative_to(ROOT))],
        "matching_contract": (
            "Automatic mapping uses exact normalized aliases only. Curated force_match decisions are stored separately "
            "in config/laboratory_identity_overrides.tsv and must reference a ModelID present in the pinned release. "
            "Multiple alias matches are never auto-resolved. preferred_only records a likely identity for review without "
            "including that model in laboratory candidate calculations. Non-human lines remain outside human DepMap triage."
        ),
    }
    MANIFEST.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")

    print("MCL Laboratory Panel v1.1")
    print(f"Laboratory lines: {qc['laboratory_lines_n']}")
    print(f"Human lines: {qc['human_lines_n']}")
    print(f"Human tumor lines: {qc['human_tumor_lines_n']}")
    print(f"Human matched to DepMap: {qc['human_matched_depmap_n']}")
    print(f"Human in CRISPR Atlas: {qc['human_in_crispr_atlas_n']}")
    print(f"Curated identity overrides: {qc['curated_identity_overrides_n']}")
    print(f"Human not found: {qc['not_found_human_n']}")
    print(f"Requires review: {qc['requires_review_n']}")
    print(f"Requires source confirmation: {qc['requires_source_confirmation_n']}")
    print(f"Wrote {OUTPUT.relative_to(ROOT)}")
    print(f"Wrote {OUTPUT_TSV.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
