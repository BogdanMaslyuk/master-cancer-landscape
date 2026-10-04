from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
PROCESSED = ROOT / "data" / "processed"
RUNTIME = ROOT / "data" / "runtime" / "explorer"
GENE_EFFECT = PROCESSED / "depmap_model_gene_effect.parquet"
ATLAS = PROCESSED / "depmap_crispr_model_atlas.parquet"
OUTPUT = RUNTIME / "gene_dependency_summary.parquet"
MANIFEST = RUNTIME / "gene_dependency_summary.json"
DEPENDENCY_THRESHOLD = -0.5
MIN_CANCER_MODELS = 3


def _dependency_type(global_fraction: float, best_fraction: float | None, enrichment: float | None) -> tuple[str, str]:
    best = best_fraction if best_fraction is not None and np.isfinite(best_fraction) else 0.0
    enrich = enrichment if enrichment is not None and np.isfinite(enrichment) else 0.0
    if global_fraction >= 0.70:
        return "broad_baseline", "Широкая базовая зависимость"
    if global_fraction >= 0.30:
        return "broad_tumor", "Широкая опухолевая зависимость"
    if best >= 0.50 and enrich >= 0.25:
        return "cancer_enriched", "Обогащена в этом типе опухоли"
    if global_fraction <= 0.15 and best >= 0.30:
        return "selective", "Селективная зависимость"
    if global_fraction < 0.05 and best < 0.25:
        return "weak", "Слабая / невыраженная"
    return "intermediate", "Промежуточная зависимость"


def _specificity_label(score: float | None, dependency_type: str) -> str:
    if dependency_type in {"broad_baseline", "broad_tumor"}:
        return "низкая"
    if score is None or not np.isfinite(score):
        return "не определена"
    if score >= 0.50:
        return "очень высокая"
    if score >= 0.30:
        return "высокая"
    if score >= 0.15:
        return "средняя"
    return "низкая"


def main() -> None:
    if not GENE_EFFECT.exists():
        raise SystemExit(f"Missing {GENE_EFFECT.relative_to(ROOT)}. Build DepMap multi-omics indexes first.")
    if not ATLAS.exists():
        raise SystemExit(f"Missing {ATLAS.relative_to(ROOT)}. Build the CRISPR Cancer Atlas first.")

    frame = pd.read_parquet(GENE_EFFECT)
    atlas = pd.read_parquet(ATLAS)
    if "model_id" not in frame.columns or "model_id" not in atlas.columns:
        raise SystemExit("Both input files must contain model_id")

    frame["model_id"] = frame["model_id"].astype(str)
    atlas["model_id"] = atlas["model_id"].astype(str)
    atlas = atlas.drop_duplicates("model_id").set_index("model_id")
    model_ids = frame["model_id"].tolist()
    meta = atlas.reindex(model_ids)

    genes = [c for c in frame.columns if c != "model_id"]
    values = frame[genes].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=np.float32, copy=False)
    valid = np.isfinite(values)
    dependent = valid & (values <= DEPENDENCY_THRESHOLD)

    available_n = valid.sum(axis=0).astype(np.int32)
    dependent_n = dependent.sum(axis=0).astype(np.int32)
    global_fraction = np.divide(
        dependent_n,
        available_n,
        out=np.full(len(genes), np.nan, dtype=np.float32),
        where=available_n > 0,
    )
    global_median = np.nanmedian(values, axis=0)

    best_fraction = np.full(len(genes), np.nan, dtype=np.float32)
    best_median = np.full(len(genes), np.nan, dtype=np.float32)
    best_models_n = np.zeros(len(genes), dtype=np.int32)
    best_dependent_n = np.zeros(len(genes), dtype=np.int32)
    best_cancer_id = np.full(len(genes), None, dtype=object)
    best_cancer_name = np.full(len(genes), None, dtype=object)
    best_organ = np.full(len(genes), None, dtype=object)

    cancer_ids = meta.get("mcl_cancer_id", pd.Series(index=meta.index, dtype=object)).fillna("").astype(str)
    for cancer_id in sorted(x for x in cancer_ids.unique() if x):
        row_mask = cancer_ids.to_numpy() == cancer_id
        if int(row_mask.sum()) < MIN_CANCER_MODELS:
            continue
        sub = values[row_mask, :]
        sub_valid = np.isfinite(sub)
        n = sub_valid.sum(axis=0).astype(np.int32)
        dep_n = (sub_valid & (sub <= DEPENDENCY_THRESHOLD)).sum(axis=0).astype(np.int32)
        frac = np.divide(dep_n, n, out=np.full(len(genes), np.nan, dtype=np.float32), where=n >= MIN_CANCER_MODELS)
        median = np.nanmedian(sub, axis=0)

        current = best_fraction
        better = np.isfinite(frac) & (
            ~np.isfinite(current)
            | (frac > current + 1e-7)
            | ((np.abs(frac - current) <= 1e-7) & (median < best_median))
        )
        if not better.any():
            continue

        meta_rows = meta.loc[row_mask]
        cancer_name = str(meta_rows.get("mcl_cancer_name", pd.Series(dtype=object)).dropna().iloc[0]) if "mcl_cancer_name" in meta_rows and meta_rows["mcl_cancer_name"].notna().any() else cancer_id
        organ = str(meta_rows.get("mcl_organ_ru", pd.Series(dtype=object)).dropna().iloc[0]) if "mcl_organ_ru" in meta_rows and meta_rows["mcl_organ_ru"].notna().any() else None

        best_fraction[better] = frac[better]
        best_median[better] = median[better]
        best_models_n[better] = n[better]
        best_dependent_n[better] = dep_n[better]
        best_cancer_id[better] = cancer_id
        best_cancer_name[better] = cancer_name
        best_organ[better] = organ

    enrichment = best_fraction - global_fraction
    rows: list[dict[str, object]] = []
    for i, gene in enumerate(genes):
        gf = float(global_fraction[i]) if np.isfinite(global_fraction[i]) else 0.0
        bf = float(best_fraction[i]) if np.isfinite(best_fraction[i]) else None
        enrich = float(enrichment[i]) if np.isfinite(enrichment[i]) else None
        dep_type, dep_label = _dependency_type(gf, bf, enrich)
        score = max(0.0, enrich) if enrich is not None else None
        rows.append(
            {
                "gene_symbol": str(gene).upper(),
                "dependency_threshold": DEPENDENCY_THRESHOLD,
                "dependency_models_n": int(dependent_n[i]),
                "dependency_models_total_n": int(available_n[i]),
                "dependency_fraction": gf,
                "global_median_gene_effect": float(global_median[i]) if np.isfinite(global_median[i]) else None,
                "dependency_type": dep_type,
                "dependency_type_ru": dep_label,
                "best_cancer_id": best_cancer_id[i],
                "best_cancer_name": best_cancer_name[i],
                "best_cancer_organ_ru": best_organ[i],
                "best_cancer_dependency_models_n": int(best_dependent_n[i]),
                "best_cancer_models_n": int(best_models_n[i]),
                "best_cancer_dependency_fraction": bf,
                "best_cancer_median_gene_effect": float(best_median[i]) if np.isfinite(best_median[i]) else None,
                "specificity_score": score,
                "specificity_label_ru": _specificity_label(score, dep_type),
            }
        )

    out = pd.DataFrame(rows)
    RUNTIME.mkdir(parents=True, exist_ok=True)
    out.to_parquet(OUTPUT, index=False, compression="zstd")
    manifest = {
        "dependency_threshold": DEPENDENCY_THRESHOLD,
        "threshold_semantics": "MCL descriptive navigation threshold; Gene Effect <= threshold is counted as a strong loss-of-fitness dependency.",
        "models_n": len(model_ids),
        "genes_n": len(genes),
        "cancer_groups_min_n": MIN_CANCER_MODELS,
        "source_gene_effect": str(GENE_EFFECT.relative_to(ROOT)),
        "source_atlas": str(ATLAS.relative_to(ROOT)),
        "output": str(OUTPUT.relative_to(ROOT)),
    }
    MANIFEST.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")

    print(f"Gene dependency summary: {len(genes)} genes × {len(model_ids)} models")
    print(f"Operational dependency threshold: Gene Effect <= {DEPENDENCY_THRESHOLD}")
    print(f"Wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
