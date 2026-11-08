from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
PROCESSED = ROOT / "data" / "processed"
RUNTIME = ROOT / "data" / "runtime" / "explorer"
GENE_EFFECT = PROCESSED / "depmap_model_gene_effect.parquet"
GENE_DEPENDENCY = PROCESSED / "depmap_model_gene_dependency.parquet"
ATLAS = PROCESSED / "depmap_crispr_model_atlas.parquet"
OUTPUT = RUNTIME / "gene_dependency_summary.parquet"
MANIFEST = RUNTIME / "gene_dependency_summary.json"
DEPENDENCY_PROBABILITY_THRESHOLD = 0.5
GENE_EFFECT_FALLBACK_THRESHOLD = -0.5
MIN_CANCER_MODELS = 5


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


def _align_probability(effect: pd.DataFrame, genes: list[str]) -> tuple[np.ndarray | None, str]:
    if not GENE_DEPENDENCY.exists():
        return None, "gene_effect_descriptive_fallback"
    probability = pd.read_parquet(GENE_DEPENDENCY)
    if "model_id" not in probability.columns:
        return None, "gene_effect_descriptive_fallback"
    probability["model_id"] = probability["model_id"].astype(str)
    probability = probability.drop_duplicates("model_id").set_index("model_id")
    available = [gene for gene in genes if gene in probability.columns]
    if not available:
        return None, "gene_effect_descriptive_fallback"
    aligned = probability.reindex(effect["model_id"].astype(str))
    matrix = np.full((len(effect), len(genes)), np.nan, dtype=np.float32)
    gene_index = {gene: i for i, gene in enumerate(genes)}
    for gene in available:
        matrix[:, gene_index[gene]] = pd.to_numeric(aligned[gene], errors="coerce").to_numpy(dtype=np.float32)
    return matrix, "probability_of_dependency"


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
    gene_effect_valid = np.isfinite(values)
    probability, call_method = _align_probability(frame, genes)

    if probability is not None:
        valid = np.isfinite(probability)
        dependent = valid & (probability > DEPENDENCY_PROBABILITY_THRESHOLD)
    else:
        # Backward-compatible descriptive fallback only. This is explicitly labelled
        # and should not be confused with the DepMap binary dependency convention.
        valid = gene_effect_valid
        dependent = valid & (values <= GENE_EFFECT_FALLBACK_THRESHOLD)

    available_n = valid.sum(axis=0).astype(np.int32)
    dependent_n = dependent.sum(axis=0).astype(np.int32)
    global_fraction = np.divide(
        dependent_n,
        available_n,
        out=np.full(len(genes), np.nan, dtype=np.float32),
        where=available_n > 0,
    )
    global_median = np.nanmedian(values, axis=0)
    global_median_probability = (
        np.nanmedian(probability, axis=0) if probability is not None else np.full(len(genes), np.nan, dtype=np.float32)
    )

    best_fraction = np.full(len(genes), np.nan, dtype=np.float32)
    best_median = np.full(len(genes), np.nan, dtype=np.float32)
    best_probability_median = np.full(len(genes), np.nan, dtype=np.float32)
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

        sub_effect = values[row_mask, :]
        sub_call_values = probability[row_mask, :] if probability is not None else sub_effect
        sub_valid = np.isfinite(sub_call_values)
        n = sub_valid.sum(axis=0).astype(np.int32)
        if probability is not None:
            dep_n = (sub_valid & (sub_call_values > DEPENDENCY_PROBABILITY_THRESHOLD)).sum(axis=0).astype(np.int32)
        else:
            dep_n = (sub_valid & (sub_call_values <= GENE_EFFECT_FALLBACK_THRESHOLD)).sum(axis=0).astype(np.int32)
        frac = np.divide(dep_n, n, out=np.full(len(genes), np.nan, dtype=np.float32), where=n >= MIN_CANCER_MODELS)
        effect_median = np.nanmedian(sub_effect, axis=0)
        probability_median = (
            np.nanmedian(probability[row_mask, :], axis=0)
            if probability is not None
            else np.full(len(genes), np.nan, dtype=np.float32)
        )

        current = best_fraction
        better = np.isfinite(frac) & (
            ~np.isfinite(current)
            | (frac > current + 1e-7)
            | ((np.abs(frac - current) <= 1e-7) & (effect_median < best_median))
        )
        if not better.any():
            continue

        meta_rows = meta.loc[row_mask]
        cancer_name = (
            str(meta_rows["mcl_cancer_name"].dropna().iloc[0])
            if "mcl_cancer_name" in meta_rows.columns and meta_rows["mcl_cancer_name"].notna().any()
            else cancer_id
        )
        organ = (
            str(meta_rows["mcl_organ_ru"].dropna().iloc[0])
            if "mcl_organ_ru" in meta_rows.columns and meta_rows["mcl_organ_ru"].notna().any()
            else None
        )

        best_fraction[better] = frac[better]
        best_median[better] = effect_median[better]
        best_probability_median[better] = probability_median[better]
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
                "dependency_call_method": call_method,
                "dependency_probability_threshold": DEPENDENCY_PROBABILITY_THRESHOLD if probability is not None else None,
                "gene_effect_fallback_threshold": GENE_EFFECT_FALLBACK_THRESHOLD if probability is None else None,
                "dependency_models_n": int(dependent_n[i]),
                "dependency_models_total_n": int(available_n[i]),
                "dependency_fraction": gf,
                "global_median_dependency_probability": float(global_median_probability[i]) if np.isfinite(global_median_probability[i]) else None,
                "global_median_gene_effect": float(global_median[i]) if np.isfinite(global_median[i]) else None,
                "dependency_type": dep_type,
                "dependency_type_ru": dep_label,
                "best_cancer_id": best_cancer_id[i],
                "best_cancer_name": best_cancer_name[i],
                "best_cancer_organ_ru": best_organ[i],
                "best_cancer_dependency_models_n": int(best_dependent_n[i]),
                "best_cancer_models_n": int(best_models_n[i]),
                "best_cancer_dependency_fraction": bf,
                "best_cancer_median_dependency_probability": (
                    float(best_probability_median[i]) if np.isfinite(best_probability_median[i]) else None
                ),
                "best_cancer_median_gene_effect": float(best_median[i]) if np.isfinite(best_median[i]) else None,
                "specificity_score": score,
                "specificity_label_ru": _specificity_label(score, dep_type),
            }
        )

    out = pd.DataFrame(rows)
    RUNTIME.mkdir(parents=True, exist_ok=True)
    out.to_parquet(OUTPUT, index=False, compression="zstd")
    manifest = {
        "dependency_call_method": call_method,
        "dependency_probability_threshold": DEPENDENCY_PROBABILITY_THRESHOLD if probability is not None else None,
        "gene_effect_fallback_threshold": GENE_EFFECT_FALLBACK_THRESHOLD if probability is None else None,
        "threshold_semantics": (
            "Binary dependent/non-dependent calls use CRISPR Probability of Dependency > 0.5 when available. "
            "If the probability matrix is missing, Gene Effect <= -0.5 is retained only as an explicitly labelled descriptive fallback."
        ),
        "models_n": len(model_ids),
        "genes_n": len(genes),
        "cancer_groups_min_n": MIN_CANCER_MODELS,
        "source_gene_effect": str(GENE_EFFECT.relative_to(ROOT)),
        "source_gene_dependency": str(GENE_DEPENDENCY.relative_to(ROOT)) if GENE_DEPENDENCY.exists() else None,
        "source_atlas": str(ATLAS.relative_to(ROOT)),
        "output": str(OUTPUT.relative_to(ROOT)),
    }
    MANIFEST.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")

    print(f"Gene dependency summary: {len(genes)} genes × {len(model_ids)} models")
    print(f"Binary dependency method: {call_method}")
    if probability is not None:
        print(f"Probability of Dependency threshold: > {DEPENDENCY_PROBABILITY_THRESHOLD}")
    else:
        print(f"WARNING: CRISPRGeneDependency missing; descriptive Gene Effect fallback <= {GENE_EFFECT_FALLBACK_THRESHOLD}")
    print(f"Minimum cancer-group size for specificity labels: {MIN_CANCER_MODELS}")
    print(f"Wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
