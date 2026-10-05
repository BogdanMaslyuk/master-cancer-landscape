from __future__ import annotations

import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RUNTIME = ROOT / "data" / "runtime" / "pharmacology"
PROCESSED = ROOT / "data" / "processed"
RESPONSES = RUNTIME / "responses.parquet"
TARGETS = RUNTIME / "target_evidence.parquet"
GENE_EFFECT = PROCESSED / "depmap_model_gene_effect.parquet"
GENE_DEPENDENCY = PROCESSED / "depmap_model_gene_dependency.parquet"
ATLAS = PROCESSED / "depmap_crispr_model_atlas.parquet"
OUTPUT = RUNTIME / "target_concordance_v2.parquet"
SUMMARY = RUNTIME / "target_concordance_v2_summary.json"

MIN_SHARED_MODELS = 20
MIN_LINEAGE_MODELS = 3
DEPENDENCY_PROBABILITY_THRESHOLD = 0.5


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value).strip()


def _safe_int(value: object, default: int = 0) -> int:
    """Convert pandas/NumPy scalar counters to int without crashing on NaN/NA."""
    try:
        if value is None or pd.isna(value):
            return default
        return int(value)
    except (TypeError, ValueError, OverflowError):
        return default


def _action_class(value: object) -> str:
    text = _text(value).lower()
    if not text:
        return "unknown"
    lof_tokens = (
        "inhibitor", "inhibition", "antagonist", "blocker", "degrader",
        "suppressor", "negative modulator", "inverse agonist",
    )
    gof_text = text.replace("inverse agonist", "")
    has_lof = any(token in text for token in lof_tokens)
    has_gof = bool(re.search(r"\bagonist\b", gof_text)) or any(
        token in gof_text for token in ("activator", "positive modulator")
    )
    if has_lof and has_gof:
        return "mixed_direction"
    if has_lof:
        return "loss_of_function_like"
    if has_gof:
        return "gain_of_function_like"
    return "unknown"


def _join_unique(values: pd.Series) -> str:
    items: list[str] = []
    for value in values:
        text = _text(value)
        if text and text not in items:
            items.append(text)
    return " | ".join(items)


def _aggregate_target_pairs(targets: pd.DataFrame) -> pd.DataFrame:
    """One row per compound × target; never cherry-pick the most favorable annotation."""
    rows: list[dict[str, object]] = []
    for (compound_id, gene), group in targets.groupby(["compound_id", "target_gene"], sort=False):
        action_classes = {
            _action_class(value) for value in group.get("action", pd.Series(dtype=object))
        } - {"unknown"}
        if not action_classes:
            action_class = "unknown"
        elif len(action_classes) == 1:
            action_class = next(iter(action_classes))
        else:
            action_class = "mixed_direction"
        rows.append(
            {
                "compound_id": str(compound_id),
                "target_gene": str(gene).upper(),
                "action": _join_unique(group.get("action", pd.Series(dtype=object))),
                "action_class": action_class,
                "evidence_type": _join_unique(group.get("evidence_type", pd.Series(dtype=object))),
                "confidence": _join_unique(group.get("confidence", pd.Series(dtype=object))),
                "directness": _join_unique(group.get("directness", pd.Series(dtype=object))),
                "target_evidence_rows_n": int(len(group)),
            }
        )
    return pd.DataFrame(rows)


def _endpoint_orientation(source: object, endpoint: object) -> str:
    if _text(source).upper() == "PRISM" and _text(endpoint).upper() == "LFC":
        return "lower_is_more_sensitive"
    return "unknown"


def _spearman(x: pd.Series, y: pd.Series) -> float | None:
    pair = pd.DataFrame(
        {"x": pd.to_numeric(x, errors="coerce"), "y": pd.to_numeric(y, errors="coerce")}
    ).dropna()
    if len(pair) < 3 or pair["x"].nunique() < 2 or pair["y"].nunique() < 2:
        return None
    rho = pair["x"].rank(method="average").corr(
        pair["y"].rank(method="average"), method="pearson"
    )
    return None if pd.isna(rho) else float(rho)


def _lineage_adjusted_spearman(frame: pd.DataFrame) -> tuple[float | None, int, int]:
    needed = frame[["gene_effect", "response_value", "lineage"]].copy()
    needed["gene_effect"] = pd.to_numeric(needed["gene_effect"], errors="coerce")
    needed["response_value"] = pd.to_numeric(needed["response_value"], errors="coerce")
    needed["lineage"] = needed["lineage"].fillna("").astype(str)
    needed = needed.dropna(subset=["gene_effect", "response_value"])
    counts = needed["lineage"].value_counts()
    valid_lineages = set(counts[counts >= MIN_LINEAGE_MODELS].index)
    needed = needed[needed["lineage"].isin(valid_lineages)].copy()
    if len(needed) < MIN_SHARED_MODELS or len(valid_lineages) < 2:
        return None, int(len(needed)), int(len(valid_lineages))

    needed["rank_ge"] = needed["gene_effect"].rank(method="average")
    needed["rank_response"] = needed["response_value"].rank(method="average")
    needed["ge_residual"] = (
        needed["rank_ge"] - needed.groupby("lineage")["rank_ge"].transform("mean")
    )
    needed["response_residual"] = (
        needed["rank_response"] - needed.groupby("lineage")["rank_response"].transform("mean")
    )
    if needed["ge_residual"].nunique() < 2 or needed["response_residual"].nunique() < 2:
        return None, int(len(needed)), int(len(valid_lineages))
    rho = needed["ge_residual"].corr(needed["response_residual"], method="pearson")
    return (
        None if pd.isna(rho) else float(rho),
        int(len(needed)),
        int(len(valid_lineages)),
    )


def _normal_cdf(value: float) -> float:
    return 0.5 * (1.0 + math.erf(value / math.sqrt(2.0)))


def _p_approx(rho: float | None, effective_n: int) -> float | None:
    """Navigation-only normal approximation; effective_n is reduced for lineage fixed effects."""
    if rho is None or effective_n < 10 or abs(rho) >= 1:
        return None
    denom = max(1e-12, 1.0 - rho * rho)
    t = abs(rho) * math.sqrt(max(0.0, (effective_n - 2) / denom))
    return float(max(0.0, min(1.0, 2.0 * (1.0 - _normal_cdf(t)))))


def _bh(values: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    out = pd.Series(np.nan, index=values.index, dtype=float)
    valid = numeric.dropna().sort_values()
    m = len(valid)
    running = 1.0
    for reverse_rank, (idx, value) in enumerate(reversed(list(valid.items())), start=1):
        rank = m - reverse_rank + 1
        candidate = min(1.0, float(value) * m / rank)
        running = min(running, candidate)
        out.loc[idx] = running
    return out


def _label(row: pd.Series) -> tuple[str, str]:
    n = _safe_int(row.get("models_n"))
    rho = row.get("primary_rho")
    q = row.get("q_value")
    delta = row.get("median_response_delta_dependent_minus_other")
    action = _text(row.get("action_class"))
    orientation = _text(row.get("response_orientation"))
    dependency_measured_n = _safe_int(row.get("dependency_probability_models_n"))

    if n < MIN_SHARED_MODELS:
        return "insufficient", "Недостаточно общих моделей для устойчивой профильной оценки."
    if action != "loss_of_function_like":
        return (
            "direction_not_resolved",
            "Направление действия не является однозначно loss-of-function-like; прямое сопоставление с CRISPR knockout запрещено.",
        )
    if orientation != "lower_is_more_sensitive":
        return "endpoint_not_resolved", "Направление чувствительности для endpoint не закреплено."
    if dependency_measured_n < MIN_SHARED_MODELS:
        return (
            "dependency_probability_missing",
            "Недостаточно Probability of Dependency для корректного dependent/non-dependent сравнения.",
        )
    if rho is None or pd.isna(rho) or delta is None or pd.isna(delta):
        return "inconclusive", "Недостаточно вариации для оценки согласованности."

    rho = float(rho)
    delta = float(delta)
    q_value = None if q is None or pd.isna(q) else float(q)
    if rho >= 0.30 and delta < 0 and q_value is not None and q_value < 0.05:
        return (
            "strong_support",
            "Фармакологический профиль согласуется с CRISPR-профилем мишени после доступной поправки на lineage.",
        )
    if rho >= 0.20 and delta < 0:
        return "supportive", "Есть профильная согласованность, но она не достигает строгого уровня MCL."
    if rho <= -0.20 or delta > 0:
        return (
            "discordant",
            "Профиль не согласуется с простой моделью ингибирование мишени → CRISPR loss-of-function фенотип.",
        )
    return "inconclusive", "Связь слабая или неоднозначная."


def main() -> None:
    for path in (RESPONSES, TARGETS, GENE_EFFECT, GENE_DEPENDENCY, ATLAS):
        if not path.exists():
            raise SystemExit(
                f"Missing {path.relative_to(ROOT)}. Candidate v2 requires CRISPRGeneDependency and the pinned CRISPR Atlas."
            )

    responses = pd.read_parquet(RESPONSES)
    targets = pd.read_parquet(TARGETS)
    ge = pd.read_parquet(GENE_EFFECT).set_index("model_id")
    dp = pd.read_parquet(GENE_DEPENDENCY).set_index("model_id")
    atlas = pd.read_parquet(ATLAS)
    atlas["model_id"] = atlas["model_id"].astype(str)
    lineage_col = "oncotree_lineage" if "oncotree_lineage" in atlas.columns else "mcl_cancer_id"
    lineage = atlas[["model_id", lineage_col]].drop_duplicates("model_id").rename(
        columns={lineage_col: "lineage"}
    )

    responses["model_id"] = responses["model_id"].astype(str)
    responses["compound_id"] = responses["compound_id"].astype(str)
    targets["compound_id"] = targets["compound_id"].astype(str)
    targets["target_gene"] = targets["target_gene"].astype(str).str.upper().str.strip()
    ge.index = ge.index.astype(str)
    dp.index = dp.index.astype(str)

    available_genes = set(ge.columns) & set(dp.columns)
    grouping = [
        c for c in (
            "compound_id", "source", "source_release", "endpoint", "assay_type",
            "dose", "dose_unit", "exposure_time_h"
        ) if c in responses.columns
    ]
    target_pairs = _aggregate_target_pairs(targets)

    rows: list[dict[str, object]] = []
    for keys, group in responses.groupby(grouping, dropna=False, sort=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        meta = dict(zip(grouping, keys))
        compound_id = str(meta["compound_id"])
        candidates = target_pairs[target_pairs["compound_id"] == compound_id]
        if candidates.empty:
            continue
        base = pd.DataFrame(
            {
                "model_id": group["model_id"].astype(str),
                "response_value": pd.to_numeric(group["value"], errors="coerce"),
            }
        ).dropna().groupby("model_id", as_index=False)["response_value"].mean()
        base = base.merge(lineage, on="model_id", how="left")
        orientation = _endpoint_orientation(meta.get("source"), meta.get("endpoint"))

        for _, target in candidates.iterrows():
            gene = _text(target.get("target_gene")).upper()
            common_meta = {
                **meta,
                "target_gene": gene,
                "action": target.get("action"),
                "action_class": target.get("action_class"),
                "evidence_type": target.get("evidence_type"),
                "confidence": target.get("confidence"),
                "directness": target.get("directness"),
                "target_evidence_rows_n": _safe_int(target.get("target_evidence_rows_n")),
                "response_orientation": orientation,
            }
            if not gene or gene not in available_genes:
                rows.append(
                    {**common_meta, "models_n": 0, "status": "target_not_in_dependency_indexes"}
                )
                continue

            effect = pd.to_numeric(ge[gene], errors="coerce").rename("gene_effect").reset_index()
            effect.columns = ["model_id", "gene_effect"]
            probability = pd.to_numeric(dp[gene], errors="coerce").rename(
                "dependency_probability"
            ).reset_index()
            probability.columns = ["model_id", "dependency_probability"]
            joined = base.merge(effect, on="model_id", how="inner").merge(
                probability, on="model_id", how="left"
            )
            joined = joined.dropna(subset=["response_value", "gene_effect"])

            raw_rho = _spearman(joined["gene_effect"], joined["response_value"])
            adjusted_rho, adjusted_n, lineages_n = _lineage_adjusted_spearman(joined)
            primary_rho = adjusted_rho if adjusted_rho is not None else raw_rho
            primary_method = (
                "lineage_fixed_effect_rank_residual"
                if adjusted_rho is not None
                else "raw_spearman"
            )
            if adjusted_rho is not None:
                effective_n = max(0, adjusted_n - (lineages_n - 1))
            else:
                effective_n = int(len(joined))
            p_value = _p_approx(primary_rho, effective_n)

            dep_measured = joined.dropna(subset=["dependency_probability"])
            dependent = dep_measured[
                dep_measured["dependency_probability"] > DEPENDENCY_PROBABILITY_THRESHOLD
            ]
            other = dep_measured[
                dep_measured["dependency_probability"] <= DEPENDENCY_PROBABILITY_THRESHOLD
            ]
            dep_median = float(dependent["response_value"].median()) if not dependent.empty else None
            other_median = float(other["response_value"].median()) if not other.empty else None
            delta = None if dep_median is None or other_median is None else dep_median - other_median

            rows.append(
                {
                    **common_meta,
                    "models_n": int(len(joined)),
                    "dependency_probability_models_n": int(len(dep_measured)),
                    "dependent_models_n": int(len(dependent)),
                    "other_models_n": int(len(other)),
                    "spearman_raw_rho": raw_rho,
                    "spearman_lineage_adjusted_rho": adjusted_rho,
                    "lineage_adjusted_models_n": adjusted_n,
                    "lineages_n": lineages_n,
                    "primary_rho": primary_rho,
                    "primary_rho_method": primary_method,
                    "p_value_effective_n": effective_n,
                    "p_value_approx": p_value,
                    "median_response_dependent": dep_median,
                    "median_response_other": other_median,
                    "median_response_delta_dependent_minus_other": delta,
                    "dependency_probability_threshold": DEPENDENCY_PROBABILITY_THRESHOLD,
                    "status": "analyzed",
                }
            )

    result = pd.DataFrame(rows)
    if result.empty:
        result.to_parquet(OUTPUT, index=False)
        SUMMARY.write_text(
            json.dumps({"status": "empty", "built_at": _now()}, indent=2),
            encoding="utf-8",
        )
        print("No pairs available for concordance v2.")
        return

    result["q_value"] = _bh(
        result.get("p_value_approx", pd.Series(np.nan, index=result.index))
    )
    labels = result.apply(_label, axis=1)
    result["concordance_label"] = [x[0] for x in labels]
    result["interpretation_ru"] = [x[1] for x in labels]
    result = result.sort_values(
        ["concordance_label", "primary_rho", "models_n"],
        ascending=[True, False, False],
        na_position="last",
    )

    RUNTIME.mkdir(parents=True, exist_ok=True)
    result.to_parquet(OUTPUT, index=False, compression="zstd")
    counts = result["concordance_label"].value_counts(dropna=False).to_dict()
    summary = {
        "contract": "mcl-pharmacology-target-concordance-v2",
        "status": "available",
        "built_at": _now(),
        "pairs_n": int(len(result)),
        "compounds_n": int(result["compound_id"].nunique()),
        "targets_n": int(result["target_gene"].replace("", np.nan).dropna().nunique()),
        "dependency_probability_threshold": DEPENDENCY_PROBABILITY_THRESHOLD,
        "min_shared_models": MIN_SHARED_MODELS,
        "min_lineage_models": MIN_LINEAGE_MODELS,
        "label_counts": {str(k): int(v) for k, v in counts.items()},
        "method_note_ru": (
            "Одна пара вещество × мишень анализируется один раз; противоречивые направления действия не cherry-pick'ятся. "
            "Gene Effect используется непрерывно. Dependent/non-dependent группы определяются по Probability of Dependency > 0.5. "
            "Основная корреляция по возможности рассчитывается после удаления средних lineage-эффектов из рангов; raw Spearman хранится отдельно. "
            "P/q являются навигационной аппроксимацией с уменьшенным effective N для lineage fixed effects и не заменяют независимую статистическую валидацию. "
            "Анализ не доказывает причинный механизм или target engagement."
        ),
    }
    SUMMARY.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")

    print("MCL pharmacology target concordance v2")
    print(f"Compound-target analyses: {summary['pairs_n']}")
    print(f"Compounds: {summary['compounds_n']}")
    print(f"Targets: {summary['targets_n']}")
    for label, count in sorted(summary["label_counts"].items()):
        print(f"  {label}: {count}")
    print(f"Wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
