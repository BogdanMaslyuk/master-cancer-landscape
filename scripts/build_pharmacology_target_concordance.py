from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RUNTIME = ROOT / "data" / "runtime" / "pharmacology"
RESPONSES = RUNTIME / "responses.parquet"
TARGETS = RUNTIME / "target_evidence.parquet"
GENE_EFFECT = ROOT / "data" / "processed" / "depmap_model_gene_effect.parquet"
OUTPUT = RUNTIME / "target_concordance.parquet"
SUMMARY = RUNTIME / "target_concordance_summary.json"

MIN_SHARED_MODELS = 20
DEPENDENCY_THRESHOLD = -0.5


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value).strip()


def _action_class(value: object) -> str:
    text = _text(value).lower()
    if not text:
        return "unknown"
    lof_tokens = (
        "inhibitor", "inhibition", "antagonist", "blocker", "degrader",
        "suppressor", "negative modulator", "inverse agonist",
    )
    gof_tokens = ("agonist", "activator", "positive modulator")
    if any(token in text for token in lof_tokens):
        return "loss_of_function_like"
    if any(token in text for token in gof_tokens):
        return "gain_of_function_like"
    return "unknown"


def _endpoint_orientation(source: object, endpoint: object) -> str:
    source_text = _text(source).upper()
    endpoint_text = _text(endpoint).upper()
    # PRISM LFC: more negative values represent stronger depletion / lower viability.
    if source_text == "PRISM" and endpoint_text == "LFC":
        return "lower_is_more_sensitive"
    # Keep future datasets conservative until each source adapter declares semantics.
    return "unknown"


def _rankdata(values: pd.Series) -> pd.Series:
    return pd.to_numeric(values, errors="coerce").rank(method="average")


def _spearman(x: pd.Series, y: pd.Series) -> float | None:
    pair = pd.DataFrame({"x": pd.to_numeric(x, errors="coerce"), "y": pd.to_numeric(y, errors="coerce")}).dropna()
    if len(pair) < 3 or pair["x"].nunique() < 2 or pair["y"].nunique() < 2:
        return None
    rx = _rankdata(pair["x"])
    ry = _rankdata(pair["y"])
    value = rx.corr(ry, method="pearson")
    return None if pd.isna(value) else float(value)


def _normal_cdf(value: float) -> float:
    return 0.5 * (1.0 + math.erf(value / math.sqrt(2.0)))


def _spearman_p_approx(rho: float | None, n: int) -> float | None:
    # Large-sample approximation. We label it explicitly as approximate rather than
    # presenting it as an exact Spearman test. It avoids a hard SciPy dependency.
    if rho is None or n < 10 or abs(rho) >= 1.0:
        return None
    denom = max(1e-12, 1.0 - rho * rho)
    t = abs(rho) * math.sqrt(max(0.0, (n - 2) / denom))
    # For the sample sizes used here, the normal approximation is conservative enough
    # for navigation. Confirmatory inference belongs in source-specific analyses.
    return float(max(0.0, min(1.0, 2.0 * (1.0 - _normal_cdf(t)))))


def _bh_adjust(values: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    result = pd.Series(np.nan, index=values.index, dtype=float)
    valid = numeric.dropna().sort_values()
    m = len(valid)
    if not m:
        return result
    adjusted = pd.Series(index=valid.index, dtype=float)
    running = 1.0
    for rank_from_end, (idx, p_value) in enumerate(reversed(list(valid.items())), start=1):
        rank = m - rank_from_end + 1
        candidate = min(1.0, float(p_value) * m / rank)
        running = min(running, candidate)
        adjusted.loc[idx] = running
    result.loc[adjusted.index] = adjusted
    return result


def _label(row: pd.Series) -> tuple[str, str]:
    n = int(row.get("models_n") or 0)
    rho = row.get("spearman_rho")
    q = row.get("q_value")
    delta = row.get("median_response_delta_dependent_minus_other")
    action = _text(row.get("action_class"))
    orientation = _text(row.get("response_orientation"))

    if n < MIN_SHARED_MODELS:
        return "insufficient", "Недостаточно общих клеточных моделей для устойчивой интерпретации."
    if action != "loss_of_function_like":
        return "direction_not_resolved", "Направление фармакологического действия нельзя напрямую сопоставить с CRISPR loss-of-function."
    if orientation != "lower_is_more_sensitive":
        return "endpoint_not_resolved", "Для этого endpoint направление чувствительности ещё не закреплено в MCL."
    if rho is None or delta is None or pd.isna(rho) or pd.isna(delta):
        return "inconclusive", "Недостаточно вариации для оценки согласованности."

    rho = float(rho)
    delta = float(delta)
    q_value = None if q is None or pd.isna(q) else float(q)

    # When lower response = greater sensitivity and lower Gene Effect = stronger
    # dependency, a positive rho plus a negative dependent-vs-other response delta
    # supports the expected inhibitor/CRISPR relationship.
    if rho >= 0.30 and delta < 0 and q_value is not None and q_value < 0.05:
        return "strong_support", "Чувствительность к препарату хорошо согласуется с CRISPR-зависимостью его аннотированной мишени."
    if rho >= 0.20 and delta < 0:
        return "supportive", "Есть согласованность между фармакологической чувствительностью и CRISPR-зависимостью, но она не достигает строгого уровня MCL."
    if rho <= -0.20 or delta > 0:
        return "discordant", "Наблюдаемая фармакологическая чувствительность не согласуется с простой моделью ингибирование мишени → фенотип CRISPR loss-of-function."
    return "inconclusive", "Связь слабая или неоднозначная; механизм нельзя подтвердить этим слоем данных."


def main() -> None:
    if not RESPONSES.exists() or not TARGETS.exists():
        raise SystemExit("Build MCL Pharmacology Layer v1 first: scripts/build_pharmacology_layer.py")
    if not GENE_EFFECT.exists():
        raise SystemExit("Missing data/processed/depmap_model_gene_effect.parquet")

    responses = pd.read_parquet(RESPONSES)
    targets = pd.read_parquet(TARGETS)
    if responses.empty or targets.empty:
        pd.DataFrame().to_parquet(OUTPUT, index=False)
        SUMMARY.write_text(json.dumps({"status": "empty", "built_at": datetime.now(timezone.utc).isoformat()}, indent=2), encoding="utf-8")
        print("No pharmacology responses or target annotations available for concordance analysis.")
        return

    responses["model_id"] = responses["model_id"].astype(str)
    responses["compound_id"] = responses["compound_id"].astype(str)
    targets["compound_id"] = targets["compound_id"].astype(str)
    targets["target_gene"] = targets["target_gene"].astype(str).str.upper().str.strip()

    ge = pd.read_parquet(GENE_EFFECT)
    ge["model_id"] = ge["model_id"].astype(str)
    ge = ge.set_index("model_id", drop=True)
    available_genes = set(ge.columns)

    grouping = [c for c in (
        "compound_id", "source", "source_release", "endpoint", "assay_type",
        "dose", "dose_unit", "exposure_time_h",
    ) if c in responses.columns]

    target_pairs = targets[[c for c in (
        "compound_id", "target_gene", "action", "evidence_type", "confidence", "directness"
    ) if c in targets.columns]].drop_duplicates()

    rows: list[dict[str, object]] = []
    for keys, group in responses.groupby(grouping, dropna=False, sort=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        meta = dict(zip(grouping, keys))
        compound_id = str(meta["compound_id"])
        candidates = target_pairs[target_pairs["compound_id"] == compound_id]
        if candidates.empty:
            continue

        response_values = pd.to_numeric(group["value"], errors="coerce")
        base = pd.DataFrame({"model_id": group["model_id"].astype(str), "response_value": response_values}).dropna()
        base = base.groupby("model_id", as_index=False)["response_value"].mean()
        orientation = _endpoint_orientation(meta.get("source"), meta.get("endpoint"))

        for _, target in candidates.iterrows():
            gene = _text(target.get("target_gene")).upper()
            if not gene or gene not in available_genes:
                rows.append({
                    **meta,
                    "target_gene": gene,
                    "action": target.get("action"),
                    "action_class": _action_class(target.get("action")),
                    "evidence_type": target.get("evidence_type"),
                    "confidence": target.get("confidence"),
                    "directness": target.get("directness"),
                    "response_orientation": orientation,
                    "models_n": 0,
                    "status": "target_not_in_gene_effect_index",
                })
                continue

            effect = pd.to_numeric(ge[gene], errors="coerce").rename("gene_effect").reset_index()
            joined = base.merge(effect, on="model_id", how="inner").dropna(subset=["response_value", "gene_effect"])
            n = int(len(joined))
            rho = _spearman(joined["gene_effect"], joined["response_value"])
            p_value = _spearman_p_approx(rho, n)
            dependent = joined[joined["gene_effect"] <= DEPENDENCY_THRESHOLD]
            other = joined[joined["gene_effect"] > DEPENDENCY_THRESHOLD]
            dep_median = float(dependent["response_value"].median()) if not dependent.empty else None
            other_median = float(other["response_value"].median()) if not other.empty else None
            delta = None if dep_median is None or other_median is None else dep_median - other_median

            rows.append({
                **meta,
                "target_gene": gene,
                "action": target.get("action"),
                "action_class": _action_class(target.get("action")),
                "evidence_type": target.get("evidence_type"),
                "confidence": target.get("confidence"),
                "directness": target.get("directness"),
                "response_orientation": orientation,
                "models_n": n,
                "dependent_models_n": int(len(dependent)),
                "other_models_n": int(len(other)),
                "spearman_rho": rho,
                "p_value_approx": p_value,
                "median_response_dependent": dep_median,
                "median_response_other": other_median,
                "median_response_delta_dependent_minus_other": delta,
                "dependency_threshold": DEPENDENCY_THRESHOLD,
                "status": "analyzed" if n else "no_overlap",
            })

    result = pd.DataFrame(rows)
    if result.empty:
        result.to_parquet(OUTPUT, index=False, compression="zstd")
        SUMMARY.write_text(json.dumps({"status": "empty", "built_at": datetime.now(timezone.utc).isoformat()}, indent=2), encoding="utf-8")
        print("No compound-target pairs could be analyzed.")
        return

    result["q_value"] = _bh_adjust(result.get("p_value_approx", pd.Series(np.nan, index=result.index)))
    labels = result.apply(_label, axis=1)
    result["concordance_label"] = [x[0] for x in labels]
    result["interpretation_ru"] = [x[1] for x in labels]

    sort_cols = [c for c in ("concordance_label", "spearman_rho", "models_n") if c in result.columns]
    if sort_cols:
        result = result.sort_values(sort_cols, ascending=[True, False, False][: len(sort_cols)], na_position="last")

    RUNTIME.mkdir(parents=True, exist_ok=True)
    result.to_parquet(OUTPUT, index=False, compression="zstd")
    counts = result["concordance_label"].value_counts(dropna=False).to_dict()
    summary = {
        "status": "available",
        "built_at": datetime.now(timezone.utc).isoformat(),
        "pairs_n": int(len(result)),
        "compounds_n": int(result["compound_id"].nunique()) if "compound_id" in result.columns else 0,
        "targets_n": int(result["target_gene"].replace("", np.nan).dropna().nunique()) if "target_gene" in result.columns else 0,
        "min_shared_models": MIN_SHARED_MODELS,
        "dependency_threshold": DEPENDENCY_THRESHOLD,
        "label_counts": {str(k): int(v) for k, v in counts.items()},
        "method_note_ru": (
            "Для PRISM LFC более отрицательный ответ означает более сильное снижение жизнеспособности; "
            "более отрицательный Chronos Gene Effect означает более сильную CRISPR-зависимость. "
            "Поэтому положительная Spearman-корреляция ожидается для ингибитора, если аннотированная "
            "мишень действительно объясняет часть профиля чувствительности. p-value является крупновыборочным "
            "приближением; q-value используется только как навигационная оценка и не заменяет подтверждающий анализ."
        ),
    }
    SUMMARY.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")

    print("MCL pharmacology target concordance")
    print(f"Compound-target analyses: {summary['pairs_n']}")
    print(f"Compounds: {summary['compounds_n']}")
    print(f"Targets: {summary['targets_n']}")
    for label, count in sorted(summary["label_counts"].items()):
        print(f"  {label}: {count}")
    print(f"Wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
