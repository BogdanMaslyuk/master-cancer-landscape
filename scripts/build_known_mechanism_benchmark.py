from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config" / "known_mechanism_benchmark.tsv"
HYPOTHESES = ROOT / "data" / "runtime" / "pharmacology" / "candidate_hypotheses.parquet"
OUTPUT = ROOT / "outputs" / "qc" / "known_mechanism_benchmark.tsv"
SUMMARY = ROOT / "outputs" / "qc" / "known_mechanism_benchmark.json"

STATUS_ORDER = {
    "priority_for_in_vitro": 0,
    "supported_hypothesis": 1,
    "exploratory_hypothesis": 2,
    "insufficient_evidence": 3,
}


def main() -> None:
    if not CONFIG.exists():
        raise SystemExit(f"Missing {CONFIG.relative_to(ROOT)}")
    if not HYPOTHESES.exists():
        raise SystemExit(
            "Missing candidate_hypotheses.parquet. Run scripts/build_candidate_hypotheses.py first."
        )

    seeds = pd.read_csv(CONFIG, sep="\t", dtype=str).fillna("")
    hypotheses = pd.read_parquet(HYPOTHESES)
    hypotheses["target_gene"] = hypotheses["target_gene"].astype(str).str.upper()
    hypotheses["preferred_name"] = hypotheses.get(
        "preferred_name", pd.Series("", index=hypotheses.index)
    ).fillna("").astype(str)
    hypotheses["_status_order"] = hypotheses.get(
        "priority_status", pd.Series("", index=hypotheses.index)
    ).map(STATUS_ORDER).fillna(99)

    rows: list[dict[str, object]] = []
    for _, seed in seeds.iterrows():
        query = str(seed["compound_name_query"]).strip().casefold()
        gene = str(seed["target_gene"]).strip().upper()
        hit = hypotheses[
            hypotheses["preferred_name"].str.casefold().str.contains(query, regex=False)
            & hypotheses["target_gene"].eq(gene)
        ].copy()

        if hit.empty:
            rows.append(
                {
                    **seed.to_dict(),
                    "status": "not_found",
                    "matched_compounds_n": 0,
                    "matched_hypotheses_n": 0,
                    "best_priority_status": None,
                    "best_compound_name": None,
                    "best_cancer_name": None,
                    "best_hypothesis_id": None,
                    "joint_support_models_n": 0,
                    "concordance_label": None,
                    "note_ru": (
                        "Пара не найдена в текущем PRISM/MCL пересечении. Это не отрицательная биологическая "
                        "валидация: вещество могло отсутствовать, иметь другое имя или не иметь нужной аннотации."
                    ),
                }
            )
            continue

        hit = hit.sort_values(
            ["_status_order", "joint_support_models_n", "models_n"],
            ascending=[True, False, False],
            na_position="last",
        )
        best = hit.iloc[0]
        best_status = str(best.get("priority_status") or "")
        if best_status in {"priority_for_in_vitro", "supported_hypothesis"}:
            benchmark_status = "mechanism_recovered"
            note = (
                "Известная вещество-мишень пара восстанавливается как одна из поддержанных гипотез. "
                "Это техническая/биологическая положительная проверка маршрута, но не доказательство всех деталей механизма."
            )
        elif best_status == "exploratory_hypothesis":
            benchmark_status = "partially_recovered"
            note = (
                "Известная пара присутствует, но текущие данные дают только исследовательскую поддержку. "
                "Нужно проверить контекст, покрытие и причины несогласованности."
            )
        else:
            benchmark_status = "present_not_recovered"
            note = (
                "Известная пара присутствует в аннотациях, но v1 не поднимает её до поддержанной гипотезы. "
                "Это важный диагностический сигнал для методики MCL, а не повод автоматически считать механизм неверным."
            )

        rows.append(
            {
                **seed.to_dict(),
                "status": benchmark_status,
                "matched_compounds_n": int(hit["compound_id"].nunique()),
                "matched_hypotheses_n": int(len(hit)),
                "best_priority_status": best.get("priority_status"),
                "best_compound_name": best.get("preferred_name"),
                "best_cancer_name": best.get("mcl_cancer_name"),
                "best_hypothesis_id": best.get("hypothesis_id"),
                "joint_support_models_n": int(best.get("joint_support_models_n") or 0),
                "concordance_label": best.get("concordance_label"),
                "note_ru": note,
            }
        )

    out = pd.DataFrame(rows)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUTPUT, sep="\t", index=False)
    counts = out["status"].value_counts().to_dict() if not out.empty else {}
    payload = {
        "benchmark_pairs_n": int(len(out)),
        "status_counts": {str(k): int(v) for k, v in counts.items()},
        "interpretation_ru": (
            "Known Mechanism Benchmark v1 проверяет, способен ли текущий маршрут MCL восстановить несколько "
            "общеизвестных вещество-мишень связей. Он пока не является строгим тестом молекулярного опухолевого "
            "контекста, потому что мутации, амплификации и фьюжены ещё не включены в статус гипотезы."
        ),
        "config": str(CONFIG.relative_to(ROOT)),
        "input": str(HYPOTHESES.relative_to(ROOT)),
        "output": str(OUTPUT.relative_to(ROOT)),
    }
    SUMMARY.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")

    print("MCL Known Mechanism Benchmark v1")
    print(f"Benchmark pairs: {len(out)}")
    for key, value in sorted(counts.items()):
        print(f"  {key}: {value}")
    print(f"Wrote {OUTPUT.relative_to(ROOT)}")
    print(f"Wrote {SUMMARY.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
