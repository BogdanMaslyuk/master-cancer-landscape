from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT / "data" / "runtime" / "pyz_target_evidence"
SHORTLIST = EVIDENCE / "pyz_target_validation_shortlist.parquet"
STRUCTURES_CONFIG = ROOT / "config" / "pyz_target_validation_structures.tsv"
OUT_TSV = EVIDENCE / "pyz_target_validation_queue.tsv"
OUT_PARQUET = EVIDENCE / "pyz_target_validation_queue.parquet"
STRUCTURES_OUT = EVIDENCE / "pyz_target_validation_structures.tsv"
MANIFEST = EVIDENCE / "pyz_target_validation_queue_manifest.json"
QC = ROOT / "outputs" / "qc" / "pyz_target_validation_queue_qc.tsv"


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _text(value: Any) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    return str(value).strip()


def _num(value: Any) -> float | None:
    try:
        if value is None or pd.isna(value):
            return None
        x = float(value)
        return x if np.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def _potency_class(value: Any) -> str:
    x = _num(value)
    if x is None:
        return "unknown"
    if x <= 100:
        return "known_ligand_le_100nM"
    if x <= 1000:
        return "known_ligand_100_1000nM"
    if x <= 10000:
        return "known_ligand_1_10uM"
    return "outside_tls_cutoff"


def _assay_hint(target: str) -> str:
    if target == "CDK4":
        return (
            "После вычислительной проверки: прямой биохимический kinase assay CDK4/cyclin-D "
            "и независимый binding/target-engagement метод; клеточный эффект отдельно."
        )
    if target == "BCL2L1":
        return (
            "После вычислительной проверки: прямой BCL-xL binding или BH3-peptide displacement assay; "
            "клеточный апоптоз отдельно и только в релевантном контексте."
        )
    return "Прямой биохимический/биофизический assay мишени после ортогональной вычислительной проверки."


def main() -> None:
    missing = [p for p in (SHORTLIST, STRUCTURES_CONFIG) if not p.exists()]
    if missing:
        raise SystemExit("Missing validation inputs:\n" + "\n".join(str(p.relative_to(ROOT)) for p in missing))

    EVIDENCE.mkdir(parents=True, exist_ok=True)
    QC.parent.mkdir(parents=True, exist_ok=True)

    shortlist = pd.read_parquet(SHORTLIST).copy()
    if shortlist.empty:
        raise SystemExit("PYZ target validation shortlist is empty.")

    shortlist["own_compound_id"] = shortlist["own_compound_id"].astype(str).str.strip()
    shortlist["target_gene"] = shortlist["target_gene"].astype(str).str.upper().str.strip()
    shortlist["best_ligand_id"] = shortlist["best_ligand_id"].astype(str).str.strip()
    shortlist["best_tanimoto_morgan_r2_2048"] = pd.to_numeric(
        shortlist["best_tanimoto_morgan_r2_2048"], errors="coerce"
    )
    potency_col = "best_lowest_reported_value_nm_across_endpoints"
    shortlist[potency_col] = pd.to_numeric(shortlist[potency_col], errors="coerce")

    reuse = shortlist["best_ligand_id"].value_counts().to_dict()
    shortlist["same_best_ligand_shortlist_pairs_n"] = shortlist["best_ligand_id"].map(reuse).fillna(0).astype(int)

    # One representative per repeatedly reused nearest-ligand cluster: highest 2D similarity.
    representative_ids: set[tuple[str, str]] = set()
    for ligand_id, group in shortlist.groupby("best_ligand_id", sort=False):
        if not ligand_id:
            continue
        best = group.sort_values(
            ["best_tanimoto_morgan_r2_2048", potency_col],
            ascending=[False, True],
            na_position="last",
        ).iloc[0]
        representative_ids.add((str(best["own_compound_id"]), str(ligand_id)))

    phases: list[str] = []
    reasons: list[str] = []
    for row in shortlist.to_dict("records"):
        own_id = str(row["own_compound_id"])
        ligand_id = str(row["best_ligand_id"])
        potency = _num(row.get(potency_col))
        reuse_n = int(row.get("same_best_ligand_shortlist_pairs_n") or 0)
        is_rep = (own_id, ligand_id) in representative_ids

        if potency is not None and potency <= 1000:
            phase = "phase_1_orthogonal_validation"
            reason = (
                "Tanimoto >=0.35 и ближайший измеренный лиганд имеет прямую активность <=1 µM; "
                "пара первой очереди для независимой structure-based проверки и прямого assay."
            )
        elif reuse_n >= 2 and is_rep:
            phase = "phase_2_representative_chemotype"
            reason = (
                "Несколько PYZ сходятся к одному ближайшему известному лиганду; сначала проверяется "
                "представитель с максимальным Tanimoto, чтобы не дублировать одно и то же химическое предположение."
            )
        elif reuse_n >= 2:
            phase = "phase_3_expand_within_chemotype_if_supported"
            reason = (
                "Член повторяющегося химического кластера; проверять после представителя только если "
                "ортогональный метод даст согласованный сигнал."
            )
        else:
            phase = "phase_2_secondary_validation"
            reason = (
                "Структурное сходство заметное, но ближайший лиганд слабее 1 µM или биологическая поддержка неполная; "
                "вторая очередь после наиболее информативных пар."
            )
        phases.append(phase)
        reasons.append(reason)

    shortlist["validation_phase"] = phases
    shortlist["validation_reason_ru"] = reasons
    shortlist["known_ligand_potency_class"] = shortlist[potency_col].map(_potency_class)
    shortlist["same_best_ligand_cluster_representative"] = [
        (str(own), str(ligand)) in representative_ids
        for own, ligand in zip(shortlist["own_compound_id"], shortlist["best_ligand_id"])
    ]
    shortlist["recommended_direct_assay_ru"] = shortlist["target_gene"].map(_assay_hint)
    shortlist["docking_interpretation_rule_ru"] = (
        "Docking используется только как ортогональная приоритизация. Низкая энергия/правдоподобная поза не доказывает "
        "связывание, Ki/IC50, селективность, target engagement или противоопухолевый механизм PYZ."
    )
    shortlist["direct_pyz_target_evidence"] = "absent_not_measured"
    shortlist["built_at"] = _now()

    phase_order = {
        "phase_1_orthogonal_validation": 0,
        "phase_2_representative_chemotype": 1,
        "phase_2_secondary_validation": 2,
        "phase_3_expand_within_chemotype_if_supported": 3,
    }
    shortlist["_phase_order"] = shortlist["validation_phase"].map(phase_order).fillna(99)
    queue = shortlist.sort_values(
        ["_phase_order", "best_tanimoto_morgan_r2_2048", potency_col],
        ascending=[True, False, True],
        na_position="last",
    ).drop(columns=["_phase_order"]).reset_index(drop=True)

    structures = pd.read_csv(STRUCTURES_CONFIG, sep="\t", dtype=str, keep_default_na=False)
    structures["target_gene"] = structures["target_gene"].astype(str).str.upper().str.strip()
    structures = structures[structures["target_gene"].isin(set(queue["target_gene"]))].copy()
    structures.to_csv(STRUCTURES_OUT, sep="\t", index=False)

    queue.to_parquet(OUT_PARQUET, index=False, compression="zstd")
    queue.to_csv(OUT_TSV, sep="\t", index=False)

    phase_counts = queue["validation_phase"].value_counts().to_dict()
    unique_best_ligands = int(queue["best_ligand_id"].nunique())
    qc_rows = [
        {"check": "shortlist_rows_preserved", "value": len(queue), "expected": len(shortlist), "status": "PASS"},
        {"check": "targets_with_structure_plan", "value": structures["target_gene"].nunique(), "expected": queue["target_gene"].nunique(), "status": "PASS" if structures["target_gene"].nunique() == queue["target_gene"].nunique() else "WARNING"},
        {"check": "phase1_pairs", "value": int((queue["validation_phase"] == "phase_1_orthogonal_validation").sum()), "expected": ">=1", "status": "PASS" if (queue["validation_phase"] == "phase_1_orthogonal_validation").any() else "WARNING"},
        {"check": "direct_pyz_target_measurements", "value": 0, "expected": 0, "status": "INFO"},
    ]
    pd.DataFrame(qc_rows).to_csv(QC, sep="\t", index=False)

    manifest = {
        "contract": "mcl-pyz-target-validation-queue-v1",
        "built_at": _now(),
        "shortlist_pairs_n": int(len(queue)),
        "targets_n": int(queue["target_gene"].nunique()),
        "unique_best_experimental_ligands_n": unique_best_ligands,
        "phase_counts": {str(k): int(v) for k, v in phase_counts.items()},
        "priority_logic": (
            "No composite 0-100 score. Phase 1 requires Tanimoto >=0.35 (inherited from shortlist) plus a nearest directly measured ligand <=1 µM. "
            "Repeated nearest-ligand clusters are de-duplicated experimentally by testing the highest-similarity representative first."
        ),
        "structure_strategy": (
            "Every target entering the queue must have a curated inhibitor-bound structure plan with mandatory co-crystal redocking before PYZ docking. "
            "CDK4 uses an ensemble to reduce single-conformation dependence; BCL2L1 uses a human BCL-xL inhibitor-bound reference structure."
        ),
        "scientific_guardrails_ru": [
            "Очередь не является доказательством мишени PYZ и не присваивает вероятности связывания.",
            "Tanimoto 0.35-0.50 остается низким структурным сходством; это только основание для проверки.",
            "Активность известного лиганда не переносится количественно на PYZ.",
            "Перед интерпретацией docking обязателен redocking сокристаллизованного лиганда и оценка воспроизводимости кармана.",
            "Даже успешный docking требует прямого биохимического/биофизического подтверждения PYZ-target interaction.",
        ],
        "outputs": [
            str(OUT_TSV.relative_to(ROOT)),
            str(OUT_PARQUET.relative_to(ROOT)),
            str(STRUCTURES_OUT.relative_to(ROOT)),
            str(QC.relative_to(ROOT)),
        ],
    }
    MANIFEST.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")

    print("MCL PYZ Target Validation Queue v1")
    print(f"Pairs: {len(queue)}; targets: {queue['target_gene'].nunique()}; unique nearest ligands: {unique_best_ligands}")
    print("Validation phases:")
    for phase, count in sorted(phase_counts.items(), key=lambda x: phase_order.get(str(x[0]), 99)):
        print(f"  {phase}: {count}")
    print("Phase 1:")
    p1 = queue[queue["validation_phase"].eq("phase_1_orthogonal_validation")]
    for row in p1.to_dict("records"):
        endpoint = _text(row.get("best_lowest_value_endpoint_type"))
        value = _num(row.get(potency_col))
        potency = f"{endpoint}={value:g} nM" if endpoint and value is not None else "direct potency unavailable"
        print(
            f"  {row['own_compound_id']} -> {row['target_gene']}: "
            f"Tanimoto={float(row['best_tanimoto_morgan_r2_2048']):.3f}; {potency}; ligand={row['best_ligand_id']}"
        )
    print("Curated docking structures:")
    for row in structures.to_dict("records"):
        print(
            f"  {row['target_gene']} [{row['validation_role']}]: {row['pdb_id']} — "
            f"{row['protein_context']} + {row['co_crystal_ligand']}"
        )
    print(f"Queue: {OUT_TSV.relative_to(ROOT)}")
    print(f"Structures: {STRUCTURES_OUT.relative_to(ROOT)}")
    print(f"QC: {QC.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
