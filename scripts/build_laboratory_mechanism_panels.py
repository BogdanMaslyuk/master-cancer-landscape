from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

try:
    import duckdb
except ImportError as exc:
    raise SystemExit(
        'DuckDB is required. Install analysis extras: .\\.venv\\Scripts\\python.exe -m pip install -e ".[analysis]"'
    ) from exc


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config" / "laboratory_mechanism_panels.tsv"
PROCESSED = ROOT / "data" / "processed"
RUNTIME = ROOT / "data" / "runtime" / "laboratory"
PANEL = PROCESSED / "laboratory_panel.parquet"
CONTEXT = PROCESSED / "depmap_model_target_context.parquet"
CANDIDATES = RUNTIME / "laboratory_candidate_hypotheses.parquet"
CANDIDATE_MODELS = RUNTIME / "laboratory_candidate_models.parquet"
OUTPUT = RUNTIME / "laboratory_mechanism_panels.parquet"
MODEL_OUTPUT = RUNTIME / "laboratory_mechanism_panel_models.parquet"
MANIFEST = RUNTIME / "laboratory_mechanism_panels_manifest.json"
QC = ROOT / "outputs" / "qc" / "laboratory_mechanism_panels_qc.tsv"


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


def _split(value: Any) -> list[str]:
    return [x.strip() for x in _text(value).split("|") if x.strip()]


def _json(values: list[Any]) -> str:
    return json.dumps(values, ensure_ascii=False)


def _clean_scalar(value: Any) -> Any:
    try:
        if value is None or pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass
    if hasattr(value, "item"):
        try:
            return value.item()
        except (TypeError, ValueError):
            pass
    return value


def _context_status(panel_id: str, rows: pd.DataFrame) -> tuple[str, str]:
    if rows.empty:
        return "context_missing", "Молекулярный контекст для этой модели в текущем индексе отсутствует."

    panel_key = panel_id.upper()
    if "MDM2" in panel_key:
        tp53 = rows[rows["target_gene"].astype(str).str.upper().eq("TP53")]
        if tp53.empty:
            return "tp53_context_missing", "TP53 не найден в текущем индексированном контексте модели."
        row = tp53.iloc[0]
        if bool(row.get("has_likely_lof", False)):
            return (
                "tp53_likely_lof_concern",
                "Для TP53 есть флаг LikelyLoF. Это механистическое предостережение для ингибирования MDM2 и требует отдельной проверки p53-функции.",
            )
        if bool(row.get("has_hotspot", False)):
            return (
                "tp53_hotspot_requires_variant_interpretation",
                "Для TP53 есть hotspot-аннотация. Её функциональное направление нельзя автоматически считать потерей функции; нужен разбор конкретного варианта.",
            )
        return (
            "tp53_no_flagged_disruptive_annotation",
            "В текущем слое нет TP53 LikelyLoF/hotspot-флага. Это не доказывает TP53 wild-type и не доказывает сохранную функцию p53.",
        )

    if "PIK3CA" in panel_key:
        notes: list[str] = []
        for gene in ("PIK3CA", "PTEN", "KRAS"):
            hit = rows[rows["target_gene"].astype(str).str.upper().eq(gene)]
            if hit.empty:
                continue
            row = hit.iloc[0]
            flags: list[str] = []
            if bool(row.get("has_hotspot", False)):
                flags.append("hotspot")
            if bool(row.get("has_likely_lof", False)):
                flags.append("LikelyLoF")
            if flags:
                notes.append(f"{gene}: {', '.join(flags)}")
        suffix = "; ".join(notes) if notes else "выделенных hotspot/LikelyLoF-флагов нет"
        return (
            "pi3k_context_available_direction_unresolved",
            "PI3K-контекст доступен; " + suffix + ". Эти признаки показываются описательно и не получают автоматического направления чувствительности/резистентности.",
        )

    return "context_available_direction_unresolved", "Молекулярный контекст доступен, но его направление не задано автоматическим правилом."


def _load_context(model_ids: list[str], genes: list[str]) -> pd.DataFrame:
    if not model_ids or not genes or not CONTEXT.exists():
        return pd.DataFrame()
    con = duckdb.connect(database=":memory:")
    model_df = pd.DataFrame({"model_id": sorted(set(model_ids))})
    gene_df = pd.DataFrame({"target_gene": sorted(set(g.upper() for g in genes))})
    con.register("selected_models", model_df)
    con.register("selected_genes", gene_df)
    path = CONTEXT.resolve().as_posix().replace("'", "''")
    columns = set(
        con.execute(f"DESCRIBE SELECT * FROM read_parquet('{path}')").df()["column_name"].astype(str)
    )
    wanted = [
        "model_id", "target_gene", "gene_effect", "dependency_probability",
        "expression_log2_tpm1", "expression_percentile_in_cancer",
        "copy_number_relative", "copy_number_percentile_in_cancer",
        "has_hotspot", "has_likely_lof", "protein_changes_json",
    ]
    expressions = [
        f"c.{col}" if col in columns else f"NULL AS {col}"
        for col in wanted
    ]
    query = f"""
        SELECT {', '.join(expressions)}
        FROM read_parquet('{path}') c
        INNER JOIN selected_models m ON CAST(c.model_id AS VARCHAR) = m.model_id
        INNER JOIN selected_genes g ON upper(CAST(c.target_gene AS VARCHAR)) = g.target_gene
    """
    frame = con.execute(query).df()
    con.close()
    if not frame.empty:
        frame["model_id"] = frame["model_id"].astype(str)
        frame["target_gene"] = frame["target_gene"].astype(str).str.upper()
    return frame


def main() -> None:
    required = [CONFIG, PANEL, CANDIDATES, CANDIDATE_MODELS]
    missing = [str(path.relative_to(ROOT)) for path in required if not path.exists()]
    if missing:
        raise SystemExit("Missing mechanism-panel inputs:\n" + "\n".join(missing))

    config = pd.read_csv(CONFIG, sep="\t", dtype=str).fillna("")
    panel = pd.read_parquet(PANEL)
    candidates = pd.read_parquet(CANDIDATES)
    model_evidence = pd.read_parquet(CANDIDATE_MODELS)

    panel["_lab_norm"] = panel["lab_name"].map(_norm)
    candidates["_compound_norm"] = candidates["preferred_name"].map(_norm)
    candidates["target_gene"] = candidates["target_gene"].astype(str).str.upper()

    selected_specs: list[dict[str, Any]] = []
    selected_model_ids: set[str] = set()
    selected_genes: set[str] = set()

    for _, spec in config.iterrows():
        panel_id = _text(spec["panel_id"])
        target_gene = _text(spec["target_gene"]).upper()
        compounds = {_norm(x) for x in _split(spec["compound_names"])}
        cancer_tokens = [x.casefold() for x in _split(spec.get("cancer_name_contains", ""))]
        context_genes = [x.upper() for x in _split(spec["context_genes"])]
        focus_names = _split(spec["focus_lab_lines"])
        focus_norms = {_norm(x) for x in focus_names}

        selected = candidates[
            candidates["target_gene"].eq(target_gene)
            & candidates["_compound_norm"].isin(compounds)
        ].copy()
        if cancer_tokens:
            cancer_series = selected["mcl_cancer_name"].fillna("").astype(str).str.casefold()
            mask = pd.Series(False, index=selected.index)
            for token in cancer_tokens:
                mask |= cancer_series.str.contains(token, regex=False)
            selected = selected[mask]

        focus = panel[panel["_lab_norm"].isin(focus_norms)].copy()
        for value in focus["model_id"].dropna().astype(str):
            if value.strip():
                selected_model_ids.add(value.strip())
        selected_genes.update(context_genes)
        selected_genes.add(target_gene)

        selected_specs.append({
            "spec": spec,
            "panel_id": panel_id,
            "target_gene": target_gene,
            "context_genes": context_genes,
            "focus": focus,
            "selected": selected,
        })

    context = _load_context(sorted(selected_model_ids), sorted(selected_genes))

    panel_rows: list[dict[str, Any]] = []
    model_rows: list[dict[str, Any]] = []

    for bundle in selected_specs:
        spec = bundle["spec"]
        panel_id = bundle["panel_id"]
        target_gene = bundle["target_gene"]
        context_genes = bundle["context_genes"]
        focus = bundle["focus"]
        selected = bundle["selected"]

        hypothesis_ids = set(selected["hypothesis_id"].astype(str)) if not selected.empty else set()
        evidence = model_evidence[
            model_evidence["hypothesis_id"].astype(str).isin(hypothesis_ids)
        ].copy() if hypothesis_ids else pd.DataFrame()

        positive_lines = sorted(set(
            evidence.loc[evidence.get("laboratory_model_role", pd.Series(index=evidence.index, dtype=object)).eq("positive"), "lab_name"].dropna().astype(str)
        )) if not evidence.empty else []
        negative_lines = sorted(set(
            evidence.loc[evidence.get("laboratory_model_role", pd.Series(index=evidence.index, dtype=object)).eq("negative_same_cancer"), "lab_name"].dropna().astype(str)
        )) if not evidence.empty else []
        discordant_lines = sorted(set(
            evidence.loc[evidence.get("laboratory_model_role", pd.Series(index=evidence.index, dtype=object)).astype(str).str.startswith("discordant"), "lab_name"].dropna().astype(str)
        )) if not evidence.empty else []

        for _, hypothesis in selected.iterrows():
            cancer_id = _text(hypothesis.get("mcl_cancer_id"))
            eligible_focus = focus[
                focus["control_role"].astype(str).eq("general_non_tumor_control")
                | focus.get("mcl_cancer_id", pd.Series("", index=focus.index)).fillna("").astype(str).eq(cancer_id)
            ].copy()
            for _, line in eligible_focus.iterrows():
                model_id = _text(line.get("model_id"))
                lab_name = _text(line.get("lab_name"))
                evidence_hit = evidence[
                    evidence["hypothesis_id"].astype(str).eq(str(hypothesis["hypothesis_id"]))
                    & evidence["lab_name"].astype(str).map(_norm).eq(_norm(lab_name))
                ] if not evidence.empty else pd.DataFrame()
                ev = evidence_hit.iloc[0] if not evidence_hit.empty else pd.Series(dtype=object)

                ctx = context[
                    context["model_id"].astype(str).eq(model_id)
                    & context["target_gene"].astype(str).str.upper().isin(context_genes)
                ].copy() if model_id and not context.empty else pd.DataFrame()
                ctx_records: list[dict[str, Any]] = []
                if not ctx.empty:
                    for _, c in ctx.iterrows():
                        ctx_records.append({
                            "gene": _text(c.get("target_gene")).upper(),
                            "dependency_probability": _clean_scalar(c.get("dependency_probability")),
                            "gene_effect": _clean_scalar(c.get("gene_effect")),
                            "expression_log2_tpm1": _clean_scalar(c.get("expression_log2_tpm1")),
                            "expression_percentile_in_cancer": _clean_scalar(c.get("expression_percentile_in_cancer")),
                            "copy_number_relative": _clean_scalar(c.get("copy_number_relative")),
                            "copy_number_percentile_in_cancer": _clean_scalar(c.get("copy_number_percentile_in_cancer")),
                            "has_hotspot": bool(c.get("has_hotspot", False)) if pd.notna(c.get("has_hotspot")) else False,
                            "has_likely_lof": bool(c.get("has_likely_lof", False)) if pd.notna(c.get("has_likely_lof")) else False,
                            "protein_changes_json": _clean_scalar(c.get("protein_changes_json")),
                        })
                status, status_ru = _context_status(panel_id, ctx)

                model_rows.append({
                    "panel_id": panel_id,
                    "panel_name_ru": _text(spec["panel_name_ru"]),
                    "hypothesis_id": _text(hypothesis.get("hypothesis_id")),
                    "preferred_name": _text(hypothesis.get("preferred_name")),
                    "compound_id": _text(hypothesis.get("compound_id")),
                    "target_gene": target_gene,
                    "mcl_cancer_id": cancer_id,
                    "mcl_cancer_name": _text(hypothesis.get("mcl_cancer_name")),
                    "lab_id": _text(line.get("lab_id")),
                    "lab_name": lab_name,
                    "model_id": model_id or None,
                    "laboratory_role": _text(line.get("laboratory_role")),
                    "laboratory_model_role": _text(ev.get("laboratory_model_role")) or (
                        "general_non_tumor_control" if _text(line.get("control_role")) == "general_non_tumor_control" else "not_evaluated_in_candidate_triage"
                    ),
                    "response_value": _clean_scalar(ev.get("response_value")),
                    "relative_sensitivity": _clean_scalar(ev.get("relative_sensitivity")),
                    "dependency_probability": _clean_scalar(ev.get("dependency_probability")),
                    "gene_effect": _clean_scalar(ev.get("gene_effect")),
                    "mechanism_context_status": status,
                    "mechanism_context_summary_ru": status_ru,
                    "mechanism_context_json": _json(ctx_records),
                    "built_at": _now(),
                })

        candidate_names = sorted(set(selected["preferred_name"].dropna().astype(str))) if not selected.empty else []
        candidate_cancers = sorted(set(selected["mcl_cancer_name"].dropna().astype(str))) if not selected.empty else []
        panel_rows.append({
            "panel_id": panel_id,
            "panel_name_ru": _text(spec["panel_name_ru"]),
            "target_gene": target_gene,
            "compound_names_json": _json(_split(spec["compound_names"])),
            "context_genes_json": _json(context_genes),
            "focus_lab_lines_json": _json(_split(spec["focus_lab_lines"])),
            "primary_question_ru": _text(spec["primary_question_ru"]),
            "guardrail_ru": _text(spec["guardrail_ru"]),
            "candidate_hypotheses_n": int(len(selected)),
            "candidate_names_json": _json(candidate_names),
            "candidate_cancers_json": _json(candidate_cancers),
            "positive_lab_lines_json": _json(positive_lines),
            "negative_lab_lines_json": _json(negative_lines),
            "discordant_lab_lines_json": _json(discordant_lines),
            "focus_lines_found_n": int(len(focus)),
            "focus_lines_mapped_n": int(focus["model_id"].notna().sum()) if not focus.empty else 0,
            "built_at": _now(),
        })

    panels_out = pd.DataFrame(panel_rows)
    models_out = pd.DataFrame(model_rows)
    RUNTIME.mkdir(parents=True, exist_ok=True)
    QC.parent.mkdir(parents=True, exist_ok=True)
    panels_out.to_parquet(OUTPUT, index=False, compression="zstd")
    models_out.to_parquet(MODEL_OUTPUT, index=False, compression="zstd")

    manifest = {
        "contract": "mcl-laboratory-mechanism-context-v1",
        "built_at": _now(),
        "panels_n": int(len(panels_out)),
        "panel_model_rows_n": int(len(models_out)),
        "context_source": str(CONTEXT.relative_to(ROOT)),
        "ranking_contract": "No composite score. Mechanism context is explanatory and can raise a concern, but does not automatically upgrade a candidate.",
        "tp53_contract": "TP53 LikelyLoF is a mechanistic concern for MDM2 inhibition. Hotspot is not automatically treated as loss-of-function. Absence of flagged variants is not proof of functional wild-type p53.",
        "pi3k_contract": "PIK3CA/PTEN/KRAS annotations are displayed as pathway context only until curated direction-aware biomarker rules are added.",
        "outputs": [str(OUTPUT.relative_to(ROOT)), str(MODEL_OUTPUT.relative_to(ROOT))],
    }
    MANIFEST.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")

    qc_rows = []
    for _, row in panels_out.iterrows():
        qc_rows.append({
            "panel_id": row["panel_id"],
            "candidate_hypotheses_n": int(row["candidate_hypotheses_n"]),
            "focus_lines_found_n": int(row["focus_lines_found_n"]),
            "focus_lines_mapped_n": int(row["focus_lines_mapped_n"]),
        })
    pd.DataFrame(qc_rows).to_csv(QC, sep="\t", index=False)

    print("MCL Laboratory Mechanism Context v1")
    print(f"Mechanism panels: {len(panels_out)}")
    for _, row in panels_out.iterrows():
        print(
            f"  {row['panel_id']}: hypotheses={int(row['candidate_hypotheses_n'])}; "
            f"focus lines={int(row['focus_lines_found_n'])}; mapped={int(row['focus_lines_mapped_n'])}"
        )
    print(f"Panel-model rows: {len(models_out)}")
    print(f"Wrote {OUTPUT.relative_to(ROOT)}")
    print(f"Wrote {MODEL_OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
