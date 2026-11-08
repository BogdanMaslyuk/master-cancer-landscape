from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
import pyarrow.parquet as pq

try:
    import duckdb
except ImportError as exc:
    raise SystemExit(
        '.\\.venv\\Scripts\\python.exe -m pip install -e ".[analysis]" is required.'
    ) from exc


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config" / "laboratory_mechanism_panels.tsv"
PROCESSED = ROOT / "data" / "processed"
RAW_DEPMAP = ROOT / "data" / "raw" / "depmap"
PHARM = ROOT / "data" / "runtime" / "pharmacology"
RUNTIME = ROOT / "data" / "runtime" / "laboratory"

PANEL = PROCESSED / "laboratory_panel.parquet"
TARGET_CONTEXT = PROCESSED / "depmap_model_target_context.parquet"
GENE_EFFECT = PROCESSED / "depmap_model_gene_effect.parquet"
GENE_DEPENDENCY = PROCESSED / "depmap_model_gene_dependency.parquet"
EXPRESSION = PROCESSED / "depmap_model_expression.parquet"
COPY_NUMBER = PROCESSED / "depmap_model_copy_number.parquet"
MULTIOMICS_MANIFEST = PROCESSED / "depmap_model_multiomics_manifest.json"
CANDIDATES = RUNTIME / "laboratory_candidate_hypotheses.parquet"
CANDIDATE_MODELS = RUNTIME / "laboratory_candidate_models.parquet"
RESPONSES = PHARM / "responses.parquet"
OUTPUT = RUNTIME / "laboratory_mechanism_panels.parquet"
MODEL_OUTPUT = RUNTIME / "laboratory_mechanism_panel_models.parquet"
MANIFEST = RUNTIME / "laboratory_mechanism_panels_manifest.json"
QC = ROOT / "outputs" / "qc" / "laboratory_mechanism_panels_qc.tsv"

PRISM_ACTIVE_LFC_THRESHOLD = -1.0
SELECTIVE_PERCENTILE = 0.75
NONRESPONSIVE_PERCENTILE = 0.25
DEPENDENCY_PROBABILITY_THRESHOLD = 0.5


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


def _safe_bool(value: Any) -> bool:
    try:
        if value is None or pd.isna(value):
            return False
    except (TypeError, ValueError):
        pass
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"true", "1", "yes", "y", "t"}


def _resolve_release() -> str | None:
    if MULTIOMICS_MANIFEST.exists():
        payload = json.loads(MULTIOMICS_MANIFEST.read_text(encoding="utf-8"))
        release = _text(payload.get("depmap_release"))
        if release:
            return release
    if RAW_DEPMAP.exists():
        releases = sorted((p.name for p in RAW_DEPMAP.iterdir() if p.is_dir()), reverse=True)
        if releases:
            return releases[0]
    return None


def _parquet_genes(path: Path) -> set[str]:
    if not path.exists():
        return set()
    return set(pq.ParquetFile(path).schema.names) - {"model_id"}


def _wide_long(path: Path, model_ids: list[str], genes: list[str], value_name: str) -> pd.DataFrame:
    if not path.exists() or not model_ids or not genes:
        return pd.DataFrame(columns=["model_id", "target_gene", value_name])
    available = _parquet_genes(path)
    keep = [gene for gene in genes if gene in available]
    if not keep:
        return pd.DataFrame(columns=["model_id", "target_gene", value_name])
    frame = pd.read_parquet(path, columns=["model_id", *keep])
    frame["model_id"] = frame["model_id"].astype(str)
    frame = frame[frame["model_id"].isin(set(model_ids))]
    if frame.empty:
        return pd.DataFrame(columns=["model_id", "target_gene", value_name])
    out = frame.melt(id_vars="model_id", var_name="target_gene", value_name=value_name)
    out["target_gene"] = out["target_gene"].astype(str).str.upper()
    out[value_name] = pd.to_numeric(out[value_name], errors="coerce")
    return out


def _pick_column(header: list[str], *candidates: str) -> str | None:
    normalized = {str(c).strip().lower().replace("_", ""): str(c) for c in header}
    for candidate in candidates:
        key = candidate.strip().lower().replace("_", "")
        if key in normalized:
            return normalized[key]
    return None


def _load_mutations(model_ids: list[str], genes: list[str]) -> pd.DataFrame:
    columns = ["model_id", "target_gene", "has_hotspot", "has_likely_lof", "protein_changes_json"]
    release = _resolve_release()
    if not release:
        return pd.DataFrame(columns=columns)
    release_dir = RAW_DEPMAP / release
    source = next(
        (
            release_dir / name
            for name in ("OmicsSomaticMutations.csv", "OmicsSomaticMutationsProfile.csv")
            if (release_dir / name).exists()
        ),
        None,
    )
    if source is None:
        return pd.DataFrame(columns=columns)

    header = pd.read_csv(source, nrows=0).columns.tolist()
    model_col = _pick_column(header, "ModelID", "DepMap_ID", "DepMapID")
    gene_col = _pick_column(header, "HugoSymbol", "Hugo_Symbol", "gene", "Gene")
    if not model_col or not gene_col:
        return pd.DataFrame(columns=columns)
    protein_col = _pick_column(header, "ProteinChange", "HGVSp", "VepHGVSp")
    hotspot_col = _pick_column(header, "Hotspot")
    lof_col = _pick_column(header, "LikelyLoF", "LikelyLof")
    usecols = list(dict.fromkeys([model_col, gene_col, *[x for x in (protein_col, hotspot_col, lof_col) if x]]))
    model_set = set(model_ids)
    gene_set = {g.upper() for g in genes}
    hits: list[pd.DataFrame] = []
    for chunk in pd.read_csv(source, usecols=usecols, chunksize=100_000, low_memory=False):
        mids = chunk[model_col].map(_text)
        symbols = chunk[gene_col].map(lambda x: _text(x).upper())
        mask = mids.isin(model_set) & symbols.isin(gene_set)
        if not mask.any():
            continue
        part = chunk.loc[mask].copy()
        part["model_id"] = mids.loc[mask]
        part["target_gene"] = symbols.loc[mask]
        hits.append(part)
    if not hits:
        return pd.DataFrame(columns=columns)

    frame = pd.concat(hits, ignore_index=True, sort=False)
    rows: list[dict[str, Any]] = []
    for (model_id, gene), group in frame.groupby(["model_id", "target_gene"], sort=False):
        changes: list[str] = []
        if protein_col:
            for value in group[protein_col]:
                text = _text(value)
                if text and text not in changes:
                    changes.append(text)
        rows.append({
            "model_id": str(model_id),
            "target_gene": str(gene).upper(),
            "has_hotspot": bool(group[hotspot_col].map(_safe_bool).any()) if hotspot_col else False,
            "has_likely_lof": bool(group[lof_col].map(_safe_bool).any()) if lof_col else False,
            "protein_changes_json": json.dumps(changes, ensure_ascii=False),
        })
    return pd.DataFrame(rows, columns=columns)


def _load_target_context_fallback(model_ids: list[str], genes: list[str]) -> pd.DataFrame:
    wanted = [
        "model_id", "target_gene", "expression_percentile_in_cancer",
        "copy_number_percentile_in_cancer", "has_hotspot", "has_likely_lof",
        "protein_changes_json",
    ]
    if not TARGET_CONTEXT.exists() or not model_ids or not genes:
        return pd.DataFrame(columns=wanted)
    con = duckdb.connect(database=":memory:")
    con.register("selected_models", pd.DataFrame({"model_id": sorted(set(model_ids))}))
    con.register("selected_genes", pd.DataFrame({"target_gene": sorted(set(g.upper() for g in genes))}))
    path = TARGET_CONTEXT.resolve().as_posix().replace("'", "''")
    available = set(
        con.execute(f"DESCRIBE SELECT * FROM read_parquet('{path}')").df()["column_name"].astype(str)
    )
    expressions = [f"c.{col}" if col in available else f"NULL AS {col}" for col in wanted]
    frame = con.execute(
        f"""
        SELECT {', '.join(expressions)}
        FROM read_parquet('{path}') c
        INNER JOIN selected_models m ON CAST(c.model_id AS VARCHAR) = m.model_id
        INNER JOIN selected_genes g ON upper(CAST(c.target_gene AS VARCHAR)) = g.target_gene
        """
    ).df()
    con.close()
    if not frame.empty:
        frame["model_id"] = frame["model_id"].astype(str)
        frame["target_gene"] = frame["target_gene"].astype(str).str.upper()
    return frame


def _coalesce_columns(frame: pd.DataFrame, primary: str, fallback: str, output: str) -> None:
    a = frame[primary] if primary in frame.columns else pd.Series(pd.NA, index=frame.index)
    b = frame[fallback] if fallback in frame.columns else pd.Series(pd.NA, index=frame.index)
    frame[output] = a.where(a.notna(), b)


def _load_context(model_ids: list[str], genes: list[str]) -> pd.DataFrame:
    model_ids = sorted(set(model_ids))
    genes = sorted(set(g.upper() for g in genes))
    if not model_ids or not genes:
        return pd.DataFrame()

    base = pd.MultiIndex.from_product(
        [model_ids, genes], names=["model_id", "target_gene"]
    ).to_frame(index=False)
    layers = [
        _wide_long(GENE_EFFECT, model_ids, genes, "gene_effect"),
        _wide_long(GENE_DEPENDENCY, model_ids, genes, "dependency_probability"),
        _wide_long(EXPRESSION, model_ids, genes, "expression_log2_tpm1"),
        _wide_long(COPY_NUMBER, model_ids, genes, "copy_number_relative"),
    ]
    out = base
    for layer in layers:
        out = out.merge(layer, on=["model_id", "target_gene"], how="left")

    fallback = _load_target_context_fallback(model_ids, genes)
    if not fallback.empty:
        fallback = fallback.rename(columns={
            "has_hotspot": "fallback_has_hotspot",
            "has_likely_lof": "fallback_has_likely_lof",
            "protein_changes_json": "fallback_protein_changes_json",
        })
        out = out.merge(fallback, on=["model_id", "target_gene"], how="left")
    else:
        out["expression_percentile_in_cancer"] = pd.NA
        out["copy_number_percentile_in_cancer"] = pd.NA
        out["fallback_has_hotspot"] = pd.NA
        out["fallback_has_likely_lof"] = pd.NA
        out["fallback_protein_changes_json"] = pd.NA

    mutations = _load_mutations(model_ids, genes)
    if not mutations.empty:
        mutations = mutations.rename(columns={
            "has_hotspot": "raw_has_hotspot",
            "has_likely_lof": "raw_has_likely_lof",
            "protein_changes_json": "raw_protein_changes_json",
        })
        out = out.merge(mutations, on=["model_id", "target_gene"], how="left")
    else:
        out["raw_has_hotspot"] = pd.NA
        out["raw_has_likely_lof"] = pd.NA
        out["raw_protein_changes_json"] = pd.NA

    _coalesce_columns(out, "raw_has_hotspot", "fallback_has_hotspot", "has_hotspot")
    _coalesce_columns(out, "raw_has_likely_lof", "fallback_has_likely_lof", "has_likely_lof")
    _coalesce_columns(out, "raw_protein_changes_json", "fallback_protein_changes_json", "protein_changes_json")
    out["has_hotspot"] = out["has_hotspot"].fillna(False).map(_safe_bool)
    out["has_likely_lof"] = out["has_likely_lof"].fillna(False).map(_safe_bool)
    out["protein_changes_json"] = out["protein_changes_json"].fillna("[]")
    return out


def _load_responses(compound_ids: list[str], model_ids: list[str]) -> pd.DataFrame:
    columns = ["model_id", "compound_id", "response_value", "relative_sensitivity"]
    if not RESPONSES.exists() or not compound_ids or not model_ids:
        return pd.DataFrame(columns=columns)
    con = duckdb.connect(database=":memory:")
    con.register("selected_compounds", pd.DataFrame({"compound_id": sorted(set(compound_ids))}))
    con.register("selected_models", pd.DataFrame({"model_id": sorted(set(model_ids))}))
    path = RESPONSES.resolve().as_posix().replace("'", "''")
    frame = con.execute(
        f"""
        WITH response_by_model AS (
          SELECT CAST(r.model_id AS VARCHAR) AS model_id,
                 CAST(r.compound_id AS VARCHAR) AS compound_id,
                 avg(CAST(r.value AS DOUBLE)) AS response_value
          FROM read_parquet('{path}') r
          INNER JOIN selected_compounds c
            ON CAST(r.compound_id AS VARCHAR) = c.compound_id
          WHERE upper(CAST(r.source AS VARCHAR)) = 'PRISM'
            AND upper(CAST(r.endpoint AS VARCHAR)) = 'LFC'
            AND r.value IS NOT NULL
          GROUP BY model_id, compound_id
        ), ranked AS (
          SELECT *,
                 1.0 - percent_rank() OVER (
                   PARTITION BY compound_id ORDER BY response_value ASC
                 ) AS relative_sensitivity
          FROM response_by_model
        )
        SELECT r.*
        FROM ranked r
        INNER JOIN selected_models m ON r.model_id = m.model_id
        """
    ).df()
    con.close()
    return frame if not frame.empty else pd.DataFrame(columns=columns)


def _classify_model_role(response: pd.Series, target_context: pd.Series, is_control: bool) -> str:
    if is_control:
        return "general_non_tumor_control"
    response_value = _clean_scalar(response.get("response_value"))
    relative = _clean_scalar(response.get("relative_sensitivity"))
    probability = _clean_scalar(target_context.get("dependency_probability"))
    if response_value is None:
        return "not_evaluated_no_prism"
    if probability is None:
        return "unclassified_missing_dependency_probability"
    response_value = float(response_value)
    relative = float(relative) if relative is not None else None
    probability = float(probability)
    active = response_value <= PRISM_ACTIVE_LFC_THRESHOLD
    priority_sensitive = bool(
        active and relative is not None and relative >= SELECTIVE_PERCENTILE
    )
    clearly_nonresponsive = bool(
        response_value > PRISM_ACTIVE_LFC_THRESHOLD
        and relative is not None
        and relative <= NONRESPONSIVE_PERCENTILE
    )
    dependent = probability > DEPENDENCY_PROBABILITY_THRESHOLD
    if priority_sensitive and dependent:
        return "positive"
    if clearly_nonresponsive and not dependent:
        return "negative_panel_comparator"
    if priority_sensitive and not dependent:
        return "discordant_sensitive_without_dependency"
    if not active and dependent:
        return "discordant_dependency_without_activity"
    return "unclassified"


def _context_status(panel_id: str, rows: pd.DataFrame) -> tuple[str, str]:
    if rows.empty:
        return "context_missing", "Молекулярный контекст для этой модели в текущем индексе отсутствует."

    panel_key = panel_id.upper()
    if "MDM2" in panel_key:
        tp53 = rows[rows["target_gene"].astype(str).str.upper().eq("TP53")]
        if tp53.empty:
            return "tp53_context_missing", "TP53 не найден в текущем индексированном контексте модели."
        row = tp53.iloc[0]
        if _safe_bool(row.get("has_likely_lof")):
            return (
                "tp53_likely_lof_concern",
                "Для TP53 есть флаг LikelyLoF. Это механистическое предостережение для ингибирования MDM2 и требует отдельной проверки p53-функции.",
            )
        if _safe_bool(row.get("has_hotspot")):
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
        missing: list[str] = []
        for gene in ("PIK3CA", "PTEN", "KRAS"):
            hit = rows[rows["target_gene"].astype(str).str.upper().eq(gene)]
            if hit.empty:
                missing.append(gene)
                continue
            row = hit.iloc[0]
            flags: list[str] = []
            if _safe_bool(row.get("has_hotspot")):
                flags.append("hotspot")
            if _safe_bool(row.get("has_likely_lof")):
                flags.append("LikelyLoF")
            if flags:
                notes.append(f"{gene}: {', '.join(flags)}")
        suffix = "; ".join(notes) if notes else "выделенных hotspot/LikelyLoF-флагов нет"
        if missing:
            suffix += "; нет данных: " + ", ".join(missing)
        return (
            "pi3k_context_available_direction_unresolved",
            "PI3K-контекст доступен; " + suffix + ". Эти признаки показываются описательно и не получают автоматического направления чувствительности/резистентности.",
        )

    return "context_available_direction_unresolved", "Молекулярный контекст доступен, но его направление не задано автоматическим правилом."


def _scope_mask(focus: pd.DataFrame, tokens: list[str]) -> pd.Series:
    if not tokens:
        return pd.Series(False, index=focus.index)
    combined = pd.Series("", index=focus.index, dtype="object")
    for column in (
        "mcl_cancer_name", "oncotree_lineage", "model_oncotree_lineage",
        "oncotree_primary_disease", "model_oncotree_primary_disease",
    ):
        if column in focus.columns:
            combined = combined + " " + focus[column].fillna("").astype(str)
    combined = combined.str.casefold()
    mask = pd.Series(False, index=focus.index)
    for token in tokens:
        mask |= combined.str.contains(token.casefold(), regex=False)
    return mask


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
    selected_compound_ids: set[str] = set()

    for _, spec in config.iterrows():
        panel_id = _text(spec["panel_id"])
        target_gene = _text(spec["target_gene"]).upper()
        compounds = {_norm(x) for x in _split(spec["compound_names"])}
        cancer_tokens = [x.casefold() for x in _split(spec.get("cancer_name_contains", ""))]
        context_genes = [x.upper() for x in _split(spec["context_genes"])]
        focus_names = _split(spec["focus_lab_lines"])
        focus_norms = {_norm(x) for x in focus_names}

        focus = panel[panel["_lab_norm"].isin(focus_norms)].copy()
        non_control_focus = focus[
            ~focus.get("control_role", pd.Series("", index=focus.index)).astype(str).eq("general_non_tumor_control")
        ].copy()
        focus_cancer_ids = set(
            non_control_focus.get("mcl_cancer_id", pd.Series(dtype=object)).dropna().astype(str)
        )

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
        elif focus_cancer_ids:
            selected = selected[selected["mcl_cancer_id"].astype(str).isin(focus_cancer_ids)]

        for value in focus["model_id"].dropna().astype(str):
            if value.strip():
                selected_model_ids.add(value.strip())
        selected_genes.update(context_genes)
        selected_genes.add(target_gene)
        selected_compound_ids.update(selected.get("compound_id", pd.Series(dtype=object)).dropna().astype(str))

        selected_specs.append({
            "spec": spec,
            "panel_id": panel_id,
            "target_gene": target_gene,
            "context_genes": context_genes,
            "cancer_tokens": cancer_tokens,
            "focus": focus,
            "selected": selected,
        })

    context = _load_context(sorted(selected_model_ids), sorted(selected_genes))
    responses = _load_responses(sorted(selected_compound_ids), sorted(selected_model_ids))

    panel_rows: list[dict[str, Any]] = []
    model_rows: list[dict[str, Any]] = []

    for bundle in selected_specs:
        spec = bundle["spec"]
        panel_id = bundle["panel_id"]
        target_gene = bundle["target_gene"]
        context_genes = bundle["context_genes"]
        cancer_tokens = bundle["cancer_tokens"]
        focus = bundle["focus"]
        selected = bundle["selected"]
        panel_model_buffer: list[dict[str, Any]] = []

        hypothesis_ids = set(selected["hypothesis_id"].astype(str)) if not selected.empty else set()
        evidence = model_evidence[
            model_evidence["hypothesis_id"].astype(str).isin(hypothesis_ids)
        ].copy() if hypothesis_ids else pd.DataFrame()

        for _, hypothesis in selected.iterrows():
            cancer_id = _text(hypothesis.get("mcl_cancer_id"))
            control_mask = focus.get("control_role", pd.Series("", index=focus.index)).astype(str).eq("general_non_tumor_control")
            if cancer_tokens:
                eligible_focus = focus[control_mask | _scope_mask(focus, cancer_tokens)].copy()
            else:
                cancer_mask = focus.get("mcl_cancer_id", pd.Series("", index=focus.index)).fillna("").astype(str).eq(cancer_id)
                eligible_focus = focus[control_mask | cancer_mask].copy()

            for _, line in eligible_focus.iterrows():
                model_id = _text(line.get("model_id"))
                lab_name = _text(line.get("lab_name"))
                is_control = _text(line.get("control_role")) == "general_non_tumor_control"
                candidate_hit = evidence[
                    evidence["hypothesis_id"].astype(str).eq(str(hypothesis["hypothesis_id"]))
                    & evidence["lab_name"].astype(str).map(_norm).eq(_norm(lab_name))
                ] if not evidence.empty else pd.DataFrame()
                candidate_role = _text(candidate_hit.iloc[0].get("laboratory_model_role")) if not candidate_hit.empty else None

                ctx = context[
                    context["model_id"].astype(str).eq(model_id)
                    & context["target_gene"].astype(str).str.upper().isin(context_genes)
                ].copy() if model_id and not context.empty else pd.DataFrame()
                target_ctx = ctx[ctx["target_gene"].astype(str).str.upper().eq(target_gene)]
                target_row = target_ctx.iloc[0] if not target_ctx.empty else pd.Series(dtype=object)
                response_hit = responses[
                    responses["model_id"].astype(str).eq(model_id)
                    & responses["compound_id"].astype(str).eq(_text(hypothesis.get("compound_id")))
                ] if model_id and not responses.empty else pd.DataFrame()
                response_row = response_hit.iloc[0] if not response_hit.empty else pd.Series(dtype=object)
                model_role = _classify_model_role(response_row, target_row, is_control)

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
                            "has_hotspot": _safe_bool(c.get("has_hotspot")),
                            "has_likely_lof": _safe_bool(c.get("has_likely_lof")),
                            "protein_changes_json": _clean_scalar(c.get("protein_changes_json")),
                        })
                status, status_ru = _context_status(panel_id, ctx)

                row = {
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
                    "laboratory_model_role": model_role,
                    "candidate_triage_role": candidate_role,
                    "response_value": _clean_scalar(response_row.get("response_value")),
                    "relative_sensitivity": _clean_scalar(response_row.get("relative_sensitivity")),
                    "dependency_probability": _clean_scalar(target_row.get("dependency_probability")),
                    "gene_effect": _clean_scalar(target_row.get("gene_effect")),
                    "mechanism_context_status": status,
                    "mechanism_context_summary_ru": status_ru,
                    "mechanism_context_json": _json(ctx_records),
                    "built_at": _now(),
                }
                model_rows.append(row)
                panel_model_buffer.append(row)

        panel_models = pd.DataFrame(panel_model_buffer)
        positive_lines = sorted(set(
            panel_models.loc[panel_models["laboratory_model_role"].eq("positive"), "lab_name"].dropna().astype(str)
        )) if not panel_models.empty else []
        negative_lines = sorted(set(
            panel_models.loc[panel_models["laboratory_model_role"].eq("negative_panel_comparator"), "lab_name"].dropna().astype(str)
        )) if not panel_models.empty else []
        discordant_lines = sorted(set(
            panel_models.loc[panel_models["laboratory_model_role"].astype(str).str.startswith("discordant"), "lab_name"].dropna().astype(str)
        )) if not panel_models.empty else []

        mapped_focus_ids = set(focus["model_id"].dropna().astype(str)) if not focus.empty else set()
        panel_context = context[context["model_id"].astype(str).isin(mapped_focus_ids)].copy() if not context.empty else pd.DataFrame()
        observed_genes: list[str] = []
        for gene in context_genes:
            hit = panel_context[panel_context["target_gene"].astype(str).str.upper().eq(gene)] if not panel_context.empty else pd.DataFrame()
            if not hit.empty:
                quantitative = [
                    "gene_effect", "dependency_probability", "expression_log2_tpm1", "copy_number_relative"
                ]
                observed = any(
                    column in hit.columns and hit[column].notna().any()
                    for column in quantitative
                ) or bool(hit.get("has_hotspot", pd.Series(False, index=hit.index)).map(_safe_bool).any()) \
                  or bool(hit.get("has_likely_lof", pd.Series(False, index=hit.index)).map(_safe_bool).any())
                if observed:
                    observed_genes.append(gene)
        missing_context_genes = [g for g in context_genes if g not in observed_genes]

        output_focus_lines = set(panel_models["lab_name"].dropna().astype(str)) if not panel_models.empty else set()
        expected_focus_lines = set(focus["lab_name"].dropna().astype(str)) if not focus.empty else set()
        omitted_focus_lines = sorted(expected_focus_lines - output_focus_lines)

        candidate_names = sorted(set(selected["preferred_name"].dropna().astype(str))) if not selected.empty else []
        candidate_cancers = sorted(set(selected["mcl_cancer_name"].dropna().astype(str))) if not selected.empty else []
        panel_rows.append({
            "panel_id": panel_id,
            "panel_name_ru": _text(spec["panel_name_ru"]),
            "target_gene": target_gene,
            "compound_names_json": _json(_split(spec["compound_names"])),
            "context_genes_json": _json(context_genes),
            "context_genes_observed_json": _json(observed_genes),
            "context_genes_missing_json": _json(missing_context_genes),
            "focus_lab_lines_json": _json(_split(spec["focus_lab_lines"])),
            "omitted_focus_lines_json": _json(omitted_focus_lines),
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
            "focus_lines_in_output_n": int(len(output_focus_lines)),
            "built_at": _now(),
        })

    panels_out = pd.DataFrame(panel_rows)
    models_out = pd.DataFrame(model_rows)
    RUNTIME.mkdir(parents=True, exist_ok=True)
    QC.parent.mkdir(parents=True, exist_ok=True)
    panels_out.to_parquet(OUTPUT, index=False, compression="zstd")
    models_out.to_parquet(MODEL_OUTPUT, index=False, compression="zstd")

    manifest = {
        "contract": "mcl-laboratory-mechanism-context-v1.1",
        "built_at": _now(),
        "panels_n": int(len(panels_out)),
        "panel_model_rows_n": int(len(models_out)),
        "context_sources": [
            str(path.relative_to(ROOT))
            for path in (GENE_EFFECT, GENE_DEPENDENCY, EXPRESSION, COPY_NUMBER, TARGET_CONTEXT)
            if path.exists()
        ],
        "response_source": str(RESPONSES.relative_to(ROOT)) if RESPONSES.exists() else None,
        "context_contract": (
            "Mechanism context genes are loaded on demand from the pinned DepMap gene-level multi-omics matrices. "
            "They are not restricted to the pharmacological target catalog. Curated panel membership can span multiple "
            "MCL cancer IDs inside one disease family; exact cancer-ID matching is therefore not used when a panel has an explicit cancer_name_contains scope."
        ),
        "ranking_contract": "No composite score. Mechanism context is explanatory and can raise a concern, but does not automatically upgrade a candidate.",
        "tp53_contract": "TP53 LikelyLoF is a mechanistic concern for MDM2 inhibition. Hotspot is not automatically treated as loss-of-function. Absence of flagged variants is not proof of functional wild-type p53.",
        "pi3k_contract": "PIK3CA/PTEN/KRAS annotations are displayed as pathway context only until curated direction-aware biomarker rules are added.",
        "outputs": [str(OUTPUT.relative_to(ROOT)), str(MODEL_OUTPUT.relative_to(ROOT))],
    }
    MANIFEST.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")

    qc_rows: list[dict[str, Any]] = []
    for _, row in panels_out.iterrows():
        qc_rows.append({
            "panel_id": row["panel_id"],
            "candidate_hypotheses_n": int(row["candidate_hypotheses_n"]),
            "focus_lines_found_n": int(row["focus_lines_found_n"]),
            "focus_lines_mapped_n": int(row["focus_lines_mapped_n"]),
            "focus_lines_in_output_n": int(row["focus_lines_in_output_n"]),
            "context_genes_observed_json": row["context_genes_observed_json"],
            "context_genes_missing_json": row["context_genes_missing_json"],
            "omitted_focus_lines_json": row["omitted_focus_lines_json"],
        })
    pd.DataFrame(qc_rows).to_csv(QC, sep="\t", index=False)

    print("MCL Laboratory Mechanism Context v1.1")
    print(f"Mechanism panels: {len(panels_out)}")
    for _, row in panels_out.iterrows():
        observed = ", ".join(json.loads(row["context_genes_observed_json"])) or "none"
        missing_genes = ", ".join(json.loads(row["context_genes_missing_json"])) or "none"
        omitted = ", ".join(json.loads(row["omitted_focus_lines_json"])) or "none"
        print(
            f"  {row['panel_id']}: hypotheses={int(row['candidate_hypotheses_n'])}; "
            f"focus lines={int(row['focus_lines_found_n'])}; mapped={int(row['focus_lines_mapped_n'])}; "
            f"in output={int(row['focus_lines_in_output_n'])}"
        )
        print(f"    context genes observed: {observed}; missing: {missing_genes}")
        print(f"    omitted focus lines: {omitted}")
    print(f"Panel-model rows: {len(models_out)}")
    print(f"Wrote {OUTPUT.relative_to(ROOT)}")
    print(f"Wrote {MODEL_OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
