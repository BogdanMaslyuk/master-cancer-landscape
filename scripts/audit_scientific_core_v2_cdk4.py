from __future__ import annotations

"""Read-only provenance audit for six laboratory cell lines and CDK4.

Reads the existing processed/runtime snapshots; does not alter them. The
result is a diagnostic, not an estimate of PYZ activity or assay readiness.
"""

import argparse
import hashlib
import json
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
GENE = "CDK4"
LAB_NAMES = ("HCT 116", "MCF-7", "T-47D", "HuH-7", "JURKAT", "U-937")
INPUTS = {
    "laboratory": "data/processed/laboratory_panel.parquet",
    "atlas": "data/processed/depmap_crispr_model_atlas.parquet",
    "crispr": "data/processed/depmap_model_target_context.parquet",
    "responses": "data/runtime/pharmacology/responses.parquet",
    "links": "data/runtime/pharmacology/model_compound_target_links.parquet",
    "compounds": "data/runtime/pharmacology/compound_catalog.parquet",
    "target_evidence": "data/runtime/pharmacology/target_evidence.parquet",
    "hypotheses": "data/runtime/pharmacology/candidate_hypotheses_v2.parquet",
    "candidate_models": "data/runtime/pharmacology/candidate_hypothesis_models_v2.parquet",
}
OPTIONAL_INPUTS = {
    "wide_probability": "data/processed/depmap_model_gene_dependency.parquet",
    "wide_gene_effect": "data/processed/depmap_model_gene_effect.parquet",
}


def _norm(value: object) -> str:
    return re.sub(r"[^A-Z0-9]", "", str(value or "").upper())


def _txt(value: object) -> str | None:
    if value is None or pd.isna(value):
        return None
    word = str(value).strip()
    return word if word else None


def _hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read(path: Path, where: str = "", args: list[object] | None = None,
          columns: str = "*") -> pd.DataFrame:
    # DuckDB performs Parquet predicate pushdown; no mutation of the input.
    import duckdb

    quoted = path.resolve().as_posix().replace("'", "''")
    sql = f"SELECT {columns} FROM read_parquet('{quoted}')"
    if where:
        sql += " WHERE " + where
    con = duckdb.connect(database=":memory:")
    try:
        return con.execute(sql, args or []).df()
    finally:
        con.close()


def _where_ids(column: str, ids: list[str]) -> tuple[str, list[str]]:
    if not ids:
        return "FALSE", []
    return f"CAST({column} AS VARCHAR) IN (" + ",".join("?" for _ in ids) + ")", ids


def identify_lines(lab: pd.DataFrame, atlas: pd.DataFrame) -> tuple[pd.DataFrame, list[dict]]:
    """Never infer physical laboratory availability from atlas-name matching."""
    issues: list[dict] = []
    rows = []
    for name in LAB_NAMES:
        hits = lab[lab["lab_name"].map(_norm).eq(_norm(name))]
        ids = sorted(set(x for x in hits.get("model_id", pd.Series(dtype=object)).map(_txt) if x))
        if len(hits) == 0:
            status, mid = "not_in_lab_registry", None
        elif len(ids) == 1:
            status, mid = "mapped_in_registry_unverified_physical_availability", ids[0]
        elif len(ids) > 1:
            status, mid = "ambiguous_mapping", None
        else:
            status, mid = "not_mapped", None
        if status != "mapped_in_registry_unverified_physical_availability":
            issues.append({"kind": "laboratory_mapping", "entity": name, "detail": status})
        exact = atlas.loc[atlas["model_id"].astype(str).eq(mid)] if mid else atlas.iloc[:0]
        if len(exact) > 1:
            issues.append({"kind": "duplicate_atlas_model", "entity": name, "detail": mid})
        a = exact.iloc[0] if len(exact) == 1 else None
        rows.append({
            "lab_name": name, "model_id": mid, "laboratory_mapping": status,
            "atlas_name": a.get("cell_line_name") if a is not None else None,
            "cancer_id": a.get("mcl_cancer_id") if a is not None else None,
            "cancer_name": a.get("mcl_cancer_name") if a is not None else None,
            "source_registry_rows": len(hits),
            "physical_stock_verified": False,  # static list is not physical QC
        })
    return pd.DataFrame(rows), issues


def analyze_rows(
    lines: pd.DataFrame, context: pd.DataFrame, responses: pd.DataFrame,
    links: pd.DataFrame, compounds: pd.DataFrame,
    target_evidence: pd.DataFrame, hypotheses: pd.DataFrame,
    candidate_models: pd.DataFrame,
) -> dict[str, pd.DataFrame]:
    """Preserve source-record grain; never assign a PRISM measurement to PYZ."""
    issues: list[dict] = []
    known = set(lines["model_id"].dropna().astype(str))
    crispr_rows = []
    for _, line in lines.iterrows():
        mid = line["model_id"]
        hits = context.loc[
            context["model_id"].astype(str).eq(str(mid))
            & context["target_gene"].astype(str).str.upper().eq(GENE)
        ] if mid else context.iloc[:0]
        if len(hits) != 1:
            issues.append({
                "kind": "crispr_context_cardinality",
                "entity": str(line["lab_name"]),
                "detail": f"matching_rows={len(hits)}",
            })
        data = hits.iloc[0] if len(hits) == 1 else None
        prob = pd.to_numeric(pd.Series([data.get("dependency_probability") if data is not None else None]), errors="coerce").iloc[0]
        ge = pd.to_numeric(pd.Series([data.get("gene_effect") if data is not None else None]), errors="coerce").iloc[0]
        if pd.notna(prob) and not 0 <= float(prob) <= 1:
            issues.append({"kind": "invalid_dependency_probability", "entity": str(mid), "detail": str(prob)})
            prob = float("nan")
        crispr_rows.append({
            **line.to_dict(), "target_gene": GENE,
            "dependency_probability": prob, "gene_effect": ge,
            "genetic_status": ("unknown" if pd.isna(prob) else
                               "dependent" if float(prob) > 0.5 else "nondependent_measured"),
            "crispr_context_rows": len(hits),
            "depmap_release": data.get("depmap_release") if data is not None else None,
        })
    models = pd.DataFrame(crispr_rows)
    links = links.loc[
        links["target_gene"].astype(str).str.upper().eq(GENE)
        & links["model_id"].astype(str).isin(known)
    ].copy()
    links["model_id"] = links["model_id"].astype(str)
    links["compound_id"] = links["compound_id"].astype(str)
    link_keys = links[["model_id", "compound_id"]].drop_duplicates()
    if len(links) > len(link_keys):
        issues.append({
            "kind": "multiple_target_links",
            "entity": GENE,
            "detail": f"source_rows={len(links)}; distinct_model_compound_pairs={len(link_keys)}",
        })
    responses = responses.loc[
        responses["model_id"].astype(str).isin(known)
        & responses["source"].astype(str).str.upper().eq("PRISM")
        & responses["endpoint"].astype(str).str.upper().eq("LFC")
    ].copy()
    responses["model_id"] = responses["model_id"].astype(str)
    responses["compound_id"] = responses["compound_id"].astype(str)
    # An inner join with DISTINCT link keys avoids manufacturing extra
    # measurements from multiple annotations for the same compound.
    prism = responses.merge(link_keys, on=["model_id", "compound_id"], how="inner", validate="many_to_one")
    if "preferred_name" in compounds.columns:
        names = compounds[["compound_id", "preferred_name"]].copy()
        names["compound_id"] = names["compound_id"].astype(str)
        if names["compound_id"].duplicated().any():
            issues.append({"kind": "duplicate_compound_identifier", "entity": "compound_catalog", "detail": "Join suppressed"})
        else:
            prism = prism.merge(names, on="compound_id", how="left", validate="many_to_one")
    if "observation_id" in prism and prism["observation_id"].dropna().duplicated().any():
        issues.append({"kind": "duplicate_prism_observation_id", "entity": GENE, "detail": "Review all source rows"})
    target_evidence = target_evidence.loc[
        target_evidence["target_gene"].astype(str).str.upper().eq(GENE)
        & target_evidence["compound_id"].astype(str).isin(set(prism["compound_id"]))
    ].copy()
    # Candidate v2 is a truncated convenience layer, never our full source.
    relevant_hyp = hypotheses.loc[hypotheses["target_gene"].astype(str).str.upper().eq(GENE)].copy()
    candidate = candidate_models.merge(
        relevant_hyp[["hypothesis_id", "compound_id", "mcl_cancer_id", "priority_status"]].drop_duplicates("hypothesis_id"),
        on="hypothesis_id", how="inner", suffixes=("", "_hypothesis"), validate="many_to_one",
    )
    candidate = candidate.loc[candidate["model_id"].astype(str).isin(known)].copy()
    role_counts = (candidate.groupby("model_id")["model_role"].nunique()
                   if not candidate.empty and "model_role" in candidate else pd.Series(dtype=int))
    for mid, count in role_counts.items():
        if count > 1:
            issues.append({"kind": "multi_role_across_reference_compounds", "entity": str(mid), "detail": f"roles={count}"})
    for _, row in models.iterrows():
        if row["genetic_status"] == "unknown":
            issues.append({"kind": "missing_crispr_measurement", "entity": str(row["lab_name"]), "detail": GENE})
    return {
        "models": models,
        "prism_measurements": prism,
        "target_annotations": target_evidence,
        "candidate_v2_rows": candidate,
        "target_links": links,
        "issues": pd.DataFrame(issues, columns=["kind", "entity", "detail"]),
    }



def trace_role_filter(candidate: pd.DataFrame, prism: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Report exactly what the legacy priority filter excludes.

    Candidate v2 is still a truncated source. This reports only observed
    candidate rows and never substitutes them for all PRISM observations.
    """
    required = {"model_id", "compound_id", "model_role", "priority_status"}
    missing = required - set(candidate.columns)
    if missing:
        raise ValueError(f"Candidate v2 role trace requires columns: {sorted(missing)}")

    trace = candidate.copy()
    for key in ("model_id", "compound_id"):
        trace[key] = trace[key].astype(str)
    trace["included_by_legacy_priority_filter"] = trace["priority_status"].astype(str).eq(
        "priority_for_in_vitro"
    )
    trace["excluded_by_legacy_priority_filter"] = ~trace["included_by_legacy_priority_filter"]
    trace["role_is_compound_specific"] = True

    if not prism.empty and {"model_id", "compound_id"}.issubset(prism.columns):
        measured = prism.copy()
        measured["model_id"] = measured["model_id"].astype(str)
        measured["compound_id"] = measured["compound_id"].astype(str)
        grp = measured.groupby(["model_id", "compound_id"], dropna=False)
        counts = grp.size().reset_index(name="normalized_prism_records_n")
        trace = trace.merge(
            counts, on=["model_id", "compound_id"], how="left", validate="many_to_one"
        )
        if "preferred_name" in measured.columns:
            names = (
                measured.groupby(["model_id", "compound_id"], dropna=False)["preferred_name"]
                .first().reset_index(name="reference_compound_name")
            )
            trace = trace.merge(
                names, on=["model_id", "compound_id"], how="left", validate="many_to_one"
            )
    else:
        trace["normalized_prism_records_n"] = 0

    trace["normalized_prism_records_n"] = (
        pd.to_numeric(trace["normalized_prism_records_n"], errors="coerce")
        .fillna(0).astype(int)
    )
    grouping = ["model_id"]
    parts = []
    for model_id, group in trace.groupby(grouping, dropna=False):
        selected = group.loc[group["included_by_legacy_priority_filter"]]
        parts.append({
            "model_id": model_id,
            "candidate_v2_rows": len(group),
            "candidate_v2_distinct_compounds": group["compound_id"].nunique(),
            "candidate_positive_rows": int(group["model_role"].eq("positive").sum()),
            "candidate_discordant_rows": int(group["model_role"].astype(str).str.startswith("discordant").sum()),
            "legacy_included_rows": len(selected),
            "legacy_included_positive_rows": int(selected["model_role"].eq("positive").sum()),
            "legacy_excluded_rows": len(group) - len(selected),
            "legacy_excluded_positive_rows": int(
                group["excluded_by_legacy_priority_filter"].mul(
                    group["model_role"].eq("positive")
                ).sum()
            ),
        })
    impact = pd.DataFrame(parts, columns=[
        "model_id", "candidate_v2_rows", "candidate_v2_distinct_compounds",
        "candidate_positive_rows", "candidate_discordant_rows", "legacy_included_rows",
        "legacy_included_positive_rows", "legacy_excluded_rows",
        "legacy_excluded_positive_rows",
    ])
    return trace, impact

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=ROOT,
                        help="Read-only scientific-data checkout; can differ from the review worktree")
    parser.add_argument("--output-root", type=Path, default=ROOT / "outputs" / "audits" / "scientific_core_v2")
    args = parser.parse_args()
    source = args.data_root.resolve()
    paths = {key: source / rel for key, rel in INPUTS.items()}
    missing = [str(p) for p in paths.values() if not p.is_file()]
    if missing:
        raise SystemExit("Missing inputs (nothing was modified):\n" + "\n".join(missing))
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    directory = args.output_root.resolve() / f"CDK4_{timestamp}"
    if directory.exists():
        raise SystemExit(f"Output folder already exists; no overwrite: {directory}")

    lab = _read(paths["laboratory"])
    atlas = _read(paths["atlas"])
    lines, initial_issues = identify_lines(lab, atlas)
    ids = lines["model_id"].dropna().astype(str).unique().tolist()
    id_sql, id_args = _where_ids("model_id", ids)
    context = _read(paths["crispr"], f"UPPER(CAST(target_gene AS VARCHAR)) = ? AND {id_sql}",
                    [GENE, *id_args])
    links = _read(paths["links"], f"UPPER(CAST(target_gene AS VARCHAR)) = ? AND {id_sql}",
                  [GENE, *id_args])
    responses = _read(paths["responses"], id_sql, id_args)
    # Optional metadata sources are filtered after loading the small lookup.
    compounds = _read(paths["compounds"])
    known_compounds = sorted(set(links.get("compound_id", pd.Series(dtype=object)).dropna().astype(str)))
    cmp_sql, cmp_args = _where_ids("compound_id", known_compounds)
    evidence = _read(paths["target_evidence"], f"UPPER(CAST(target_gene AS VARCHAR)) = ? AND {cmp_sql}",
                     [GENE, *cmp_args])
    hyp = _read(paths["hypotheses"], "UPPER(CAST(target_gene AS VARCHAR)) = ?", [GENE])
    hyp_sql, hyp_args = _where_ids("hypothesis_id", hyp["hypothesis_id"].dropna().astype(str).tolist())
    candidates = _read(paths["candidate_models"], hyp_sql, hyp_args)
    result = analyze_rows(lines, context, responses, links, compounds, evidence, hyp, candidates)
    result["candidate_role_trace"], result["legacy_priority_filter_impact"] = trace_role_filter(
        result["candidate_v2_rows"], result["prism_measurements"]
    )
    if initial_issues:
        result["issues"] = pd.concat(
            [result["issues"], pd.DataFrame(initial_issues)], ignore_index=True
        )
    # Check the independent processed wide matrices, without overwriting context values.
    for source_key, value_column in (("wide_probability", "dependency_probability"),
                                     ("wide_gene_effect", "gene_effect")):
        path = source / OPTIONAL_INPUTS[source_key]
        if not path.exists():
            result["issues"] = pd.concat([result["issues"], pd.DataFrame([{
                "kind": "wide_matrix_missing", "entity": source_key, "detail": str(path),
            }])], ignore_index=True)
            continue
        try:
            wide = _read(path, id_sql, id_args, columns='model_id, "CDK4"')
            if GENE not in wide.columns:
                raise ValueError(f"{GENE} column not present")
            check = wide[["model_id", GENE]].rename(columns={GENE: f"{value_column}_wide"})
            result["models"] = result["models"].merge(check, on="model_id", how="left", validate="one_to_one")
            left = pd.to_numeric(result["models"][value_column], errors="coerce")
            right = pd.to_numeric(result["models"][f"{value_column}_wide"], errors="coerce")
            divergent = (left.notna() & right.notna() & (left - right).abs().gt(1e-5))
            for _, row in result["models"].loc[divergent].iterrows():
                result["issues"] = pd.concat([result["issues"], pd.DataFrame([{
                    "kind": "wide_context_disagreement", "entity": str(row["model_id"]),
                    "detail": value_column,
                }])], ignore_index=True)
        except (ValueError, KeyError) as exc:
            result["issues"] = pd.concat([result["issues"], pd.DataFrame([{
                "kind": "wide_matrix_unverifiable", "entity": source_key, "detail": str(exc),
            }])], ignore_index=True)

    inventory = []
    for label, path in {**paths, **{
        k: source / v for k, v in OPTIONAL_INPUTS.items()
    }}.items():
        inventory.append({
            "label": label, "relative_path": str(path.relative_to(source)),
            "exists": path.exists(), "bytes": path.stat().st_size if path.exists() else None,
            "sha256": _hash(path) if path.exists() else None,
        })
    git = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"],
                         capture_output=True, text=True, check=False)
    manifest = {
        "contract": "mcl-scientific-core-v2-cdk4-provenance-audit-v1",
        "audit_time_utc": datetime.now(timezone.utc).isoformat(),
        "code_commit": git.stdout.strip() if git.returncode == 0 else None,
        "data_root": str(source),
        "target": GENE, "laboratory_lines": list(LAB_NAMES),
        "note": "Audit of normalized source records, not measured PYZ response, not complete provenance to original raw PRISM.",
        "counts": {k: len(v) for k, v in result.items()},
        "files": inventory,
    }
    directory.mkdir(parents=True, exist_ok=False)
    for label, frame in result.items():
        frame.to_csv(directory / f"{label}.tsv", sep="\t", index=False, encoding="utf-8-sig")
    pd.DataFrame(inventory).to_csv(directory / "input_inventory.tsv", sep="\t", index=False)
    (directory / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print("MCL Scientific Core v2 CDK4 audit")
    print(f"Report: {directory}")
    print(f"Models: {len(result['models'])}; linked PRISM observations: {len(result['prism_measurements'])}")
    print(f"Candidate v2 records: {len(result['candidate_v2_rows'])}; issues: {len(result['issues'])}")
    print("Legacy priority filter impact:")
    print(result["legacy_priority_filter_impact"].to_string(index=False))
    print(result["models"][["lab_name", "model_id", "cancer_id",
                             "dependency_probability", "gene_effect",
                             "genetic_status"]].to_string(index=False))
    print("Audit only: no scientific source tables or existing analysis files changed.")


if __name__ == "__main__":
    main()
