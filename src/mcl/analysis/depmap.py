from __future__ import annotations

import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import yaml
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

from mcl.models.depmap import (
    DepMapContextAuditRow,
    DepMapEvidence,
    DepMapKrasSensitivityAuditRow,
    DepMapKrasSensitivityEvidence,
)
from mcl.models.provenance import ProvenanceRecord
from mcl.models.qc import QCRecord
from mcl.sources.depmap import (
    depmap_file_inventory,
    load_gene_matrix,
    load_model_metadata,
    load_relevant_mutations,
    load_sequenced_model_ids,
    resolve_depmap_files,
    row_has_driver_signal,
)
from mcl.utils.hash import sha256_file


PARSER_VERSION = "depmap-m3-v0.2"
KRAS_SENSITIVITY_PARSER_VERSION = "depmap-m3.1-v0.1"
OFFICIAL_SOURCE = "https://depmap.org/portal/data_page/"


def _rel_or_abs(path: Path, root: Path) -> str:
    try:
        return str(path.resolve().relative_to(root.resolve()))
    except ValueError:
        return str(path.resolve())


def _norm(value: object) -> str:
    return re.sub(r"\s+", " ", str(value or "").strip()).casefold()


def load_context_config(path: Path) -> dict[str, dict]:
    cfg = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    contexts = cfg.get("contexts", {})
    return {str(k): v for k, v in contexts.items()}


def _metadata_matches(row: pd.Series, depmap_cfg: dict) -> bool:
    """Match a model to the disease component of a context using explicit exact aliases."""
    checks: list[bool] = []
    mappings = {
        "oncotree_codes": "OncotreeCode",
        "depmap_model_types": "DepmapModelType",
        "oncotree_subtypes": "OncotreeSubtype",
        "oncotree_primary_diseases": "OncotreePrimaryDisease",
    }
    for cfg_key, col in mappings.items():
        allowed = depmap_cfg.get("metadata", {}).get(cfg_key, []) or []
        if allowed and col in row.index:
            checks.append(_norm(row.get(col)) in {_norm(x) for x in allowed})
    # Lineage is a guardrail, not sufficient by itself to define the disease.
    lineages = depmap_cfg.get("metadata", {}).get("oncotree_lineages", []) or []
    lineage_ok = True
    if lineages and "OncotreeLineage" in row.index:
        lineage_ok = _norm(row.get("OncotreeLineage")) in {_norm(x) for x in lineages}
    return bool(checks) and any(checks) and lineage_ok


def _variant_payload(row: pd.Series) -> dict:
    keep = [
        "GeneSymbol",
        "ProteinChange",
        "VariantClassification",
        "Hotspot",
        "HessDriver",
        "LikelyDriver",
        "Oncogenic",
        "CosmicHotspot",
        "VariantInfo",
        "Description",
        "DNAChange",
        "CDSChange",
    ]
    return {k: str(row.get(k, "")) for k in keep if k in row.index and str(row.get(k, ""))}


def _regex_match(value: object, patterns: Iterable[str]) -> bool:
    text = str(value or "")
    return any(re.search(pattern, text, flags=re.IGNORECASE) for pattern in patterns)


def _combined_text(row: pd.Series, columns: Iterable[str]) -> str:
    return " | ".join(str(row.get(c, "")) for c in columns if str(row.get(c, "")))


def _classify_alteration(
    model_id: str,
    mutations_by_model: dict[str, pd.DataFrame],
    mutation_cfg: dict,
    driver_columns: list[str],
    text_columns: list[str],
) -> tuple[str, list[dict], list[dict]]:
    """Return alteration status and auditable variant lists.

    Statuses:
      defining_present, defining_absent, ambiguous
    Sequencing availability is handled separately before this function.
    """
    rows = mutations_by_model.get(model_id)
    if rows is None:
        rows = pd.DataFrame()
    mode = str(mutation_cfg.get("mode", "")).strip()

    if mode == "exact_protein_change_present":
        gene = str(mutation_cfg["gene"]).upper()
        patterns = mutation_cfg.get("protein_change_regexes", [])
        gene_rows = rows[rows.get("GeneSymbol", pd.Series(dtype=str)).astype(str).str.upper() == gene] if not rows.empty else rows
        qualifying, other = [], []
        for _, row in gene_rows.iterrows():
            payload = _variant_payload(row)
            if _regex_match(row.get("ProteinChange", ""), patterns):
                qualifying.append(payload)
            else:
                other.append(payload)
        if qualifying:
            return "defining_present", qualifying, other
        if other and bool(mutation_cfg.get("other_variants_ambiguous", False)):
            return "ambiguous", [], other
        return "defining_absent", [], other

    if mode == "absence_of_defining_variants":
        defining_rules = mutation_cfg.get("defining_variants", []) or []
        relevant_genes = {str(x["gene"]).upper() for x in defining_rules}
        relevant_rows = rows[
            rows.get("GeneSymbol", pd.Series(dtype=str)).astype(str).str.upper().isin(relevant_genes)
        ] if not rows.empty else rows
        qualifying, other = [], []
        for _, row in relevant_rows.iterrows():
            gene = str(row.get("GeneSymbol", "")).upper()
            matched = False
            for rule in defining_rules:
                if gene == str(rule["gene"]).upper() and _regex_match(
                    row.get("ProteinChange", ""), rule.get("protein_change_regexes", [])
                ):
                    matched = True
                    break
            (qualifying if matched else other).append(_variant_payload(row))
        if qualifying:
            return "defining_present", qualifying, other
        if other and bool(mutation_cfg.get("other_relevant_variants_ambiguous", True)):
            return "ambiguous", [], other
        return "defining_absent", [], other

    if mode == "driver_or_hotspot_gene_variant":
        gene = str(mutation_cfg["gene"]).upper()
        protein_patterns = mutation_cfg.get("protein_change_regexes", []) or []
        text_patterns = mutation_cfg.get("text_regexes", []) or []
        gene_rows = rows[rows.get("GeneSymbol", pd.Series(dtype=str)).astype(str).str.upper() == gene] if not rows.empty else rows
        qualifying, other = [], []
        for _, row in gene_rows.iterrows():
            by_driver = row_has_driver_signal(row, driver_columns)
            by_protein = _regex_match(row.get("ProteinChange", ""), protein_patterns)
            by_text = _regex_match(_combined_text(row, text_columns), text_patterns)
            payload = _variant_payload(row)
            (qualifying if (by_driver or by_protein or by_text) else other).append(payload)
        if qualifying:
            return "defining_present", qualifying, other
        if other:
            return "ambiguous", [], other
        return "defining_absent", [], []

    raise ValueError(f"Unsupported DepMap mutation mode: {mode!r}")


def build_context_membership(
    models: pd.DataFrame,
    mutations: pd.DataFrame,
    sequenced_models: set[str],
    contexts: dict[str, dict],
    mutation_schema: dict[str, object],
) -> tuple[dict[str, dict[str, list[str]]], list[DepMapContextAuditRow]]:
    if "ModelID" not in models.columns:
        raise ValueError("Model metadata must contain ModelID")
    models = models.copy()
    models["ModelID"] = models["ModelID"].astype(str)
    if models["ModelID"].duplicated().any():
        dups = models.loc[models["ModelID"].duplicated(keep=False), "ModelID"].unique().tolist()
        raise ValueError(f"Model metadata contains duplicate ModelID values: {dups[:10]}")

    mutations_by_model = {
        str(mid): group.copy()
        for mid, group in mutations.groupby("ModelID", sort=False)
    } if not mutations.empty else {}
    driver_cols = [str(x) for x in mutation_schema.get("driver_flags", [])]
    text_cols = [str(x) for x in mutation_schema.get("text_columns", [])]

    memberships: dict[str, dict[str, list[str]]] = {}
    audit: list[DepMapContextAuditRow] = []
    for cancer_id, cfg in contexts.items():
        depmap_cfg = cfg.get("depmap", {})
        if not depmap_cfg:
            continue
        mutation_cfg = depmap_cfg.get("mutation", {})
        context_mode = depmap_cfg.get("context_assignment", "defining_present")
        comparator_mode = depmap_cfg.get("comparator_assignment", "defining_absent")
        groups = {"context": [], "comparator": [], "excluded": []}
        for _, row in models.iterrows():
            model_id = str(row["ModelID"])
            disease_match = _metadata_matches(row, depmap_cfg)
            if not disease_match:
                continue
            sequenced = model_id in sequenced_models
            if not sequenced:
                status = "unsequenced"
                assigned, reason = "excluded", "No default WES/WGS profile in OmicsProfiles.csv"
                qual, other = [], []
            else:
                status, qual, other = _classify_alteration(
                    model_id,
                    mutations_by_model,
                    mutation_cfg,
                    driver_columns=driver_cols,
                    text_columns=text_cols,
                )
                if status == context_mode:
                    assigned, reason = "context", f"Alteration status = {status}"
                elif status == comparator_mode:
                    assigned, reason = "comparator", f"Alteration status = {status}"
                else:
                    assigned, reason = "excluded", f"Alteration status = {status}; not safe for context/comparator"
            groups[assigned].append(model_id)
            audit.append(
                DepMapContextAuditRow(
                    cancer_id=cancer_id,
                    model_id=model_id,
                    cell_line_name=str(row.get("CellLineName", "")) or None,
                    depmap_model_type=str(row.get("DepmapModelType", "")) or None,
                    oncotree_lineage=str(row.get("OncotreeLineage", "")) or None,
                    oncotree_primary_disease=str(row.get("OncotreePrimaryDisease", "")) or None,
                    oncotree_subtype=str(row.get("OncotreeSubtype", "")) or None,
                    oncotree_code=str(row.get("OncotreeCode", "")) or None,
                    disease_metadata_match=disease_match,
                    sequencing_available=sequenced,
                    alteration_status=status,
                    assigned_group=assigned,
                    assignment_reason=reason,
                    qualifying_variants_json=json.dumps(qual, ensure_ascii=False, sort_keys=True),
                    other_relevant_variants_json=json.dumps(other, ensure_ascii=False, sort_keys=True),
                )
            )
        memberships[cancer_id] = groups
    return memberships, audit


def cliffs_delta(x: Iterable[float], y: Iterable[float]) -> float | None:
    a = np.asarray(list(x), dtype=float)
    b = np.asarray(list(y), dtype=float)
    a = a[~np.isnan(a)]
    b = b[~np.isnan(b)]
    if len(a) == 0 or len(b) == 0:
        return None
    greater = 0
    less = 0
    for value in a:
        greater += int(np.sum(value > b))
        less += int(np.sum(value < b))
    return float((greater - less) / (len(a) * len(b)))


def _median(series: pd.Series) -> float | None:
    s = pd.to_numeric(series, errors="coerce").dropna()
    return float(s.median()) if len(s) else None


def _mean(series: pd.Series) -> float | None:
    s = pd.to_numeric(series, errors="coerce").dropna()
    return float(s.mean()) if len(s) else None


def _iqr(series: pd.Series) -> float | None:
    s = pd.to_numeric(series, errors="coerce").dropna()
    if len(s) == 0:
        return None
    return float(s.quantile(0.75) - s.quantile(0.25))


def _min(series: pd.Series) -> float | None:
    s = pd.to_numeric(series, errors="coerce").dropna()
    return float(s.min()) if len(s) else None


def _max(series: pd.Series) -> float | None:
    s = pd.to_numeric(series, errors="coerce").dropna()
    return float(s.max()) if len(s) else None


def _subset(matrix: pd.DataFrame, model_ids: list[str], gene: str) -> pd.Series:
    ids = [x for x in model_ids if x in matrix.index]
    return pd.to_numeric(matrix.loc[ids, gene], errors="coerce") if ids else pd.Series(dtype=float)



def validate_depmap_inputs(
    root: Path,
    release: str,
    depmap_dir: Path | None = None,
) -> tuple[list[dict], list[QCRecord], dict]:
    files = resolve_depmap_files(root, release, depmap_dir)
    inventory = depmap_file_inventory(files, release)
    qc: list[QCRecord] = []
    pairs = pd.read_csv(root / "data/input/cancer_target_pairs_wave1.tsv", sep="\t", dtype=str).fillna("")
    symbols = list(dict.fromkeys(pairs["HGNC_symbol"].tolist()))

    try:
        models = load_model_metadata(files["models"])
        if models["ModelID"].duplicated().any():
            qc.append(QCRecord(
                severity="ERROR", check="model_id_unique",
                observed=str(int(models["ModelID"].duplicated().sum())), expected="0",
                message="Model metadata contains duplicate ModelID values."
            ))
        else:
            qc.append(QCRecord(
                severity="INFO", check="model_id_unique", observed=str(len(models)), expected="unique",
                message="Model metadata has unique ModelID values."
            ))
        metadata_cols = {"OncotreeCode", "DepmapModelType", "OncotreeSubtype", "OncotreePrimaryDisease"}
        if not metadata_cols.intersection(models.columns):
            qc.append(QCRecord(
                severity="ERROR", check="model_disease_metadata",
                observed=str(list(models.columns)), expected="OncoTree/DepMap disease metadata",
                message="Model metadata lacks columns required for explicit disease-context filtering."
            ))
    except Exception as exc:
        qc.append(QCRecord(severity="ERROR", check="model_file_schema", message=str(exc)))
        models = pd.DataFrame()

    try:
        seq_ids, seq_schema = load_sequenced_model_ids(files["omics_profiles"])
        if not seq_ids:
            qc.append(QCRecord(
                severity="ERROR", check="sequencing_profiles", observed="0", expected=">0",
                message="No WES/WGS model profiles found in OmicsProfiles.csv."
            ))
        else:
            qc.append(QCRecord(
                severity="INFO", check="sequencing_profiles", observed=str(len(seq_ids)), expected=">0",
                message="WES/WGS model coverage is available for comparator eligibility."
            ))
        if seq_schema.get("default_entry_column") is None:
            qc.append(QCRecord(
                severity="WARNING", check="omics_default_profile_flag",
                observed="not found", expected="default-entry flag if supplied by release",
                message="OmicsProfiles has no recognized default-entry flag; all WES/WGS profiles will count as sequencing coverage."
            ))
    except Exception as exc:
        seq_schema = {}
        qc.append(QCRecord(severity="ERROR", check="omics_profiles_schema", message=str(exc)))

    for logical in ("gene_effect", "gene_dependency"):
        try:
            matrix = load_gene_matrix(files[logical], symbols)
            qc.append(QCRecord(
                severity="INFO", check=f"{logical}_targets",
                observed=f"{len(matrix)} models; {len(symbols)} Wave 1 targets",
                expected=f"all {len(symbols)} targets",
                message=f"{files[logical].name} contains all Wave 1 target symbols with unique ModelIDs."
            ))
        except Exception as exc:
            qc.append(QCRecord(severity="ERROR", check=f"{logical}_schema", message=str(exc)))

    contexts = load_context_config(root / "config/cancer_contexts.yaml")
    thresholds_cfg = yaml.safe_load((root / "config/thresholds.yaml").read_text(encoding="utf-8")) or {}
    depmap_thresholds = thresholds_cfg.get("depmap", thresholds_cfg.get("future_depmap", {}))
    global_min_context_n = int(depmap_thresholds.get("minimum_context_n_warning", 5))
    global_min_comparator_n = int(depmap_thresholds.get("minimum_comparator_n_warning", 5))
    mutation_genes: set[str] = set()
    needs_protein_change = False
    for cfg in contexts.values():
        mcfg = cfg.get("depmap", {}).get("mutation", {})
        if "gene" in mcfg:
            mutation_genes.add(str(mcfg["gene"]))
        for rule in mcfg.get("defining_variants", []) or []:
            mutation_genes.add(str(rule["gene"]))
        if mcfg.get("protein_change_regexes") or mcfg.get("defining_variants"):
            needs_protein_change = True
    try:
        muts, mut_schema = load_relevant_mutations(files["mutations"], sorted(mutation_genes))
        if needs_protein_change and mut_schema.get("protein_change") is None:
            qc.append(QCRecord(
                severity="ERROR", check="mutation_protein_change",
                observed="not found", expected="ProteinChange/HGVSp column",
                message="Exact KRAS/IDH molecular-context rules require a protein-change annotation."
            ))
        else:
            qc.append(QCRecord(
                severity="INFO", check="mutation_schema",
                observed=f"{len(muts)} relevant mutation rows",
                expected="model/gene/protein annotations",
                message="Mutation schema supports the configured Wave 1 context rules."
            ))
        if not mut_schema.get("driver_flags"):
            qc.append(QCRecord(
                severity="WARNING", check="mutation_driver_flags",
                observed="none recognized", expected="driver/hotspot annotation if available",
                message="FLT3 context will rely on explicit ITD/TKD pattern rules because no recognized driver/hotspot flag was found."
            ))
    except Exception as exc:
        mut_schema = {}
        muts = pd.DataFrame()
        qc.append(QCRecord(severity="ERROR", check="mutation_file_schema", message=str(exc)))

    # Before any target statistics, verify that the pre-specified context rules
    # actually produce auditable groups in this exact DepMap release.
    if not models.empty and 'seq_ids' in locals() and mut_schema:
        try:
            memberships, audit_rows = build_context_membership(
                models=models,
                mutations=muts,
                sequenced_models=seq_ids,
                contexts=contexts,
                mutation_schema=mut_schema,
            )
            for cancer_id, groups in memberships.items():
                context_n = len(groups.get("context", []))
                comparator_n = len(groups.get("comparator", []))
                context_audit = [x for x in audit_rows if x.cancer_id == cancer_id]
                disease_models_n = len(context_audit)
                sequenced_n = sum(1 for x in context_audit if x.sequencing_available)
                status_counts: dict[str, int] = {}
                for x in context_audit:
                    status_counts[x.alteration_status] = status_counts.get(x.alteration_status, 0) + 1
                qc.append(QCRecord(
                    severity="INFO", check="context_rule_diagnostics", entity_id=cancer_id,
                    observed=(
                        f"disease_metadata_matches={disease_models_n}; sequenced={sequenced_n}; "
                        f"alteration_statuses={json.dumps(status_counts, sort_keys=True)}; "
                        f"groups=context:{context_n},comparator:{comparator_n},excluded:{len(groups.get('excluded', []))}"
                    ),
                    expected="auditable nonzero disease and context counts",
                    message="Diagnostic breakdown of disease matching, sequencing coverage, mutation classification, and final cohort assignment."
                ))
                context_min = int(contexts[cancer_id].get("minimum_context_n_warning", global_min_context_n))
                comparator_min = global_min_comparator_n
                if context_n == 0:
                    qc.append(QCRecord(
                        severity="ERROR", check="context_rule_feasibility", entity_id=cancer_id,
                        observed=f"context=0; comparator={comparator_n}", expected="context > 0",
                        message="Configured molecular-context rule produced no models in this release; do not proceed by guessing a broader cohort."
                    ))
                elif comparator_n == 0:
                    qc.append(QCRecord(
                        severity="WARNING", check="context_rule_feasibility", entity_id=cancer_id,
                        observed=f"context={context_n}; comparator=0", expected="comparator > 0 when context-specific testing is intended",
                        message="Context exists but no same-disease comparator was formed; descriptive dependency can still be calculated, but context-specific inference is unavailable."
                    ))
                else:
                    qc.append(QCRecord(
                        severity="INFO", check="context_rule_feasibility", entity_id=cancer_id,
                        observed=f"context={context_n}; comparator={comparator_n}", expected="both groups > 0",
                        message="Configured disease/mutation rules form both context and comparator groups."
                    ))

                if 0 < context_n < context_min:
                    qc.append(QCRecord(
                        severity="WARNING", check="context_sample_size", entity_id=cancer_id,
                        observed=str(context_n), expected=f">={context_min}",
                        message="Molecular-context cohort is below the pre-specified sample-size warning threshold; quantitative results may be calculated but must remain exploratory."
                    ))
                if 0 < comparator_n < comparator_min:
                    qc.append(QCRecord(
                        severity="WARNING", check="comparator_sample_size", entity_id=cancer_id,
                        observed=str(comparator_n), expected=f">={comparator_min}",
                        message="Comparator cohort is below the pre-specified sample-size warning threshold; context-selective inference must remain exploratory."
                    ))
        except Exception as exc:
            qc.append(QCRecord(severity="ERROR", check="context_rule_feasibility", message=str(exc)))

    meta = {
        "release": release,
        "files": {k: str(v) for k, v in files.items()},
        "omics_profiles_schema": seq_schema,
        "mutation_schema": {k: v for k, v in mut_schema.items() if k != "all_columns"},
    }
    return inventory, qc, meta


def analyze_depmap_wave1(
    root: Path,
    release: str,
    depmap_dir: Path | None = None,
) -> tuple[
    list[DepMapEvidence],
    list[DepMapContextAuditRow],
    list[ProvenanceRecord],
    list[dict],
    dict,
]:
    files = resolve_depmap_files(root, release, depmap_dir)
    inventory = depmap_file_inventory(files, release)
    contexts = load_context_config(root / "config/cancer_contexts.yaml")
    pairs = pd.read_csv(root / "data/input/cancer_target_pairs_wave1.tsv", sep="\t", dtype=str).fillna("")
    symbols = list(dict.fromkeys(pairs["HGNC_symbol"].tolist()))

    models = load_model_metadata(files["models"])
    sequenced_models, profile_schema = load_sequenced_model_ids(files["omics_profiles"])

    mutation_genes: set[str] = set()
    for cfg in contexts.values():
        mcfg = cfg.get("depmap", {}).get("mutation", {})
        if "gene" in mcfg:
            mutation_genes.add(str(mcfg["gene"]))
        for rule in mcfg.get("defining_variants", []) or []:
            mutation_genes.add(str(rule["gene"]))
    mutations, mut_schema = load_relevant_mutations(files["mutations"], sorted(mutation_genes))

    memberships, audit = build_context_membership(
        models=models,
        mutations=mutations,
        sequenced_models=sequenced_models,
        contexts=contexts,
        mutation_schema=mut_schema,
    )
    gene_effect = load_gene_matrix(files["gene_effect"], symbols)
    dependency = load_gene_matrix(files["gene_dependency"], symbols)

    thresholds = yaml.safe_load((root / "config/thresholds.yaml").read_text(encoding="utf-8")) or {}
    dcfg = thresholds.get("depmap", thresholds.get("future_depmap", {}))
    min_context_n = int(dcfg.get("minimum_context_n_warning", 3))
    min_comparator_n = int(dcfg.get("minimum_comparator_n_warning", 3))
    dep_threshold = float(dcfg.get("dependency_probability_threshold", 0.5))
    broad_threshold = float(dcfg.get("broad_dependency_fraction_warning", 0.8))
    min_test_n = int(dcfg.get("minimum_n_for_statistical_test", 2))

    retrieved_at = datetime.now(timezone.utc).isoformat()
    evidences: list[DepMapEvidence] = []
    provenance: list[ProvenanceRecord] = []
    p_positions: list[int] = []
    p_values: list[float] = []

    for _, pair in pairs.iterrows():
        cancer_id = pair["Cancer_ID"]
        gene = pair["HGNC_symbol"]
        groups = memberships.get(cancer_id, {"context": [], "comparator": []})
        context_ids = groups.get("context", [])
        comparator_ids = groups.get("comparator", [])
        c_ge = _subset(gene_effect, context_ids, gene)
        r_ge = _subset(gene_effect, comparator_ids, gene)
        c_dp = _subset(dependency, context_ids, gene)
        r_dp = _subset(dependency, comparator_ids, gene)
        all_ge = pd.to_numeric(gene_effect[gene], errors="coerce").dropna()
        all_dp = pd.to_numeric(dependency[gene], errors="coerce").dropna()

        c_ge_valid = c_ge.dropna()
        r_ge_valid = r_ge.dropna()
        c_dp_valid = c_dp.dropna()
        r_dp_valid = r_dp.dropna()
        low_context_sample = len(c_ge_valid) < min_context_n
        low_comparator_sample = bool(comparator_ids) and len(r_ge_valid) < min_comparator_n
        low_sample = low_context_sample or low_comparator_sample

        p_value = None
        delta = None
        stat_test = None
        performed = False
        if len(c_ge_valid) >= min_test_n and len(r_ge_valid) >= min_test_n:
            result = mannwhitneyu(c_ge_valid, r_ge_valid, alternative="two-sided", method="auto")
            p_value = float(result.pvalue)
            delta = cliffs_delta(c_ge_valid, r_ge_valid)
            stat_test = "Mann–Whitney U, two-sided"
            performed = True

        c_med = _median(c_ge)
        r_med = _median(r_ge)
        delta_ge = None if c_med is None or r_med is None else float(c_med - r_med)
        dep_frac = float((c_dp_valid > dep_threshold).mean()) if len(c_dp_valid) else None
        r_dep_frac = float((r_dp_valid > dep_threshold).mean()) if len(r_dp_valid) else None
        broad_frac = float((all_dp > dep_threshold).mean()) if len(all_dp) else None
        context_cfg = contexts[cancer_id].get("depmap", {})
        context_def = str(context_cfg.get("context_definition", contexts[cancer_id].get("name", cancer_id)))
        comparator_def = str(context_cfg.get("comparator_definition", contexts[cancer_id].get("comparator", "")))

        ev = DepMapEvidence(
            evidence_id=pair["Evidence_ID"],
            cancer_id=cancer_id,
            target_id=pair["Target_ID"],
            target_symbol=gene,
            depmap_release=release,
            context_definition=context_def,
            comparator_definition=comparator_def,
            context_models_n=len(context_ids),
            comparator_models_n=len(comparator_ids),
            context_gene_effect_n=len(c_ge_valid),
            comparator_gene_effect_n=len(r_ge_valid),
            context_dependency_probability_n=len(c_dp_valid),
            comparator_dependency_probability_n=len(r_dp_valid),
            dependent_cell_lines_n=int((c_dp_valid > dep_threshold).sum()),
            dependency_fraction=dep_frac,
            median_gene_effect=c_med,
            mean_gene_effect=_mean(c_ge),
            gene_effect_iqr=_iqr(c_ge),
            median_dependency_probability=_median(c_dp),
            min_gene_effect=_min(c_ge),
            max_gene_effect=_max(c_ge),
            missing_context_gene_effect_n=max(0, len(context_ids) - len(c_ge_valid)),
            comparator_dependent_cell_lines_n=int((r_dp_valid > dep_threshold).sum()),
            comparator_dependency_fraction=r_dep_frac,
            comparator_median_gene_effect=r_med,
            comparator_median_dependency_probability=_median(r_dp),
            delta_gene_effect=delta_ge,
            cliffs_delta=delta,
            statistical_test=stat_test,
            p_value=p_value,
            q_value=None,
            low_sample_size=low_sample,
            statistical_test_performed=performed,
            all_models_gene_effect_n=len(all_ge),
            all_models_dependency_probability_n=len(all_dp),
            all_models_median_gene_effect=_median(all_ge),
            broad_dependency_fraction=broad_frac,
            broad_dependency_warning=bool(broad_frac is not None and broad_frac >= broad_threshold),
            gene_effect_file=_rel_or_abs(files["gene_effect"], root),
            dependency_probability_file=_rel_or_abs(files["gene_dependency"], root),
            model_file=_rel_or_abs(files["models"], root),
            mutation_file=_rel_or_abs(files["mutations"], root),
            omics_profiles_file=_rel_or_abs(files["omics_profiles"], root),
            retrieved_at=retrieved_at,
            notes=(
                "Chronos Gene Effect is treated as a continuous functional-dependency measure. "
                "Dependent/non-dependent calls use CRISPRGeneDependency probability > "
                f"{dep_threshold}. A negative delta_gene_effect means stronger dependency in the "
                "molecular context than in the same-disease comparator. Cliff's delta is context "
                "vs comparator; negative values likewise indicate lower (stronger-dependency) "
                "Gene Effect in context. broad_dependency_warning is a project screening flag, "
                "not an official DepMap common-essential classification and not a safety verdict."
            ),
        )
        evidences.append(ev)
        if p_value is not None:
            p_positions.append(len(evidences) - 1)
            p_values.append(p_value)

    if p_values:
        qvals = multipletests(p_values, alpha=0.05, method="fdr_bh")[1]
        for pos, q in zip(p_positions, qvals):
            evidences[pos].q_value = float(q)

    file_hashes = {k: sha256_file(v) for k, v in files.items()}
    for ev in evidences:
        entity_id = f"{ev.cancer_id}|{ev.target_id}"
        provenance_specs = [
            ("median_gene_effect", ev.median_gene_effect, "DepMap CRISPR Gene Effect", files["gene_effect"], "Observed"),
            ("dependency_fraction", ev.dependency_fraction, "DepMap CRISPR Gene Dependency", files["gene_dependency"], "Observed"),
            ("median_dependency_probability", ev.median_dependency_probability, "DepMap CRISPR Gene Dependency", files["gene_dependency"], "Observed"),
            ("context_models_n", ev.context_models_n, "DepMap Model + Omics", files["models"], "Derived"),
            ("delta_gene_effect", ev.delta_gene_effect, "MCL DepMap analysis", files["gene_effect"], "Derived"),
            ("cliffs_delta", ev.cliffs_delta, "MCL DepMap analysis", files["gene_effect"], "Derived"),
            ("p_value", ev.p_value, "MCL DepMap analysis", files["gene_effect"], "Derived"),
            ("q_value", ev.q_value, "MCL DepMap analysis", files["gene_effect"], "Derived"),
            ("broad_dependency_fraction", ev.broad_dependency_fraction, "DepMap CRISPR Gene Dependency", files["gene_dependency"], "Derived"),
        ]
        for field_name, value, source_name, raw_path, nature in provenance_specs:
            provenance.append(
                ProvenanceRecord(
                    entity_type="CancerTargetEvidence",
                    entity_id=entity_id,
                    field_name=field_name,
                    value=None if value is None else str(value),
                    source_name=source_name,
                    source_record_id=ev.evidence_id,
                    source_release=release,
                    source_url_or_endpoint=OFFICIAL_SOURCE,
                    retrieved_at=retrieved_at,
                    evidence_nature=nature,
                    parser_version=PARSER_VERSION,
                    raw_file=_rel_or_abs(raw_path, root),
                    raw_record_hash=file_hashes[
                        "gene_effect" if raw_path == files["gene_effect"] else
                        "gene_dependency" if raw_path == files["gene_dependency"] else
                        "models"
                    ],
                )
            )

    meta = {
        "release": release,
        "parser_version": PARSER_VERSION,
        "source_url": OFFICIAL_SOURCE,
        "retrieved_at": retrieved_at,
        "files": {
            key: {
                "path": _rel_or_abs(path, root),
                "sha256": file_hashes[key],
                "size_bytes": path.stat().st_size,
            }
            for key, path in files.items()
        },
        "omics_profiles_schema": profile_schema,
        "mutation_schema": {
            k: v for k, v in mut_schema.items() if k != "all_columns"
        },
        "multiplicity_correction": "Benjamini–Hochberg across all Wave 1 comparisons with calculable p-values",
        "dependency_probability_threshold": dep_threshold,
        "minimum_context_n_warning": min_context_n,
        "minimum_comparator_n_warning": min_comparator_n,
        "minimum_n_for_statistical_test": min_test_n,
    }
    return evidences, audit, provenance, inventory, meta


def _target_variant_label(patterns: list[str]) -> str:
    """Return a human-readable KRAS protein-change label from an exact configured regex."""
    for pattern in patterns:
        match = re.fullmatch(r"\^p\\\.([A-Za-z0-9_*]+)\$", str(pattern))
        if match:
            return f"p.{match.group(1)}"
    return "configured KRAS variant"


def build_kras_sensitivity_membership(
    models: pd.DataFrame,
    mutations: pd.DataFrame,
    sequenced_models: set[str],
    contexts: dict[str, dict],
    mutation_schema: dict[str, object],
) -> tuple[dict[str, dict[str, list[str]]], list[DepMapKrasSensitivityAuditRow]]:
    """Split same-disease KRAS contexts into target variant, WT proxy, and other driver/hotspot KRAS.

    WT proxy is deliberately conservative: the model must have default WES/WGS coverage and
    no KRAS variant row in OmicsSomaticMutations. A non-target KRAS variant without a recognized
    driver/hotspot annotation is excluded as ambiguous rather than called WT.
    """
    if "ModelID" not in models.columns:
        raise ValueError("Model metadata must contain ModelID")
    models = models.copy()
    models["ModelID"] = models["ModelID"].astype(str)
    mutations_by_model = {
        str(mid): group.copy()
        for mid, group in mutations.groupby("ModelID", sort=False)
    } if not mutations.empty else {}
    driver_cols = [str(x) for x in mutation_schema.get("driver_flags", [])]

    memberships: dict[str, dict[str, list[str]]] = {}
    audit: list[DepMapKrasSensitivityAuditRow] = []

    for cancer_id, cfg in contexts.items():
        depmap_cfg = cfg.get("depmap", {})
        mutation_cfg = depmap_cfg.get("mutation", {})
        if (
            str(mutation_cfg.get("mode", "")) != "exact_protein_change_present"
            or str(mutation_cfg.get("gene", "")).upper() != "KRAS"
        ):
            continue

        patterns = list(mutation_cfg.get("protein_change_regexes", []) or [])
        target_label = _target_variant_label(patterns)
        groups = {
            "target_variant": [],
            "kras_wildtype_proxy": [],
            "other_kras_driver_hotspot": [],
            "excluded": [],
        }

        for _, row in models.iterrows():
            if not _metadata_matches(row, depmap_cfg):
                continue
            model_id = str(row["ModelID"])
            sequenced = model_id in sequenced_models
            target_rows: list[dict] = []
            other_rows: list[dict] = []

            if not sequenced:
                status = "unsequenced"
                assigned = "excluded"
                reason = "No default WES/WGS profile in OmicsProfiles.csv"
            else:
                rows = mutations_by_model.get(model_id, pd.DataFrame())
                if rows.empty:
                    kras_rows = rows
                else:
                    kras_rows = rows[
                        rows.get("GeneSymbol", pd.Series(dtype=str)).astype(str).str.upper() == "KRAS"
                    ]

                target_mask: list[bool] = []
                if not kras_rows.empty:
                    for _, mut in kras_rows.iterrows():
                        is_target = _regex_match(mut.get("ProteinChange", ""), patterns)
                        target_mask.append(is_target)
                        payload = _variant_payload(mut)
                        if is_target:
                            target_rows.append(payload)
                        else:
                            other_rows.append(payload)

                if target_rows:
                    status = "target_variant"
                    assigned = "target_variant"
                    reason = f"Exact configured KRAS target variant {target_label} detected"
                elif kras_rows.empty:
                    status = "kras_wildtype_proxy"
                    assigned = "kras_wildtype_proxy"
                    reason = "Sequenced model with no KRAS variant record; used as KRAS WT proxy"
                else:
                    other_driver = any(
                        row_has_driver_signal(mut, driver_cols)
                        for _, mut in kras_rows.iterrows()
                    )
                    if other_driver:
                        status = "other_kras_driver_hotspot"
                        assigned = "other_kras_driver_hotspot"
                        reason = "Non-target KRAS variant with DepMap driver/hotspot annotation"
                    else:
                        status = "ambiguous_kras_variant"
                        assigned = "excluded"
                        reason = "Non-target KRAS variant lacks recognized driver/hotspot annotation"

            groups[assigned].append(model_id)
            audit.append(
                DepMapKrasSensitivityAuditRow(
                    cancer_id=cancer_id,
                    target_kras_variant=target_label,
                    model_id=model_id,
                    cell_line_name=str(row.get("CellLineName", "")) or None,
                    depmap_model_type=str(row.get("DepmapModelType", "")) or None,
                    oncotree_lineage=str(row.get("OncotreeLineage", "")) or None,
                    oncotree_primary_disease=str(row.get("OncotreePrimaryDisease", "")) or None,
                    oncotree_subtype=str(row.get("OncotreeSubtype", "")) or None,
                    oncotree_code=str(row.get("OncotreeCode", "")) or None,
                    sequencing_available=sequenced,
                    kras_status=status,
                    assigned_group=assigned,
                    assignment_reason=reason,
                    target_variants_json=json.dumps(target_rows, ensure_ascii=False, sort_keys=True),
                    other_kras_variants_json=json.dumps(other_rows, ensure_ascii=False, sort_keys=True),
                )
            )
        memberships[cancer_id] = groups
    return memberships, audit


def analyze_depmap_kras_sensitivity(
    root: Path,
    release: str,
    depmap_dir: Path | None = None,
) -> tuple[list[DepMapKrasSensitivityEvidence], list[DepMapKrasSensitivityAuditRow], dict]:
    """Run M3.1 sensitivity analysis for KRAS G12C LUAD and KRAS G12D pancreatic cancer."""
    files = resolve_depmap_files(root, release, depmap_dir)
    contexts = load_context_config(root / "config/cancer_contexts.yaml")
    pairs = pd.read_csv(root / "data/input/cancer_target_pairs_wave1.tsv", sep="\t", dtype=str).fillna("")

    eligible_contexts = {
        cancer_id: cfg
        for cancer_id, cfg in contexts.items()
        if str(cfg.get("depmap", {}).get("mutation", {}).get("mode", "")) == "exact_protein_change_present"
        and str(cfg.get("depmap", {}).get("mutation", {}).get("gene", "")).upper() == "KRAS"
    }
    if not eligible_contexts:
        raise ValueError("No exact KRAS molecular contexts are configured for sensitivity analysis")

    kras_pairs = pairs[pairs["Cancer_ID"].isin(eligible_contexts)].copy()
    symbols = list(dict.fromkeys(kras_pairs["HGNC_symbol"].tolist()))
    models = load_model_metadata(files["models"])
    sequenced_models, profile_schema = load_sequenced_model_ids(files["omics_profiles"])
    mutations, mut_schema = load_relevant_mutations(files["mutations"], ["KRAS"])
    memberships, audit = build_kras_sensitivity_membership(
        models=models,
        mutations=mutations,
        sequenced_models=sequenced_models,
        contexts=eligible_contexts,
        mutation_schema=mut_schema,
    )
    gene_effect = load_gene_matrix(files["gene_effect"], symbols)
    dependency = load_gene_matrix(files["gene_dependency"], symbols)

    thresholds = yaml.safe_load((root / "config/thresholds.yaml").read_text(encoding="utf-8")) or {}
    dcfg = thresholds.get("depmap", thresholds.get("future_depmap", {}))
    min_context_n = int(dcfg.get("minimum_context_n_warning", 5))
    min_comparator_n = int(dcfg.get("minimum_comparator_n_warning", 5))
    min_test_n = int(dcfg.get("minimum_n_for_statistical_test", 2))
    dep_threshold = float(dcfg.get("dependency_probability_threshold", 0.5))
    retrieved_at = datetime.now(timezone.utc).isoformat()

    rows: list[DepMapKrasSensitivityEvidence] = []
    p_positions: list[int] = []
    p_values: list[float] = []

    comparison_specs = [
        (
            "vs_kras_wildtype_proxy",
            "kras_wildtype_proxy",
            "same-disease models with WES/WGS coverage and no detected KRAS variant",
        ),
        (
            "vs_other_kras_driver_hotspot",
            "other_kras_driver_hotspot",
            "same-disease models carrying a non-target KRAS variant annotated as driver/hotspot",
        ),
    ]

    for _, pair in kras_pairs.iterrows():
        cancer_id = str(pair["Cancer_ID"])
        gene = str(pair["HGNC_symbol"])
        groups = memberships[cancer_id]
        context_ids = groups["target_variant"]
        target_patterns = eligible_contexts[cancer_id]["depmap"]["mutation"].get("protein_change_regexes", []) or []
        target_label = _target_variant_label(list(target_patterns))

        for comparison_type, comparator_key, comparator_definition in comparison_specs:
            comparator_ids = groups[comparator_key]
            c_ge = _subset(gene_effect, context_ids, gene)
            r_ge = _subset(gene_effect, comparator_ids, gene)
            c_dp = _subset(dependency, context_ids, gene)
            r_dp = _subset(dependency, comparator_ids, gene)
            c_ge_valid, r_ge_valid = c_ge.dropna(), r_ge.dropna()
            c_dp_valid, r_dp_valid = c_dp.dropna(), r_dp.dropna()

            low_sample = (
                len(c_ge_valid) < min_context_n
                or (bool(comparator_ids) and len(r_ge_valid) < min_comparator_n)
            )
            p_value = None
            delta = None
            stat_test = None
            performed = False
            if len(c_ge_valid) >= min_test_n and len(r_ge_valid) >= min_test_n:
                result = mannwhitneyu(c_ge_valid, r_ge_valid, alternative="two-sided", method="auto")
                p_value = float(result.pvalue)
                delta = cliffs_delta(c_ge_valid, r_ge_valid)
                stat_test = "Mann–Whitney U, two-sided"
                performed = True

            c_med = _median(c_ge)
            r_med = _median(r_ge)
            delta_ge = None if c_med is None or r_med is None else float(c_med - r_med)
            dep_frac = float((c_dp_valid > dep_threshold).mean()) if len(c_dp_valid) else None
            r_dep_frac = float((r_dp_valid > dep_threshold).mean()) if len(r_dp_valid) else None

            row = DepMapKrasSensitivityEvidence(
                sensitivity_id=f"{pair['Evidence_ID']}|{comparison_type}",
                evidence_id=str(pair["Evidence_ID"]),
                cancer_id=cancer_id,
                target_id=str(pair["Target_ID"]),
                target_symbol=gene,
                depmap_release=release,
                target_kras_variant=f"KRAS {target_label}",
                comparison_type=comparison_type,
                context_definition=f"same-disease models with KRAS {target_label}",
                comparator_definition=comparator_definition,
                context_models_n=len(context_ids),
                comparator_models_n=len(comparator_ids),
                context_gene_effect_n=len(c_ge_valid),
                comparator_gene_effect_n=len(r_ge_valid),
                context_dependency_probability_n=len(c_dp_valid),
                comparator_dependency_probability_n=len(r_dp_valid),
                dependent_cell_lines_n=int((c_dp_valid > dep_threshold).sum()),
                dependency_fraction=dep_frac,
                median_gene_effect=c_med,
                median_dependency_probability=_median(c_dp),
                comparator_dependent_cell_lines_n=int((r_dp_valid > dep_threshold).sum()),
                comparator_dependency_fraction=r_dep_frac,
                comparator_median_gene_effect=r_med,
                comparator_median_dependency_probability=_median(r_dp),
                delta_gene_effect=delta_ge,
                cliffs_delta=delta,
                statistical_test=stat_test,
                p_value=p_value,
                q_value=None,
                low_sample_size=low_sample,
                statistical_test_performed=performed,
                retrieved_at=retrieved_at,
                notes=(
                    "M3.1 sensitivity analysis. The primary M3 comparator remains unchanged. "
                    "KRAS WT is a sequencing-backed proxy defined as no KRAS variant record in "
                    "OmicsSomaticMutations among models with a default WES/WGS profile. Other KRAS "
                    "mutant comparator requires a non-target KRAS variant with a recognized DepMap "
                    "driver/hotspot annotation. Ambiguous KRAS variants are excluded. Negative "
                    "delta_gene_effect and Cliff's delta indicate stronger dependency in the exact "
                    "KRAS-variant context."
                ),
            )
            rows.append(row)
            if p_value is not None:
                p_positions.append(len(rows) - 1)
                p_values.append(p_value)

    if p_values:
        qvals = multipletests(p_values, alpha=0.05, method="fdr_bh")[1]
        for pos, q in zip(p_positions, qvals):
            rows[pos].q_value = float(q)

    file_hashes = {k: sha256_file(v) for k, v in files.items()}
    meta = {
        "release": release,
        "parser_version": KRAS_SENSITIVITY_PARSER_VERSION,
        "analysis": "M3.1 KRAS comparator sensitivity analysis",
        "source_url": OFFICIAL_SOURCE,
        "retrieved_at": retrieved_at,
        "eligible_contexts": sorted(eligible_contexts),
        "comparison_types": [x[0] for x in comparison_specs],
        "classification": {
            "target_variant": "exact configured KRAS protein change",
            "kras_wildtype_proxy": "WES/WGS-covered same-disease model with no KRAS variant record",
            "other_kras_driver_hotspot": "non-target KRAS variant with DepMap driver/hotspot annotation",
            "ambiguous_kras_variant": "KRAS variant present but not target and no recognized driver/hotspot annotation; excluded",
        },
        "multiplicity_correction": "Benjamini–Hochberg across all calculable M3.1 KRAS sensitivity comparisons",
        "dependency_probability_threshold": dep_threshold,
        "minimum_context_n_warning": min_context_n,
        "minimum_comparator_n_warning": min_comparator_n,
        "minimum_n_for_statistical_test": min_test_n,
        "omics_profiles_schema": profile_schema,
        "mutation_schema": {k: v for k, v in mut_schema.items() if k != "all_columns"},
        "files": {
            key: {
                "path": _rel_or_abs(path, root),
                "sha256": file_hashes[key],
                "size_bytes": path.stat().st_size,
            }
            for key, path in files.items()
        },
    }
    return rows, audit, meta
