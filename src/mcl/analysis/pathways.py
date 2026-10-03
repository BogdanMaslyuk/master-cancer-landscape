from __future__ import annotations

import json
from collections import Counter
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from mcl.sources.gprofiler import (
    DEFAULT_CONVERT_ENDPOINT,
    DEFAULT_PROFILE_ENDPOINT,
    DEFAULT_VERSIONS_ENDPOINT,
    fetch_gprofiler_versions,
    run_gprofiler_convert,
    run_gprofiler_profile,
)


PATHWAY_ANALYSIS_VERSION = "m3.3-pathways-v0.3"



def _normalize_entrez_id(value: Any) -> str | None:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    text = str(value).strip()
    if not text or text.lower() in {"nan", "none"}:
        return None
    if text.endswith(".0") and text[:-2].isdigit():
        text = text[:-2]
    return text if text.isdigit() else None


def build_entrez_identifier_map(frame: pd.DataFrame) -> dict[str, str]:
    if "gene_symbol" not in frame.columns or "entrez_gene_id" not in frame.columns:
        raise ValueError("Genome-wide table lacks gene_symbol or entrez_gene_id")
    mapping: dict[str, str] = {}
    for symbol, entrez in frame[["gene_symbol", "entrez_gene_id"]].itertuples(index=False, name=None):
        symbol_text = str(symbol).strip() if pd.notna(symbol) else ""
        entrez_text = _normalize_entrez_id(entrez)
        if not symbol_text or entrez_text is None:
            continue
        previous = mapping.get(symbol_text)
        if previous is not None and previous != entrez_text:
            raise ValueError(
                f"Conflicting Entrez identifiers for {symbol_text}: {previous} vs {entrez_text}"
            )
        mapping[symbol_text] = entrez_text
    return mapping




def resolve_ensembl_identifiers(
    symbol_to_entrez: dict[str, str],
    convert_rows: list[dict[str, Any]],
) -> pd.DataFrame:
    """Resolve Entrez IDs to one canonical Ensembl gene per MCL symbol.

    g:Convert can legitimately return multiple Ensembl genes for one numeric ID
    (for example a canonical gene plus a readthrough transcript locus). MCL uses
    the HGNC-approved symbol already attached to the DepMap gene to disambiguate
    such cases whenever exactly one conversion row has a matching gene name.
    """
    by_incoming: dict[str, list[dict[str, Any]]] = {}
    for row in convert_rows:
        incoming = str(row.get("incoming") or "").strip()
        if incoming:
            by_incoming.setdefault(incoming, []).append(row)

    records: list[dict[str, Any]] = []
    for symbol, entrez in sorted(symbol_to_entrez.items()):
        rows = by_incoming.get(str(entrez), [])
        candidates: list[tuple[str, str]] = []
        for row in rows:
            converted = str(row.get("converted") or "").strip()
            name = str(row.get("name") or "").strip()
            if not converted or converted.lower() in {"nan", "none", "n/a"}:
                continue
            if not converted.startswith("ENSG"):
                continue
            pair = (converted, name)
            if pair not in candidates:
                candidates.append(pair)

        exact = sorted({ensg for ensg, name in candidates if name.casefold() == symbol.casefold()})
        all_ids = sorted({ensg for ensg, _ in candidates})
        if len(exact) == 1:
            selected = exact[0]
            status = "resolved_symbol_match"
        elif len(all_ids) == 1:
            selected = all_ids[0]
            status = "resolved_single"
        elif not all_ids:
            selected = None
            status = "unmapped"
        else:
            selected = None
            status = "ambiguous"

        records.append({
            "gene_symbol": symbol,
            "entrez_gene_id": entrez,
            "ensembl_gene_id": selected,
            "mapping_status": status,
            "candidate_ensembl_ids_json": json.dumps(all_ids, ensure_ascii=False),
            "candidate_names_json": json.dumps(
                [{"ensembl_gene_id": ensg, "name": name} for ensg, name in candidates],
                ensure_ascii=False,
            ),
        })
    return pd.DataFrame(records)


def _convert_in_chunks(
    identifiers: list[str],
    *,
    organism: str,
    endpoint: str,
    numeric_ns: str,
    chunk_size: int,
) -> tuple[list[dict[str, Any]], dict[str, dict[str, Any]]]:
    rows: list[dict[str, Any]] = []
    raw: dict[str, dict[str, Any]] = {}
    for start in range(0, len(identifiers), chunk_size):
        chunk = identifiers[start:start + chunk_size]
        response = run_gprofiler_convert(
            chunk,
            organism=organism,
            target="ENSG",
            numeric_ns=numeric_ns,
            endpoint=endpoint,
        )
        rows.extend(response.get("result") or [])
        raw[f"identifier_conversion_{start // chunk_size + 1:03d}"] = response
    return rows, raw


def _intersecting_input_ids(response: dict[str, Any], row: dict[str, Any]) -> list[str]:
    intersections = row.get("intersections") or []
    if not intersections:
        return []
    query_name = str(row.get("query") or "query_1")
    qmeta = (
        ((response.get("meta") or {}).get("genes_metadata") or {})
        .get("query", {})
        .get(query_name, {})
    )
    ensgs = list(qmeta.get("ensgs") or [])
    mapping = qmeta.get("mapping") or {}
    if not ensgs:
        return []

    reverse: dict[str, list[str]] = {}
    if isinstance(mapping, dict):
        for incoming, mapped in mapping.items():
            values = mapped if isinstance(mapped, list) else [mapped]
            for ensg in values:
                if ensg:
                    reverse.setdefault(str(ensg), []).append(str(incoming))

    hits: list[str] = []
    for idx, evidence in enumerate(intersections):
        if idx >= len(ensgs) or not evidence:
            continue
        ensg = str(ensgs[idx])
        incoming = reverse.get(ensg) or [ensg]
        for value in incoming:
            if value not in hits:
                hits.append(value)
    return hits

def load_pathway_config(path: Path) -> dict[str, Any]:
    data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if not isinstance(data, dict):
        raise ValueError("Pathway configuration must be a YAML mapping")
    comparisons = data.get("comparisons") or []
    if not comparisons:
        raise ValueError("Pathway configuration contains no comparisons")
    return data


def _bool_series(frame: pd.DataFrame, column: str, default: bool = False) -> pd.Series:
    if column not in frame.columns:
        return pd.Series(default, index=frame.index, dtype=bool)
    values = frame[column]
    if pd.api.types.is_bool_dtype(values):
        return values.fillna(default).astype(bool)
    normalized = values.astype(str).str.strip().str.lower()
    mapping = {"true": True, "1": True, "yes": True, "false": False, "0": False, "no": False}
    mapped = normalized.map(mapping)
    return mapped.fillna(default).astype(bool)


def eligible_gene_set(
    frame: pd.DataFrame,
    *,
    exclude_broad_dependency: bool = True,
    exclude_low_sample_size: bool = True,
) -> set[str]:
    """Genes that were technically eligible to enter a ranked candidate list.

    Direction of effect is intentionally *not* part of the background definition.
    The universe should represent genes that could have been selected, not genes
    that happened to show the selected direction in the observed data.
    """
    if "gene_symbol" not in frame.columns or "delta_gene_effect" not in frame.columns:
        raise ValueError("Genome-wide table lacks gene_symbol or delta_gene_effect")

    mask = pd.to_numeric(frame["delta_gene_effect"], errors="coerce").notna()
    if exclude_broad_dependency:
        mask &= ~_bool_series(frame, "broad_dependency_warning")
    if exclude_low_sample_size:
        mask &= ~_bool_series(frame, "low_sample_size")

    genes = frame.loc[mask, "gene_symbol"].dropna().astype(str).str.strip()
    return {x for x in genes if x}


def select_top_candidates(
    frame: pd.DataFrame,
    *,
    comparison_label: str,
    top_n: int = 100,
    require_negative_delta: bool = True,
    exclude_broad_dependency: bool = True,
    exclude_low_sample_size: bool = True,
) -> pd.DataFrame:
    if top_n <= 0:
        raise ValueError("top_n must be positive")
    eligible = eligible_gene_set(
        frame,
        exclude_broad_dependency=exclude_broad_dependency,
        exclude_low_sample_size=exclude_low_sample_size,
    )
    work = frame[frame["gene_symbol"].astype(str).isin(eligible)].copy()
    work["delta_gene_effect"] = pd.to_numeric(work["delta_gene_effect"], errors="coerce")
    if require_negative_delta:
        work = work[work["delta_gene_effect"] < 0].copy()
    work = work.sort_values(["delta_gene_effect", "gene_symbol"], ascending=[True, True]).head(top_n)
    work.insert(0, "rank", range(1, len(work) + 1))
    work.insert(0, "comparison", comparison_label)
    return work.reset_index(drop=True)


def build_candidate_recurrence(candidate_long: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "gene_symbol",
        "comparisons_n",
        "comparisons",
        "best_rank",
        "mean_rank",
        "mean_delta_gene_effect",
    ]
    if candidate_long.empty:
        return pd.DataFrame(columns=columns)

    grouped = []
    for gene, rows in candidate_long.groupby("gene_symbol", sort=False):
        ranks = pd.to_numeric(rows["rank"], errors="coerce")
        deltas = pd.to_numeric(rows["delta_gene_effect"], errors="coerce")
        grouped.append(
            {
                "gene_symbol": gene,
                "comparisons_n": int(rows["comparison"].nunique()),
                "comparisons": "; ".join(rows["comparison"].astype(str).tolist()),
                "best_rank": int(ranks.min()),
                "mean_rank": round(float(ranks.mean()), 2),
                "mean_delta_gene_effect": round(float(deltas.mean()), 6),
            }
        )
    result = pd.DataFrame(grouped)
    return result.sort_values(
        ["comparisons_n", "mean_rank", "gene_symbol"],
        ascending=[False, True, True],
    ).reset_index(drop=True)


def build_background_sets(
    eligible_by_comparison: dict[str, set[str]],
    thresholds: dict[str, int],
) -> dict[str, set[str]]:
    counts: Counter[str] = Counter()
    for genes in eligible_by_comparison.values():
        counts.update(set(genes))
    return {
        name: {gene for gene, n in counts.items() if n >= int(minimum)}
        for name, minimum in thresholds.items()
    }


def normalize_gprofiler_result(
    response: dict[str, Any],
    analysis_set: str,
    *,
    excluded_term_ids: set[str] | None = None,
    input_id_to_symbol: dict[str, str] | None = None,
) -> pd.DataFrame:
    rows = response.get("result") or []
    excluded_term_ids = excluded_term_ids or set()
    input_id_to_symbol = input_id_to_symbol or {}
    normalized = []
    for row in rows:
        term_id = str(row.get("native") or "")
        if term_id in excluded_term_ids:
            continue
        input_ids = _intersecting_input_ids(response, row)
        symbols = [input_id_to_symbol.get(x, x) for x in input_ids]
        normalized.append(
            {
                "analysis_set": analysis_set,
                "source": row.get("source"),
                "term_id": term_id,
                "term_name": row.get("name"),
                "description": row.get("description"),
                "p_value_adjusted": row.get("p_value"),
                "significant": bool(row.get("significant", False)),
                "query_size": row.get("query_size"),
                "term_size": row.get("term_size"),
                "intersection_size": row.get("intersection_size"),
                "effective_domain_size": row.get("effective_domain_size"),
                "precision": row.get("precision"),
                "recall": row.get("recall"),
                "parents_json": json.dumps(row.get("parents") or [], ensure_ascii=False),
                "intersecting_input_ids_json": json.dumps(input_ids, ensure_ascii=False),
                "intersecting_gene_symbols_json": json.dumps(symbols, ensure_ascii=False),
                "intersections_json": json.dumps(
                    row.get("intersections") or [], ensure_ascii=False
                ),
            }
        )
    return pd.DataFrame(normalized)


def build_pathway_inputs(
    root: Path,
    config: dict[str, Any],
) -> tuple[
    pd.DataFrame,
    pd.DataFrame,
    dict[str, set[str]],
    dict[str, set[str]],
    dict[str, str],
    dict[str, Any],
]:
    selection = config.get("candidate_selection") or {}
    top_n = int(selection.get("top_n_per_comparison", 100))
    require_negative = bool(selection.get("require_negative_delta", True))
    exclude_broad = bool(selection.get("exclude_broad_dependency", True))
    exclude_low = bool(selection.get("exclude_low_sample_size", True))

    candidate_frames: list[pd.DataFrame] = []
    eligible_by_comparison: dict[str, set[str]] = {}
    identifier_map: dict[str, str] = {}
    comparison_meta = []

    for spec in config.get("comparisons") or []:
        label = str(spec["label"])
        cancer_id = str(spec["cancer_id"])
        comparison = str(spec["comparison"])
        key = f"{cancer_id}__{comparison}"
        path = root / "data/processed/depmap_genomewide" / f"{key}_genes.tsv"
        if not path.exists():
            raise FileNotFoundError(
                f"Missing genome-wide input for pathway analysis: {path}"
            )
        frame = pd.read_csv(path, sep="\t", low_memory=False)
        local_identifier_map = build_entrez_identifier_map(frame)
        for symbol, entrez in local_identifier_map.items():
            previous = identifier_map.get(symbol)
            if previous is not None and previous != entrez:
                raise ValueError(
                    f"Conflicting Entrez identifiers across comparisons for {symbol}: "
                    f"{previous} vs {entrez}"
                )
            identifier_map[symbol] = entrez
        eligible = eligible_gene_set(
            frame,
            exclude_broad_dependency=exclude_broad,
            exclude_low_sample_size=exclude_low,
        )
        selected = select_top_candidates(
            frame,
            comparison_label=label,
            top_n=top_n,
            require_negative_delta=require_negative,
            exclude_broad_dependency=exclude_broad,
            exclude_low_sample_size=exclude_low,
        )
        candidate_frames.append(selected)
        eligible_by_comparison[label] = eligible
        comparison_meta.append(
            {
                "label": label,
                "cancer_id": cancer_id,
                "comparison": comparison,
                "input": str(path),
                "genes_total_n": int(len(frame)),
                "eligible_background_n": int(len(eligible)),
                "selected_top_n": int(len(selected)),
            }
        )

    candidate_long = pd.concat(candidate_frames, ignore_index=True) if candidate_frames else pd.DataFrame()
    recurrence = build_candidate_recurrence(candidate_long)

    set_cfg = config.get("analysis_sets") or {}
    thresholds = {
        name: int((details or {}).get("minimum_comparisons", 1))
        for name, details in set_cfg.items()
    }
    backgrounds = build_background_sets(eligible_by_comparison, thresholds)
    candidate_sets = {
        name: set(
            recurrence.loc[
                recurrence["comparisons_n"] >= minimum,
                "gene_symbol",
            ].astype(str)
        )
        for name, minimum in thresholds.items()
    }

    meta = {
        "analysis_version": PATHWAY_ANALYSIS_VERSION,
        "candidate_selection": selection,
        "comparisons": comparison_meta,
        "analysis_sets": {
            name: {
                "minimum_comparisons": thresholds[name],
                "query_genes_n": len(candidate_sets.get(name, set())),
                "background_genes_n": len(backgrounds.get(name, set())),
            }
            for name in thresholds
        },
    }
    return candidate_long, recurrence, candidate_sets, backgrounds, identifier_map, meta


def analyze_pathways(
    root: Path,
    config_path: Path | None = None,
) -> tuple[
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    dict[str, dict[str, Any]],
    dict[str, Any],
]:
    config_path = config_path or (root / "config/pathways.yaml")
    config = load_pathway_config(config_path)
    candidate_long, recurrence, candidate_sets, backgrounds, identifier_map, meta = build_pathway_inputs(
        root, config
    )

    enrichment_cfg = config.get("enrichment") or {}
    resolution_cfg = enrichment_cfg.get("identifier_resolution") or {}
    organism = str(config.get("organism") or "hsapiens")
    source_versions_path = root / "config/source_versions.yaml"
    source_versions = (
        yaml.safe_load(source_versions_path.read_text(encoding="utf-8")) or {}
        if source_versions_path.exists()
        else {}
    )
    gprofiler_cfg = source_versions.get("gprofiler") or {}
    profile_endpoint = str(gprofiler_cfg.get("profile_endpoint") or DEFAULT_PROFILE_ENDPOINT)
    convert_endpoint = str(gprofiler_cfg.get("convert_endpoint") or DEFAULT_CONVERT_ENDPOINT)
    versions_endpoint = str(
        gprofiler_cfg.get("data_versions_endpoint") or DEFAULT_VERSIONS_ENDPOINT
    )
    sources = [str(x) for x in enrichment_cfg.get("sources") or ["GO:BP", "REAC", "KEGG", "CORUM"]]
    method = str(enrichment_cfg.get("significance_threshold_method") or "g_SCS")
    user_threshold = float(enrichment_cfg.get("user_threshold", 0.05))
    convert_numeric_ns = str(resolution_cfg.get("numeric_ns") or "ENTREZGENE_ACC")
    convert_chunk_size = max(100, int(resolution_cfg.get("chunk_size", 1000)))
    excluded_term_ids = {str(x) for x in enrichment_cfg.get("excluded_term_ids") or []}

    # Convert only genes that can enter one of the current analysis backgrounds.
    needed_symbols = set().union(*backgrounds.values()) if backgrounds else set()
    needed_symbols.update(set().union(*candidate_sets.values()) if candidate_sets else set())
    source_missing = sorted(g for g in needed_symbols if g not in identifier_map)
    if source_missing:
        raise ValueError(
            "Missing Entrez identifiers before canonical Ensembl resolution: "
            + ", ".join(source_missing[:20])
        )
    symbol_to_entrez = {g: identifier_map[g] for g in sorted(needed_symbols)}
    convert_ids = sorted(set(symbol_to_entrez.values()), key=int)
    convert_rows, conversion_raw = _convert_in_chunks(
        convert_ids,
        organism=organism,
        endpoint=convert_endpoint,
        numeric_ns=convert_numeric_ns,
        chunk_size=convert_chunk_size,
    )
    identifier_resolution = resolve_ensembl_identifiers(symbol_to_entrez, convert_rows)
    resolution_by_symbol = identifier_resolution.set_index("gene_symbol", drop=False)

    raw_responses: dict[str, dict[str, Any]] = dict(conversion_raw)
    enrichment_frames: list[pd.DataFrame] = []
    background_rows: list[dict[str, Any]] = []

    for set_name, query_genes in candidate_sets.items():
        background = backgrounds.get(set_name, set())
        query_rows = identifier_resolution[
            identifier_resolution["gene_symbol"].isin(query_genes)
        ].copy()
        background_resolved = identifier_resolution[
            identifier_resolution["gene_symbol"].isin(background)
        ].copy()

        query_unresolved = query_rows[query_rows["ensembl_gene_id"].isna()]
        background_unresolved = background_resolved[background_resolved["ensembl_gene_id"].isna()]
        query_ensg = sorted(
            set(query_rows["ensembl_gene_id"].dropna().astype(str))
        )
        background_ensg = sorted(
            set(background_resolved["ensembl_gene_id"].dropna().astype(str))
        )
        if not query_ensg:
            continue
        if not query_unresolved.empty:
            details_text = ", ".join(
                f"{r.gene_symbol}({r.mapping_status})"
                for r in query_unresolved.itertuples(index=False)
            )
            raise ValueError(
                f"Unresolved canonical Ensembl identifiers for pathway query {set_name}: "
                + details_text[:1000]
            )
        missing_from_background = sorted(set(query_ensg) - set(background_ensg))
        if missing_from_background:
            raise ValueError(
                f"Resolved pathway query {set_name} is not contained in its background: "
                + ", ".join(missing_from_background[:20])
            )

        ensg_to_symbol = {
            str(r.ensembl_gene_id): str(r.gene_symbol)
            for r in background_resolved.dropna(subset=["ensembl_gene_id"]).itertuples(index=False)
        }
        details = meta.setdefault("analysis_sets", {}).setdefault(set_name, {})
        details.update(
            {
                "source_identifier_type": "ENTREZGENE",
                "profile_identifier_type": "ENSG",
                "query_identifiers_n": len(query_ensg),
                "background_identifiers_n": len(background_ensg),
                "query_unresolved_symbols": query_unresolved["gene_symbol"].astype(str).tolist(),
                "background_unresolved_symbols_n": int(len(background_unresolved)),
                "background_unresolved_symbols_sample": background_unresolved["gene_symbol"].astype(str).head(20).tolist(),
            }
        )

        for gene in sorted(background):
            row = resolution_by_symbol.loc[gene]
            background_rows.append(
                {
                    "analysis_set": set_name,
                    "gene_symbol": gene,
                    "entrez_gene_id": row["entrez_gene_id"],
                    "ensembl_gene_id": row["ensembl_gene_id"],
                    "mapping_status": row["mapping_status"],
                }
            )

        response = run_gprofiler_profile(
            query_ensg,
            background_ensg,
            organism=organism,
            sources=sources,
            significance_threshold_method=method,
            user_threshold=user_threshold,
            all_results=bool(enrichment_cfg.get("all_results", True)),
            domain_scope=str(enrichment_cfg.get("domain_scope") or "custom"),
            no_evidences=bool(enrichment_cfg.get("no_evidences", False)),
            numeric_ns="",
            endpoint=profile_endpoint,
        )
        raw_responses[set_name] = response
        normalized = normalize_gprofiler_result(
            response,
            set_name,
            excluded_term_ids=excluded_term_ids,
            input_id_to_symbol=ensg_to_symbol,
        )
        if not normalized.empty:
            enrichment_frames.append(normalized)
            details["gprofiler_query_sizes"] = sorted(
                {int(x) for x in normalized["query_size"].dropna().tolist()}
            )
            details["gprofiler_effective_domain_sizes"] = sorted(
                {int(x) for x in normalized["effective_domain_size"].dropna().tolist()}
            )

    enrichment = (
        pd.concat(enrichment_frames, ignore_index=True)
        if enrichment_frames
        else pd.DataFrame(
            columns=[
                "analysis_set",
                "source",
                "term_id",
                "term_name",
                "description",
                "p_value_adjusted",
                "significant",
                "query_size",
                "term_size",
                "intersection_size",
                "effective_domain_size",
                "precision",
                "recall",
                "parents_json",
                "intersecting_input_ids_json",
                "intersecting_gene_symbols_json",
                "intersections_json",
            ]
        )
    )
    backgrounds_long = pd.DataFrame(background_rows)
    versions = fetch_gprofiler_versions(
        organism=organism, endpoint=versions_endpoint
    )
    resolution_counts = {
        str(k): int(v)
        for k, v in identifier_resolution["mapping_status"].value_counts().to_dict().items()
    }
    meta.update(
        {
            "config_path": str(config_path),
            "organism": organism,
            "sources": sources,
            "significance_threshold_method": method,
            "user_threshold": user_threshold,
            "source_identifier_type": "ENTREZGENE",
            "profile_identifier_type": "ENSG",
            "identifier_resolution": {
                "service": "g:Profiler g:Convert",
                "target_namespace": "ENSG",
                "numeric_ns": convert_numeric_ns,
                "convert_endpoint": convert_endpoint,
                "chunk_size": convert_chunk_size,
                "input_symbols_n": len(symbol_to_entrez),
                "input_entrez_ids_n": len(convert_ids),
                "resolution_status_counts": resolution_counts,
            },
            "excluded_term_ids": sorted(excluded_term_ids),
            "gprofiler_profile_endpoint": profile_endpoint,
            "gprofiler_versions_endpoint": versions_endpoint,
            "gprofiler_versions": versions,
        }
    )
    return (
        candidate_long,
        recurrence,
        backgrounds_long,
        identifier_resolution,
        enrichment,
        raw_responses,
        meta,
    )
