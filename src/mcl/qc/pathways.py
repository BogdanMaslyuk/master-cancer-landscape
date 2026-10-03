from __future__ import annotations

import pandas as pd

from mcl.models.qc import QCRecord


def pathway_enrichment_qc(
    recurrence: pd.DataFrame,
    backgrounds: pd.DataFrame,
    enrichment: pd.DataFrame,
    meta: dict,
) -> list[QCRecord]:
    qc: list[QCRecord] = []
    entity = "M3.3"

    if recurrence.empty:
        qc.append(QCRecord(
            severity="ERROR",
            check="candidate_recurrence_nonempty",
            entity_id=entity,
            observed="0",
            expected=">0",
            message="No recurrent candidate table was produced.",
        ))
        return qc

    qc.append(QCRecord(
        severity="INFO",
        check="recurrent_gene_universe",
        entity_id=entity,
        observed=str(len(recurrence)),
        expected=">0",
        message="Number of unique genes observed in at least one per-comparison top candidate list.",
    ))

    analysis_sets = meta.get("analysis_sets") or {}
    for set_name, details in analysis_sets.items():
        query_n = int(details.get("query_genes_n", 0))
        background_n = int(details.get("background_genes_n", 0))
        if query_n == 0:
            qc.append(QCRecord(
                severity="ERROR",
                check="pathway_query_nonempty",
                entity_id=set_name,
                observed="0",
                expected=">0",
                message="Pathway-analysis query set is empty.",
            ))
        if background_n < query_n:
            qc.append(QCRecord(
                severity="ERROR",
                check="background_contains_query",
                entity_id=set_name,
                observed=f"background={background_n}; query={query_n}",
                expected="background >= query",
                message="Custom enrichment background cannot be smaller than its query set.",
            ))
        else:
            qc.append(QCRecord(
                severity="INFO",
                check="pathway_set_sizes",
                entity_id=set_name,
                observed=f"query={query_n}; background={background_n}",
                expected="eligibility-matched custom background",
                message="Pathway query and custom background sizes.",
            ))

        query_ids_n = int(details.get("query_identifiers_n", query_n))
        background_ids_n = int(details.get("background_identifiers_n", background_n))
        unresolved_query = list(details.get("query_unresolved_symbols") or [])
        unresolved_background_n = int(details.get("background_unresolved_symbols_n", 0))
        if unresolved_query:
            qc.append(QCRecord(
                severity="ERROR",
                check="pathway_query_identifier_mapping",
                entity_id=set_name,
                observed=", ".join(unresolved_query[:20]),
                expected="all query genes resolved to one canonical Ensembl gene",
                message="One or more pathway-query genes could not be resolved unambiguously to ENSG.",
            ))
        elif query_ids_n != query_n:
            qc.append(QCRecord(
                severity="WARNING",
                check="pathway_query_identifier_mapping",
                entity_id=set_name,
                observed=f"genes={query_n}; unique_ENSG={query_ids_n}",
                expected="one unique ENSG per query gene",
                message="Canonical mapping produced fewer unique ENSG identifiers than query genes.",
            ))
        else:
            qc.append(QCRecord(
                severity="INFO",
                check="pathway_query_identifier_mapping",
                entity_id=set_name,
                observed=str(query_ids_n),
                expected=str(query_n),
                message="All pathway-query genes resolved to one canonical Ensembl gene.",
            ))

        if unresolved_background_n:
            qc.append(QCRecord(
                severity="WARNING",
                check="pathway_background_identifier_mapping",
                entity_id=set_name,
                observed=str(unresolved_background_n),
                expected="0 unresolved background genes",
                message="Some eligible background genes were excluded because canonical ENSG mapping was unresolved.",
            ))
        elif background_ids_n != background_n:
            qc.append(QCRecord(
                severity="WARNING",
                check="pathway_background_identifier_mapping",
                entity_id=set_name,
                observed=f"genes={background_n}; unique_ENSG={background_ids_n}",
                expected="one unique ENSG per eligible background gene",
                message="Canonical mapping collapsed two or more eligible background symbols onto the same ENSG.",
            ))
        else:
            qc.append(QCRecord(
                severity="INFO",
                check="pathway_background_identifier_mapping",
                entity_id=set_name,
                observed=str(background_ids_n),
                expected=str(background_n),
                message="All eligible background genes resolved one-to-one to canonical Ensembl identifiers.",
            ))

        observed_query_sizes = [int(x) for x in details.get("gprofiler_query_sizes") or []]
        if observed_query_sizes and observed_query_sizes != [query_ids_n]:
            qc.append(QCRecord(
                severity="WARNING",
                check="gprofiler_query_size_match",
                entity_id=set_name,
                observed=str(observed_query_sizes),
                expected=str(query_ids_n),
                message="g:Profiler query size differs from the canonical ENSG query submitted by MCL.",
            ))
        elif observed_query_sizes:
            qc.append(QCRecord(
                severity="INFO",
                check="gprofiler_query_size_match",
                entity_id=set_name,
                observed=str(observed_query_sizes[0]),
                expected=str(query_ids_n),
                message="g:Profiler query size matches the canonical ENSG query submitted by MCL.",
            ))

        observed_domains = [int(x) for x in details.get("gprofiler_effective_domain_sizes") or []]
        if observed_domains and observed_domains != [background_ids_n]:
            qc.append(QCRecord(
                severity="WARNING",
                check="gprofiler_background_size_match",
                entity_id=set_name,
                observed=str(observed_domains),
                expected=str(background_ids_n),
                message="g:Profiler effective domain differs from the canonical ENSG custom background submitted by MCL.",
            ))
        elif observed_domains:
            qc.append(QCRecord(
                severity="INFO",
                check="gprofiler_background_size_match",
                entity_id=set_name,
                observed=str(observed_domains[0]),
                expected=str(background_ids_n),
                message="g:Profiler effective domain matches the canonical ENSG custom background.",
            ))

    if backgrounds.empty:
        qc.append(QCRecord(
            severity="ERROR",
            check="pathway_background_nonempty",
            entity_id=entity,
            observed="0",
            expected=">0",
            message="No custom statistical background genes were produced.",
        ))

    if enrichment.empty:
        qc.append(QCRecord(
            severity="WARNING",
            check="enrichment_results_nonempty",
            entity_id=entity,
            observed="0",
            expected=">0 terms if the external service returns annotated overlaps",
            message="g:Profiler returned no enrichment terms; this is not automatically an analysis failure.",
        ))
        return qc

    p = pd.to_numeric(enrichment.get("p_value_adjusted"), errors="coerce")
    if p.isna().any():
        qc.append(QCRecord(
            severity="ERROR",
            check="enrichment_pvalue_complete",
            entity_id=entity,
            observed=str(int(p.isna().sum())),
            expected="0 missing adjusted p-values",
            message="One or more enrichment rows lack a corrected p-value.",
        ))

    for set_name, rows in enrichment.groupby("analysis_set"):
        significant = rows["significant"].fillna(False).astype(bool)
        qc.append(QCRecord(
            severity="INFO",
            check="significant_pathway_terms",
            entity_id=str(set_name),
            observed=str(int(significant.sum())),
            expected="descriptive count",
            message="Significant functional terms after the configured g:Profiler correction method.",
        ))
        for source, source_rows in rows.groupby("source"):
            source_sig = source_rows["significant"].fillna(False).astype(bool)
            qc.append(QCRecord(
                severity="INFO",
                check="pathway_source_terms",
                entity_id=f"{set_name}|{source}",
                observed=f"all={len(source_rows)}; significant={int(source_sig.sum())}",
                expected="descriptive count",
                message="Returned and significant terms for this functional datasource.",
            ))

    return qc
