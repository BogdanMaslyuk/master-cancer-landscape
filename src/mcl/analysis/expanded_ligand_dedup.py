"""P0.3: conservative duplicate classification without losing original records."""
from __future__ import annotations

import hashlib
import json
from collections import defaultdict

from .expanded_ligand_measurements import clean, mirror_group_key


def _id(prefix: str, parts: tuple) -> str:
    value = "|".join(str(part) for part in parts)
    return prefix + hashlib.sha256(value.encode("utf-8")).hexdigest()[:18]


def classify_duplicates(rows: list[dict]) -> list[dict]:
    """Classify copies and ambiguous pairs; do not collapse independent assays.

    A shared DOI alone never proves duplication. This function is deliberately
    structure-string conservative until RDKit exact-identity curation is added.
    """
    outputs = [{
        **row, "duplicate_class": "unmatched_independent_candidate",
        "duplicate_cluster_id": "", "count_as_independent": "pending_curation",
        "duplicate_reason": "",
    } for row in rows]
    id_groups: dict[tuple, list[int]] = defaultdict(list)
    for i, row in enumerate(rows):
        key = (clean(row.get("source")), clean(row.get("source_record_id")),
               clean(row.get("endpoint")))
        if key[0] and key[1] and key[2]:
            id_groups[key].append(i)
    for key, indices in id_groups.items():
        if len(indices) < 2:
            continue
        # Exact IDs may still be paired with conflicting values in a damaged extract.
        signatures = {
            (clean(rows[i].get("target_gene")), clean(rows[i].get("original_smiles")),
             clean(rows[i].get("relation")), clean(rows[i].get("numeric_value_nm")))
            for i in indices
        }
        consistent = len(signatures) == 1
        cluster = _id("SRC-", key)
        for j, index in enumerate(sorted(indices)):
            outputs[index].update({
                "duplicate_class": "duplicate_same_source_record" if consistent else "conflicting_source_record",
                "duplicate_cluster_id": cluster,
                "duplicate_reason": "same provider, record ID and endpoint",
                "count_as_independent": ("pending_curation" if not consistent
                                         else "primary_candidate" if j == 0 else "no"),
            })

    mirrors: dict[str, list[int]] = defaultdict(list)
    for i, row in enumerate(rows):
        key = mirror_group_key(row)
        if key:
            mirrors[key].append(i)
    for key, indices in mirrors.items():
        if len({rows[i].get("source") for i in indices}) < 2:
            continue
        signatures = {
            (clean(rows[i].get("target_gene")), clean(rows[i].get("original_smiles")),
             clean(rows[i].get("endpoint")).upper(), clean(rows[i].get("relation")),
             clean(rows[i].get("numeric_value_nm")))
            for i in indices
        }
        consistent = len(signatures) == 1
        # The ChEMBL source record is the canonical pointer for explicit mirrors.
        canonical = sorted(indices, key=lambda i: (rows[i].get("source") != "ChEMBL", i))[0]
        for i in indices:
            outputs[i].update({
                "duplicate_class": "confirmed_source_mirror" if consistent else "conflicting_source_mirror",
                "duplicate_cluster_id": key,
                "duplicate_reason": "explicit shared ChEMBL original activity ID",
                "count_as_independent": (
                    "pending_curation" if not consistent
                    else "primary_candidate" if i == canonical else "no"
                ),
            })

    # Flag publication-level resemblance, but never merge merely on DOI.
    likely: dict[tuple, list[int]] = defaultdict(list)
    for i, row in enumerate(rows):
        if outputs[i]["duplicate_class"] not in {"unmatched_independent_candidate"}:
            continue
        doi = clean(row.get("doi")).lower()
        if not doi:
            continue
        signature = (
            doi, clean(row.get("target_gene")), clean(row.get("original_smiles")),
            clean(row.get("endpoint")).upper(), clean(row.get("relation")),
            clean(row.get("numeric_value_nm")),
        )
        likely[signature].append(i)
    for signature, indices in likely.items():
        if len(indices) < 2:
            continue
        # Identical values in one paper could still represent independent repeats.
        cluster = _id("PROB-", signature)
        for i in indices:
            outputs[i].update({
                "duplicate_class": "probable_publication_duplicate_requires_review",
                "duplicate_cluster_id": cluster,
                "duplicate_reason": "shared DOI/structure/target/endpoint/value; independent assay not established",
                "count_as_independent": "pending_curation",
            })
    return outputs


def summarize_duplicates(rows: list[dict]) -> dict:
    counts: dict[str, int] = {}
    for row in rows:
        name = row["duplicate_class"]
        counts[name] = counts.get(name, 0) + 1
    return {
        "total_measurement_rows": len(rows),
        "duplicate_class_counts": dict(sorted(counts.items())),
        "confirmed_no_double_count_n": sum(r["count_as_independent"] == "no" for r in rows),
        "unresolved_review_n": sum(r["count_as_independent"] == "pending_curation" for r in rows),
        "independent_experiment_count_verified": False,
        "note": "A primary candidate is not a final independent experiment until assay provenance is reviewed.",
    }
