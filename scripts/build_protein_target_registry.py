from __future__ import annotations

import argparse
import csv
import http.client
import io
import json
import ssl
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
TARGET_CATALOG = ROOT / "data" / "runtime" / "pharmacology" / "target_catalog.parquet"
GENE_REFERENCE = ROOT / "data" / "processed" / "gene_explorer" / "gene_reference.parquet"
OUTPUT_DIR = ROOT / "data" / "processed" / "target_registry"
OUTPUT = OUTPUT_DIR / "protein_targets.parquet"
MANIFEST = OUTPUT_DIR / "protein_target_registry_manifest.json"
CACHE = OUTPUT_DIR / "uniprot_target_cache.json"
QC = ROOT / "outputs" / "qc" / "protein_target_registry_unresolved.tsv"

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_FIELDS = "accession,id,protein_name,gene_primary,gene_names,protein_families,length"


class UniProtReleaseChanged(RuntimeError):
    pass


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


def _json_list(value: Any) -> list[str]:
    if isinstance(value, list):
        return [str(x).strip() for x in value if str(x).strip()]
    text = _text(value)
    if not text:
        return []
    try:
        parsed = json.loads(text)
    except (json.JSONDecodeError, TypeError, ValueError):
        return []
    if not isinstance(parsed, list):
        return []
    out: list[str] = []
    for item in parsed:
        item_text = str(item).strip()
        if item_text and item_text not in out:
            out.append(item_text)
    return out


def _preferred_protein_name(raw_name: str) -> str | None:
    text = _text(raw_name)
    if not text:
        return None
    # UniProt TSV renders the recommended name first, followed by alternative/short
    # names in parentheses. The complete raw field is retained separately.
    preferred = text.split(" (", 1)[0].strip()
    return preferred or text


def _pick(row: dict[str, str], *names: str) -> str:
    for name in names:
        value = _text(row.get(name))
        if value:
            return value
    return ""


def _load_cache() -> tuple[dict[str, dict[str, str]], set[str], str | None]:
    if not CACHE.exists():
        return {}, set(), None
    try:
        payload = json.loads(CACHE.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        print(f"WARNING: cannot read UniProt cache {CACHE.relative_to(ROOT)}; starting with an empty cache")
        return {}, set(), None
    raw_records = payload.get("records") or {}
    records = {
        str(accession): {str(k): str(v) for k, v in row.items()}
        for accession, row in raw_records.items()
        if isinstance(row, dict)
    }
    no_result = {str(x) for x in (payload.get("no_result_accessions") or [])}
    release = _text(payload.get("uniprot_release")) or None
    return records, no_result, release


def _save_cache(
    records: dict[str, dict[str, str]],
    no_result: set[str],
    release: str | None,
) -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    payload = {
        "cache_contract": "mcl-uniprot-target-cache-v1",
        "updated_at": _now(),
        "uniprot_release": release,
        "records": records,
        "no_result_accessions": sorted(no_result),
    }
    temp = CACHE.with_suffix(".json.tmp")
    temp.write_text(json.dumps(payload, ensure_ascii=False), encoding="utf-8")
    temp.replace(CACHE)


def _merge_release(current: str | None, incoming: str | None) -> str | None:
    incoming = _text(incoming) or None
    if current and incoming and current != incoming:
        raise UniProtReleaseChanged(
            f"UniProt release changed during a resumable build: cached={current}, current={incoming}. "
            f"Delete {CACHE.relative_to(ROOT)} and rerun to create a single-release snapshot."
        )
    return current or incoming


def _fetch_batch(
    accessions: list[str],
    timeout: int,
    retries: int,
) -> tuple[list[dict[str, str]], str | None]:
    if not accessions:
        return [], None
    clauses = " OR ".join(f"accession:{accession}" for accession in accessions)
    query = f"({clauses}) AND organism_id:9606 AND reviewed:true"
    url = f"{UNIPROT_SEARCH}?{urlencode({'query': query, 'format': 'tsv', 'fields': UNIPROT_FIELDS, 'size': 500})}"
    last_error: Exception | None = None
    attempts = max(1, int(retries))

    for attempt in range(attempts):
        try:
            request = Request(
                url,
                headers={
                    "User-Agent": "MasterCancerLandscape/0.1 protein-target-registry",
                    "Accept": "text/tab-separated-values",
                    "Connection": "close",
                },
            )
            with urlopen(request, timeout=max(10, int(timeout))) as response:
                payload = response.read().decode("utf-8")
                release = response.headers.get("X-UniProt-Release") or response.headers.get("x-uniprot-release")
            reader = csv.DictReader(io.StringIO(payload), delimiter="\t")
            return [dict(row) for row in reader], release
        except (
            HTTPError,
            URLError,
            TimeoutError,
            ConnectionError,
            OSError,
            ssl.SSLError,
            http.client.HTTPException,
            UnicodeDecodeError,
        ) as exc:
            last_error = exc
            if attempt + 1 < attempts:
                delay = min(30, 2 ** attempt)
                print(
                    f"WARNING: UniProt request failed ({type(exc).__name__}: {exc}); "
                    f"retry {attempt + 2}/{attempts} in {delay}s"
                )
                time.sleep(delay)

    raise RuntimeError(f"UniProtKB request failed after {attempts} attempts: {last_error}")


def _fetch_accessions(
    accessions: list[str],
    batch_size: int,
    timeout: int,
    retries: int,
) -> tuple[dict[str, dict[str, str]], str | None]:
    records, no_result, release = _load_cache()
    wanted = set(accessions)
    records = {k: v for k, v in records.items() if k in wanted}
    no_result &= wanted

    completed = set(records) | no_result
    pending = [accession for accession in accessions if accession not in completed]
    if completed:
        print(
            f"UniProt cache: {len(completed)}/{len(accessions)} accessions already completed; "
            f"{len(pending)} remain"
        )
    if not pending:
        return records, release

    batch_size = max(1, min(int(batch_size), 80))
    network_failures: list[str] = []

    for start in range(0, len(pending), batch_size):
        batch = pending[start : start + batch_size]
        batch_completed: set[str] = set()
        try:
            rows, batch_release = _fetch_batch(batch, timeout, retries)
            release = _merge_release(release, batch_release)
            returned: set[str] = set()
            for row in rows:
                accession = _pick(row, "Entry", "accession")
                if accession:
                    records[accession] = row
                    returned.add(accession)
            for accession in batch:
                if accession not in returned:
                    no_result.add(accession)
                batch_completed.add(accession)
            _save_cache(records, no_result, release)
        except UniProtReleaseChanged:
            _save_cache(records, no_result, release)
            raise
        except RuntimeError as batch_error:
            print(
                f"WARNING: batch request failed after retries ({batch_error}). "
                "Falling back to one accession at a time."
            )
            for accession in batch:
                try:
                    one_rows, one_release = _fetch_batch([accession], timeout, retries)
                    release = _merge_release(release, one_release)
                    returned = False
                    for row in one_rows:
                        row_accession = _pick(row, "Entry", "accession")
                        if row_accession:
                            records[row_accession] = row
                            if row_accession == accession:
                                returned = True
                    if not returned:
                        no_result.add(accession)
                    batch_completed.add(accession)
                    _save_cache(records, no_result, release)
                except UniProtReleaseChanged:
                    _save_cache(records, no_result, release)
                    raise
                except RuntimeError as accession_error:
                    network_failures.append(accession)
                    print(
                        f"WARNING: could not resolve {accession} after retries: {accession_error}. "
                        "Progress is saved; this accession will be retried on the next run."
                    )

        current_completed = len((set(records) | no_result) & wanted)
        print(f"UniProt completed {current_completed}/{len(accessions)} accessions")

    _save_cache(records, no_result, release)
    still_pending = [a for a in accessions if a not in records and a not in no_result]
    if still_pending:
        examples = ", ".join(still_pending[:10])
        raise RuntimeError(
            f"UniProt network retrieval remains incomplete for {len(still_pending)} accessions "
            f"({examples}{' ...' if len(still_pending) > 10 else ''}). "
            "All successful progress has been cached. Rerun the same command; only unfinished accessions will be requested."
        )

    return records, release


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Build the MCL protein target registry. Pharmacology target annotations remain gene-level; "
            "UniProtKB/Swiss-Prot metadata is added only as an explicit gene-to-protein mapping layer."
        )
    )
    parser.add_argument("--batch-size", type=int, default=40)
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--retries", type=int, default=6)
    parser.add_argument("--limit", type=int, default=None, help="Optional target limit for smoke testing.")
    parser.add_argument(
        "--reset-cache",
        action="store_true",
        help="Delete the resumable UniProt cache before downloading a fresh single-release snapshot.",
    )
    args = parser.parse_args()

    if args.reset_cache and CACHE.exists():
        CACHE.unlink()
        print(f"Removed UniProt cache: {CACHE.relative_to(ROOT)}")

    if not TARGET_CATALOG.exists():
        raise SystemExit("Missing target_catalog.parquet. Run scripts/build_pharmacology_layer.py first.")
    if not GENE_REFERENCE.exists():
        raise SystemExit(
            "Missing data/processed/gene_explorer/gene_reference.parquet. "
            "Run scripts/build_gene_reference_snapshot.py first."
        )

    targets = pd.read_parquet(TARGET_CATALOG, columns=["target_gene"])
    target_genes = sorted(set(targets["target_gene"].dropna().astype(str).str.upper().str.strip()))
    if args.limit is not None:
        target_genes = target_genes[: max(1, int(args.limit))]

    reference = pd.read_parquet(GENE_REFERENCE)
    if "gene_symbol" not in reference.columns or "uniprot_swissprot_ids_json" not in reference.columns:
        raise SystemExit("gene_reference.parquet does not contain the required Swiss-Prot mapping columns.")
    reference["gene_symbol"] = reference["gene_symbol"].astype(str).str.upper().str.strip()
    reference = reference[reference["gene_symbol"].isin(target_genes)].copy()
    reference_lookup = {str(row["gene_symbol"]): row for _, row in reference.iterrows()}

    gene_accessions: dict[str, list[str]] = {}
    all_accessions: list[str] = []
    for gene in target_genes:
        row = reference_lookup.get(gene)
        ids = _json_list(row.get("uniprot_swissprot_ids_json")) if row is not None else []
        gene_accessions[gene] = ids
        for accession in ids:
            if accession not in all_accessions:
                all_accessions.append(accession)

    print(f"Pharmacology target genes: {len(target_genes)}")
    print(f"Unique Swiss-Prot accessions to resolve: {len(all_accessions)}")
    try:
        records_by_accession, uniprot_release = _fetch_accessions(
            all_accessions,
            batch_size=args.batch_size,
            timeout=args.timeout,
            retries=args.retries,
        )
    except UniProtReleaseChanged as exc:
        raise SystemExit(str(exc)) from exc
    except RuntimeError as exc:
        raise SystemExit(str(exc)) from exc

    retrieved_at = _now()
    rows: list[dict[str, Any]] = []
    unresolved: list[dict[str, Any]] = []
    for gene in target_genes:
        accessions = gene_accessions.get(gene, [])
        protein_records = [records_by_accession[a] for a in accessions if a in records_by_accession]
        names = [_preferred_protein_name(_pick(row, "Protein names", "Protein Names")) for row in protein_records]
        names = [name for name in names if name]

        if not accessions:
            status = "no_swissprot_mapping"
        elif not protein_records:
            status = "uniprot_unresolved"
        elif len(accessions) == 1 and len(protein_records) == 1:
            status = "unique_swissprot"
        else:
            status = "multiple_swissprot"

        unique = status == "unique_swissprot"
        record = protein_records[0] if unique else {}
        primary_accession = _pick(record, "Entry", "accession") if unique else ""
        raw_name = _pick(record, "Protein names", "Protein Names") if unique else ""
        preferred_name = _preferred_protein_name(raw_name) if unique else None
        rows.append(
            {
                "target_gene": gene,
                "source_resolution": "gene_mapped",
                "protein_mapping_status": status,
                "protein_preferred_name": preferred_name,
                "protein_name_raw": raw_name or None,
                "protein_names_json": json.dumps(names, ensure_ascii=False),
                "uniprot_primary_accession": primary_accession or None,
                "uniprot_accessions_json": json.dumps(accessions, ensure_ascii=False),
                "uniprot_entry_name": (_pick(record, "Entry Name", "id") or None) if unique else None,
                "uniprot_gene_primary": (_pick(record, "Gene Names (primary)", "Gene Names (Primary)") or None) if unique else None,
                "uniprot_gene_names": (_pick(record, "Gene Names") or None) if unique else None,
                "protein_families": (_pick(record, "Protein families", "Protein Families") or None) if unique else None,
                "protein_length": pd.to_numeric(_pick(record, "Length"), errors="coerce") if unique else None,
                "mapping_source": "MyGene.info Swiss-Prot mapping + UniProtKB reviewed entry",
                "mapping_retrieved_at": retrieved_at,
                "uniprot_release": uniprot_release,
            }
        )
        if status != "unique_swissprot":
            unresolved.append(
                {
                    "target_gene": gene,
                    "protein_mapping_status": status,
                    "uniprot_accessions_json": json.dumps(accessions, ensure_ascii=False),
                    "resolved_accessions_json": json.dumps(
                        [a for a in accessions if a in records_by_accession], ensure_ascii=False
                    ),
                }
            )

    frame = pd.DataFrame(rows).sort_values("target_gene").reset_index(drop=True)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    QC.parent.mkdir(parents=True, exist_ok=True)
    frame.to_parquet(OUTPUT, index=False, compression="zstd")
    pd.DataFrame(unresolved).to_csv(QC, sep="\t", index=False)

    counts = frame["protein_mapping_status"].value_counts().to_dict() if not frame.empty else {}
    manifest = {
        "contract": "mcl-protein-target-registry-v1",
        "built_at": _now(),
        "mapping_retrieved_at": retrieved_at,
        "source_resolution": "gene_mapped",
        "target_genes_n": len(target_genes),
        "unique_swissprot_accessions_n": len(all_accessions),
        "mapping_status_counts": {str(k): int(v) for k, v in counts.items()},
        "uniprot_release": uniprot_release,
        "uniprot_endpoint": UNIPROT_SEARCH,
        "cache": str(CACHE.relative_to(ROOT)),
        "output": str(OUTPUT.relative_to(ROOT)),
        "qc": str(QC.relative_to(ROOT)),
        "interpretation_ru": (
            "Фармакологический источник по-прежнему разрешён на уровне гена. Название белка и UniProt ID "
            "добавляются как справочное сопоставление gene -> reviewed protein и не доказывают, что источник "
            "измерял конкретную изоформу или протеоформу."
        ),
    }
    MANIFEST.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")

    print("MCL Protein Target Registry v1")
    print(f"Target genes: {manifest['target_genes_n']}")
    for status, count in sorted(manifest["mapping_status_counts"].items()):
        print(f"  {status}: {count}")
    print(f"UniProt release: {uniprot_release or 'not reported'}")
    print(f"Wrote {OUTPUT.relative_to(ROOT)}")
    print(f"QC {QC.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
