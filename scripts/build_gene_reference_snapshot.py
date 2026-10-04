from __future__ import annotations

import argparse
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
API_DIR = ROOT / "apps" / "api"
if str(API_DIR) not in sys.path:
    sys.path.insert(0, str(API_DIR))

from mcl_api.gene_explorer import GeneExplorerStore  # noqa: E402
from mcl_api.store import MCLDataStore  # noqa: E402


MYGENE_QUERY_URL = "https://mygene.info/v3/query"
FIELDS = (
    "symbol,name,alias,entrezgene,ensembl.gene,uniprot.Swiss-Prot,type_of_gene,"
    "go,pathway,pdb,interpro,pfam"
)


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _as_list(value: Any) -> list[Any]:
    if value is None:
        return []
    return value if isinstance(value, list) else [value]


def _strings(value: Any) -> list[str]:
    out: list[str] = []
    for item in _as_list(value):
        if isinstance(item, dict):
            continue
        text = str(item).strip()
        if text and text not in out:
            out.append(text)
    return out


def _ensembl_ids(value: Any) -> list[str]:
    out: list[str] = []
    for item in _as_list(value):
        if isinstance(item, dict):
            candidate = item.get("gene")
        else:
            candidate = item
        text = str(candidate or "").strip()
        if text and text not in out:
            out.append(text)
    return out


def _uniprot_ids(value: Any) -> list[str]:
    if isinstance(value, dict):
        value = value.get("Swiss-Prot")
    return _strings(value)


def _term_rows(query_symbol: str, payload: dict[str, Any], retrieved_at: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    go = payload.get("go") or {}
    if isinstance(go, dict):
        for aspect, source in (("BP", "GO:BP"), ("MF", "GO:MF"), ("CC", "GO:CC")):
            for term in _as_list(go.get(aspect)):
                if not isinstance(term, dict):
                    continue
                term_id = str(term.get("id") or "").strip()
                term_name = str(term.get("term") or "").strip()
                if not term_id and not term_name:
                    continue
                rows.append(
                    {
                        "gene_symbol": query_symbol,
                        "source": source,
                        "term_id": term_id or None,
                        "term_name": term_name or term_id or None,
                        "source_version": None,
                        "retrieved_at": retrieved_at,
                        "aggregator": "MyGene.info v3",
                        "aggregator_endpoint": MYGENE_QUERY_URL,
                    }
                )

    pathways = payload.get("pathway") or {}
    if isinstance(pathways, dict):
        source_map = {
            "reactome": "REAC",
            "kegg": "KEGG",
            "wikipathways": "WikiPathways",
            "biocarta": "BioCarta",
            "netpath": "NetPath",
            "pid": "PID",
        }
        for raw_source, entries in pathways.items():
            source = source_map.get(str(raw_source).lower(), str(raw_source))
            for term in _as_list(entries):
                if isinstance(term, dict):
                    term_id = str(term.get("id") or term.get("term") or "").strip()
                    term_name = str(term.get("name") or term.get("term") or "").strip()
                else:
                    term_id = ""
                    term_name = str(term or "").strip()
                if not term_id and not term_name:
                    continue
                rows.append(
                    {
                        "gene_symbol": query_symbol,
                        "source": source,
                        "term_id": term_id or None,
                        "term_name": term_name or term_id or None,
                        "source_version": None,
                        "retrieved_at": retrieved_at,
                        "aggregator": "MyGene.info v3",
                        "aggregator_endpoint": MYGENE_QUERY_URL,
                    }
                )
    return rows


def _post_batch(symbols: list[str], timeout: int, retries: int) -> list[dict[str, Any]]:
    data = urlencode(
        {
            "q": ",".join(symbols),
            "scopes": "symbol",
            "species": "human",
            "fields": FIELDS,
        }
    ).encode("utf-8")
    request = Request(
        MYGENE_QUERY_URL,
        data=data,
        method="POST",
        headers={
            "Content-Type": "application/x-www-form-urlencoded",
            "User-Agent": "MasterCancerLandscape/0.1 gene-reference-snapshot",
        },
    )
    last_error: Exception | None = None
    for attempt in range(max(1, retries)):
        try:
            with urlopen(request, timeout=timeout) as response:
                payload = json.loads(response.read().decode("utf-8"))
            if not isinstance(payload, list):
                raise RuntimeError("Unexpected MyGene response: batch query did not return a list")
            return [row for row in payload if isinstance(row, dict)]
        except (HTTPError, URLError, TimeoutError, json.JSONDecodeError, RuntimeError) as exc:
            last_error = exc
            if attempt + 1 < retries:
                time.sleep(2 ** attempt)
    raise RuntimeError(f"MyGene batch query failed after {retries} attempts: {last_error}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Build a local, timestamped human-gene reference snapshot for the current MCL gene universe using MyGene.info v3. "
            "The snapshot stores source terms and aggregator provenance; it does not silently replace direct GO/Reactome/KEGG releases."
        )
    )
    parser.add_argument("--batch-size", type=int, default=400, help="Gene symbols per MyGene POST request (default: 400).")
    parser.add_argument("--timeout", type=int, default=120, help="HTTP timeout per request in seconds.")
    parser.add_argument("--retries", type=int, default=3, help="Retries per batch.")
    parser.add_argument("--limit", type=int, default=None, help="Optional gene limit for smoke tests.")
    args = parser.parse_args()

    store = MCLDataStore(ROOT)
    explorer = GeneExplorerStore(ROOT, store)
    catalog = explorer.catalog_frame()
    if catalog.empty or "gene_symbol" not in catalog.columns:
        raise SystemExit("Current MCL gene universe is empty. Build genome-wide/Gene Explorer inputs first.")
    symbols = sorted(set(catalog["gene_symbol"].dropna().astype(str).str.upper()))
    if args.limit is not None:
        symbols = symbols[: max(1, int(args.limit))]

    output_dir = ROOT / "data" / "processed" / "gene_explorer"
    output_dir.mkdir(parents=True, exist_ok=True)
    retrieved_at = _utc_now()
    references: list[dict[str, Any]] = []
    terms: list[dict[str, Any]] = []
    unresolved: list[str] = []

    batch_size = max(1, min(int(args.batch_size), 1000))
    print(f"MCL gene universe to annotate: {len(symbols)}")
    print(f"Source: {MYGENE_QUERY_URL}")

    for start in range(0, len(symbols), batch_size):
        batch = symbols[start : start + batch_size]
        payload = _post_batch(batch, timeout=max(10, int(args.timeout)), retries=max(1, int(args.retries)))
        by_query: dict[str, list[dict[str, Any]]] = {}
        for row in payload:
            query = str(row.get("query") or "").strip().upper()
            if query:
                by_query.setdefault(query, []).append(row)

        for symbol in batch:
            candidates = [x for x in by_query.get(symbol, []) if not x.get("notfound")]
            exact = [
                x for x in candidates
                if int(x.get("taxid") or 9606) == 9606 and str(x.get("symbol") or "").strip().upper() == symbol
            ]
            hit = exact[0] if exact else (candidates[0] if candidates else None)
            if hit is None:
                unresolved.append(symbol)
                continue
            aliases = _strings(hit.get("alias"))
            ensembl = _ensembl_ids((hit.get("ensembl") or {}).get("gene") if isinstance(hit.get("ensembl"), dict) else hit.get("ensembl"))
            uniprot = _uniprot_ids(hit.get("uniprot"))
            references.append(
                {
                    "gene_symbol": symbol,
                    "matched_symbol": str(hit.get("symbol") or symbol),
                    "gene_name": hit.get("name"),
                    "aliases_json": json.dumps(aliases, ensure_ascii=False),
                    "entrez_gene_id": hit.get("entrezgene"),
                    "ensembl_gene_ids_json": json.dumps(ensembl, ensure_ascii=False),
                    "uniprot_swissprot_ids_json": json.dumps(uniprot, ensure_ascii=False),
                    "type_of_gene": hit.get("type_of_gene"),
                    "pdb_ids_json": json.dumps(_strings(hit.get("pdb")), ensure_ascii=False),
                    "interpro_json": json.dumps(_as_list(hit.get("interpro")), ensure_ascii=False),
                    "pfam_json": json.dumps(_as_list(hit.get("pfam")), ensure_ascii=False),
                    "retrieved_at": retrieved_at,
                    "aggregator": "MyGene.info v3",
                    "aggregator_endpoint": MYGENE_QUERY_URL,
                }
            )
            terms.extend(_term_rows(symbol, hit, retrieved_at))

        done = min(start + batch_size, len(symbols))
        print(f"Annotated {done}/{len(symbols)} genes")

    reference_frame = pd.DataFrame(references)
    term_frame = pd.DataFrame(terms)
    if not reference_frame.empty:
        reference_frame = reference_frame.drop_duplicates("gene_symbol", keep="first").sort_values("gene_symbol")
    if not term_frame.empty:
        term_frame = term_frame.drop_duplicates(["gene_symbol", "source", "term_id", "term_name"]).sort_values(["gene_symbol", "source", "term_name"])

    reference_path = output_dir / "gene_reference.parquet"
    terms_path = output_dir / "gene_reference_terms.parquet"
    reference_frame.to_parquet(reference_path, index=False, compression="zstd")
    term_frame.to_parquet(terms_path, index=False, compression="zstd")

    manifest = {
        "built_at": _utc_now(),
        "retrieved_at": retrieved_at,
        "aggregator": "MyGene.info v3",
        "endpoint": MYGENE_QUERY_URL,
        "documentation": "https://docs.mygene.info/en/latest/doc/query_service.html",
        "query_scope": "symbol",
        "species": "human",
        "fields": FIELDS,
        "gene_universe_n": len(symbols),
        "resolved_genes_n": int(reference_frame["gene_symbol"].nunique()) if not reference_frame.empty else 0,
        "unresolved_genes_n": len(unresolved),
        "unresolved_genes": unresolved,
        "term_rows_n": int(len(term_frame)),
        "files": [reference_path.name, terms_path.name],
        "provenance_note": (
            "This is a timestamped secondary aggregation snapshot. GO/pathway term IDs and names are retained, but direct source-release versions may be unavailable in MyGene payloads and therefore remain null rather than being inferred."
        ),
    }
    (output_dir / "gene_reference_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    print(f"Resolved genes: {manifest['resolved_genes_n']}")
    print(f"Unresolved genes: {manifest['unresolved_genes_n']}")
    print(f"Formal GO/pathway rows: {manifest['term_rows_n']}")
    print("Wrote data/processed/gene_explorer/gene_reference.parquet")
    print("Wrote data/processed/gene_explorer/gene_reference_terms.parquet")
    print("Wrote data/processed/gene_explorer/gene_reference_manifest.json")


if __name__ == "__main__":
    main()
