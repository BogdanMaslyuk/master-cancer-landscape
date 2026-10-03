from __future__ import annotations

from typing import Any

import httpx


DEFAULT_PROFILE_ENDPOINT = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
DEFAULT_CONVERT_ENDPOINT = "https://biit.cs.ut.ee/gprofiler/api/convert/convert/"
DEFAULT_VERSIONS_ENDPOINT = "https://biit.cs.ut.ee/gprofiler/api/util/data_versions"
USER_AGENT = "master-cancer-landscape/0.5.0-dev"


def fetch_gprofiler_versions(
    organism: str = "hsapiens",
    endpoint: str = DEFAULT_VERSIONS_ENDPOINT,
    timeout: float = 60.0,
) -> dict[str, Any]:
    """Fetch g:Profiler datasource-version metadata for provenance."""
    with httpx.Client(
        timeout=timeout,
        follow_redirects=True,
        headers={"User-Agent": USER_AGENT},
    ) as client:
        response = client.get(endpoint, params={"organism": organism})
        response.raise_for_status()
        data = response.json()
    if not isinstance(data, dict):
        raise RuntimeError("g:Profiler data_versions endpoint returned a non-object response")
    return data


def run_gprofiler_convert(
    query: list[str],
    *,
    organism: str = "hsapiens",
    target: str = "ENSG",
    numeric_ns: str = "ENTREZGENE_ACC",
    endpoint: str = DEFAULT_CONVERT_ENDPOINT,
    timeout: float = 120.0,
) -> dict[str, Any]:
    """Resolve identifiers explicitly with g:Convert before enrichment.

    MCL uses this as a deterministic identifier-normalization gate so that
    g:GOSt receives canonical Ensembl gene IDs rather than re-interpreting
    numeric identifiers internally during the enrichment request.
    """
    query_unique = list(dict.fromkeys(str(x).strip() for x in query if str(x).strip()))
    if not query_unique:
        raise ValueError("g:Profiler conversion query must contain at least one identifier")

    payload: dict[str, Any] = {
        "organism": organism,
        "query": query_unique,
        "target": target,
    }
    if numeric_ns:
        payload["numeric_ns"] = numeric_ns

    with httpx.Client(
        timeout=timeout,
        follow_redirects=True,
        headers={"User-Agent": USER_AGENT},
    ) as client:
        response = client.post(endpoint, json=payload)
        response.raise_for_status()
        data = response.json()

    if not isinstance(data, dict) or "result" not in data:
        raise RuntimeError("g:Profiler convert endpoint returned an unexpected response")
    return data


def run_gprofiler_profile(
    query: list[str],
    background: list[str],
    *,
    organism: str = "hsapiens",
    sources: list[str] | None = None,
    significance_threshold_method: str = "g_SCS",
    user_threshold: float = 0.05,
    all_results: bool = True,
    domain_scope: str = "custom",
    no_evidences: bool = True,
    numeric_ns: str = "",
    endpoint: str = DEFAULT_PROFILE_ENDPOINT,
    timeout: float = 120.0,
) -> dict[str, Any]:
    """Run g:GOSt enrichment with an explicit custom statistical background."""
    if not query:
        raise ValueError("g:Profiler query must contain at least one gene")
    if not background:
        raise ValueError("g:Profiler background must contain at least one gene")

    query_unique = list(dict.fromkeys(str(x).strip() for x in query if str(x).strip()))
    background_unique = list(
        dict.fromkeys(str(x).strip() for x in background if str(x).strip())
    )
    missing = sorted(set(query_unique) - set(background_unique))
    if missing:
        raise ValueError(
            "All query genes must belong to the custom background; missing: "
            + ", ".join(missing[:20])
        )

    payload: dict[str, Any] = {
        "organism": organism,
        "query": query_unique,
        "background": background_unique,
        "sources": sources or ["GO:BP", "REAC", "KEGG", "CORUM"],
        "significance_threshold_method": significance_threshold_method,
        "user_threshold": float(user_threshold),
        "all_results": bool(all_results),
        "domain_scope": domain_scope,
        "no_evidences": bool(no_evidences),
        "ordered": False,
    }
    if numeric_ns:
        payload["numeric_ns"] = numeric_ns

    with httpx.Client(
        timeout=timeout,
        follow_redirects=True,
        headers={"User-Agent": USER_AGENT},
    ) as client:
        response = client.post(endpoint, json=payload)
        response.raise_for_status()
        data = response.json()

    if not isinstance(data, dict) or "result" not in data:
        raise RuntimeError("g:Profiler profile endpoint returned an unexpected response")
    return data
