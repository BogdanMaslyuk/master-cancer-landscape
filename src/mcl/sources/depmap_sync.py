from __future__ import annotations

import csv
import io
import json
import os
import re
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable
from urllib.parse import urlparse

import httpx
import pandas as pd

from mcl.sources.depmap import (
    FILE_CANDIDATES,
    MODEL_ID_ALIASES,
    GENE_ALIASES,
)
from mcl.utils.hash import sha256_file

PUBLIC_CATALOG_URL = "https://depmap.org/portal/api/no-captcha/download/files"
PROTECTED_CATALOG_URL = "https://depmap.org/portal/api/download/files"

RELEASE_RE = re.compile(r"(?<!\d)(?P<yy>\d{2})Q(?P<q>[1-4])(?!\d)", re.IGNORECASE)


@dataclass
class SyncFileStatus:
    logical_name: str
    expected_names: str
    local_filename: str | None
    local_path: str
    status: str
    catalog_found: bool | None
    catalog_filename: str | None
    downloaded: bool
    size_bytes: int | None
    sha256: str | None
    header_check: str | None
    note: str = ""


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def parse_catalog_csv(text: str) -> list[dict[str, str]]:
    """Parse the DepMap download catalog without assuming a fixed column order.

    DepMap deliberately exposes a no-CAPTCHA metadata endpoint but may evolve
    the catalog schema. We therefore preserve all columns and use value-based
    matching for releases and filenames.
    """
    if not text.strip():
        return []
    reader = csv.DictReader(io.StringIO(text))
    if not reader.fieldnames:
        return []
    rows: list[dict[str, str]] = []
    for row in reader:
        rows.append({str(k): "" if v is None else str(v) for k, v in row.items() if k is not None})
    return rows


def fetch_catalog(
    url: str = PUBLIC_CATALOG_URL,
    bearer_token: str | None = None,
    client: httpx.Client | None = None,
) -> tuple[list[dict[str, str]], dict[str, str]]:
    own = client is None
    headers = {"User-Agent": "Master-Cancer-Landscape/0.3.6"}
    if bearer_token:
        headers["Authorization"] = f"Bearer {bearer_token}"
    c = client or httpx.Client(timeout=60.0, follow_redirects=True, headers=headers)
    try:
        response = c.get(url, headers=headers)
        response.raise_for_status()
        rows = parse_catalog_csv(response.text)
        return rows, {
            "url": url,
            "retrieved_at": _utc_now(),
            "content_type": response.headers.get("content-type", ""),
            "row_count": str(len(rows)),
        }
    finally:
        if own:
            c.close()


def _release_key(value: str) -> tuple[int, int] | None:
    m = RELEASE_RE.search(value.strip())
    if not m:
        return None
    return int(m.group("yy")), int(m.group("q"))


def catalog_releases(rows: Iterable[dict[str, str]]) -> list[str]:
    found: set[str] = set()
    for row in rows:
        for value in row.values():
            text = str(value).strip().upper()
            match = RELEASE_RE.search(text)
            if match:
                found.add(f"{match.group('yy')}Q{match.group('q')}".upper())
    return sorted(found, key=lambda x: _release_key(x) or (-1, -1))


def latest_catalog_release(rows: Iterable[dict[str, str]]) -> str | None:
    releases = catalog_releases(rows)
    return releases[-1] if releases else None


def _row_contains_release(row: dict[str, str], release: str) -> bool:
    target = release.strip().upper()
    return any(target in str(v).strip().upper() for v in row.values())


def _basename(value: str) -> str:
    text = str(value).strip().replace("\\", "/")
    return text.rsplit("/", 1)[-1]


def _find_catalog_row(
    rows: Iterable[dict[str, str]], release: str, candidates: Iterable[str]
) -> tuple[dict[str, str] | None, str | None]:
    names = {x.lower(): x for x in candidates}
    rows = list(rows)
    release_rows = [r for r in rows if _row_contains_release(r, release)]
    pool = release_rows if release_rows else rows
    for row in pool:
        for value in row.values():
            base = _basename(str(value))
            if base.lower() in names:
                return row, base
    return None, None


def _extract_url(row: dict[str, str] | None) -> str | None:
    if not row:
        return None
    preferred = ("url", "download_url", "downloadurl", "file_url", "fileurl")
    lowered = {k.lower(): v for k, v in row.items()}
    for key in preferred:
        value = str(lowered.get(key, "")).strip()
        if value.startswith(("https://", "http://")):
            return value
    for value in row.values():
        text = str(value).strip()
        if text.startswith(("https://", "http://")):
            return text
    return None


def _header_check(logical_name: str, path: Path) -> tuple[bool, str]:
    try:
        cols = pd.read_csv(path, nrows=0).columns.tolist()
    except Exception as exc:  # pragma: no cover - pandas gives many parser subclasses
        return False, f"cannot read CSV header: {exc}"
    if not cols:
        return False, "empty CSV header"
    if logical_name == "models":
        ok = any(c in cols for c in MODEL_ID_ALIASES)
        return ok, "ModelID column found" if ok else "no recognized ModelID column"
    if logical_name == "mutations":
        model_ok = any(c in cols for c in MODEL_ID_ALIASES)
        gene_ok = any(c in cols for c in GENE_ALIASES)
        ok = model_ok and gene_ok
        return ok, "ModelID and gene columns found" if ok else "missing ModelID or gene column"
    if logical_name == "omics_profiles":
        model_ok = any(c in cols for c in MODEL_ID_ALIASES)
        datatype_ok = any(c in cols for c in ("Datatype", "DataType", "datatype", "data_type"))
        ok = model_ok and datatype_ok
        return ok, "ModelID and Datatype columns found" if ok else "missing ModelID or Datatype column"
    if logical_name in {"gene_effect", "gene_dependency"}:
        ok = len(cols) >= 2
        return ok, f"matrix header with {len(cols)} columns" if ok else "matrix has fewer than 2 columns"
    return True, f"header read ({len(cols)} columns)"


def _download_file(url: str, destination: Path, token: str, client: httpx.Client | None = None) -> None:
    own = client is None
    headers = {
        "Authorization": f"Bearer {token}",
        "User-Agent": "Master-Cancer-Landscape/0.3.6",
    }
    c = client or httpx.Client(timeout=300.0, follow_redirects=True, headers=headers)
    part = destination.with_suffix(destination.suffix + ".part")
    try:
        with c.stream("GET", url, headers=headers) as response:
            response.raise_for_status()
            ctype = response.headers.get("content-type", "").lower()
            if "text/html" in ctype:
                raise RuntimeError("download endpoint returned HTML instead of a data file")
            destination.parent.mkdir(parents=True, exist_ok=True)
            with part.open("wb") as fh:
                for chunk in response.iter_bytes():
                    fh.write(chunk)
        part.replace(destination)
    finally:
        if part.exists():
            part.unlink(missing_ok=True)
        if own:
            c.close()


def sync_depmap(
    root: Path,
    release: str,
    depmap_dir: Path | None = None,
    check_catalog: bool = True,
    bearer_token: str | None = None,
    client: httpx.Client | None = None,
) -> dict:
    """Synchronize the pinned DepMap release conservatively.

    Current public DepMap metadata are discoverable programmatically, while raw
    file downloads are CAPTCHA-protected. If DepMap exposes a URL-bearing catalog
    to an authorized bearer token in the future, the same command can download
    missing files. It never fabricates or scrapes a CAPTCHA-protected URL.
    """
    root = root.resolve()
    base = (depmap_dir or (root / "data/raw/depmap" / release)).resolve()
    base.mkdir(parents=True, exist_ok=True)
    token = bearer_token or os.getenv("DEPMAP_BEARER_TOKEN") or None

    catalog: list[dict[str, str]] = []
    catalog_meta: dict[str, str] = {}
    catalog_error: str | None = None
    if check_catalog:
        try:
            catalog, catalog_meta = fetch_catalog(PUBLIC_CATALOG_URL, client=client)
        except Exception as exc:
            catalog_error = f"{type(exc).__name__}: {exc}"

    latest = latest_catalog_release(catalog) if catalog else None
    available_releases = catalog_releases(catalog) if catalog else []

    protected_catalog: list[dict[str, str]] | None = None
    protected_error: str | None = None
    statuses: list[SyncFileStatus] = []

    for logical_name, candidates in FILE_CANDIDATES.items():
        existing = next((base / name for name in candidates if (base / name).exists()), None)
        public_row, public_name = _find_catalog_row(catalog, release, candidates) if catalog else (None, None)
        catalog_found = (public_row is not None) if catalog else None
        downloaded = False
        note = ""

        if existing is None and token:
            if protected_catalog is None:
                try:
                    protected_catalog, _ = fetch_catalog(
                        PROTECTED_CATALOG_URL, bearer_token=token, client=client
                    )
                except Exception as exc:
                    protected_catalog = []
                    protected_error = f"{type(exc).__name__}: {exc}"
            prow, pname = _find_catalog_row(protected_catalog or [], release, candidates)
            download_url = _extract_url(prow)
            if download_url:
                target_name = pname or candidates[0]
                target = base / target_name
                _download_file(download_url, target, token=token, client=client)
                existing = target
                downloaded = True
                note = "downloaded through authorized DepMap URL"
            elif protected_error:
                note = "authorized download unavailable: " + protected_error
            else:
                note = "authorized catalog returned no downloadable URL; manual portal download still required"

        if existing is not None:
            ok, header_note = _header_check(logical_name, existing)
            status = "ready" if ok else "invalid_header"
            statuses.append(
                SyncFileStatus(
                    logical_name=logical_name,
                    expected_names=" or ".join(candidates),
                    local_filename=existing.name,
                    local_path=str(existing),
                    status=status,
                    catalog_found=catalog_found,
                    catalog_filename=public_name,
                    downloaded=downloaded,
                    size_bytes=existing.stat().st_size,
                    sha256=sha256_file(existing),
                    header_check=header_note,
                    note=note,
                )
            )
        else:
            if not note:
                if catalog_found:
                    note = "listed by DepMap; download via official portal because public file URLs are CAPTCHA-protected"
                elif catalog_error:
                    note = "catalog unavailable; local file is missing"
                else:
                    note = "file is missing; confirm the release in the official DepMap portal"
            statuses.append(
                SyncFileStatus(
                    logical_name=logical_name,
                    expected_names=" or ".join(candidates),
                    local_filename=None,
                    local_path=str(base / candidates[0]),
                    status="manual_required" if not token else "download_unavailable",
                    catalog_found=catalog_found,
                    catalog_filename=public_name,
                    downloaded=False,
                    size_bytes=None,
                    sha256=None,
                    header_check=None,
                    note=note,
                )
            )

    ready = all(s.status == "ready" for s in statuses)
    return {
        "depmap_release": release,
        "raw_dir": str(base),
        "ready": ready,
        "checked_at": _utc_now(),
        "public_catalog_url": PUBLIC_CATALOG_URL,
        "catalog_retrieved": bool(catalog),
        "catalog_error": catalog_error,
        "catalog_meta": catalog_meta,
        "catalog_releases": available_releases,
        "latest_catalog_release": latest,
        "newer_release_available": bool(latest and _release_key(latest) and _release_key(release) and _release_key(latest) > _release_key(release)),
        "bearer_token_configured": bool(token),
        "protected_catalog_error": protected_error,
        "files": [asdict(x) for x in statuses],
    }
