from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import time
import zipfile
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import urljoin

import httpx

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "external" / "preclinical_db"
DOWNLOADS_PAGE = "https://www.preclinicaldata.org/downloads"


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


class _LinkParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.hrefs: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag.lower() != "a":
            return
        for key, value in attrs:
            if key.lower() == "href" and value:
                self.hrefs.append(value)


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _retry_after_seconds(response: httpx.Response, fallback: float) -> float:
    value = response.headers.get("Retry-After", "").strip()
    if not value:
        return fallback
    try:
        return max(float(value), 0.0)
    except ValueError:
        try:
            target = parsedate_to_datetime(value)
            if target.tzinfo is None:
                target = target.replace(tzinfo=timezone.utc)
            return max((target - datetime.now(timezone.utc)).total_seconds(), 0.0)
        except (TypeError, ValueError, OverflowError):
            return fallback


def _get_with_backoff(
    client: httpx.Client,
    url: str,
    *,
    attempts: int = 5,
    stream: bool = False,
) -> httpx.Response:
    """GET with conservative retry for transient 429/5xx responses."""
    last_response: httpx.Response | None = None
    for attempt in range(1, attempts + 1):
        if stream:
            request = client.build_request("GET", url)
            response = client.send(request, stream=True)
        else:
            response = client.get(url)
        last_response = response
        if response.status_code not in {429, 500, 502, 503, 504}:
            return response
        if attempt == attempts:
            return response
        wait = _retry_after_seconds(response, fallback=min(3 * (2 ** (attempt - 1)), 30))
        print(
            f"Transient HTTP {response.status_code} from {url}; "
            f"retry {attempt}/{attempts - 1} after {wait:.0f}s."
        )
        response.close()
        time.sleep(wait)
    assert last_response is not None
    return last_response


def _discover_download_url(client: httpx.Client, release: str) -> str:
    response = _get_with_backoff(client, DOWNLOADS_PAGE)
    if response.status_code == 429:
        raise SystemExit(
            "Preclinical Database downloads page is rate-limiting automated requests (HTTP 429).\n"
            "The pipeline will use a local archive automatically if available. Download the official "
            f"preclinical-db-{release}.zip in a browser and either place it in your Downloads folder or run "
            f"with --archive <path>. Official page: {DOWNLOADS_PAGE}"
        )
    response.raise_for_status()
    parser = _LinkParser()
    parser.feed(response.text)
    expected = f"preclinical-db-{release}.zip".lower()
    candidates = [urljoin(DOWNLOADS_PAGE, href) for href in parser.hrefs if expected in href.lower()]
    if len(candidates) == 1:
        return candidates[0]
    if not candidates:
        match = re.search(
            r'href=["\']([^"\']*preclinical-db-' + re.escape(release) + r'\.zip[^"\']*)',
            response.text,
            re.I,
        )
        if match:
            return urljoin(DOWNLOADS_PAGE, match.group(1))
        raise SystemExit(
            f"Could not discover preclinical-db-{release}.zip from {DOWNLOADS_PAGE}. "
            "Use --url with the exact official download URL or --archive with a local ZIP."
        )
    raise SystemExit(f"Ambiguous download links for release {release}: {candidates}")


def _safe_extract(archive: Path, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    root = destination.resolve()
    with zipfile.ZipFile(archive) as zf:
        for member in zf.infolist():
            target = (destination / member.filename).resolve()
            if root != target and root not in target.parents:
                raise SystemExit(f"Unsafe ZIP member path: {member.filename}")
        zf.extractall(destination)


def _find_local_archive(release: str, explicit: str) -> Path | None:
    if explicit.strip():
        path = Path(explicit).expanduser().resolve()
        if not path.exists():
            raise SystemExit(f"Local PCDB archive not found: {path}")
        return path

    name = f"preclinical-db-{release}.zip"
    candidates = [
        BASE / release / name,
        Path.home() / "Downloads" / name,
        Path.home() / "Desktop" / name,
    ]
    for path in candidates:
        if path.exists():
            return path.resolve()
    return None


def _copy_local_archive(source: Path, destination: Path) -> None:
    if source.resolve() == destination.resolve():
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    tmp = destination.with_suffix(".zip.part")
    shutil.copy2(source, tmp)
    tmp.replace(destination)


def _download_archive(client: httpx.Client, url: str, destination: Path) -> None:
    tmp = destination.with_suffix(".zip.part")
    response = _get_with_backoff(client, url, stream=True)
    try:
        if response.status_code == 429:
            raise SystemExit(
                "Preclinical Database archive download is rate-limited (HTTP 429). "
                "Download the official ZIP in a browser and rerun; the script checks your Downloads folder automatically."
            )
        response.raise_for_status()
        with tmp.open("wb") as handle:
            for chunk in response.iter_bytes(1024 * 1024):
                handle.write(chunk)
    finally:
        response.close()
    tmp.replace(destination)


def main() -> None:
    parser = argparse.ArgumentParser(description="Download and extract an official Preclinical Database bulk release.")
    parser.add_argument("--release", default="v1.0")
    parser.add_argument("--url", default="", help="Optional exact official ZIP URL; otherwise discovered from downloads page.")
    parser.add_argument(
        "--archive",
        default="",
        help="Optional local official ZIP. If omitted, the script also checks ~/Downloads and ~/Desktop.",
    )
    parser.add_argument("--force", action="store_true", help="Redownload or recopy even if a local repository ZIP already exists.")
    args = parser.parse_args()

    release_dir = BASE / args.release
    archive = release_dir / f"preclinical-db-{args.release}.zip"
    extracted = release_dir / "extracted"
    release_dir.mkdir(parents=True, exist_ok=True)

    local_source = _find_local_archive(args.release, args.archive)
    source_kind = "existing_repository_archive" if archive.exists() else ""
    source_url = ""
    source_local_path = ""

    # Important: if the archive already exists locally, do not touch the website just to rediscover its URL.
    if args.force or not archive.exists():
        if local_source is not None:
            _copy_local_archive(local_source, archive)
            source_kind = "local_archive"
            source_local_path = str(local_source)
        else:
            headers = {
                "User-Agent": (
                    "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
                    "AppleWebKit/537.36 (KHTML, like Gecko) Chrome/150.0 Safari/537.36"
                ),
                "Accept": "text/html,application/xhtml+xml,application/zip;q=0.9,*/*;q=0.8",
                "Accept-Language": "en-US,en;q=0.9",
            }
            with httpx.Client(timeout=120, follow_redirects=True, headers=headers) as client:
                url = args.url.strip() or _discover_download_url(client, args.release)
                _download_archive(client, url, archive)
                source_kind = "official_http_download"
                source_url = url

    if not archive.exists():
        raise SystemExit(f"PCDB archive was not created: {archive}")
    if not zipfile.is_zipfile(archive):
        raise SystemExit(f"Downloaded/copied PCDB file is not a valid ZIP archive: {archive}")

    digest = _sha256(archive)
    if extracted.exists() and args.force:
        shutil.rmtree(extracted)
    if not extracted.exists() or not any(extracted.rglob("pcdb_complete.tsv")):
        _safe_extract(archive, extracted)

    complete = list(extracted.rglob("pcdb_complete.tsv"))
    if len(complete) != 1:
        raise SystemExit(f"Expected exactly one pcdb_complete.tsv after extraction; found {len(complete)}")

    manifest = {
        "contract": "mcl-preclinical-db-source-v1.1",
        "fetched_at": _now(),
        "release": args.release,
        "source_page": DOWNLOADS_PAGE,
        "source_kind": source_kind,
        "download_url": source_url,
        "local_source_path": source_local_path,
        "archive": str(archive.relative_to(ROOT)),
        "archive_sha256": digest,
        "pcdb_complete": str(complete[0].relative_to(ROOT)),
        "note": (
            "Archive SHA256 is recorded locally. Compare it with the official checksum displayed on the PCDB downloads page "
            "when freezing a release for publication. HTTP 429 is handled with backoff; a manually downloaded official ZIP "
            "can be supplied without changing downstream provenance."
        ),
    }
    (release_dir / "source_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    print("MCL Preclinical Database source v1.1")
    print(f"Release: {args.release}")
    print(f"Source kind: {source_kind}")
    if source_url:
        print(f"URL: {source_url}")
    if source_local_path:
        print(f"Local source: {source_local_path}")
    print(f"Archive SHA256: {digest}")
    print(f"pcdb_complete.tsv: {complete[0].relative_to(ROOT)}")


if __name__ == "__main__":
    main()
