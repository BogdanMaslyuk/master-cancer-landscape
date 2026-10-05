from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import zipfile
from datetime import datetime, timezone
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


def _discover_download_url(client: httpx.Client, release: str) -> str:
    response = client.get(DOWNLOADS_PAGE)
    response.raise_for_status()
    parser = _LinkParser()
    parser.feed(response.text)
    expected = f"preclinical-db-{release}.zip".lower()
    candidates = [urljoin(DOWNLOADS_PAGE, href) for href in parser.hrefs if expected in href.lower()]
    if len(candidates) == 1:
        return candidates[0]
    if not candidates:
        # Defensive fallback for pages whose link is rendered outside a conventional <a> tag.
        match = re.search(r'href=["\']([^"\']*preclinical-db-' + re.escape(release) + r'\.zip[^"\']*)', response.text, re.I)
        if match:
            return urljoin(DOWNLOADS_PAGE, match.group(1))
        raise SystemExit(
            f"Could not discover preclinical-db-{release}.zip from {DOWNLOADS_PAGE}. "
            "Use --url with the exact official download URL."
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


def main() -> None:
    parser = argparse.ArgumentParser(description="Download and extract an official Preclinical Database bulk release.")
    parser.add_argument("--release", default="v1.0")
    parser.add_argument("--url", default="", help="Optional exact official ZIP URL; otherwise discovered from downloads page.")
    parser.add_argument("--force", action="store_true", help="Redownload even if a local ZIP already exists.")
    args = parser.parse_args()

    release_dir = BASE / args.release
    archive = release_dir / f"preclinical-db-{args.release}.zip"
    extracted = release_dir / "extracted"
    release_dir.mkdir(parents=True, exist_ok=True)

    with httpx.Client(timeout=120, follow_redirects=True, headers={"User-Agent": "MCL-research-pipeline/1.0"}) as client:
        url = args.url.strip() or _discover_download_url(client, args.release)
        if args.force or not archive.exists():
            tmp = archive.with_suffix(".zip.part")
            with client.stream("GET", url) as response:
                response.raise_for_status()
                with tmp.open("wb") as handle:
                    for chunk in response.iter_bytes(1024 * 1024):
                        handle.write(chunk)
            tmp.replace(archive)

    digest = _sha256(archive)
    if extracted.exists() and args.force:
        shutil.rmtree(extracted)
    if not extracted.exists() or not any(extracted.rglob("pcdb_complete.tsv")):
        _safe_extract(archive, extracted)

    complete = list(extracted.rglob("pcdb_complete.tsv"))
    if len(complete) != 1:
        raise SystemExit(f"Expected exactly one pcdb_complete.tsv after extraction; found {len(complete)}")

    manifest = {
        "contract": "mcl-preclinical-db-source-v1",
        "fetched_at": _now(),
        "release": args.release,
        "source_page": DOWNLOADS_PAGE,
        "download_url": url,
        "archive": str(archive.relative_to(ROOT)),
        "archive_sha256": digest,
        "pcdb_complete": str(complete[0].relative_to(ROOT)),
        "note": "Archive SHA256 is recorded locally. Compare it with the official checksum displayed on the PCDB downloads page when freezing a release for publication.",
    }
    (release_dir / "source_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")

    print("MCL Preclinical Database source")
    print(f"Release: {args.release}")
    print(f"URL: {url}")
    print(f"Archive SHA256: {digest}")
    print(f"pcdb_complete.tsv: {complete[0].relative_to(ROOT)}")


if __name__ == "__main__":
    main()
