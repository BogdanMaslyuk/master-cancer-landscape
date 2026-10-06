from __future__ import annotations

import argparse
import shutil
import subprocess
import time
from pathlib import Path
from typing import Any

import httpx


FIGSHARE_ARTICLE_ID = 25917643
API_URL = f"https://api.figshare.com/v2/articles/{FIGSHARE_ARTICLE_ID}"
ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "data" / "raw" / "pharmacology" / "prism_24q2"
VALID_SUFFIXES = (".csv", ".txt", ".tsv")
DOWNLOAD_TIMEOUT = httpx.Timeout(connect=30.0, read=20.0, write=60.0, pool=30.0)
MAX_RETRIES = 8
CHUNK_SIZE = 1024 * 1024


def _key(row: dict[str, Any]) -> str:
    return str(row.get("name") or "").lower().replace("-", "_")


def _select_files(all_files: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], str]:
    files = [row for row in all_files if _key(row).endswith(VALID_SUFFIXES)]
    lfc = [row for row in files if "lfc" in _key(row)]
    treatment = [row for row in files if "treatment_info" in _key(row)]
    cell_info = [row for row in files if "cell_line_info" in _key(row)]

    if lfc and treatment:
        selected = [*lfc, *treatment]
        if cell_info:
            selected.extend(cell_info[:1])
        return list({str(row.get("id") or row.get("name")): row for row in selected}.values()), "long LFC + treatment metadata"

    matrices = [row for row in files if "extended_primary_data_matrix" in _key(row)]
    compounds = [row for row in files if "extended_primary_compound_list" in _key(row)]
    if matrices and compounds:
        return [matrices[0], compounds[0]], "extended primary matrix fallback"

    raise SystemExit(
        "No supported PRISM 24Q2 source pair was found. Expected either LFC + Treatment_Info "
        "or Extended_Primary_Data_Matrix + Extended_Primary_Compound_List."
    )


def _mb(value: int) -> float:
    return value / 1024 / 1024


def _prepare_partial(target: Path, expected_size: int) -> Path:
    part = target.with_name(target.name + ".part")

    if target.exists() and expected_size and target.stat().st_size == expected_size:
        return target

    # Preserve progress created by older versions that wrote directly to the final file.
    if target.exists():
        target_size = target.stat().st_size
        part_size = part.stat().st_size if part.exists() else -1
        if target_size > part_size:
            if part.exists():
                part.unlink()
            target.replace(part)
        else:
            target.unlink()

    return part


def _validate_and_finalize(working: Path, target: Path, expected_size: int) -> None:
    actual_size = working.stat().st_size if working.exists() else 0
    if expected_size and actual_size != expected_size:
        raise IOError(
            f"Downloaded size mismatch: got {actual_size} bytes, expected {expected_size} bytes"
        )
    working.replace(target)
    print(f"  wrote {target.relative_to(ROOT)} ({_mb(actual_size):.1f} MB)")


def _curl_executable() -> str | None:
    # Modern Windows ships curl.exe. Prefer it for large Figshare downloads because
    # its Range resume/retry handling is more robust than a long-lived Python stream.
    return shutil.which("curl.exe") or shutil.which("curl")


def _download_with_curl(url: str, target: Path, expected_size: int) -> bool:
    curl = _curl_executable()
    if not curl:
        return False

    working = _prepare_partial(target, expected_size)
    if working == target:
        print(f"Already present: {target.name}")
        return True

    offset = working.stat().st_size if working.exists() else 0
    if offset:
        print(f"  curl: continuing from {_mb(offset):.1f} MB")
    else:
        print("  curl: robust resumable download")

    command = [
        curl,
        "-L",
        "--fail",
        "--show-error",
        "--progress-bar",
        "--retry", "20",
        "--retry-delay", "2",
        "--retry-all-errors",
        "--connect-timeout", "20",
        # Abort a genuinely stalled transfer instead of hanging forever. curl will
        # then retry and resume from the bytes already written.
        "--speed-limit", "1024",
        "--speed-time", "20",
        "-C", "-",
        "-o", str(working),
        url,
    ]

    result = subprocess.run(command, check=False)
    if result.returncode != 0:
        partial_size = working.stat().st_size if working.exists() else 0
        raise SystemExit(
            f"curl download failed with exit code {result.returncode}. "
            f"Partial file preserved at {working} ({_mb(partial_size):.1f} MB). "
            "Re-run the same command to continue."
        )

    try:
        _validate_and_finalize(working, target, expected_size)
    except OSError as exc:
        raise SystemExit(
            f"Download finished but validation failed: {exc}. Partial file preserved at {working}."
        ) from exc
    return True


def _download_with_httpx(client: httpx.Client, url: str, target: Path, expected_size: int) -> None:
    working = _prepare_partial(target, expected_size)
    if working == target:
        print(f"Already present: {target.name}")
        return

    for attempt in range(1, MAX_RETRIES + 1):
        offset = working.stat().st_size if working.exists() else 0
        headers = {"Range": f"bytes={offset}-"} if offset > 0 else {}
        if offset:
            print(f"  resume from {_mb(offset):.1f} MB (attempt {attempt}/{MAX_RETRIES})")
        elif attempt > 1:
            print(f"  retry from start (attempt {attempt}/{MAX_RETRIES})")

        try:
            with client.stream("GET", url, headers=headers, timeout=DOWNLOAD_TIMEOUT) as response:
                if offset > 0 and response.status_code == 200:
                    offset = 0
                    if working.exists():
                        working.unlink()
                response.raise_for_status()

                mode = "ab" if offset > 0 and response.status_code == 206 else "wb"
                written = offset if mode == "ab" else 0
                total = expected_size or 0
                with working.open(mode) as handle:
                    for chunk in response.iter_bytes(chunk_size=CHUNK_SIZE):
                        if not chunk:
                            continue
                        handle.write(chunk)
                        written += len(chunk)
                        if total:
                            pct = min(100.0, written * 100.0 / total)
                            print(
                                f"  {_mb(written):.1f}/{_mb(total):.1f} MB  {pct:5.1f}%",
                                end="\r",
                                flush=True,
                            )
                        else:
                            print(f"  {_mb(written):.1f} MB", end="\r", flush=True)

            _validate_and_finalize(working, target, expected_size)
            return

        except (httpx.HTTPError, OSError) as exc:
            partial_size = working.stat().st_size if working.exists() else 0
            print(
                f"\n  connection interrupted after {_mb(partial_size):.1f} MB: "
                f"{type(exc).__name__}: {exc}"
            )
            if attempt >= MAX_RETRIES:
                raise SystemExit(
                    f"Download failed after {MAX_RETRIES} attempts. Partial file preserved at {working}"
                ) from exc
            delay = min(30, 2 ** (attempt - 1))
            print(f"  retrying in {delay} s ...")
            time.sleep(delay)


def _download(client: httpx.Client, url: str, target: Path, expected_size: int) -> None:
    if _download_with_curl(url, target, expected_size):
        return
    print("  curl.exe not found; falling back to Python resumable downloader")
    _download_with_httpx(client, url, target, expected_size)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Download the minimal PRISM Repurposing Public 24Q2 files needed by MCL."
    )
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--list-only", action="store_true")
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)

    with httpx.Client(follow_redirects=True, timeout=DOWNLOAD_TIMEOUT) as client:
        article = client.get(API_URL)
        article.raise_for_status()
        files, mode = _select_files(article.json().get("files", []))
        print(f"PRISM Repurposing Public 24Q2 source mode: {mode}")
        print("Files selected:")
        total_size = 0
        for row in files:
            size = int(row.get("size") or 0)
            total_size += size
            print(f"- {row.get('name')} ({_mb(size):.1f} MB)")
        print(f"Total selected: {_mb(total_size):.1f} MB")
        if args.list_only:
            return

        for row in files:
            name = str(row["name"])
            target = out / name
            expected_size = int(row.get("size") or 0)
            if target.exists() and expected_size and target.stat().st_size == expected_size:
                print(f"Already present: {target.name}")
                continue
            url = row.get("download_url")
            if not url:
                print(f"Skip {name}: Figshare did not provide download_url")
                continue
            print(f"Downloading {name} ...")
            _download(client, str(url), target, expected_size)


if __name__ == "__main__":
    main()
