from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import httpx


FIGSHARE_ARTICLE_ID = 25917643
API_URL = f"https://api.figshare.com/v2/articles/{FIGSHARE_ARTICLE_ID}"
ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "data" / "raw" / "pharmacology" / "prism_24q2"
VALID_SUFFIXES = (".csv", ".txt", ".tsv")


def _key(row: dict[str, Any]) -> str:
    return str(row.get("name") or "").lower().replace("-", "_")


def _select_files(all_files: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], str]:
    files = [row for row in all_files if _key(row).endswith(VALID_SUFFIXES)]
    lfc = [row for row in files if "lfc" in _key(row)]
    treatment = [row for row in files if "treatment_info" in _key(row)]
    cell_info = [row for row in files if "cell_line_info" in _key(row)]

    # Prefer the long-format release because it is the smallest sufficient input
    # and already carries DepMap/ACH model identifiers. Cell-line metadata is useful
    # for QC but not required to build the response table.
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


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Download the minimal PRISM Repurposing Public 24Q2 files needed by MCL."
    )
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--list-only", action="store_true")
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)

    with httpx.Client(follow_redirects=True, timeout=120.0) as client:
        article = client.get(API_URL)
        article.raise_for_status()
        files, mode = _select_files(article.json().get("files", []))
        print(f"PRISM Repurposing Public 24Q2 source mode: {mode}")
        print("Files selected:")
        total_size = 0
        for row in files:
            size = int(row.get("size") or 0)
            total_size += size
            print(f"- {row.get('name')} ({size / 1024 / 1024:.1f} MB)")
        print(f"Total selected: {total_size / 1024 / 1024:.1f} MB")
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
            with client.stream("GET", url, timeout=None) as response:
                response.raise_for_status()
                total = int(response.headers.get("content-length") or expected_size or 0)
                written = 0
                with target.open("wb") as handle:
                    for chunk in response.iter_bytes(chunk_size=1024 * 1024):
                        handle.write(chunk)
                        written += len(chunk)
                        if total:
                            print(
                                f"  {written / 1024 / 1024:.0f}/{total / 1024 / 1024:.0f} MB",
                                end="\r",
                            )
            print(f"  wrote {target.relative_to(ROOT)}{' ' * 24}")


if __name__ == "__main__":
    main()
