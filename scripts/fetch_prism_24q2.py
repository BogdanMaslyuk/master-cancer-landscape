from __future__ import annotations

import argparse
from pathlib import Path

import httpx


FIGSHARE_ARTICLE_ID = 25917643
API_URL = f"https://api.figshare.com/v2/articles/{FIGSHARE_ARTICLE_ID}"
ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "data" / "raw" / "pharmacology" / "prism_24q2"


def _wanted(name: str) -> bool:
    key = name.lower()
    return any(
        token in key
        for token in (
            "lfc",
            "treatment_info",
            "treatment-info",
            "cell_line_info",
            "cell-line-info",
            "extended_primary_data_matrix",
            "extended_primary_compound_list",
        )
    ) and key.endswith((".csv", ".txt", ".tsv"))


def main() -> None:
    parser = argparse.ArgumentParser(description="Download the minimal PRISM Repurposing Public 24Q2 files needed by MCL.")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--list-only", action="store_true")
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)

    with httpx.Client(follow_redirects=True, timeout=120.0) as client:
        article = client.get(API_URL)
        article.raise_for_status()
        files = [row for row in article.json().get("files", []) if _wanted(str(row.get("name") or ""))]
        if not files:
            raise SystemExit("No expected PRISM 24Q2 files were found in the Figshare article metadata.")
        print("PRISM Repurposing Public 24Q2 files selected:")
        for row in files:
            print(f"- {row.get('name')} ({row.get('size', 0) / 1024 / 1024:.1f} MB)")
        if args.list_only:
            return

        for row in files:
            name = str(row["name"])
            target = out / name
            if target.exists() and target.stat().st_size == int(row.get("size") or -1):
                print(f"Already present: {target.name}")
                continue
            url = row.get("download_url")
            if not url:
                print(f"Skip {name}: Figshare did not provide download_url")
                continue
            print(f"Downloading {name} ...")
            with client.stream("GET", url, timeout=None) as response:
                response.raise_for_status()
                total = int(response.headers.get("content-length") or row.get("size") or 0)
                written = 0
                with target.open("wb") as handle:
                    for chunk in response.iter_bytes(chunk_size=1024 * 1024):
                        handle.write(chunk)
                        written += len(chunk)
                        if total:
                            print(f"  {written / 1024 / 1024:.0f}/{total / 1024 / 1024:.0f} MB", end="\r")
            print(f"  wrote {target.relative_to(ROOT)}{' ' * 24}")


if __name__ == "__main__":
    main()
