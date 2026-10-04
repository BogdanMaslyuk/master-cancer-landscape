from __future__ import annotations

import sys
from pathlib import Path
from time import perf_counter

from fastapi.testclient import TestClient


ROOT = Path(__file__).resolve().parents[1]
API_DIR = ROOT / "apps" / "api"
if str(API_DIR) not in sys.path:
    sys.path.insert(0, str(API_DIR))

from mcl_api.main_fast import app  # noqa: E402


CHECKS = (
    ("health", "/health"),
    ("genes", "/api/genes/search?page_size=3"),
    ("gene facets", "/api/genes/facets"),
    ("gene matrix", "/api/gene-matrix?limit=3"),
    ("pathways", "/api/pathways?significant_only=true&limit=5"),
)


def main() -> None:
    client = TestClient(app)
    slow: list[tuple[str, float]] = []
    for label, path in CHECKS:
        start = perf_counter()
        response = client.get(path)
        elapsed = perf_counter() - start
        if response.status_code != 200:
            raise SystemExit(f"ERROR: {label} returned HTTP {response.status_code}: {response.text[:500]}")
        print(f"PASS  {label:<12} {elapsed:7.3f}s  {path}")
        if elapsed > 5.0:
            slow.append((label, elapsed))

    if slow:
        print("WARNING: slow first-load endpoints detected:")
        for label, elapsed in slow:
            print(f"  {label}: {elapsed:.3f}s")
    print("Explorer API smoke verification: PASS")


if __name__ == "__main__":
    main()
