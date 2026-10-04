from __future__ import annotations

import ast
from pathlib import Path


API_PACKAGE = Path(__file__).resolve().parents[1] / "mcl_api"
ROUTER_DIR = API_PACKAGE / "routers"
FORBIDDEN_ROUTER_MODULES = {"numpy", "pandas", "pyarrow", "scipy", "statsmodels"}


def _state_imports(path: Path) -> set[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    imported: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.module in {"state", "..state"}:
            imported.update(alias.name for alias in node.names)
    return imported


def _top_level_imports(path: Path) -> set[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    imported: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.update(alias.name.split(".", 1)[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.add(node.module.lstrip(".").split(".", 1)[0])
    return imported


def test_http_layer_does_not_import_data_stores_directly():
    files = [API_PACKAGE / "main.py", *sorted(ROUTER_DIR.glob("*.py"))]
    offenders: dict[str, list[str]] = {}

    for path in files:
        names = _state_imports(path)
        forbidden = sorted(name for name in names if name == "store" or name.endswith("_store"))
        if forbidden:
            offenders[str(path.relative_to(API_PACKAGE))] = forbidden

    assert offenders == {}, f"HTTP layer must use services, not stores: {offenders}"


def test_http_layer_does_not_own_lru_caches():
    files = [API_PACKAGE / "main.py", *sorted(ROUTER_DIR.glob("*.py"))]
    offenders = []
    for path in files:
        text = path.read_text(encoding="utf-8")
        if "lru_cache" in text:
            offenders.append(str(path.relative_to(API_PACKAGE)))

    assert offenders == [], f"Request-level caching belongs in services: {offenders}"


def test_routers_do_not_import_scientific_dataframes_or_statistics_libraries():
    offenders: dict[str, list[str]] = {}
    for path in sorted(ROUTER_DIR.glob("*.py")):
        forbidden = sorted(_top_level_imports(path) & FORBIDDEN_ROUTER_MODULES)
        if forbidden:
            offenders[str(path.relative_to(API_PACKAGE))] = forbidden

    assert offenders == {}, f"Routers must remain HTTP-only: {offenders}"


def test_api_package_does_not_reference_raw_or_interim_data_roots():
    offenders: list[str] = []
    for path in sorted(API_PACKAGE.rglob("*.py")):
        text = path.read_text(encoding="utf-8").replace("\\", "/").lower()
        if "data/raw" in text or "data/interim" in text:
            offenders.append(str(path.relative_to(API_PACKAGE)))

    assert offenders == [], (
        "Interactive API must not read raw/interim source layers directly: "
        f"{offenders}"
    )
