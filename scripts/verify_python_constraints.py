from __future__ import annotations

import re
import tomllib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CONSTRAINTS = ROOT / "constraints" / "python-3.13.txt"
PROJECT_FILES = [ROOT / "pyproject.toml", ROOT / "apps" / "api" / "pyproject.toml"]


def normalize(name: str) -> str:
    return re.sub(r"[-_.]+", "-", name).lower()


def requirement_name(requirement: str) -> str:
    head = requirement.split(";", 1)[0].strip()
    match = re.match(r"^([A-Za-z0-9_.-]+)", head)
    if not match:
        raise ValueError(f"Cannot parse requirement: {requirement}")
    return normalize(match.group(1))


def load_constraints() -> dict[str, str]:
    if not CONSTRAINTS.exists():
        raise SystemExit(f"Missing Python constraints file: {CONSTRAINTS.relative_to(ROOT)}")

    pinned: dict[str, str] = {}
    for raw in CONSTRAINTS.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        requirement = line.split(";", 1)[0].strip()
        if "==" not in requirement:
            raise SystemExit(f"Constraint is not exactly pinned: {line}")
        name, version = requirement.split("==", 1)
        normalized = normalize(name.strip())
        if normalized in pinned:
            raise SystemExit(f"Duplicate constraint for {normalized}")
        if not version.strip():
            raise SystemExit(f"Empty pinned version for {normalized}")
        pinned[normalized] = version.strip()
    return pinned


def declared_dependencies(path: Path) -> set[str]:
    payload = tomllib.loads(path.read_text(encoding="utf-8"))
    project = payload.get("project") or {}
    requirements = list(project.get("dependencies") or [])
    for values in (project.get("optional-dependencies") or {}).values():
        requirements.extend(values or [])
    return {requirement_name(str(item)) for item in requirements}


def main() -> None:
    pinned = load_constraints()
    declared: set[str] = set()
    for path in PROJECT_FILES:
        declared.update(declared_dependencies(path))

    missing = sorted(name for name in declared if name not in pinned)
    if missing:
        joined = ", ".join(missing)
        raise SystemExit(
            "Python dependency constraints are stale. Missing direct project dependencies: " + joined
        )

    print("Python dependency constraints: PASS")
    print(f"Pinned packages: {len(pinned)}")
    print(f"Declared dependencies covered: {len(declared)}")


if __name__ == "__main__":
    main()
