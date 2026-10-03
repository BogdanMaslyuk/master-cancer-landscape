from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml


def project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def load_yaml(path: str | Path) -> dict[str, Any]:
    with Path(path).open("r", encoding="utf-8") as fh:
        data = yaml.safe_load(fh)
    return data or {}


def load_project_config(root: Path | None = None) -> dict[str, Any]:
    root = root or project_root()
    return load_yaml(root / "config" / "project.yaml")
