from __future__ import annotations

import os
from pathlib import Path


def resolve_mcl_root() -> Path:
    configured = os.environ.get("MCL_ROOT")
    if configured:
        return Path(configured).expanduser().resolve()
    # <repo>/apps/api/mcl_api/settings.py
    return Path(__file__).resolve().parents[3]


MCL_ROOT = resolve_mcl_root()
