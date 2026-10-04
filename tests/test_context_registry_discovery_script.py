from __future__ import annotations

import importlib.util
from pathlib import Path


def _load_discovery_module():
    root = Path(__file__).resolve().parents[1]
    path = root / "scripts" / "discover_cancer_context_registry.py"
    spec = importlib.util.spec_from_file_location("discover_cancer_context_registry", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_load_crispr_model_ids_accepts_header_only_probe(tmp_path: Path):
    path = tmp_path / "CRISPRGeneEffect.csv"
    path.write_text(
        "Unnamed: 0,KRAS (3845),TP53 (7157)\n"
        "ACH-000001,-1.2,-0.4\n"
        "ACH-000002,-0.7,-0.9\n",
        encoding="utf-8",
    )

    module = _load_discovery_module()
    assert module._load_crispr_model_ids(path) == {"ACH-000001", "ACH-000002"}
