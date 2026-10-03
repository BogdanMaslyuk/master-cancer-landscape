from __future__ import annotations

from pathlib import Path
import yaml

from mcl.models.opentargets import (
    DiseaseMappingScope,
    MappingConfidence,
    OpenTargetsDiseaseMapping,
)


def load_opentargets_disease_mappings(config_path: str | Path) -> dict[str, OpenTargetsDiseaseMapping]:
    data = yaml.safe_load(Path(config_path).read_text(encoding="utf-8")) or {}
    contexts = data.get("contexts", {})
    mappings: dict[str, OpenTargetsDiseaseMapping] = {}
    for cancer_id, row in contexts.items():
        ot = row.get("opentargets")
        if not ot:
            continue
        mappings[cancer_id] = OpenTargetsDiseaseMapping(
            cancer_id=cancer_id,
            cancer_name=row["name"],
            molecular_context=ot["molecular_context"],
            disease_id=ot["disease_id"],
            expected_disease_name=ot["expected_disease_name"],
            mapping_scope=DiseaseMappingScope(ot["mapping_scope"]),
            mapping_confidence=MappingConfidence(ot["mapping_confidence"]),
            molecular_context_encoded_in_ot=bool(ot.get("molecular_context_encoded_in_ot", False)),
            rationale=ot["rationale"],
            mapping_source_url=ot["mapping_source_url"],
        )
    return mappings
