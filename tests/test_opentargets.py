import json
from pathlib import Path

from mcl.models.opentargets import (
    DiseaseMappingScope,
    MappingConfidence,
    OpenTargetsAssociation,
)
from mcl.normalize.disease_mapping import load_opentargets_disease_mappings
from mcl.qc.opentargets import opentargets_qc
from mcl.sources.opentargets import ASSOCIATED_TARGETS_QUERY, TRACTABILITY_QUERY

ROOT = Path(__file__).resolve().parents[1]


def make_row(**updates):
    base = dict(
        evidence_id="EVID-TEST",
        cancer_id="CANCER-001",
        target_id="TGT-KRAS",
        target_ensembl_id="ENSG00000133703",
        target_symbol="KRAS",
        disease_id="EFO_0000571",
        disease_name="lung adenocarcinoma",
        expected_disease_name="lung adenocarcinoma",
        molecular_context="KRAS p.G12C",
        mapping_scope=DiseaseMappingScope.EXACT_DISEASE,
        mapping_confidence=MappingConfidence.HIGH,
        molecular_context_encoded_in_ot=False,
        association_found=True,
        association_score=0.8,
        direct_association=True,
        genetic_association_score=0.4,
        somatic_mutation_score=0.7,
        known_drug_score=0.9,
        literature_score=0.5,
        datatype_scores_json=json.dumps({"known_drug": 0.9}),
        datasource_scores_json=json.dumps({"chembl": 0.9}),
        tractability_json=json.dumps([{"label": "Small Molecule", "modality": "SM", "value": True}]),
        source_release="26.09",
        api_endpoint="https://api.platform.opentargets.org/api/v4/graphql",
        retrieved_at="2026-09-26T00:00:00+00:00",
        raw_association_file="data/raw/opentargets/26.09/test.json",
        raw_tractability_file="data/raw/opentargets/26.09/tractability.json",
        raw_record_hash="abc",
        notes="test",
    )
    base.update(updates)
    return OpenTargetsAssociation(**base)


def test_wave1_disease_mappings_are_explicit_and_molecular_context_not_encoded():
    mappings = load_opentargets_disease_mappings(ROOT / "config/cancer_contexts.yaml")
    assert set(mappings) == {"CANCER-001", "CANCER-004", "CANCER-011", "CANCER-013"}
    assert mappings["CANCER-001"].disease_id == "MONDO_0005061"
    assert mappings["CANCER-013"].disease_id == "MONDO_0018874"
    assert mappings["CANCER-004"].disease_id == "MONDO_0005184"
    assert mappings["CANCER-011"].disease_id == "MONDO_0018177"
    assert mappings["CANCER-004"].mapping_scope == DiseaseMappingScope.EXACT_DISEASE
    assert mappings["CANCER-011"].mapping_scope == DiseaseMappingScope.EXACT_DISEASE
    assert all(not x.molecular_context_encoded_in_ot for x in mappings.values())


def test_graphql_query_requests_direct_associations_only():
    assert "enableIndirect: false" in ASSOCIATED_TARGETS_QUERY
    assert "datatypeScores" in ASSOCIATED_TARGETS_QUERY
    assert "datasourceScores" in ASSOCIATED_TARGETS_QUERY
    assert "tractability" in TRACTABILITY_QUERY


def test_opentargets_qc_passes_valid_row():
    qc = opentargets_qc([make_row()], expected_n=1)
    assert not [x for x in qc if x.severity == "ERROR"]
    assert not [x for x in qc if x.severity == "WARNING"]


def test_absent_association_is_information_not_zero_or_error():
    qc = opentargets_qc([make_row(association_found=False, association_score=None)], expected_n=1)
    assert not [x for x in qc if x.severity == "ERROR"]
    assert [x for x in qc if x.check == "association_absent" and x.severity == "INFO"]


def test_duplicate_cancer_target_is_error():
    qc = opentargets_qc([make_row(), make_row(evidence_id="EVID-TEST2")], expected_n=2)
    assert [x for x in qc if x.check == "unique_cancer_target" and x.severity == "ERROR"]


def test_score_outside_zero_one_is_error():
    qc = opentargets_qc([make_row(association_score=1.5)], expected_n=1)
    assert [x for x in qc if x.check == "association_score_range" and x.severity == "ERROR"]


def test_collect_opentargets_wave1_offline_roundtrip(tmp_path, monkeypatch):
    import pandas as pd
    import yaml
    from mcl.analysis import opentargets as analysis_ot

    (tmp_path / "data/input").mkdir(parents=True)
    (tmp_path / "data/processed").mkdir(parents=True)
    (tmp_path / "config").mkdir(parents=True)

    (tmp_path / "data/input/cancer_target_pairs_wave1.tsv").write_text(
        "Evidence_ID\tCancer_ID\tTarget_ID\tCancer_name\tMolecular_context\tHGNC_symbol\n"
        "EVID-1\tCANCER-001\tTGT-KRAS\tLung adenocarcinoma — KRAS G12C\tKRAS p.G12C\tKRAS\n",
        encoding="utf-8",
    )
    (tmp_path / "data/processed/target_identifiers.tsv").write_text(
        "target_id\thgnc_symbol\tensembl_gene_id\nTGT-KRAS\tKRAS\tENSG00000133703\n",
        encoding="utf-8",
    )
    cfg = {
        "contexts": {
            "CANCER-001": {
                "name": "Lung adenocarcinoma — KRAS G12C",
                "opentargets": {
                    "disease_id": "EFO_0000571",
                    "expected_disease_name": "lung adenocarcinoma",
                    "mapping_scope": "exact_disease",
                    "mapping_confidence": "high",
                    "molecular_context": "KRAS p.G12C",
                    "molecular_context_encoded_in_ot": False,
                    "rationale": "test",
                    "mapping_source_url": "https://platform.opentargets.org/disease/EFO_0000571",
                },
            }
        }
    }
    (tmp_path / "config/cancer_contexts.yaml").write_text(
        yaml.safe_dump(cfg, sort_keys=False, allow_unicode=True), encoding="utf-8"
    )

    class FakeClient:
        def __init__(self, endpoint):
            self.endpoint = endpoint

        def fetch_meta(self):
            return {
                "data": {
                    "meta": {
                        "apiVersion": {"x": "26", "y": "9", "z": "0"},
                        "dataVersion": {"year": "26", "month": "09"},
                    }
                }
            }

        def fetch_disease(self, disease_id):
            assert disease_id == "EFO_0000571"
            return {"data": {"disease": {"id": "EFO_0000571", "name": "lung adenocarcinoma"}}}

        def resolve_disease_names(self, names):
            raise AssertionError("Resolver should not be called when configured ID is current")

        def fetch_target_tractability(self, ensembl_id):
            assert ensembl_id == "ENSG00000133703"
            return {
                "data": {
                    "target": {
                        "id": "ENSG00000133703",
                        "approvedSymbol": "KRAS",
                        "tractability": [{"label": "Small Molecule", "modality": "SM", "value": True}],
                    }
                }
            }

        def fetch_associated_targets_page(self, disease_id, index, size):
            assert disease_id == "EFO_0000571"
            assert index == 0
            return {
                "data": {
                    "disease": {
                        "id": "EFO_0000571",
                        "name": "lung adenocarcinoma",
                        "associatedTargets": {
                            "count": 1,
                            "rows": [
                                {
                                    "target": {"id": "ENSG00000133703", "approvedSymbol": "KRAS"},
                                    "score": 0.91,
                                    "datatypeScores": [
                                        {"id": "somatic_mutation", "score": 0.88},
                                        {"id": "known_drug", "score": 0.93},
                                    ],
                                    "datasourceScores": [{"id": "chembl", "score": 0.93}],
                                }
                            ],
                        },
                    }
                }
            }

    monkeypatch.setattr(analysis_ot, "OpenTargetsClient", FakeClient)
    rows, provenance, meta = analysis_ot.collect_opentargets_wave1(
        root=tmp_path,
        source_release="26.09-test",
        endpoint="https://example.test/graphql",
        page_size=500,
    )
    assert len(rows) == 1
    row = rows[0]
    assert row.association_found is True
    assert row.association_score == 0.91
    assert row.somatic_mutation_score == 0.88
    assert row.known_drug_score == 0.93
    assert row.molecular_context_encoded_in_ot is False
    assert len(provenance) >= 10
    assert meta["diseases"]["CANCER-001"]["returned_id"] == "EFO_0000571"
    assert (tmp_path / "data/raw/opentargets/26.09-test/target_ENSG00000133703_tractability.json").exists()
    assert (tmp_path / "data/raw/opentargets/26.09-test/platform_meta.json").exists()
    assert list((tmp_path / "data/raw/opentargets/26.09-test").glob("*associated_targets_page_0000.json"))


def test_runtime_disease_resolver_recovers_identifier_drift(tmp_path):
    from types import SimpleNamespace
    from mcl.analysis.opentargets import resolve_disease_for_release

    class FakeClient:
        def fetch_disease(self, disease_id):
            if disease_id == "EFO_OLD":
                return {"data": {"disease": None}}
            assert disease_id == "MONDO_NEW"
            return {"data": {"disease": {"id": "MONDO_NEW", "name": "lung adenocarcinoma"}}}

        def resolve_disease_names(self, names):
            assert names == ["lung adenocarcinoma"]
            return {
                "data": {
                    "mapIds": {
                        "mappings": [
                            {
                                "term": "lung adenocarcinoma",
                                "hits": [
                                    {"id": "MONDO_WRONG", "name": "lung adenocarcinoma measurement"},
                                    {"id": "MONDO_NEW", "name": "lung adenocarcinoma"},
                                ],
                            }
                        ]
                    }
                }
            }

    root = tmp_path
    raw_root = root / "data/raw/opentargets/test"
    raw_root.mkdir(parents=True)
    mapping = SimpleNamespace(disease_id="EFO_OLD", expected_disease_name="lung adenocarcinoma")
    out = resolve_disease_for_release(
        client=FakeClient(), root=root, raw_root=raw_root, cancer_id="CANCER-001", mapping=mapping
    )
    assert out["resolved_id"] == "MONDO_NEW"
    assert out["resolved_name"] == "lung adenocarcinoma"
    assert out["resolution_method"] == "mapIds_exact_name"
    assert out["id_changed"] is True


def test_runtime_disease_resolver_refuses_ambiguous_exact_names(tmp_path):
    from types import SimpleNamespace
    from mcl.analysis.opentargets import resolve_disease_for_release

    class FakeClient:
        def fetch_disease(self, disease_id):
            return {"data": {"disease": None}}

        def resolve_disease_names(self, names):
            return {
                "data": {
                    "mapIds": {
                        "mappings": [
                            {
                                "term": names[0],
                                "hits": [
                                    {"id": "MONDO_A", "name": names[0]},
                                    {"id": "EFO_B", "name": names[0]},
                                ],
                            }
                        ]
                    }
                }
            }

    root = tmp_path
    raw_root = root / "data/raw/opentargets/test"
    raw_root.mkdir(parents=True)
    mapping = SimpleNamespace(disease_id="EFO_OLD", expected_disease_name="ambiguous disease")
    out = resolve_disease_for_release(
        client=FakeClient(), root=root, raw_root=raw_root, cancer_id="CANCER-X", mapping=mapping
    )
    assert out["resolved_id"] is None
    assert out["resolution_method"] == "unresolved"
    assert "found 2" in out["resolution_error"]
