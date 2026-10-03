from pathlib import Path

from mcl.normalize.targets import load_target_seeds, normalize_targets
from mcl.qc.targets import target_normalization_qc

FIXTURES = Path(__file__).parent / "fixtures"


def test_fixture_normalization_passes_without_errors():
    seeds = load_target_seeds(FIXTURES / "targets_seed_fixture.tsv")
    identities, _ = normalize_targets(
        FIXTURES / "targets_seed_fixture.tsv",
        FIXTURES / "hgnc_complete_set_minimal.tsv",
        FIXTURES / "withdrawn_minimal.tsv",
        source_release="test-fixture",
    )
    qc = target_normalization_qc(seeds, identities)
    assert not [x for x in qc if x.severity == "ERROR"]
    assert {x.hgnc_symbol for x in identities} == {"KRAS", "PTGS2", "BRAF"}
