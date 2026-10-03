from pathlib import Path

from mcl.models import MappingStatus, TargetSeed
from mcl.sources.hgnc import HGNCLocalResolver

FIXTURES = Path(__file__).parent / "fixtures"


def resolver() -> HGNCLocalResolver:
    return HGNCLocalResolver(
        FIXTURES / "hgnc_complete_set_minimal.tsv",
        FIXTURES / "withdrawn_minimal.tsv",
    )


def test_exact_approved_symbol():
    r = resolver().resolve("KRAS")
    assert r.status == MappingStatus.EXACT_APPROVED_SYMBOL
    assert r.record["hgnc_id"] == "HGNC:6407"


def test_previous_symbol():
    r = resolver().resolve("KRAS2")
    assert r.status == MappingStatus.PREVIOUS_SYMBOL
    assert r.record["symbol"] == "KRAS"


def test_cox2_regression_resolves_to_ptgs2():
    r = resolver().resolve("COX2")
    assert r.status == MappingStatus.ALIAS_SYMBOL
    assert r.record["symbol"] == "PTGS2"


def test_withdrawn_merged():
    r = resolver().resolve("OLDKRAS")
    assert r.status == MappingStatus.WITHDRAWN_MERGED
    assert r.record["symbol"] == "KRAS"


def test_ambiguous_alias_never_guessed():
    r = resolver().resolve("DUPALIAS")
    assert r.status == MappingStatus.MANUAL_REVIEW_REQUIRED
    assert r.record is None


def test_withdrawn_split_never_guessed():
    r = resolver().resolve("OLDSPLIT")
    assert r.status == MappingStatus.MANUAL_REVIEW_REQUIRED
    assert r.record is None


def test_target_identity_has_crossrefs():
    seed = TargetSeed(Target_ID="TGT-KRAS", input_symbol="KRAS")
    x = resolver().normalize(seed)
    assert x.hgnc_symbol == "KRAS"
    assert x.uniprot_id == "P01116"
    assert x.ensembl_gene_id == "ENSG00000133703"
    assert x.ncbi_gene_id == "3845"
