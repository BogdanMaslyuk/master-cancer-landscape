from __future__ import annotations

from mcl_api.repositories.pathway_repository import PathwayRepository


class DummyPathwayStore:
    def pathway_stability(self):
        return [
            {
                "term_id": "A",
                "significant_thresholds_n": 3,
                "significant_all_thresholds": True,
            },
            {
                "term_id": "B",
                "significant_thresholds_n": 2,
                "significant_all_thresholds": False,
            },
            {
                "term_id": "C",
                "significant_thresholds_n": 1,
                "significant_all_thresholds": False,
            },
            {
                "term_id": "D",
                "significant_thresholds_n": None,
                "significant_all_thresholds": None,
            },
        ]


def test_pathway_stability_can_return_only_fully_stable_terms():
    repository = PathwayRepository(DummyPathwayStore())

    rows = repository.stability(stable_only=True, limit=50)

    assert [row["term_id"] for row in rows] == ["A"]


def test_pathway_stability_can_filter_by_recurrence_threshold():
    repository = PathwayRepository(DummyPathwayStore())

    rows = repository.stability(min_significant_thresholds=2, limit=50)

    assert [row["term_id"] for row in rows] == ["A", "B"]


def test_pathway_stability_applies_response_limit_after_filtering():
    repository = PathwayRepository(DummyPathwayStore())

    rows = repository.stability(min_significant_thresholds=1, limit=2)

    assert [row["term_id"] for row in rows] == ["A", "B"]
