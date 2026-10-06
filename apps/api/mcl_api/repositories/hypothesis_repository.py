from __future__ import annotations

from ..candidate_hypotheses import CandidateHypothesisStore


class HypothesisRepository:
    def __init__(self, store: CandidateHypothesisStore):
        self.store = store

    def summary(self):
        return self.store.summary()

    def search(
        self,
        *,
        search: str | None,
        status: str | None,
        target_gene: str | None,
        cancer: str | None,
        limit: int,
        offset: int,
    ):
        return self.store.search(
            search=search,
            status=status,
            target_gene=target_gene,
            cancer=cancer,
            limit=limit,
            offset=offset,
        )

    def detail(self, hypothesis_id: str):
        return self.store.detail(hypothesis_id)
