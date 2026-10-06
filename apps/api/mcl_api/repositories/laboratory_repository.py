from __future__ import annotations

from ..laboratory import LaboratoryPanelStore


class LaboratoryRepository:
    def __init__(self, store: LaboratoryPanelStore):
        self.store = store

    def summary(self):
        return self.store.summary()

    def lines(self, **kwargs):
        return self.store.lines(**kwargs)

    def candidates(self, **kwargs):
        return self.store.candidate_list(**kwargs)

    def mechanism_panels(self):
        return self.store.mechanism_panels()

    def candidate_detail(self, hypothesis_id: str):
        return self.store.candidate_detail(hypothesis_id)
