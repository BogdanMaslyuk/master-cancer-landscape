from __future__ import annotations

from ..phenotype import PhenotypeStore


class PhenotypeRepository:
    def __init__(self, store: PhenotypeStore):
        self.store = store

    def models(self):
        return self.store.models()

    def compounds(self):
        return self.store.compounds()

    def create(self, payload: dict):
        return self.store.create(payload)

    def detail(self, experiment_id: str):
        return self.store.detail(experiment_id)

    def import_csv(
        self, experiment_id: str, csv_text: str | None, filename: str,
        *, raw_bytes: bytes | None = None,
    ):
        return self.store.import_csv(
            experiment_id, csv_text, filename, raw_bytes=raw_bytes
        )
