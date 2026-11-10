from __future__ import annotations

from ..repositories.phenotype_repository import PhenotypeRepository


class PhenotypeService:
    def __init__(self, repository: PhenotypeRepository):
        self.repository = repository

    def models(self):
        return self.repository.models()

    def compounds(self):
        return self.repository.compounds()

    def create(self, payload: dict):
        return self.repository.create(payload)

    def detail(self, experiment_id: str):
        return self.repository.detail(experiment_id)

    def import_csv(self, experiment_id: str, csv_text: str, filename: str):
        return self.repository.import_csv(experiment_id, csv_text, filename)
