from __future__ import annotations

import base64
import binascii

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

    def import_csv(
        self, experiment_id: str, csv_text: str | None, filename: str,
        file_base64: str | None = None,
    ):
        if (csv_text is None) == (file_base64 is None):
            raise ValueError("Supply exactly one of csv_text or file_base64")
        raw_bytes = None
        if file_base64 is not None:
            try:
                raw_bytes = base64.b64decode(file_base64, validate=True)
            except (binascii.Error, ValueError) as exc:
                raise ValueError("Invalid base64 encoding of raw CSV") from exc
        return self.repository.import_csv(
            experiment_id, csv_text, filename, raw_bytes=raw_bytes
        )
