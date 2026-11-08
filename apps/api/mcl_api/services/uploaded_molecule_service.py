from __future__ import annotations

from ..repositories.uploaded_molecule_repository import UploadedMoleculeRepository


class UploadedMoleculeService:
    def __init__(self, repository: UploadedMoleculeRepository):
        self.repository = repository

    def analyze(
        self,
        smiles: str,
        name: str | None = None,
        force_recompute: bool = False,
    ):
        return self.repository.analyze(smiles, name, force_recompute)

    def detail(self, analysis_id: str):
        return self.repository.detail(analysis_id)

    def targets(self, analysis_id: str):
        return self.repository.targets(analysis_id)

    def cell_lines(self, analysis_id: str, target_gene: str | None, laboratory_only: bool):
        return self.repository.cell_lines(analysis_id, target_gene, laboratory_only)

    def experiment(self, analysis_id: str):
        return self.repository.experiment(analysis_id)

    def smiles(self, analysis_id: str):
        return self.repository.smiles(analysis_id)
