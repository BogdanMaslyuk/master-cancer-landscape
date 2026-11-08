from __future__ import annotations

from ..uploaded_molecules import UploadedMoleculeStore


class UploadedMoleculeRepository:
    def __init__(self, store: UploadedMoleculeStore):
        self.store = store

    def analyze(
        self,
        smiles: str,
        name: str | None = None,
        force_recompute: bool = False,
    ):
        return self.store.analyze(smiles, name, force_recompute)

    def detail(self, analysis_id: str):
        return self.store.detail(analysis_id)

    def targets(self, analysis_id: str):
        return self.store.targets(analysis_id)

    def cell_lines(self, analysis_id: str, target_gene: str | None, laboratory_only: bool):
        return self.store.cell_lines(analysis_id, target_gene, laboratory_only)

    def experiment(self, analysis_id: str):
        return self.store.experiment(analysis_id)

    def smiles(self, analysis_id: str):
        return self.store.smiles(analysis_id)
