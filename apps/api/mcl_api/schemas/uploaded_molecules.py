from __future__ import annotations

from pydantic import BaseModel, Field

from .common import ApiResponseModel


class MoleculeAnalyzeRequest(BaseModel):
    smiles: str = Field(min_length=1, max_length=5000)
    name: str | None = Field(default=None, max_length=200)
    force_recompute: bool = False


class MoleculeAnalysisResponse(ApiResponseModel):
    pass


class MoleculeDetailResponse(ApiResponseModel):
    pass


class MoleculeTargetsResponse(ApiResponseModel):
    pass


class MoleculeCellLinesResponse(ApiResponseModel):
    pass


class MoleculeExperimentResponse(ApiResponseModel):
    pass
