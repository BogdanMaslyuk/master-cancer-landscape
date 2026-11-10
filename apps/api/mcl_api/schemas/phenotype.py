from __future__ import annotations

from pydantic import BaseModel, Field

from .common import ApiResponseModel


class PhenotypeModelsResponse(ApiResponseModel):
    pass


class PhenotypeCompoundsResponse(ApiResponseModel):
    pass


class PhenotypeExperimentResponse(ApiResponseModel):
    pass


class PhenotypeImportResponse(ApiResponseModel):
    pass


class PhenotypeExperimentCreateRequest(BaseModel):
    title: str = Field(min_length=1, max_length=180)
    protocol_id: str = Field(min_length=1, max_length=100)
    cell_line_ids: list[str] = Field(min_length=1, max_length=20)
    compound_ids: list[str] = Field(min_length=1, max_length=65)
    exposure_h: float = Field(gt=0, le=720)
    assay_type: str = "MTT"


class PhenotypeImportRequest(BaseModel):
    filename: str = Field(min_length=5, max_length=256)
    csv_text: str | None = Field(default=None, max_length=8_000_000)
    file_base64: str | None = Field(default=None, max_length=11_000_000)
