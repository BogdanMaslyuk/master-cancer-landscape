from __future__ import annotations

from fastapi import APIRouter

from ..api_utils import guard
from ..schemas.phenotype import (
    PhenotypeModelsResponse,
    PhenotypeCompoundsResponse,
    PhenotypeExperimentResponse,
    PhenotypeExperimentCreateRequest,
    PhenotypeImportRequest,
    PhenotypeImportResponse,
)
from ..state import phenotype_store


router = APIRouter(prefix="/api/phenotype", tags=["phenotype"])


@router.get("/models", response_model=PhenotypeModelsResponse)
def list_phenotype_models():
    return guard(phenotype_store.models)


@router.get("/compounds", response_model=PhenotypeCompoundsResponse)
def list_phenotype_compounds():
    return guard(phenotype_store.compounds)


@router.post("/experiments", response_model=PhenotypeExperimentResponse)
def create_phenotype_experiment(request: PhenotypeExperimentCreateRequest):
    return guard(lambda: phenotype_store.create(request.model_dump()))


@router.get("/experiments/{experiment_id}", response_model=PhenotypeExperimentResponse)
def phenotype_experiment(experiment_id: str):
    return guard(lambda: phenotype_store.detail(experiment_id))


@router.post("/experiments/{experiment_id}/import", response_model=PhenotypeImportResponse)
def import_phenotype_csv(experiment_id: str, request: PhenotypeImportRequest):
    return guard(lambda: phenotype_store.import_csv(experiment_id, request.csv_text, request.filename))
