from __future__ import annotations

from mcl_api.main import app


EXPECTED_RESPONSE_MODELS = {
    "/api/summary": "OverviewResponse",
    "/api/genes/search": "GeneSearchResponse",
    "/api/genes/facets": "GeneFacetResponse",
    "/api/gene-matrix": "GeneMatrixResponse",
    "/api/genes/{gene_symbol}": "GeneDetailResponse",
    "/api/atlas": "AtlasResponse",
    "/api/multiomics": "MultiomicsResponse",
    "/api/atlas/{cancer_id}/cohort": "CohortResponse",
    "/api/atlas/{cancer_id}/multiomics": "MultiomicsResponse",
    "/api/atlas/{cancer_id}": "CancerContextResponse",
    "/api/models": "ModelsResponse",
    "/api/models/{model_id}/multiomics": "ModelMultiomicsResponse",
    "/api/models/{model_id}": "ModelDetailResponse",
    "/api/comparisons": "ComparisonListResponse",
    "/api/comparisons/{comparison_id}": "ComparisonSummary",
    "/api/comparisons/{comparison_id}/genes": "ComparisonGenesResponse",
    "/api/pathways": "PathwayListResponse",
    "/api/pathways/stability": "PathwayStabilityResponse",
    "/api/network": "NetworkResponse",
    "/api/qc": "QCResponse",
}


GENE_RESPONSE_MODELS = {
    "GeneSearchResponse",
    "GeneFacetResponse",
    "GeneMatrixResponse",
    "GeneDetailResponse",
}


def test_explorer_routes_publish_named_response_models():
    schema = app.openapi()

    for path, model_name in EXPECTED_RESPONSE_MODELS.items():
        response_schema = schema["paths"][path]["get"]["responses"]["200"]["content"][
            "application/json"
        ]["schema"]
        assert response_schema == {"$ref": f"#/components/schemas/{model_name}"}
        assert model_name in schema["components"]["schemas"]


def test_gene_response_models_preserve_untyped_legacy_fields():
    schema = app.openapi()

    for model_name in GENE_RESPONSE_MODELS:
        model_schema = schema["components"]["schemas"][model_name]
        assert model_schema.get("additionalProperties") is True
