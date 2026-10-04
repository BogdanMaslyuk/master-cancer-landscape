from __future__ import annotations

from mcl_api.main import app


EXPECTED_RESPONSE_MODELS = {
    "/api/summary": "OverviewResponse",
    "/api/genes": "GeneListResponse",
    "/api/genes/suggest": "GeneSuggestResponse",
    "/api/genes/stable": "StableGenesResponse",
    "/api/genes/search": "GeneSearchResponse",
    "/api/genes/facets": "GeneFacetResponse",
    "/api/gene-matrix": "GeneMatrixResponse",
    "/api/genes/{gene_symbol}/contexts": "GeneContextsResponse",
    "/api/genes/{gene_symbol}/models": "GeneModelsResponse",
    "/api/genes/{gene_symbol}/annotations": "GeneAnnotationsResponse",
    "/api/genes/{gene_symbol}/mutation-associations": "GeneMutationAssociationsResponse",
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


def test_every_api_get_route_has_explicit_named_response_schema():
    schema = app.openapi()
    missing: list[str] = []

    for path, operations in schema["paths"].items():
        if not path.startswith("/api/") or "get" not in operations:
            continue
        response_schema = operations["get"]["responses"]["200"]["content"]["application/json"][
            "schema"
        ]
        ref = response_schema.get("$ref") if isinstance(response_schema, dict) else None
        if not ref or not ref.startswith("#/components/schemas/"):
            missing.append(path)

    assert missing == [], f"API GET routes must expose named response models: {missing}"


def test_gene_response_models_preserve_untyped_legacy_fields():
    schema = app.openapi()

    for model_name in GENE_RESPONSE_MODELS:
        model_schema = schema["components"]["schemas"][model_name]
        assert model_schema.get("additionalProperties") is True
