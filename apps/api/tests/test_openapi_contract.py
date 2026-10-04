from __future__ import annotations

from mcl_api.main import app


EXPECTED_GENE_RESPONSE_MODELS = {
    "/api/genes/search": "GeneSearchResponse",
    "/api/genes/facets": "GeneFacetResponse",
    "/api/gene-matrix": "GeneMatrixResponse",
    "/api/genes/{gene_symbol}": "GeneDetailResponse",
}


def test_gene_explorer_routes_publish_named_response_models():
    schema = app.openapi()

    for path, model_name in EXPECTED_GENE_RESPONSE_MODELS.items():
        response_schema = schema["paths"][path]["get"]["responses"]["200"]["content"][
            "application/json"
        ]["schema"]
        assert response_schema == {"$ref": f"#/components/schemas/{model_name}"}
        assert model_name in schema["components"]["schemas"]


def test_gene_response_models_preserve_untyped_legacy_fields():
    schema = app.openapi()

    for model_name in EXPECTED_GENE_RESPONSE_MODELS.values():
        model_schema = schema["components"]["schemas"][model_name]
        assert model_schema.get("additionalProperties") is True
