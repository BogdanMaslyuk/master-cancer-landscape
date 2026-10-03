from __future__ import annotations

from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware

from .atlas import MCLAtlas
from .settings import MCL_ROOT
from .store import MCLDataError, MCLDataStore


app = FastAPI(
    title="MCL Explorer API",
    version="0.2.0",
    description="Read-only API over Master Cancer Landscape processed outputs.",
)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:3000", "http://127.0.0.1:3000"],
    allow_credentials=True,
    allow_methods=["GET"],
    allow_headers=["*"],
)
store = MCLDataStore(MCL_ROOT)
atlas_store = MCLAtlas(MCL_ROOT, store)


def _guard(call):
    try:
        return call()
    except MCLDataError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    except (ValueError, KeyError) as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.get("/health")
def health():
    return {"status": "ok", "mcl_root": str(MCL_ROOT)}


@app.get("/api/summary")
def summary():
    return _guard(store.summary)


@app.get("/api/atlas")
def atlas():
    return _guard(atlas_store.atlas)


@app.get("/api/atlas/{cancer_id}")
def cancer_context(cancer_id: str):
    return _guard(lambda: atlas_store.context(cancer_id))


@app.get("/api/models")
def models(
    cancer_id: str | None = None,
    group: str | None = None,
    search: str | None = None,
    sequencing_only: bool = False,
    limit: int = Query(1000, ge=1, le=5000),
):
    return _guard(
        lambda: atlas_store.models(
            cancer_id=cancer_id,
            group=group,
            search=search,
            sequencing_only=sequencing_only,
            limit=limit,
        )
    )


@app.get("/api/models/{model_id}")
def model(model_id: str):
    return _guard(lambda: atlas_store.model(model_id))


@app.get("/api/comparisons")
def comparisons():
    return _guard(store.comparison_specs)


@app.get("/api/comparisons/{comparison_id}")
def comparison(comparison_id: str):
    return _guard(lambda: store.comparison(comparison_id))


@app.get("/api/comparisons/{comparison_id}/genes")
def comparison_genes(
    comparison_id: str,
    page: int = Query(1, ge=1),
    page_size: int = Query(100, ge=1, le=2000),
    search: str | None = None,
    negative_delta_only: bool = False,
    exclude_broad: bool = False,
    exclude_low_sample: bool = False,
    fdr_only: bool = False,
    stable_only: bool = False,
    sort_by: str = "delta_gene_effect",
    sort_order: str = "asc",
):
    return _guard(
        lambda: store.comparison_genes(
            comparison_id,
            page=page,
            page_size=page_size,
            search=search,
            negative_delta_only=negative_delta_only,
            exclude_broad=exclude_broad,
            exclude_low_sample=exclude_low_sample,
            fdr_only=fdr_only,
            stable_only=stable_only,
            sort_by=sort_by,
            sort_order=sort_order,
        )
    )


@app.get("/api/genes")
def genes(search: str | None = None, limit: int = Query(200, ge=1, le=2000)):
    return _guard(lambda: store.genes(search=search, limit=limit))


@app.get("/api/genes/stable")
def stable_genes():
    return _guard(store.stable_genes)


@app.get("/api/genes/{gene_symbol}")
def gene(gene_symbol: str):
    return _guard(lambda: store.gene(gene_symbol))


@app.get("/api/pathways")
def pathways(
    top_n: int | None = None,
    source: str | None = None,
    stable_only: bool = False,
    significant_only: bool = True,
    search: str | None = None,
    limit: int = Query(1000, ge=1, le=5000),
):
    return _guard(
        lambda: store.pathways(
            top_n=top_n,
            source=source,
            stable_only=stable_only,
            significant_only=significant_only,
            search=search,
            limit=limit,
        )
    )


@app.get("/api/pathways/stability")
def pathway_stability():
    return _guard(store.pathway_stability)


@app.get("/api/network")
def network(stable_only: bool = True, limit_terms: int = Query(100, ge=1, le=500)):
    return _guard(lambda: store.network(stable_only=stable_only, limit_terms=limit_terms))


@app.get("/api/qc")
def qc():
    return _guard(store.qc)
