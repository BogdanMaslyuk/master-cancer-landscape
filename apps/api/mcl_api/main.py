from __future__ import annotations

from functools import lru_cache
from time import perf_counter

from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.middleware.cors import CORSMiddleware

from .routers.genes import router as genes_router
from .settings import MCL_ROOT
from .state import atlas_store, cohort_store, multiomics_store, store
from .store import MCLDataError


app = FastAPI(
    title="MCL Explorer API",
    version="0.7.0",
    description="Read-only API over Master Cancer Landscape processed outputs.",
)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:3000", "http://127.0.0.1:3000"],
    allow_credentials=True,
    allow_methods=["GET"],
    allow_headers=["*"],
)
app.include_router(genes_router)


@app.middleware("http")
async def add_mcl_timing(request: Request, call_next):
    start = perf_counter()
    response = await call_next(request)
    elapsed_ms = (perf_counter() - start) * 1000.0
    response.headers["X-MCL-Process-Time-Ms"] = f"{elapsed_ms:.1f}"
    if request.method == "GET" and request.url.path.startswith("/api/"):
        response.headers.setdefault("Cache-Control", "public, max-age=300")
    return response


def _guard(call):
    try:
        return call()
    except MCLDataError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    except (ValueError, KeyError) as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


def _split_genes(value: str | None) -> tuple[str, ...]:
    if not value:
        return ()
    genes: list[str] = []
    for part in value.replace(";", ",").split(","):
        gene = part.strip().upper()
        if gene and gene not in genes:
            genes.append(gene)
    return tuple(genes[:100])


@lru_cache(maxsize=1)
def _summary_cached():
    return store.summary()


@lru_cache(maxsize=1)
def _atlas_cached():
    return atlas_store.atlas()


@lru_cache(maxsize=32)
def _cohort_cached(cancer_id: str):
    return cohort_store.summary(cancer_id)


@lru_cache(maxsize=32)
def _context_cached(cancer_id: str):
    return atlas_store.context(cancer_id)


@lru_cache(maxsize=256)
def _models_cached(
    cancer_id: str | None,
    group: str | None,
    search: str | None,
    sequencing_only: bool,
    limit: int,
):
    return atlas_store.models(
        cancer_id=cancer_id,
        group=group,
        search=search,
        sequencing_only=sequencing_only,
        limit=limit,
    )


@lru_cache(maxsize=512)
def _model_cached(model_id: str):
    return atlas_store.model(model_id)


@lru_cache(maxsize=1)
def _multiomics_availability_cached():
    return multiomics_store.availability()


@lru_cache(maxsize=128)
def _context_multiomics_cached(cancer_id: str, genes: tuple[str, ...], limit: int):
    return multiomics_store.context(cancer_id, list(genes) or None, limit=limit)


@lru_cache(maxsize=512)
def _model_multiomics_cached(model_id: str, genes: tuple[str, ...], limit: int):
    return multiomics_store.model(model_id, list(genes) or None, limit=limit)


@lru_cache(maxsize=1)
def _comparisons_cached():
    return store.comparison_specs()


@lru_cache(maxsize=64)
def _comparison_cached(comparison_id: str):
    return store.comparison(comparison_id)


@lru_cache(maxsize=512)
def _comparison_genes_cached(
    comparison_id: str,
    page: int,
    page_size: int,
    search: str | None,
    negative_delta_only: bool,
    exclude_broad: bool,
    exclude_low_sample: bool,
    fdr_only: bool,
    stable_only: bool,
    sort_by: str,
    sort_order: str,
):
    return store.comparison_genes(
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


@lru_cache(maxsize=256)
def _pathways_cached(
    top_n: int | None,
    source: str | None,
    stable_only: bool,
    significant_only: bool,
    search: str | None,
    limit: int,
):
    return store.pathways(
        top_n=top_n,
        source=source,
        stable_only=stable_only,
        significant_only=significant_only,
        search=search,
        limit=limit,
    )


@lru_cache(maxsize=1)
def _pathway_stability_cached():
    return store.pathway_stability()


@lru_cache(maxsize=32)
def _network_cached(stable_only: bool, limit_terms: int):
    return store.network(stable_only=stable_only, limit_terms=limit_terms)


@lru_cache(maxsize=1)
def _qc_cached():
    return store.qc()


@app.get("/health")
def health():
    return {"status": "ok", "mcl_root": str(MCL_ROOT)}


@app.get("/api/summary")
def summary():
    return _guard(_summary_cached)


@app.get("/api/atlas")
def atlas():
    return _guard(_atlas_cached)


@app.get("/api/multiomics")
def multiomics_availability():
    return _guard(_multiomics_availability_cached)


@app.get("/api/atlas/{cancer_id}/cohort")
def cancer_model_cohort(cancer_id: str):
    return _guard(lambda: _cohort_cached(cancer_id))


@app.get("/api/atlas/{cancer_id}/multiomics")
def cancer_multiomics(
    cancer_id: str,
    genes: str | None = None,
    limit: int = Query(12, ge=1, le=50),
):
    gene_tuple = _split_genes(genes)
    return _guard(lambda: _context_multiomics_cached(cancer_id, gene_tuple, limit))


@app.get("/api/atlas/{cancer_id}")
def cancer_context(cancer_id: str):
    return _guard(lambda: _context_cached(cancer_id))


@app.get("/api/models")
def models(
    cancer_id: str | None = None,
    group: str | None = None,
    search: str | None = None,
    sequencing_only: bool = False,
    limit: int = Query(1000, ge=1, le=5000),
):
    return _guard(lambda: _models_cached(cancer_id, group, search, sequencing_only, limit))


@app.get("/api/models/{model_id}/multiomics")
def model_multiomics(
    model_id: str,
    genes: str | None = None,
    limit: int = Query(30, ge=1, le=100),
):
    gene_tuple = _split_genes(genes)
    return _guard(lambda: _model_multiomics_cached(model_id, gene_tuple, limit))


@app.get("/api/models/{model_id}")
def model(model_id: str):
    return _guard(lambda: _model_cached(model_id))


@app.get("/api/comparisons")
def comparisons():
    return _guard(_comparisons_cached)


@app.get("/api/comparisons/{comparison_id}")
def comparison(comparison_id: str):
    return _guard(lambda: _comparison_cached(comparison_id))


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
        lambda: _comparison_genes_cached(
            comparison_id,
            page,
            page_size,
            search,
            negative_delta_only,
            exclude_broad,
            exclude_low_sample,
            fdr_only,
            stable_only,
            sort_by,
            sort_order,
        )
    )


@app.get("/api/pathways")
def pathways(
    top_n: int | None = None,
    source: str | None = None,
    stable_only: bool = False,
    significant_only: bool = True,
    search: str | None = None,
    limit: int = Query(1000, ge=1, le=5000),
):
    return _guard(lambda: _pathways_cached(top_n, source, stable_only, significant_only, search, limit))


@app.get("/api/pathways/stability")
def pathway_stability():
    return _guard(_pathway_stability_cached)


@app.get("/api/network")
def network(stable_only: bool = True, limit_terms: int = Query(100, ge=1, le=500)):
    return _guard(lambda: _network_cached(stable_only, limit_terms))


@app.get("/api/qc")
def qc():
    return _guard(_qc_cached)
