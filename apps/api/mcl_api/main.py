from __future__ import annotations

from time import perf_counter

from fastapi import FastAPI, Request
from fastapi.middleware.cors import CORSMiddleware

from .api_utils import guard
from .routers.atlas import router as atlas_router
from .routers.comparisons import router as comparisons_router
from .routers.crispr_catalog import router as crispr_catalog_router
from .routers.genes import router as genes_router
from .routers.hypotheses import router as hypotheses_router
from .routers.models import router as models_router
from .routers.pathways import router as pathways_router
from .routers.pharmacology import router as pharmacology_router
from .routers.qc import router as qc_router
from .schemas.overview import OverviewResponse
from .settings import MCL_ROOT
from .state import overview_service


app = FastAPI(
    title="MCL Explorer API",
    version="0.10.0",
    description="Read-only API over Master Cancer Landscape processed outputs.",
)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:3000", "http://127.0.0.1:3000"],
    allow_credentials=True,
    allow_methods=["GET"],
    allow_headers=["*"],
)

for router in (
    atlas_router,
    crispr_catalog_router,
    models_router,
    pharmacology_router,
    hypotheses_router,
    comparisons_router,
    genes_router,
    pathways_router,
    qc_router,
):
    app.include_router(router)


@app.middleware("http")
async def add_mcl_timing(request: Request, call_next):
    start = perf_counter()
    response = await call_next(request)
    elapsed_ms = (perf_counter() - start) * 1000.0
    response.headers["X-MCL-Process-Time-Ms"] = f"{elapsed_ms:.1f}"
    if request.method == "GET" and request.url.path.startswith("/api/"):
        response.headers.setdefault("Cache-Control", "public, max-age=300")
    return response


@app.get("/health")
def health():
    return {"status":"ok", "mcl_root": str(MCL_ROOT)}


@app.get("/api/summary", response_model=OverviewResponse)
def summary():
    return guard(overview_service.summary)
