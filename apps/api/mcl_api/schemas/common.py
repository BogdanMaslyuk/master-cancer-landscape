from __future__ import annotations

from pydantic import BaseModel, ConfigDict


class ApiResponseModel(BaseModel):
    """Base model for incremental API typing without dropping legacy fields.

    Architecture v1 is introducing explicit response contracts gradually. Existing
    payloads can contain scientifically useful fields that are not typed yet, so the
    contract keeps unknown fields instead of silently stripping them during FastAPI
    response serialization.
    """

    model_config = ConfigDict(extra="allow")
