from __future__ import annotations

from typing import Any

from pydantic import RootModel

from .common import ApiResponseModel


class PathwayListResponse(RootModel[list[dict[str, Any]]]):
    pass


class PathwayStabilityResponse(RootModel[list[dict[str, Any]]]):
    pass


class NetworkResponse(ApiResponseModel):
    pass
