from __future__ import annotations

from typing import Literal
from pydantic import BaseModel, ConfigDict


class QCRecord(BaseModel):
    model_config = ConfigDict(extra="forbid")

    severity: Literal["INFO", "WARNING", "ERROR"]
    check: str
    entity_id: str | None = None
    field: str | None = None
    observed: str | None = None
    expected: str | None = None
    message: str
