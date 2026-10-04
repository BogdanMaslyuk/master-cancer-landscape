from __future__ import annotations

from typing import Any


class QCRepository:
    """Read-only data access for Explorer quality-control summaries."""

    def __init__(self, store: Any):
        self.store = store

    def qc(self):
        return self.store.qc()
