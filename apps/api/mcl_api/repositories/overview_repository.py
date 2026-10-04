from __future__ import annotations

from typing import Any


class OverviewRepository:
    """Read-only data access for the Explorer-wide summary payload."""

    def __init__(self, store: Any):
        self.store = store

    def summary(self):
        return self.store.summary()
