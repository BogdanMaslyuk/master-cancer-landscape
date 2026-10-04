from __future__ import annotations

from typing import Any


class PathwayRepository:
    """Read-only data access for pathway enrichment, stability and networks."""

    def __init__(self, store: Any):
        self.store = store

    def pathways(
        self,
        *,
        top_n: int | None,
        source: str | None,
        stable_only: bool,
        significant_only: bool,
        search: str | None,
        limit: int,
    ):
        return self.store.pathways(
            top_n=top_n,
            source=source,
            stable_only=stable_only,
            significant_only=significant_only,
            search=search,
            limit=limit,
        )

    def stability(
        self,
        *,
        stable_only: bool = False,
        min_significant_thresholds: int | None = None,
        limit: int = 5000,
    ):
        """Return only the stability rows needed by the interactive view.

        The complete term-stability table is >2 MB in the current dataset. Sending
        that full JSON payload to Next.js for every page transition created avoidable
        latency and made unrelated pages appear intermittently unavailable. Keep the
        scientific table intact, but filter the serving payload before serialization.
        """
        rows = self.store.pathway_stability()
        if stable_only:
            rows = [
                row
                for row in rows
                if str(row.get("significant_all_thresholds", "")).strip().lower()
                in {"true", "1", "yes"}
            ]
        if min_significant_thresholds is not None:
            threshold = max(0, int(min_significant_thresholds))
            rows = [
                row
                for row in rows
                if _as_int(row.get("significant_thresholds_n")) >= threshold
            ]
        return rows[: max(1, min(int(limit), 5000))]

    def network(self, *, stable_only: bool, limit_terms: int):
        return self.store.network(stable_only=stable_only, limit_terms=limit_terms)


def _as_int(value: Any) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return 0
