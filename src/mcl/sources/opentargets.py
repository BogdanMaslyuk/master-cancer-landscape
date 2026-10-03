from __future__ import annotations

import json
import logging
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import httpx

from mcl.utils.hash import sha256_bytes
from mcl.utils.io import write_json

LOGGER = logging.getLogger(__name__)

DEFAULT_ENDPOINT = "https://api.platform.opentargets.org/api/v4/graphql"
PARSER_VERSION = "opentargets-v0.2.2"

META_QUERY = """
query PlatformMeta {
  meta {
    apiVersion { x y z }
    dataVersion { year month }
  }
}
"""

DISEASE_QUERY = """
query DiseaseInfo($efoId: String!) {
  disease(efoId: $efoId) {
    id
    name
  }
}
"""

MAP_DISEASE_NAME_QUERY = """
query ResolveDiseaseNames($queryTerms: [String!]!) {
  mapIds(queryTerms: $queryTerms, entityNames: ["disease"]) {
    mappings {
      term
      hits {
        id
        name
      }
    }
  }
}
"""

ASSOCIATED_TARGETS_QUERY = """
query DiseaseAssociatedTargets($efoId: String!, $index: Int!, $size: Int!) {
  disease(efoId: $efoId) {
    id
    name
    associatedTargets(enableIndirect: false, page: {index: $index, size: $size}) {
      count
      rows {
        target { id approvedSymbol }
        score
        datasourceScores { id score }
        datatypeScores { id score }
      }
    }
  }
}
"""

TARGET_TRACTABILITY_QUERY = """
query TargetTractability($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    id
    approvedSymbol
    tractability { label modality value }
  }
}
"""


# Backwards-compatible constant name used by existing tests/docs.
TRACTABILITY_QUERY = TARGET_TRACTABILITY_QUERY


class OpenTargetsAPIError(RuntimeError):
    pass


class OpenTargetsClient:
    def __init__(
        self,
        endpoint: str = DEFAULT_ENDPOINT,
        timeout_seconds: float = 60.0,
        retries: int = 3,
        retry_backoff_seconds: float = 1.0,
    ) -> None:
        self.endpoint = endpoint
        self.timeout_seconds = timeout_seconds
        self.retries = retries
        self.retry_backoff_seconds = retry_backoff_seconds

    def _post(self, query: str, variables: dict[str, Any] | None = None) -> dict[str, Any]:
        last_error: Exception | None = None
        for attempt in range(1, self.retries + 1):
            try:
                with httpx.Client(timeout=self.timeout_seconds) as client:
                    response = client.post(
                        self.endpoint,
                        json={"query": query, "variables": variables or {}},
                        headers={
                            "User-Agent": "master-cancer-landscape/0.2.2",
                            "Content-Type": "application/json",
                        },
                    )
                response.raise_for_status()
                payload = response.json()
                # GraphQL can return HTTP 200 while the requested entity is null.
                # Errors must therefore be inspected before parsing `data`.
                if payload.get("errors"):
                    raise OpenTargetsAPIError(json.dumps(payload["errors"], ensure_ascii=False))
                if "data" not in payload:
                    raise OpenTargetsAPIError("Open Targets response has no `data` field")
                return payload
            except (httpx.HTTPError, ValueError, OpenTargetsAPIError) as exc:
                last_error = exc
                LOGGER.warning("Open Targets request attempt %s/%s failed: %s", attempt, self.retries, exc)
                if attempt < self.retries:
                    time.sleep(self.retry_backoff_seconds * attempt)
        raise OpenTargetsAPIError(f"Open Targets request failed after {self.retries} attempts: {last_error}")

    def fetch_meta(self) -> dict[str, Any]:
        return self._post(META_QUERY)

    def fetch_disease(self, disease_id: str) -> dict[str, Any]:
        return self._post(DISEASE_QUERY, {"efoId": disease_id})

    def resolve_disease_names(self, names: list[str]) -> dict[str, Any]:
        return self._post(MAP_DISEASE_NAME_QUERY, {"queryTerms": names})

    def fetch_associated_targets_page(self, disease_id: str, index: int, size: int) -> dict[str, Any]:
        return self._post(
            ASSOCIATED_TARGETS_QUERY,
            {"efoId": disease_id, "index": index, "size": size},
        )

    def fetch_target_tractability(self, ensembl_id: str) -> dict[str, Any]:
        return self._post(TARGET_TRACTABILITY_QUERY, {"ensemblId": ensembl_id})


def save_raw_response(path: str | Path, payload: dict[str, Any]) -> str:
    path = Path(path)
    write_json(path, payload)
    return sha256_bytes(path.read_bytes())


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()
