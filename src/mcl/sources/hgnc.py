from __future__ import annotations

import csv
import json
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import httpx
import time

from mcl.models import MappingStatus, TargetIdentity, TargetSeed
from mcl.utils.hash import sha256_file, sha256_text

HGNC_APPROVED_URL = (
    "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
)
HGNC_WITHDRAWN_URL = (
    "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/withdrawn.txt"
)
PARSER_VERSION = "hgnc-tsv-v1"


def _split_pipe(value: str | None) -> list[str]:
    if not value:
        return []
    return [x.strip() for x in value.split("|") if x.strip()]


def _first_present(row: dict[str, str], *names: str) -> str | None:
    for name in names:
        value = (row.get(name) or "").strip()
        if value:
            return value.split("|")[0].strip()
    return None


def _canonical_record_hash(row: dict[str, str]) -> str:
    payload = json.dumps(row, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    return sha256_text(payload)


@dataclass(frozen=True)
class Resolution:
    status: MappingStatus
    basis: str
    record: dict[str, str] | None
    note: str | None = None


class HGNCLocalResolver:
    """Resolve targets against immutable HGNC approved + withdrawn TSV snapshots."""

    def __init__(self, approved_tsv: str | Path, withdrawn_tsv: str | Path | None = None):
        self.approved_tsv = Path(approved_tsv)
        self.withdrawn_tsv = Path(withdrawn_tsv) if withdrawn_tsv else None
        self.by_symbol: dict[str, dict[str, str]] = {}
        self.by_hgnc_id: dict[str, dict[str, str]] = {}
        self.by_prev: dict[str, list[dict[str, str]]] = defaultdict(list)
        self.by_alias: dict[str, list[dict[str, str]]] = defaultdict(list)
        self.withdrawn: dict[str, list[str]] = defaultdict(list)
        self._load()

    def _load(self) -> None:
        with self.approved_tsv.open("r", encoding="utf-8-sig", newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                symbol = (row.get("symbol") or "").strip()
                hgnc_id = (row.get("hgnc_id") or "").strip()
                if not symbol:
                    continue
                self.by_symbol[symbol.upper()] = row
                if hgnc_id:
                    self.by_hgnc_id[hgnc_id.upper()] = row
                for item in _split_pipe(row.get("prev_symbol")):
                    self.by_prev[item.upper()].append(row)
                for item in _split_pipe(row.get("alias_symbol")):
                    self.by_alias[item.upper()].append(row)

        if self.withdrawn_tsv and self.withdrawn_tsv.exists():
            with self.withdrawn_tsv.open("r", encoding="utf-8-sig", newline="") as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    symbol = (row.get("WITHDRAWN_SYMBOL") or "").strip().upper()
                    merged = (row.get("MERGED_INTO_REPORT(S) (i.e HGNC_ID|SYMBOL|STATUS)") or "").strip()
                    if not symbol or not merged:
                        continue
                    replacements: list[str] = []
                    for part in merged.split(","):
                        bits = [x.strip() for x in part.split("|")]
                        if len(bits) >= 2 and bits[1]:
                            replacements.append(bits[1].upper())
                    self.withdrawn[symbol].extend(replacements)

    def resolve(self, input_symbol: str) -> Resolution:
        key = input_symbol.strip().upper()
        if not key:
            return Resolution(MappingStatus.NOT_FOUND, "empty_input", None, "Empty symbol")

        if key in self.by_symbol:
            return Resolution(MappingStatus.EXACT_APPROVED_SYMBOL, "approved_symbol", self.by_symbol[key])

        prev = self.by_prev.get(key, [])
        if len(prev) == 1:
            return Resolution(MappingStatus.PREVIOUS_SYMBOL, "previous_symbol", prev[0])
        if len(prev) > 1:
            return Resolution(
                MappingStatus.MANUAL_REVIEW_REQUIRED,
                "previous_symbol_ambiguous",
                None,
                f"Previous symbol maps to {len(prev)} approved records",
            )

        aliases = self.by_alias.get(key, [])
        if len(aliases) == 1:
            return Resolution(MappingStatus.ALIAS_SYMBOL, "alias_symbol", aliases[0])
        if len(aliases) > 1:
            return Resolution(
                MappingStatus.MANUAL_REVIEW_REQUIRED,
                "alias_symbol_ambiguous",
                None,
                f"Alias maps to {len(aliases)} approved records",
            )

        replacements = sorted(set(self.withdrawn.get(key, [])))
        if len(replacements) == 1 and replacements[0] in self.by_symbol:
            return Resolution(
                MappingStatus.WITHDRAWN_MERGED,
                "withdrawn_merged_into",
                self.by_symbol[replacements[0]],
            )
        if replacements:
            return Resolution(
                MappingStatus.MANUAL_REVIEW_REQUIRED,
                "withdrawn_split_or_unresolved",
                None,
                f"Withdrawn symbol maps to: {', '.join(replacements)}",
            )
        return Resolution(MappingStatus.NOT_FOUND, "not_found", None, "No HGNC mapping found")

    def normalize(self, seed: TargetSeed) -> TargetIdentity:
        resolution = self.resolve(seed.input_symbol)
        if resolution.record is None:
            return TargetIdentity(
                target_id=seed.target_id,
                input_symbol=seed.input_symbol,
                mapping_status=resolution.status,
                mapping_basis=resolution.basis,
                source_snapshot=str(self.approved_tsv),
                notes=resolution.note,
            )

        row = resolution.record
        return TargetIdentity(
            target_id=seed.target_id,
            input_symbol=seed.input_symbol,
            mapping_status=resolution.status,
            mapping_basis=resolution.basis,
            hgnc_id=_first_present(row, "hgnc_id"),
            hgnc_symbol=_first_present(row, "symbol"),
            hgnc_approved_name=_first_present(row, "name"),
            hgnc_status=_first_present(row, "status"),
            hgnc_date_modified=_first_present(row, "date_modified"),
            hgnc_previous_symbols=_split_pipe(row.get("prev_symbol")),
            hgnc_aliases=_split_pipe(row.get("alias_symbol")),
            ensembl_gene_id=_first_present(row, "ensembl_gene_id"),
            ncbi_gene_id=_first_present(row, "entrez_id"),
            uniprot_id=_first_present(row, "uniprot_ids", "uniprot_id"),
            raw_record_hash=_canonical_record_hash(row),
            source_snapshot=str(self.approved_tsv),
            notes=resolution.note,
        )


def _download(client: httpx.Client, url: str, destination: Path, attempts: int = 3) -> None:
    last_error: Exception | None = None
    for attempt in range(1, attempts + 1):
        try:
            with client.stream("GET", url, follow_redirects=True) as response:
                response.raise_for_status()
                destination.parent.mkdir(parents=True, exist_ok=True)
                tmp = destination.with_suffix(destination.suffix + ".part")
                with tmp.open("wb") as fh:
                    for chunk in response.iter_bytes():
                        fh.write(chunk)
                tmp.replace(destination)
                return
        except Exception as exc:
            last_error = exc
            destination.with_suffix(destination.suffix + ".part").unlink(missing_ok=True)
            if attempt < attempts:
                time.sleep(min(2 ** (attempt - 1), 8))
    assert last_error is not None
    raise last_error


def fetch_current_snapshot(raw_root: str | Path, snapshot_date: str | None = None) -> dict[str, str]:
    snapshot_date = snapshot_date or datetime.now(timezone.utc).date().isoformat()
    out_dir = Path(raw_root) / "hgnc" / snapshot_date
    approved = out_dir / "hgnc_complete_set.txt"
    withdrawn = out_dir / "withdrawn.txt"
    if approved.exists() or withdrawn.exists():
        raise FileExistsError(
            f"Snapshot directory already contains files: {out_dir}. Raw snapshots are immutable."
        )

    with httpx.Client(timeout=120, headers={"User-Agent": "MasterCancerLandscape/0.1"}) as client:
        _download(client, HGNC_APPROVED_URL, approved)
        _download(client, HGNC_WITHDRAWN_URL, withdrawn)

    meta = {
        "source": "HGNC",
        "snapshot_date": snapshot_date,
        "approved_url": HGNC_APPROVED_URL,
        "withdrawn_url": HGNC_WITHDRAWN_URL,
        "approved_file": str(approved),
        "withdrawn_file": str(withdrawn),
        "approved_sha256": sha256_file(approved),
        "withdrawn_sha256": sha256_file(withdrawn),
        "approved_bytes": str(approved.stat().st_size),
        "withdrawn_bytes": str(withdrawn.stat().st_size),
        "retrieved_at": datetime.now(timezone.utc).isoformat(),
        "parser_version": PARSER_VERSION,
    }
    (out_dir / "snapshot_metadata.json").write_text(
        json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    return meta
