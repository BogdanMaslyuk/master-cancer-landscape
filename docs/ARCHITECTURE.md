# Architecture

## Data flow

1. **Official sources** — authoritative downloadable/API data.
2. **`data/raw/`** — immutable release/date-stamped snapshots.
3. **Normalization** — typed records with explicit mapping status.
4. **`data/processed/`** — canonical Parquet tables plus human-readable TSV mirrors.
5. **QC + provenance** — every normalized field remains traceable to its raw source.
6. **Excel export** — later milestone; Excel is a presentation/curation layer, not raw storage.

## Identity rule

- Canonical target key: HGNC approved symbol.
- Stable technical identifier: HGNC ID.
- Exact approved symbol wins.
- Previous symbol, alias and withdrawn mappings are explicit mapping classes.
- Ambiguous mappings are never guessed.

## Milestone 1 boundaries

This milestone normalizes target identity only. It does not calculate biological priority,
clinical validation, tumour selectivity, druggability scores, or composite rankings.
