from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable

import pandas as pd

from mcl.utils.hash import sha256_file


FILE_CANDIDATES: dict[str, tuple[str, ...]] = {
    "models": ("Model.csv", "Models.csv"),
    "gene_effect": ("CRISPRGeneEffect.csv",),
    "gene_dependency": ("CRISPRGeneDependency.csv",),
    "mutations": ("OmicsSomaticMutations.csv",),
    "omics_profiles": ("OmicsProfiles.csv",),
}

MODEL_ID_ALIASES = ("ModelID", "DepMap_ID", "DepMapID")
GENE_ALIASES = ("HugoSymbol", "Hugo_Symbol", "GeneSymbol", "Gene", "gene")
PROTEIN_CHANGE_ALIASES = (
    "ProteinChange",
    "Protein_Change",
    "HGVSp_Short",
    "HGVSp",
    "protein_change",
)
VARIANT_CLASS_ALIASES = (
    "VariantClassification",
    "Variant_Classification",
    "Consequence",
    "VariantType",
    "Variant_Type",
)
DRIVER_FLAG_CANDIDATES = (
    "Hotspot",
    "hotspot",
    "HessDriver",
    "LikelyDriver",
    "likely_driver",
    "Oncogenic",
    "oncogenic",
    "CosmicHotspot",
    "is_driver",
    "Driver",
)
TEXT_CANDIDATES = (
    "VariantInfo",
    "Description",
    "RescueReason",
    "DNAChange",
    "CDSChange",
    "MutationDescription",
)


def resolve_depmap_files(
    root: Path,
    release: str,
    depmap_dir: Path | None = None,
) -> dict[str, Path]:
    base = (depmap_dir or (root / "data/raw/depmap" / release)).resolve()
    out: dict[str, Path] = {}
    missing: list[str] = []
    for logical_name, candidates in FILE_CANDIDATES.items():
        found = next((base / name for name in candidates if (base / name).exists()), None)
        if found is None:
            missing.append(f"{logical_name}: {' or '.join(candidates)}")
        else:
            out[logical_name] = found
    if missing:
        raise FileNotFoundError(
            "Missing required DepMap files in " + str(base) + ": " + "; ".join(missing)
        )
    return out


def depmap_file_inventory(files: dict[str, Path], release: str) -> list[dict]:
    rows = []
    for logical_name, path in files.items():
        header = pd.read_csv(path, nrows=0).columns.tolist()
        rows.append(
            {
                "depmap_release": release,
                "logical_name": logical_name,
                "filename": path.name,
                "path": str(path),
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
                "columns_n": len(header),
                "columns_json": __import__("json").dumps(header, ensure_ascii=False),
            }
        )
    return rows


def _first_existing(columns: Iterable[str], aliases: Iterable[str]) -> str | None:
    cols = set(columns)
    for alias in aliases:
        if alias in cols:
            return alias
    return None


def load_model_metadata(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, dtype=str, low_memory=False).fillna("")
    model_col = _first_existing(df.columns, MODEL_ID_ALIASES)
    if model_col is None:
        raise ValueError(f"{path.name}: no ModelID column found")
    if model_col != "ModelID":
        df = df.rename(columns={model_col: "ModelID"})
    return df


def _truthy(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin({"1", "true", "yes", "y", "t"})


def load_sequenced_model_ids(path: Path) -> tuple[set[str], dict[str, str | None]]:
    """Return models with a WES/WGS profile and schema information.

    Since 25Q2 DepMap places the default-profile flag in OmicsProfiles.csv. If a
    default flag is available, only default WES/WGS profiles are used. If the
    flag is absent, all WES/WGS profiles are accepted and that fact is surfaced
    in schema metadata rather than guessed silently.
    """

    df = pd.read_csv(path, dtype=str, low_memory=False).fillna("")
    model_col = _first_existing(df.columns, MODEL_ID_ALIASES)
    datatype_col = _first_existing(df.columns, ("Datatype", "DataType", "datatype", "data_type"))
    default_col = _first_existing(
        df.columns,
        ("IsDefaultEntryForModel", "is_default_entry", "IsDefaultEntryModel", "IsDefaultEntry", "isDefaultEntry"),
    )
    if model_col is None or datatype_col is None:
        raise ValueError(
            f"{path.name}: expected ModelID and Datatype columns; got {list(df.columns)}"
        )
    seq = df[df[datatype_col].astype(str).str.lower().isin({"wes", "wgs"})].copy()
    if default_col is not None:
        seq = seq[_truthy(seq[default_col])]
    models = set(seq[model_col].astype(str))
    return models, {
        "model_id_column": model_col,
        "datatype_column": datatype_col,
        "default_entry_column": default_col,
    }


def parse_gene_label(label: str) -> tuple[str, str | None]:
    text = str(label).strip()
    match = re.match(r"^(.*?)\s*\((\d+)\)\s*$", text)
    if match:
        return match.group(1).strip(), match.group(2)
    return text, None


def gene_column_map(path: Path) -> dict[str, str]:
    cols = pd.read_csv(path, nrows=0).columns.tolist()
    if not cols:
        raise ValueError(f"{path.name}: empty header")
    mapping: dict[str, str] = {}
    duplicates: set[str] = set()
    for col in cols[1:]:
        symbol, _ = parse_gene_label(col)
        if symbol in mapping and mapping[symbol] != col:
            duplicates.add(symbol)
        mapping[symbol] = col
    if duplicates:
        raise ValueError(f"{path.name}: duplicate gene symbols after normalization: {sorted(duplicates)}")
    return mapping


def load_gene_matrix(path: Path, target_symbols: Iterable[str]) -> pd.DataFrame:
    header = pd.read_csv(path, nrows=0).columns.tolist()
    if not header:
        raise ValueError(f"{path.name}: empty header")
    id_col = header[0]
    mapping = gene_column_map(path)
    missing = [gene for gene in target_symbols if gene not in mapping]
    if missing:
        raise ValueError(f"{path.name}: missing target genes: {missing}")
    usecols = [id_col] + [mapping[g] for g in target_symbols]
    df = pd.read_csv(path, usecols=usecols, low_memory=False)
    df = df.rename(columns={id_col: "ModelID", **{mapping[g]: g for g in target_symbols}})
    df["ModelID"] = df["ModelID"].astype(str)
    if df["ModelID"].duplicated().any():
        dup = df.loc[df["ModelID"].duplicated(keep=False), "ModelID"].unique().tolist()
        raise ValueError(f"{path.name}: duplicate ModelID values: {dup[:10]}")
    for gene in target_symbols:
        df[gene] = pd.to_numeric(df[gene], errors="coerce")
    return df.set_index("ModelID")


def mutation_schema(path: Path) -> dict[str, object]:
    cols = pd.read_csv(path, nrows=0).columns.tolist()
    schema = {
        "model_id": _first_existing(cols, MODEL_ID_ALIASES),
        "gene": _first_existing(cols, GENE_ALIASES),
        "protein_change": _first_existing(cols, PROTEIN_CHANGE_ALIASES),
        "variant_classification": _first_existing(cols, VARIANT_CLASS_ALIASES),
        "driver_flags": [c for c in DRIVER_FLAG_CANDIDATES if c in cols],
        "text_columns": [c for c in TEXT_CANDIDATES if c in cols],
        "all_columns": cols,
    }
    return schema


def load_relevant_mutations(path: Path, genes: Iterable[str]) -> tuple[pd.DataFrame, dict[str, object]]:
    schema = mutation_schema(path)
    if schema["model_id"] is None or schema["gene"] is None:
        raise ValueError(
            f"{path.name}: cannot identify ModelID/gene columns. Columns: {schema['all_columns']}"
        )
    needed = [schema["model_id"], schema["gene"]]
    for key in ("protein_change", "variant_classification"):
        if schema[key] is not None:
            needed.append(schema[key])
    needed.extend(schema["driver_flags"])
    needed.extend(schema["text_columns"])
    needed = list(dict.fromkeys(str(x) for x in needed if x is not None))
    wanted = {str(g).upper() for g in genes}

    # Mutation files can be large; chunked loading prevents unnecessary memory use.
    parts: list[pd.DataFrame] = []
    for chunk in pd.read_csv(path, usecols=needed, dtype=str, chunksize=200_000, low_memory=False):
        gene_col = str(schema["gene"])
        mask = chunk[gene_col].fillna("").astype(str).str.upper().isin(wanted)
        if mask.any():
            parts.append(chunk.loc[mask].copy())
    df = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=needed)
    rename = {str(schema["model_id"]): "ModelID", str(schema["gene"]): "GeneSymbol"}
    if schema["protein_change"] is not None:
        rename[str(schema["protein_change"])] = "ProteinChange"
    if schema["variant_classification"] is not None:
        rename[str(schema["variant_classification"])] = "VariantClassification"
    df = df.rename(columns=rename).fillna("")
    if "ProteinChange" not in df:
        df["ProteinChange"] = ""
    if "VariantClassification" not in df:
        df["VariantClassification"] = ""
    return df, schema


def row_has_driver_signal(row: pd.Series, driver_columns: Iterable[str]) -> bool:
    truthy = {"1", "true", "yes", "y", "t", "likely", "oncogenic", "hotspot"}
    for col in driver_columns:
        value = str(row.get(col, "")).strip().lower()
        if value in truthy:
            return True
        if any(token in value for token in ("oncogenic", "driver", "hotspot", "likely")) and value not in {
            "false", "no", "0", "none", "na", "nan"
        }:
            return True
    return False


def genomewide_gene_schema(path: Path) -> tuple[str, list[str], dict[str, str]]:
    """Return matrix ID column, normalized gene symbols, and raw->symbol rename map.

    DepMap CRISPR matrices encode genes as labels such as ``KRAS (3845)``.  The
    genome-wide explorer works on HGNC-style symbols while preserving the raw
    source file unchanged.
    """
    header = pd.read_csv(path, nrows=0).columns.tolist()
    if not header:
        raise ValueError(f"{path.name}: empty header")
    id_col = header[0]
    symbols: list[str] = []
    rename: dict[str, str] = {}
    seen: set[str] = set()
    duplicates: set[str] = set()
    for raw in header[1:]:
        symbol, _ = parse_gene_label(raw)
        if symbol in seen:
            duplicates.add(symbol)
        seen.add(symbol)
        symbols.append(symbol)
        rename[raw] = symbol
    if duplicates:
        raise ValueError(
            f"{path.name}: duplicate gene symbols after normalization: {sorted(duplicates)[:20]}"
        )
    return id_col, symbols, rename


def load_genomewide_matrix_subset(
    path: Path,
    model_ids: Iterable[str],
    *,
    chunksize: int = 128,
) -> pd.DataFrame:
    """Load all genes but retain only requested DepMap models.

    The public CRISPR matrices are very wide (~18.5k genes).  Reading the entire
    matrix into memory is unnecessary for a context analysis, because only a
    small set of rows (context + comparator) is needed.  This function scans the
    CSV once in row chunks and keeps only those rows.
    """
    wanted = {str(x) for x in model_ids}
    id_col, symbols, rename = genomewide_gene_schema(path)
    parts: list[pd.DataFrame] = []
    for chunk in pd.read_csv(path, chunksize=chunksize, low_memory=False):
        ids = chunk[id_col].astype(str)
        mask = ids.isin(wanted)
        if mask.any():
            parts.append(chunk.loc[mask].copy())
    if not parts:
        return pd.DataFrame(index=pd.Index([], name="ModelID"), columns=symbols, dtype=float)
    df = pd.concat(parts, ignore_index=True)
    df = df.rename(columns={id_col: "ModelID", **rename})
    df["ModelID"] = df["ModelID"].astype(str)
    if df["ModelID"].duplicated().any():
        dup = df.loc[df["ModelID"].duplicated(keep=False), "ModelID"].unique().tolist()
        raise ValueError(f"{path.name}: duplicate ModelID values in selected subset: {dup[:10]}")
    numeric = df[symbols].apply(pd.to_numeric, errors="coerce")
    numeric.insert(0, "ModelID", df["ModelID"].values)
    return numeric.set_index("ModelID")


def load_genomewide_dependency_subset_with_background(
    path: Path,
    model_ids: Iterable[str],
    *,
    dependency_threshold: float = 0.5,
    chunksize: int = 128,
) -> tuple[pd.DataFrame, pd.Series, pd.Series]:
    """Load requested dependency rows and calculate all-model dependency fractions.

    Returns ``(selected_matrix, broad_fraction, broad_n)``.  Broad fractions are
    computed from every model in the matrix without retaining the full matrix in
    memory.  They are a screening descriptor only, not a safety classification.
    """
    wanted = {str(x) for x in model_ids}
    id_col, symbols, rename = genomewide_gene_schema(path)
    dependent_counts = pd.Series(0, index=symbols, dtype="int64")
    observed_counts = pd.Series(0, index=symbols, dtype="int64")
    parts: list[pd.DataFrame] = []

    raw_gene_cols = list(rename.keys())
    for chunk in pd.read_csv(path, chunksize=chunksize, low_memory=False):
        ids = chunk[id_col].astype(str)
        numeric = chunk[raw_gene_cols].apply(pd.to_numeric, errors="coerce")
        numeric.columns = symbols
        observed_counts = observed_counts.add(numeric.notna().sum(axis=0), fill_value=0).astype("int64")
        dependent_counts = dependent_counts.add(
            (numeric > dependency_threshold).sum(axis=0), fill_value=0
        ).astype("int64")
        mask = ids.isin(wanted)
        if mask.any():
            selected = numeric.loc[mask].copy()
            selected.insert(0, "ModelID", ids.loc[mask].values)
            parts.append(selected)

    if parts:
        selected_df = pd.concat(parts, ignore_index=True)
        if selected_df["ModelID"].duplicated().any():
            dup = selected_df.loc[
                selected_df["ModelID"].duplicated(keep=False), "ModelID"
            ].unique().tolist()
            raise ValueError(f"{path.name}: duplicate ModelID values in selected subset: {dup[:10]}")
        selected_df = selected_df.set_index("ModelID")
    else:
        selected_df = pd.DataFrame(
            index=pd.Index([], name="ModelID"), columns=symbols, dtype=float
        )

    broad_fraction = dependent_counts.astype(float).div(
        observed_counts.where(observed_counts > 0).astype(float)
    )
    broad_fraction.name = "broad_dependency_fraction"
    observed_counts.name = "broad_dependency_n"
    return selected_df, broad_fraction, observed_counts
