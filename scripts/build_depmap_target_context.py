from __future__ import annotations

import argparse
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


ROOT = Path(__file__).resolve().parents[1]
PROCESSED = ROOT / "data" / "processed"
RAW_DEPMAP = ROOT / "data" / "raw" / "depmap"
PHARM = ROOT / "data" / "runtime" / "pharmacology"
ATLAS = PROCESSED / "depmap_crispr_model_atlas.parquet"
TARGET_CATALOG = PHARM / "target_catalog.parquet"
GENE_EFFECT = PROCESSED / "depmap_model_gene_effect.parquet"
GENE_DEPENDENCY = PROCESSED / "depmap_model_gene_dependency.parquet"
EXPRESSION = PROCESSED / "depmap_model_expression.parquet"
COPY_NUMBER = PROCESSED / "depmap_model_copy_number.parquet"
OUTPUT = PROCESSED / "depmap_model_target_context.parquet"
MUTATION_OUTPUT = PROCESSED / "depmap_model_target_mutations.parquet"
MANIFEST = PROCESSED / "depmap_model_target_context_manifest.json"
QC = ROOT / "outputs" / "qc" / "depmap_target_context_qc.tsv"

DEPENDENCY_PROBABILITY_THRESHOLD = 0.5
MIN_CONTEXT_MODELS = 5
GENE_LABEL_RE = re.compile(r"^(.*?)\s*\((\d+)\)\s*$")

MUTATION_STRING_COLUMNS = (
    "mutation_resolution",
    "protein_changes_json",
    "dna_changes_json",
    "variant_info_json",
    "vep_impacts_json",
    "rescue_reasons_json",
)
MUTATION_INT_COLUMNS = (
    "selected_variant_count",
    "hotspot_variant_count",
    "likely_lof_variant_count",
)
MUTATION_BOOL_COLUMNS = (
    "has_hotspot",
    "has_likely_lof",
    "has_hess_driver",
    "has_oncogene_high_impact",
    "has_tsg_high_impact",
    "has_rescued_variant",
)
MUTATION_FLOAT_COLUMNS = (
    "max_variant_allele_fraction",
    "max_read_depth",
    "hotspot_matrix_call",
    "damaging_matrix_call",
)


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _text(value: Any) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    return str(value).strip()


def _bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return _text(value).lower() in {"true", "1", "yes", "y", "t"}


def _normalize_gene(label: Any) -> str:
    text = _text(label)
    match = GENE_LABEL_RE.match(text)
    return (match.group(1).strip() if match else text).upper()


def _ensure_mutation_schema(frame: pd.DataFrame) -> pd.DataFrame:
    """Keep one stable Arrow schema across mutation sources and context chunks."""
    out = frame.copy()
    for key in ("model_id", "target_gene"):
        if key not in out.columns:
            out[key] = pd.Series(dtype="string")
        out[key] = out[key].astype("string")

    for column in MUTATION_STRING_COLUMNS:
        if column not in out.columns:
            out[column] = pd.Series(pd.NA, index=out.index, dtype="string")
        else:
            out[column] = out[column].astype("string")
    for column in (
        "protein_changes_json", "dna_changes_json", "variant_info_json",
        "vep_impacts_json", "rescue_reasons_json",
    ):
        out[column] = out[column].fillna("[]")

    for column in MUTATION_INT_COLUMNS:
        if column not in out.columns:
            out[column] = pd.Series(pd.NA, index=out.index, dtype="Int64")
        else:
            out[column] = pd.to_numeric(out[column], errors="coerce").astype("Int64")

    for column in MUTATION_BOOL_COLUMNS:
        if column not in out.columns:
            out[column] = False
        out[column] = out[column].fillna(False).astype(bool)

    for column in MUTATION_FLOAT_COLUMNS:
        if column not in out.columns:
            out[column] = np.nan
        out[column] = pd.to_numeric(out[column], errors="coerce").astype("float64")
    return out


def _resolve_release(explicit: str | None) -> str:
    if explicit:
        return explicit
    manifest = PROCESSED / "depmap_model_multiomics_manifest.json"
    if manifest.exists():
        payload = json.loads(manifest.read_text(encoding="utf-8"))
        value = _text(payload.get("depmap_release"))
        if value:
            return value
    inventory = PROCESSED / "depmap_input_inventory.tsv"
    if inventory.exists():
        frame = pd.read_csv(inventory, sep="\t", nrows=1)
        if "depmap_release" in frame.columns and not frame.empty:
            value = _text(frame.iloc[0]["depmap_release"])
            if value:
                return value
    if RAW_DEPMAP.exists():
        candidates = sorted((p.name for p in RAW_DEPMAP.iterdir() if p.is_dir()), reverse=True)
        if candidates:
            return candidates[0]
    raise SystemExit("Cannot infer DepMap release. Pass --release 26Q1.")


def _matrix_genes(path: Path) -> set[str]:
    if not path.exists():
        return set()
    return set(pq.ParquetFile(path).schema.names) - {"model_id"}


def _read_long(path: Path, genes: list[str], value_name: str) -> pd.DataFrame:
    if not path.exists() or not genes:
        return pd.DataFrame(columns=["model_id", "target_gene", value_name])
    available = _matrix_genes(path)
    keep = [gene for gene in genes if gene in available]
    if not keep:
        return pd.DataFrame(columns=["model_id", "target_gene", value_name])
    frame = pd.read_parquet(path, columns=["model_id", *keep])
    out = frame.melt(id_vars="model_id", var_name="target_gene", value_name=value_name)
    out["model_id"] = out["model_id"].astype(str)
    out["target_gene"] = out["target_gene"].astype(str).str.upper()
    out[value_name] = pd.to_numeric(out[value_name], errors="coerce")
    return out


def _first_existing(release_dir: Path, names: tuple[str, ...]) -> Path | None:
    for name in names:
        path = release_dir / name
        if path.exists():
            return path
    return None


def _column(header: list[str], *names: str) -> str | None:
    lookup = {str(c).strip().lower().replace("_", ""): str(c) for c in header}
    for name in names:
        key = name.strip().lower().replace("_", "")
        if key in lookup:
            return lookup[key]
    return None


def _build_variant_level_mutations(
    path: Path,
    model_ids: set[str],
    target_genes: set[str],
    chunksize: int,
) -> pd.DataFrame:
    header = pd.read_csv(path, nrows=0).columns.tolist()
    model_col = _column(header, "ModelID", "DepMap_ID", "DepMapID", "ModelConditionID")
    gene_col = _column(header, "HugoSymbol", "Hugo_Symbol", "gene", "Gene")
    if not model_col or not gene_col:
        return _ensure_mutation_schema(pd.DataFrame())

    optional = {
        "protein_change": _column(header, "ProteinChange", "HGVSp", "VepHGVSp"),
        "dna_change": _column(header, "DNAChange", "HGVSc", "VepHGVSc"),
        "variant_info": _column(header, "VariantInfo", "Consequence"),
        "vep_impact": _column(header, "VepImpact", "IMPACT"),
        "af": _column(header, "AF", "AlleleFraction"),
        "dp": _column(header, "DP", "Depth"),
        "hotspot": _column(header, "Hotspot"),
        "likely_lof": _column(header, "LikelyLoF", "LikelyLof"),
        "hess_driver": _column(header, "HessDriver"),
        "oncogene_high_impact": _column(header, "OncogeneHighImpact"),
        "tsg_high_impact": _column(header, "TumorSuppressorHighImpact"),
        "rescue": _column(header, "Rescue"),
        "rescue_reason": _column(header, "RescueReason"),
    }
    usecols = list(dict.fromkeys([model_col, gene_col, *[c for c in optional.values() if c]]))
    hits: list[pd.DataFrame] = []
    for chunk in pd.read_csv(path, usecols=usecols, chunksize=max(10_000, chunksize), low_memory=False):
        chunk["_model_id"] = chunk[model_col].map(_text)
        chunk["_gene"] = chunk[gene_col].map(_normalize_gene)
        # A ModelConditionID is not assumed to be equivalent to an ACH ModelID. If the
        # detailed table cannot resolve directly, _build_mutations falls back to matrices.
        mask = chunk["_model_id"].isin(model_ids) & chunk["_gene"].isin(target_genes)
        if mask.any():
            hits.append(chunk.loc[mask].copy())
    if not hits:
        return _ensure_mutation_schema(pd.DataFrame())

    frame = pd.concat(hits, ignore_index=True, sort=False)
    frame["model_id"] = frame["_model_id"]
    frame["target_gene"] = frame["_gene"]
    for field in ("hotspot", "likely_lof", "hess_driver", "oncogene_high_impact", "tsg_high_impact", "rescue"):
        col = optional[field]
        frame[f"_{field}"] = frame[col].map(_bool) if col else False
    af_col = optional["af"]
    dp_col = optional["dp"]
    frame["_af"] = pd.to_numeric(frame[af_col], errors="coerce") if af_col else np.nan
    frame["_dp"] = pd.to_numeric(frame[dp_col], errors="coerce") if dp_col else np.nan

    def unique_json(group: pd.DataFrame, field: str) -> str:
        col = optional[field]
        if not col:
            return "[]"
        values: list[str] = []
        for value in group[col]:
            text = _text(value)
            if text and text not in values:
                values.append(text)
        return json.dumps(values, ensure_ascii=False)

    rows: list[dict[str, Any]] = []
    for (model_id, gene), group in frame.groupby(["model_id", "target_gene"], sort=False):
        rows.append(
            {
                "model_id": str(model_id),
                "target_gene": str(gene),
                "mutation_resolution": "variant_level",
                "selected_variant_count": int(len(group)),
                "hotspot_variant_count": int(group["_hotspot"].sum()),
                "likely_lof_variant_count": int(group["_likely_lof"].sum()),
                "has_hotspot": bool(group["_hotspot"].any()),
                "has_likely_lof": bool(group["_likely_lof"].any()),
                "has_hess_driver": bool(group["_hess_driver"].any()),
                "has_oncogene_high_impact": bool(group["_oncogene_high_impact"].any()),
                "has_tsg_high_impact": bool(group["_tsg_high_impact"].any()),
                "has_rescued_variant": bool(group["_rescue"].any()),
                "max_variant_allele_fraction": float(group["_af"].max()) if group["_af"].notna().any() else None,
                "max_read_depth": float(group["_dp"].max()) if group["_dp"].notna().any() else None,
                "protein_changes_json": unique_json(group, "protein_change"),
                "dna_changes_json": unique_json(group, "dna_change"),
                "variant_info_json": unique_json(group, "variant_info"),
                "vep_impacts_json": unique_json(group, "vep_impact"),
                "rescue_reasons_json": unique_json(group, "rescue_reason"),
            }
        )
    return _ensure_mutation_schema(pd.DataFrame(rows))


def _read_mutation_matrix(path: Path, model_ids: set[str], target_genes: set[str], value_name: str) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=["model_id", "target_gene", value_name])
    header = pd.read_csv(path, nrows=0).columns.tolist()
    id_col = header[0]
    rename: dict[str, str] = {}
    usecols = [id_col]
    for raw in header[1:]:
        gene = _normalize_gene(raw)
        if gene in target_genes and gene not in rename.values():
            rename[raw] = gene
            usecols.append(raw)
    if len(usecols) == 1:
        return pd.DataFrame(columns=["model_id", "target_gene", value_name])
    frame = pd.read_csv(path, usecols=usecols, low_memory=False).rename(columns={id_col: "model_id", **rename})
    frame["model_id"] = frame["model_id"].map(_text)
    frame = frame[frame["model_id"].isin(model_ids)]
    out = frame.melt(id_vars="model_id", var_name="target_gene", value_name=value_name)
    out[value_name] = pd.to_numeric(out[value_name], errors="coerce").fillna(0)
    return out[out[value_name] > 0].copy()


def _build_matrix_mutations(release_dir: Path, model_ids: set[str], target_genes: set[str]) -> pd.DataFrame:
    hotspot_path = _first_existing(release_dir, ("OmicsSomaticMutationsMatrixHotspot.csv",))
    damaging_path = _first_existing(release_dir, ("OmicsSomaticMutationsMatrixDamaging.csv",))
    hotspot = _read_mutation_matrix(hotspot_path, model_ids, target_genes, "hotspot_matrix_call") if hotspot_path else pd.DataFrame()
    damaging = _read_mutation_matrix(damaging_path, model_ids, target_genes, "damaging_matrix_call") if damaging_path else pd.DataFrame()
    if hotspot.empty and damaging.empty:
        return _ensure_mutation_schema(pd.DataFrame())
    keys = ["model_id", "target_gene"]
    if hotspot.empty:
        out = damaging.copy()
        out["hotspot_matrix_call"] = 0
    elif damaging.empty:
        out = hotspot.copy()
        out["damaging_matrix_call"] = 0
    else:
        out = hotspot.merge(damaging, on=keys, how="outer")
    out[["hotspot_matrix_call", "damaging_matrix_call"]] = out[["hotspot_matrix_call", "damaging_matrix_call"]].fillna(0)
    out["mutation_resolution"] = "matrix_only"
    out["has_hotspot"] = out["hotspot_matrix_call"] > 0
    out["has_likely_lof"] = out["damaging_matrix_call"] > 0
    return _ensure_mutation_schema(out)


def _build_mutations(release_dir: Path, model_ids: set[str], target_genes: set[str], chunksize: int) -> tuple[pd.DataFrame, str | None]:
    detailed = _first_existing(release_dir, ("OmicsSomaticMutations.csv", "OmicsSomaticMutationsProfile.csv"))
    if detailed:
        detailed_frame = _build_variant_level_mutations(detailed, model_ids, target_genes, chunksize)
        if not detailed_frame.empty:
            return detailed_frame, detailed.name
        print(f"Mutation detail file {detailed.name} did not resolve direct ACH model IDs; trying mutation matrices.")
    matrix = _build_matrix_mutations(release_dir, model_ids, target_genes)
    if not matrix.empty:
        return matrix, "OmicsSomaticMutationsMatrixHotspot/Damaging"
    return _ensure_mutation_schema(pd.DataFrame()), None


def _context_percentile(frame: pd.DataFrame, value_col: str, out_col: str, n_col: str) -> None:
    group_cols = ["target_gene", "mcl_cancer_id"]
    n = frame.groupby(group_cols, dropna=False)[value_col].transform(lambda s: int(s.notna().sum()))
    pct = frame.groupby(group_cols, dropna=False)[value_col].rank(method="average", pct=True)
    frame[n_col] = n.astype("int32")
    frame[out_col] = pct.where(n >= MIN_CONTEXT_MODELS)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Build a model × pharmacological-target molecular-context index from pinned DepMap omics. "
            "Mutation, expression and relative copy number remain separate evidence axes."
        )
    )
    parser.add_argument("--release", default=None)
    parser.add_argument("--gene-chunk-size", type=int, default=64)
    parser.add_argument("--mutation-chunksize", type=int, default=100_000)
    args = parser.parse_args()

    for path in (ATLAS, TARGET_CATALOG, GENE_EFFECT):
        if not path.exists():
            raise SystemExit(f"Missing required input: {path.relative_to(ROOT)}")

    release = _resolve_release(args.release)
    release_dir = RAW_DEPMAP / release
    atlas_columns = [
        c for c in ["model_id", "mcl_cancer_id", "mcl_cancer_name", "mcl_organ_ru"]
        if c in pq.ParquetFile(ATLAS).schema.names
    ]
    atlas = pd.read_parquet(ATLAS, columns=atlas_columns)
    atlas["model_id"] = atlas["model_id"].astype(str)
    atlas = atlas.drop_duplicates("model_id")
    model_ids = set(atlas["model_id"])

    targets = pd.read_parquet(TARGET_CATALOG, columns=["target_gene"])
    target_genes = sorted(set(targets["target_gene"].dropna().astype(str).str.upper().str.strip()) - {""})
    target_gene_set = set(target_genes)

    mutations, mutation_source = _build_mutations(release_dir, model_ids, target_gene_set, args.mutation_chunksize)
    mutations = _ensure_mutation_schema(mutations)
    PROCESSED.mkdir(parents=True, exist_ok=True)
    mutations.to_parquet(MUTATION_OUTPUT, index=False, compression="zstd")

    layer_paths = {
        "gene_effect": GENE_EFFECT,
        "dependency_probability": GENE_DEPENDENCY,
        "expression_log2_tpm1": EXPRESSION,
        "copy_number_relative": COPY_NUMBER,
    }
    available_layers = {name: path.exists() for name, path in layer_paths.items()}

    if OUTPUT.exists():
        OUTPUT.unlink()
    writer: pq.ParquetWriter | None = None
    rows_written = 0
    chunks_written = 0
    chunk_size = max(8, int(args.gene_chunk_size))

    for start in range(0, len(target_genes), chunk_size):
        genes = target_genes[start:start + chunk_size]
        base = pd.MultiIndex.from_product(
            [atlas["model_id"].tolist(), genes], names=["model_id", "target_gene"]
        ).to_frame(index=False)
        base = base.merge(atlas, on="model_id", how="left")
        for value_name, path in layer_paths.items():
            long = _read_long(path, genes, value_name)
            if not long.empty:
                base = base.merge(long, on=["model_id", "target_gene"], how="left")
            elif value_name not in base.columns:
                base[value_name] = np.nan

        base["dependency_call"] = pd.Series(pd.NA, index=base.index, dtype="boolean")
        has_probability = base["dependency_probability"].notna()
        base.loc[has_probability, "dependency_call"] = (
            base.loc[has_probability, "dependency_probability"] > DEPENDENCY_PROBABILITY_THRESHOLD
        ).astype("boolean")
        base["gene_effect_depletion_ge05"] = base["gene_effect"] <= -0.5
        base["gene_effect_strong_ge1"] = base["gene_effect"] <= -1.0
        base["expression_detected"] = base["expression_log2_tpm1"].fillna(0) > 0
        _context_percentile(base, "expression_log2_tpm1", "expression_percentile_in_cancer", "expression_context_models_n")
        _context_percentile(base, "copy_number_relative", "copy_number_percentile_in_cancer", "copy_number_context_models_n")

        mut_sub = mutations[mutations["target_gene"].isin(genes)].copy()
        if not mut_sub.empty:
            base = base.merge(mut_sub, on=["model_id", "target_gene"], how="left")
        base = _ensure_mutation_schema(base)
        base["depmap_release"] = release
        base["context_built_at"] = _now()

        table = pa.Table.from_pandas(base, preserve_index=False)
        if writer is None:
            writer = pq.ParquetWriter(OUTPUT, table.schema, compression="zstd")
        elif table.schema != writer.schema:
            table = table.cast(writer.schema, safe=False)
        writer.write_table(table)
        rows_written += len(base)
        chunks_written += 1
        print(f"Context genes {min(start + chunk_size, len(target_genes))}/{len(target_genes)}; rows {rows_written}")

    if writer is not None:
        writer.close()

    QC.parent.mkdir(parents=True, exist_ok=True)
    qc_rows = [
        {"metric": "models_n", "value": len(model_ids)},
        {"metric": "target_genes_n", "value": len(target_genes)},
        {"metric": "context_rows_n", "value": rows_written},
        {"metric": "mutation_rows_n", "value": len(mutations)},
        {"metric": "gene_dependency_available", "value": available_layers["dependency_probability"]},
        {"metric": "expression_available", "value": available_layers["expression_log2_tpm1"]},
        {"metric": "copy_number_available", "value": available_layers["copy_number_relative"]},
    ]
    pd.DataFrame(qc_rows).to_csv(QC, sep="\t", index=False)

    manifest = {
        "contract": "mcl-depmap-model-target-context-v1",
        "built_at": _now(),
        "depmap_release": release,
        "models_n": len(model_ids),
        "target_genes_n": len(target_genes),
        "context_rows_n": rows_written,
        "chunks_n": chunks_written,
        "dependency_probability_threshold": DEPENDENCY_PROBABILITY_THRESHOLD,
        "minimum_context_models_for_percentile": MIN_CONTEXT_MODELS,
        "available_layers": available_layers,
        "mutation_source": mutation_source,
        "outputs": [str(OUTPUT.relative_to(ROOT)), str(MUTATION_OUTPUT.relative_to(ROOT))],
        "methodology_notes_ru": [
            "Бинарная CRISPR-зависимость определяется только по Probability of Dependency > 0.5, если этот слой доступен.",
            "Chronos Gene Effect хранится как непрерывная мера силы loss-of-function фенотипа; пороги -0.5/-1 показаны описательно, но не заменяют Probability of Dependency.",
            "Экспрессия хранится как log2(TPM+1); относительный перцентиль рассчитывается внутри типа опухоли только при n>=5.",
            "OmicsCNGeneWGS трактуется как линейное относительное copy number. MCL не называет высокий/низкий перцентиль амплификацией/делецией.",
            "Hotspot и LikelyLoF сохраняются отдельно. Hotspot не интерпретируется автоматически как gain-of-function.",
            "Молекулярный контекст сам по себе не доказывает причинный механизм ответа на препарат.",
        ],
    }
    MANIFEST.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")

    print("MCL DepMap target molecular context v1")
    print(f"DepMap release: {release}")
    print(f"Models: {len(model_ids)}")
    print(f"Targets: {len(target_genes)}")
    print(f"Context rows: {rows_written}")
    print(f"Dependency probability: {'available' if available_layers['dependency_probability'] else 'MISSING'}")
    print(f"Expression: {'available' if available_layers['expression_log2_tpm1'] else 'MISSING'}")
    print(f"Relative copy number: {'available' if available_layers['copy_number_relative'] else 'MISSING'}")
    print(f"Mutation source: {mutation_source or 'MISSING'}")
    print(f"Wrote {OUTPUT.relative_to(ROOT)}")
    print(f"Wrote {MUTATION_OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
