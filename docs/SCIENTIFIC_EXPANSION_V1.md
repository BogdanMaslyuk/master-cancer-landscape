# MCL Scientific Expansion v1

## Goal

Expand Master Cancer Landscape from a small set of hand-configured molecular contexts into a reproducible, scalable cancer dependency landscape.

Architecture v1 remains frozen as the platform baseline. Scientific Expansion v1 changes scientific coverage and analysis inputs, not the core API layering.

## Phase 1 — Cancer Context Registry v1

The first task is to discover molecular contexts that can support a scientifically interpretable DepMap comparison.

A context is not accepted merely because a mutation occurs in a cell line. It needs:

1. a defined disease scope (preferably an OncoTree subtype/code);
2. a defined molecular event;
3. CRISPR Gene Effect coverage for the context models;
4. mutation profiling for context/comparator assignment;
5. an explicit comparator design;
6. minimum sample-size checks;
7. manual biological review before promotion into `config/cancer_contexts.yaml`.

### Discovery modes

The discovery pipeline currently proposes two event types:

- `exact_protein_change` — recurrent protein-level event, e.g. BRAF p.V600E;
- `gene_altered` — recurrent qualifying alteration in a gene, useful for contexts where a gene-level altered state is biologically meaningful.

Qualifying variants must have at least one of the following DepMap annotations:

- Hess driver;
- hotspot;
- oncogene high impact;
- tumor suppressor high impact;
- likely loss-of-function;
- VEP HIGH impact.

This intentionally excludes the majority of unprioritized passenger variants from automatic discovery.

### Comparator designs

For an exact recurrent variant the registry can propose:

- `vs_gene_wildtype`: exact variant versus mutation-profiled models with no qualifying alteration in the same gene;
- `vs_other_gene_variants`: exact variant versus models carrying another qualifying alteration in the same gene.

For gene-level altered contexts the comparator is:

- `vs_gene_wildtype`.

The two exact-variant comparator designs answer different biological questions and must not be silently merged.

### Default thresholds

- discovery threshold: context n >= 3;
- primary context threshold: context n >= 5;
- primary comparator threshold: comparator n >= 5.

Statuses:

- `READY` — both primary thresholds are met;
- `EXPLORATORY_LOW_N` — context group is below the primary threshold;
- `EXPLORATORY_LOW_COMPARATOR_N` — comparator is below the primary threshold;
- `NO_COMPARATOR` — differential dependency cannot be claimed.

Priority tiers are transparent labels rather than a hidden composite score:

- Tier A — READY exact protein-change context;
- Tier B — READY gene-level altered context;
- Tier C — exploratory or no-comparator candidate.

## Discovery command

From the repository root:

```powershell
.\.venv\Scripts\python.exe .\scripts\discover_cancer_context_registry.py
```

The script reads the pinned local DepMap release from `data/raw/depmap/<release>/` and scans:

- `Model.csv`;
- `CRISPRGeneEffect.csv`;
- `OmicsSomaticMutations.csv`.

The mutation file is scanned once in chunks. Only CRISPR-covered models and qualifying high-confidence variants are retained in memory.

Outputs:

```text
data/processed/cancer_context_registry_candidates.tsv
data/processed/cancer_context_registry_candidates.parquet
outputs/reports/cancer_context_registry_candidates.json
outputs/reports/cancer_context_registry_v1.md
```

The generated registry is a discovery artifact. It does **not** automatically modify `config/cancer_contexts.yaml`.

## Promotion gate

Before a candidate becomes an official MCL context, review:

1. OncoTree disease scope and whether the grouped models really represent one interpretable disease entity;
2. exact molecular definition and whether gene-level aggregation is biologically justified;
3. context and comparator cell-line lists;
4. confounding by other variants in the same gene;
5. context/comparator sample size;
6. mutation-profile coverage;
7. CRISPR coverage;
8. expected direction and plausible mechanism;
9. disease mapping for patient-level evidence (Open Targets / later cBioPortal or TCGA);
10. whether the context should be primary discovery or exploratory only.

Only after this review should a context be promoted into `config/cancer_contexts.yaml` and then into the genome-wide comparison configuration.

## Next phases

After Registry v1:

1. promote the first reviewed batch of READY contexts;
2. automate genome-wide context-vs-comparator runs;
3. build the expanded Gene × Cancer matrix;
4. add co-dependency / synthetic-lethal analysis;
5. connect patient evidence and druggability layers.
