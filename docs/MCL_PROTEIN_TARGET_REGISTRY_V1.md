# MCL Protein Target Registry v1

## Зачем нужен отдельный реестр

Фармакологические источники часто указывают мишень как **символ гена** (`HTR2A`, `EGFR`, `BRAF`), хотя в интерфейсе исследователь ожидает видеть **белок**.

MCL теперь разделяет эти сущности:

```text
вещество
  ↓
фармакологическая target_gene-аннотация
  ↓
кодирующий ген
  ↓
справочное gene → protein сопоставление
  ↓
reviewed-белок UniProtKB/Swiss-Prot
```

Пример:

```text
target_gene: HTR2A
protein: 5-hydroxytryptamine receptor 2A
UniProt: P28223
source_resolution: gene_mapped
```

Важно: если PRISM/Repurposing Hub указал только `HTR2A`, MCL **не повышает** исходное доказательство до protein-resolved или isoform-resolved. Название белка и UniProt ID являются отдельным справочным слоем.

## Источники

1. `data/processed/gene_explorer/gene_reference.parquet` — локальный snapshot MCL с сопоставлением gene symbol → Swiss-Prot accession, построенный через MyGene.info.
2. UniProtKB reviewed / Swiss-Prot — официальное название белка, accession, семейство и длина.

## Выход

```text
data/processed/target_registry/protein_targets.parquet
```

Основные поля:

- `target_gene` — символ кодирующего гена и исходный ключ фармакологической аннотации;
- `source_resolution` — для текущего PRISM слоя `gene_mapped`;
- `protein_mapping_status`:
  - `unique_swissprot` — однозначная reviewed-запись;
  - `multiple_swissprot` — несколько reviewed-записей, одна белковая сущность автоматически не назначается;
  - `no_swissprot_mapping` — в локальном gene reference нет Swiss-Prot ID;
  - `uniprot_unresolved` — accession был, но не разрешился в текущем запросе UniProt;
- `protein_preferred_name` — официальное рекомендуемое название белка, только при однозначном сопоставлении;
- `uniprot_primary_accession`;
- `uniprot_accessions_json`;
- `protein_families`;
- `protein_length`;
- provenance-поля.

QC:

```text
outputs/qc/protein_target_registry_unresolved.tsv
```

## Построение

Сначала должен существовать Gene Reference Snapshot:

```powershell
.\.venv\Scripts\python.exe .\scripts\build_gene_reference_snapshot.py
```

После построения фармакологического каталога:

```powershell
.\.venv\Scripts\python.exe .\scripts\build_protein_target_registry.py
```

Для быстрой проверки можно ограничить число мишеней:

```powershell
.\.venv\Scripts\python.exe .\scripts\build_protein_target_registry.py --limit 50
```

После полного построения backend нужно перезапустить, потому что реестр кэшируется как компактный read-only reference layer.

## Интерпретация

CRISPR остаётся геновым измерением:

```text
5-HT2A receptor
  ↓ кодируется
HTR2A
  ↓ CRISPR knockout
Chronos Gene Effect
```

Поэтому корректная формулировка в MCL:

> «Вещество аннотировано к мишени через ген HTR2A; MCL сопоставляет HTR2A с reviewed-белком 5-hydroxytryptamine receptor 2A. CRISPR-согласованность рассчитывается по кодирующему гену HTR2A.»

Некорректно автоматически утверждать:

> «PRISM доказал действие вещества на конкретную изоформу 5-HT2A-рецептора».

Для `protein_resolved` и `isoform_resolved` потребуются отдельные источники с явной белковой/изоформной идентификацией.
