# MCL Pharmacology Layer v1

## Цель

Связать три независимых типа данных без подмены одного другим:

1. **Клеточная модель × вещество × экспериментальный ответ** — что реально произошло с моделью после воздействия вещества.
2. **Вещество × белковая мишень × доказательство** — на какие белки вещество известно или предполагается воздействовать и на каком основании.
3. **Клеточная модель × ген × CRISPR Gene Effect** — насколько потеря функции гена влияет на рост/жизнеспособность той же модели.

Совпадение этих трёх слоёв поддерживает механистическую гипотезу, но не является доказательством причинного механизма.

## Главный принцип

Не объединять IC50, AUC, GI50, viability, PRISM LFC и другие показатели в один универсальный показатель чувствительности. Каждый источник и endpoint сохраняются в исходной семантике и единицах. Сравнение допускается только внутри совместимых экспериментов или после отдельно описанной нормализации.

## Контракт данных

### `compound`

Одна строка — одно нормализованное вещество.

Основные поля:

- `compound_id` — стабильный идентификатор MCL/источника;
- `preferred_name`;
- `canonical_smiles`;
- `inchikey`;
- `pubchem_cid`;
- `chembl_id`;
- `broad_id`;
- `gdsc_id`.

### `drug_response_observation`

Одна строка — одно экспериментальное наблюдение.

Обязательная логика:

- `model_id` должен однозначно разрешаться в модель CRISPR-Атласа;
- `compound_id` должен существовать в реестре веществ;
- обязательно сохраняются `source`, `source_release`, `assay_type`, `endpoint`, `value`, `unit`;
- дополнительные поля (`dose`, `exposure_time_h`, `AUC`, `IC50`, `GI50` и т. д.) не заменяют исходный endpoint.

### `compound_target_evidence`

Одна строка — одно доказательство связи вещества с белком.

- `target_gene`;
- `action` — ингибитор/агонист/деградер и т. п., только если источник это действительно определяет;
- `evidence_type` — binding, biochemical, curated MOA, cellular, computational и т. п.;
- `activity_type/value/unit`;
- `source`, `source_assay_id`, `publication`;
- `confidence`, `directness`.

Аннотация мишени препарата не означает, что именно эта мишень объясняет ответ конкретной клеточной модели.

## CRISPR-согласованность

После materialization MCL добавляет к паре `model × compound × annotated target` Chronos Gene Effect соответствующего гена в той же модели.

Навигационные уровни v1:

- `strong_dependency`: Gene Effect <= -1.0;
- `dependency`: -1.0 < Gene Effect <= -0.5;
- `weak_or_none`: Gene Effect > -0.5;
- `not_available`: CRISPR-значение отсутствует.

Это описательные категории MCL. Они не доказывают фармакологический механизм и не заменяют probability of dependency.

## Уровни доказательности механизма

Предлагаемая шкала для будущего интерфейса:

- **A — прямое механистическое подтверждение:** target engagement, biochemical activity и генетическая rescue/perturbation validation в релевантной системе.
- **B — сильная согласованность:** известная мишень + чувствительность модели + CRISPR-зависимость + согласованный молекулярный контекст.
- **C — ассоциативная поддержка:** корреляции с экспрессией, мутациями, CNV или pharmacogenomic profile.
- **D — вычислительная гипотеза:** target prediction, similarity, docking и т. п.

Автоматическое повышение до A без прямых экспериментальных данных запрещено.

## Источники v1

Приоритет подключения:

1. PRISM Repurposing Public 24Q2 — первая реализация, так как модели используют DepMap/ACH-идентификаторы.
2. GDSC1/GDSC2.
3. CTRP/CTD².
4. ChEMBL cell-based assays + target bioactivity.
5. PubChem BioAssay.
6. NCI-60 для пересекающихся моделей.

Каждый источник должен сначала преобразовываться в нормализованный контракт, а только затем попадать в runtime API.

## Каталоги

Исходные файлы источников:

`data/raw/pharmacology/<source>/`

Нормализованные таблицы:

`data/raw/pharmacology/normalized/<source>/`

Runtime для сайта:

`data/runtime/pharmacology/`

QC:

`outputs/qc/pharmacology_*.tsv`

API никогда не читает `data/raw` напрямую.

## PRISM 24Q2: первый маршрут

```powershell
.\.venv\Scripts\python.exe .\scripts\fetch_prism_24q2.py --list-only
.\.venv\Scripts\python.exe .\scripts\fetch_prism_24q2.py
.\.venv\Scripts\python.exe .\scripts\ingest_prism_24q2.py --output-dir .\data\raw\pharmacology\normalized\prism_24q2
.\.venv\Scripts\python.exe .\scripts\build_pharmacology_layer.py
```

Перед полной загрузкой рекомендуется сначала выполнить `--list-only` и проверить имена/размеры файлов текущего Figshare release.

## QC gates

Слой считается пригодным для Explorer только если:

- неизвестные ModelID не подставляются по сходству имени;
- неизвестные compound IDs не теряются тихо;
- unresolved rows экспортируются в QC;
- ModelID связывается с CRISPR Atlas;
- endpoint и unit сохранены;
- target evidence сохраняет источник и тип доказательства;
- CRISPR-согласованность не маркируется как causal mechanism;
- source release фиксируется.

## Следующий этап

После PRISM:

1. добавить GDSC и CTRP как отдельные source adapters;
2. построить cross-source compound registry по InChIKey/структуре;
3. добавить карточку вещества;
4. добавить в карточку гена раздел «Какие препараты уже воздействовали на эту мишень?»;
5. построить drug-response × CRISPR-profile concordance для восстановления известных механизмов и последующей валидации target inference наших собственных SMILES.
