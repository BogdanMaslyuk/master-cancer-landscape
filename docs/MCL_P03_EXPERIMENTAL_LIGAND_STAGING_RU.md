# MCL Target Ligand Space v2 — P0.3: первая реализация, статус

**Это промежуточный этап P0.3.1–P0.3.2 и часть правил P0.3.4, а не законченная P0.3.**

Код размещён в feature/expanded-target-ligand-space-v2, черновой PR #3.
Не изменены v1 ligand_catalog.parquet, 65 PYZ, DepMap, PRISM, экспериментальные данные и API.

## Реализовано в коде

- scripts/ingest_expanded_target_ligands.py: контроль исходных TSV по SHA-256, версии и полям, потоковое чтение ChEMBL/BindingDB в отдельный неактивный staging Parquet.
- src/mcl/analysis/expanded_ligand_measurements.py: отдельная Ki/Kd/IC50, конверсия M/mM/µM/nM/pM, сохранение =/< /> и верхних/нижних границ, предварительные классы.
- src/mcl/analysis/expanded_ligand_dedup.py: проверяемые правила распознавания зеркал по исходному ChEMBL Activity ID, одинаковых ID, возможных публикационных повторов; одни DOI не доказывают дублирование.
- tests/test_expanded_target_ligands.py и tests/test_expanded_target_dedup.py: контрольные синтетические данные без заявления об их экспериментальном происхождении.

Текущие классы candidate_A/candidate_B — НЕ опубликованные A/B: необходима верификация структур RDKit, ассай-протокола и происхождения. BindingDB без подтверждённого метода не включается автоматически в прямое связывание; остаётся класс C. Строки ChEMBL без однозначной UniProt и молекулярного формата также не допускаются как A/B.

## Важное ограничение входа ChEMBL

Скрипт принимает **подготовленный TSV с явными соединёнными полями** activity + assay + target + target_components + molecule.
Это ещё НЕ загрузчик произвольного полного SQL-архива ChEMBL. Требуется отдельное создание воспроизводимого экспорта из закреплённого выпуска с задокументированным запросом.

Обязательные имена колонок:
activity_id, assay_chembl_id, target_chembl_id, uniprot_accession,
organism_id, target_type, confidence_score, relationship_type,
assay_type, bao_label, standard_type, standard_value, standard_units,
standard_relation, canonical_smiles, molecule_chembl_id.

BindingDB: официальный TSV с колонками BindingDB Reactant_set_id,
BindingDB MonomerID, Ligand SMILES, UniProt (SwissProt) Primary ID of Target Chain,
Target Source Organism According to Curator or DataSource, Ki (nM), Kd (nM), IC50 (nM).
Дополнительные поля первоисточника сохраняются, если присутствуют.

## Пример файла-манифеста источников

Это только пример формы, НЕ реальные скачанные релизы и НЕ валидные хеши:

~~~json
{
  "contract": "mcl-tls-v2-pinned-extracts-v1",
  "sources": [
    {
      "source": "ChEMBL",
      "release": "REPLACE_WITH_ACTUAL_RELEASE",
      "path": "your_pinned_chembl_joined.tsv",
      "sha256": "REPLACE_WITH_REAL_SHA256",
      "url": "REPLACE_WITH_OFFICIAL_SOURCE_URL",
      "license": "REPLACE_WITH_VERIFIED_TERMS",
      "format": "joined_activity_tsv"
    },
    {
      "source": "BindingDB",
      "release": "REPLACE_WITH_ACTUAL_RELEASE",
      "path": "your_pinned_bindingdb.tsv",
      "sha256": "REPLACE_WITH_REAL_SHA256",
      "url": "REPLACE_WITH_OFFICIAL_SOURCE_URL",
      "license": "REPLACE_WITH_VERIFIED_TERMS",
      "format": "bindingdb_bulk_tsv"
    }
  ]
}
~~~

Ни одно значение release/hash нельзя заполнять придуманным или текущим по умолчанию.

## Проверка на Windows

В локальной копии нужной ветки (изолировать от текущего исходного каталога):
~~~powershell
$Review = "$env:USERPROFILE\Desktop\mcl-vnext-review"
$Source = "$env:USERPROFILE\Desktop\master-cancer-landscape"
$Python = "$Source\.venv\Scripts\python.exe"
Set-Location $Review
git status --short
git fetch origin
git switch feature/expanded-target-ligand-space-v2
git pull --ff-only
$env:PYTHONPATH = "$Review\src;$Review\apps\api"
& $Python -m pytest .\tests\test_expanded_target_registry.py .\tests\test_expanded_target_ligands.py .\tests\test_expanded_target_dedup.py -q
~~~

Если существуют зафиксированные реальные входы и валидный реестр P0.2:
~~~powershell
& $Python ".\scripts\ingest_expanded_target_ligands.py" --source-manifest "C:\PATH\source_manifest.json" --registry "$Review\data\processed\target_ligand_space_v2\protein_target_universe.tsv"
~~~

Отсутствие реальных файлов / неверный SHA-256 => ошибка; не создаётся пустой «научный результат». Повторный запуск с той же версией не перезаписывает staging.

## Что остаётся до приёмки P0.3

1. Получить закреплённые официальные выпуски и построить воспроизводимый ChEMBL join с исходными идентификаторами и документированным SQL.
2. Полноценно разобрать BindingDB_Assays и Reaction Set ID; подтвердить конструкции белков и исходные публикации.
3. Реализовать RDKit exact_structure_id / parent_connectivity_id и разделение стереоизомеров и смесей.
4. Запустить контролируемую межбазовую дедупликацию на исходных данных, подготовить реестр ручных решений; текущие правила покрывают лишь проверяемые случаи.
5. Сформировать утверждённые A/B curated measurements, каталог уникальных лигандов, версионированный индекс Morgan и QC-отчёт.
6. Проверить P0.1/P0.2 на реальном локальном архиве, полные тесты GitHub и откат. До этого PR остаётся черновым, P1.1/P1.2/P1.3 не запускаются.

**Не интерпретировать предварительные каталоги как список проверенных эффектов собственных PYZ.**
