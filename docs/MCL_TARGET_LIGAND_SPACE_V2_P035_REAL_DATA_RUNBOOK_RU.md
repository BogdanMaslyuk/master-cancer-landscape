# MCL Target Ligand Space v2 — P0.3.5 / реальный прогон

## Статус

Код сборки и валидации готов к тестированию на настоящих снимках. **Полного запуска по реальному ChEMBL/BindingDB пока не было:** официальные большие файлы, рассчитанный P0.2 universe и подтверждения экспериментальных протоколов не предоставлены этому рабочему окружению. Нельзя подменять реальные записи тестовыми.

## Официальные исходники для закрепления

- ChEMBL 37 (май 2026): https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_37/
- BindingDB TSV 202611 (опубликован 31.10.2026): https://bindingdb.org/rwd/bind/chemsearch/marvin/Download.jsp

Избегать ссылок latest. Для ChEMBL нужны архив с SQLite и официальная версия, для BindingDB — полный TSV и, по возможности, BindingDB_rsid_eaids + BindingDB_Assays. Оригинальные ZIP/tar.gz сохранять в отдельном каталоге. ChEMBL SQLite при полной распаковке требует десятки гигабайт пространства. Проверять версию, размеры, хеши и действительные условия источников.

## Архитектура

1. P0.2 protein_target_universe.tsv — только однозначные reviewed-human записи.
2. scripts/export_tls_v2_pinned_sources.py — воспроизводимый SQL-join ChEMBL SQLite: activity + assay + target + UniProt component + molecular structure + publication + BAO.
3. scripts/ingest_expanded_target_ligands.py — потоковая нормализация и сохранение исходных операторов Ki/Kd/IC50.
4. config/expanded_target_assay_review_template.tsv — пустой образец документа проверки методики (НЕ настоящие approvals).
5. scripts/build_expanded_target_index.py — RDKit точная структура; отдельно parent; исключение неоднозначных форм; проверка A/B; зеркала и возможные дубликаты; индекс Morgan r2/2048.
6. scripts/validate_expanded_target_index.py — независимое повторное вычисление химической идентичности, отпечатков, смещений и SHA-256.
7. Всё сохраняется в staging и не включается в API до P0.3.6.

## Подтверждение методики

Для каждой одобренной методики нужен заполненный TSV со столбцами:

source, source_release, assay_id, method_class, reviewer, evidence_url, verification_status

Значения method_class: binding_direct или biochemical_function. Поле verification_status: approved — только после чтения реального первичного экспериментального метода. Без утверждённых методик реальный A/B каталог автоматически не создаётся. Записи BindingDB без проверенной методики остаются в C.

## Windows: обновить и проверить код

~~~powershell
$Root = "$env:USERPROFILE\Desktop\mcl-vnext-review"
Set-Location $Root
git fetch origin
git switch feature/expanded-target-ligand-space-v2
git pull --ff-only
$Python = ".\.venv\Scripts\python.exe"
$env:PYTHONPATH = "$Root\src;$Root\apps\api"
& $Python -m pytest .\tests\test_expanded_target_index.py .\tests\test_export_tls_v2_sqlite.py -q
~~~

Папка .venv должна реально существовать в выбранном рабочем каталоге; в противном случае использовать проверенное окружение другой локальной копии MCL.

## Windows: когда настоящие архивы распакованы

Пути — образцы расположения, а не подтверждение существования файлов.

~~~powershell
$Science = "D:\MCL-scientific"
$ChemblDb = "$Science\chembl_37\chembl_37.db"
$BindingDB = "$Science\bindingdb_202611\BindingDB_All_202611.tsv"
$Registry = "$Root\data\processed\target_ligand_space_v2\protein_target_universe.tsv"
$SourceDir = "$Science\tls_v2_sources"
$Review = "$Science\assay_reviews_verified.tsv"
foreach ($p in @($ChemblDb,$BindingDB,$Registry)) { if (-not (Test-Path $p)) { throw "Missing: $p" } }
$ChemblHash = (Get-FileHash $ChemblDb -Algorithm SHA256).Hash.ToLowerInvariant()

$argsSource = @(
  ".\scripts\export_tls_v2_pinned_sources.py",
  "--chembl-sqlite",$ChemblDb,"--chembl-sqlite-sha256",$ChemblHash,
  "--chembl-release","ChEMBL_37",
  "--chembl-url","https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_37/",
  "--bindingdb-tsv",$BindingDB,
  "--bindingdb-release","202611",
  "--bindingdb-url","https://bindingdb.org/rwd/bind/chemsearch/marvin/Download.jsp",
  "--bindingdb-license","REPLACE_WITH_ACTUAL_VERIFIED_LICENSE",
  "--registry",$Registry,"--output-dir",$SourceDir
)
& $Python @argsSource

& $Python ".\scripts\ingest_expanded_target_ligands.py" --source-manifest "$SourceDir\source_manifest.json" --registry $Registry --out-dir "$Science\tls_v2_normalized"

$IngestID = "REPLACE_WITH_ACTUAL_ID_PRINTED_BY_INGEST"
$IngestStage = "$Science\tls_v2_normalized\staging\$IngestID"
if (-not (Test-Path "$IngestStage\ingest_manifest.json")) { throw "Real ingest stage missing" }
if (-not (Test-Path $Review)) { throw "Verified assay review TSV missing" }

& $Python ".\scripts\build_expanded_target_index.py" --ingest-stage $IngestStage --registry $Registry --reviewed-assays $Review --output-dir "$Science\tls_v2_index"

$IndexID = "REPLACE_WITH_ACTUAL_ID_PRINTED_BY_INDEX"
& $Python ".\scripts\validate_expanded_target_index.py" --stage "$Science\tls_v2_index\staging\$IndexID"
~~~

License placeholder is deliberately NOT a valid source citation. Check and insert the actual terms for the downloaded archive. P0.2 registry, reviewed assays, ChEMBL SQLite and BindingDB TSV are NOT bundled in GitHub.

## Обязательные ограничения

- Exact structure identity сохраняет стереохимию, соль и состав; parent connectivity служит навигации и не объединяет измерения разных веществ.
- Неоднозначная стереохимия, фрагменты, радикалы и металлы направляются на ручную проверку.
- Morgan 2048 radius 2 без хиральности используется только как совместимый с прежним MCL поиск; сходство не доказывает мишень или клеточный эффект.
- Один DOI не означает один эксперимент, а зеркальные записи не увеличивают число независимых опытов.
- Значения с отношениями < и > остаются цензурированными; Ki/Kd/IC50 сохраняются как разные конечные точки.
- Нет вывода эффективности для PYZ, нет переносов из PRISM/DepMap и нет переключения /molecule-analysis.
- Публикация настоящего релиза требует P0.3.6 с полным QC, проверкой качества исходников и откатом.

## Файлы результата P0.3.5

Путь: data/runtime/target_ligand_space_v2/staging/<build_id>/

- measurements_curated.parquet
- duplicate_clusters.parquet
- ligand_catalog.parquet
- ligand_measurement_links.parquet
- fp_index/fingerprints.bin, fp_index/entries.tsv
- target_coverage.tsv
- build_manifest.json

Последний фиксирует версии и SHA-256 каждого артефакта, количество записей и экспериментальных мишеней, используемые RDKit и Morgan. Ни одно численное значение покрытия без запуска на настоящих файлах не следует публиковать как реальное.
