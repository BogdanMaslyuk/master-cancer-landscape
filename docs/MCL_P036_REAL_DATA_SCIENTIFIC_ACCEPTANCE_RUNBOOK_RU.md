# MCL Target Ligand Space v2 — P0.3.6: запуск, аудит и научная приёмка

**Статус:** программные проверки и тесты на СИНТЕТИЧЕСКИХ данных реализованы. Полные архивы ChEMBL/BindingDB, официальная выгрузка UniProt и настоящий экспертный журнал опытов сейчас не предоставлены для обработки. Поэтому реального научного PASS пока нет. PR #3 остаётся черновым.

## Порядок без опасного обновления v1

1. P0.1: запустить аудит имеющегося научного каталога. Сверить 13 core, 65 PYZ, хеши оригинальных файлов. Отчёт сохранять вне исходного архива.
2. P0.2: собрать независимый reviewed-human UniProt-реестр по реальным зафиксированным исходникам. Проверить хеш источника и версию.
3. Сформировать pilot_panel.tsv: все 13 core + от 5 до 10 дополнительных reviewed-human белков. Дополнительные выбираются по имеющимся экспериментальным данным или по явно названному исследовательскому критерию. Их включение в пилот само по себе НЕ подтверждает наличие лиганда.
4. Закрепить официальные ChEMBL SQLite и BindingDB TSV с оригинальными сжатыми архивами, версией, хешами, ссылками и проверенными условиями источников.
5. Выгрузить ChEMBL только по пилотной панели (новый --pilot-panel в export_tls_v2_pinned_sources.py). BindingDB читать потоково и отфильтровать по UniProt (новый --pilot-panel в ingest_expanded_target_ligands.py).
6. Провести настоящую проверку описаний опытов, исходных публикаций и атрибуции A/B. Пустой шаблон config/expanded_target_assay_review_template.tsv не является одобренным экспериментом.
7. Дважды независимо построить P0.3.5 из одних и тех же закреплённых входов в двух разных выходных каталогах.
8. Независимо сверить полный индекс, контрольные суммы, точные химические структуры, связи измерений и списки подозрительных повторов.
9. Запустить P0.3.6 научный валидатор. Он требует совпадающий ingest, источники, 2 сборки, неизменность исходного v1 и настоящий ручной протокол. При ошибке статус BLOCKED, активный сайт НЕ меняется.
10. После приёмки пилота сделать **отдельный полный** прогон по всему однозначно сопоставленному реестру и представить отчёт по источникам, белкам, измерениям, отказам, затратам ресурсов и полноценный план отката. P1.1/P1.2/P1.3 только после итоговой научной приёмки.

## Новые файлы

- scripts/preflight_tls_v2_real_data.py — хеши P0.1/P0.2, 13/65, снимки источников, журнал методик, свободный диск и версии окружения
- scripts/build_tls_v2_pilot_panel.py — 13 core и 5–10 вне-core, честные статусы доступности лигандов, хеш пилота
- scripts/validate_tls_v2_scientific_pilot.py — приёмка без вмешательства в рабочий v1, отчёт по каждой мишени, повторы, 2 сборки, имитационная проверка отката
- tests/test_p036_preflight_pilot.py, tests/test_p036_scientific_validation.py — проверки контрольных сценариев
- Изменены: scripts/export_tls_v2_pinned_sources.py, scripts/ingest_expanded_target_ligands.py для точного ограничения обрабатываемых белков

## Пример PowerShell: начать с доступного локального архива

~~~powershell
$Review = "$env:USERPROFILE\Desktop\mcl-vnext-review"
Set-Location $Review
git fetch origin
git switch feature/expanded-target-ligand-space-v2
git pull --ff-only
$Python = "$Review\.venv\Scripts\python.exe"
if (-not (Test-Path $Python)) { throw "Укажите существующее Python-окружение MCL" }

$Science = "D:\MCL-scientific"
$Run = "D:\MCL-runs\tls-v2-pilot-001"
$ExistingRegistry = "$Science\data\processed\target_registry\protein_targets.parquet"
$UniProt = "D:\MCL-sources\uniprot\official_reviewed_human.tsv"

& $Python ".\scripts\audit_target_space_v2_baseline.py" --data-root $Science --output "$Run\p01\baseline_audit.json" --strict
& $Python ".\scripts\build_expanded_target_registry.py" --v1-registry $ExistingRegistry --uniprot-tsv $UniProt --uniprot-release ACTUAL_UNIPROT_RELEASE --uniprot-source-url "https://rest.uniprot.org/uniprotkb/stream" --output-dir "$Run\p02"

$Baseline = "$Run\p01\baseline_audit.json"
$Registry = "$Run\p02\protein_target_universe.tsv"
& $Python ".\scripts\build_tls_v2_pilot_panel.py" --baseline $Baseline --registry $Registry --additional-n 7 --output-dir "$Run\pilot"
$Pilot = "$Run\pilot\pilot_panel.tsv"
~~~

Все приведённые буквы дисков, имена архивов и названия RELEASE — **образцы**, не реальные пути и не утверждение о наличии файлов. Не подставлять фиктивные версии. Ссылка UniProt должна соответствовать фактически скачанному и закреплённому запросу, а не просто домену ресурса.

## Пробная загрузка официальных источников

~~~powershell
$ChemblDb = "D:\MCL-sources\chembl_37\chembl_37.db"
$BindingDbTsv = "D:\MCL-sources\bindingdb_202611\BindingDB_All_202611.tsv"
$Hash = (Get-FileHash $ChemblDb -Algorithm SHA256).Hash.ToLowerInvariant()
$SourceArgs = @(".\scripts\export_tls_v2_pinned_sources.py",
 "--chembl-sqlite",$ChemblDb,"--chembl-sqlite-sha256",$Hash,
 "--chembl-release","ChEMBL_37",
 "--chembl-url","https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_37/",
 "--bindingdb-tsv",$BindingDbTsv,"--bindingdb-release","202611",
 "--bindingdb-url","https://bindingdb.org/rwd/bind/chemsearch/marvin/Download.jsp",
 "--bindingdb-license","ACTUAL_VERIFIED_LICENSE_TERMS",
 "--registry",$Registry,"--pilot-panel",$Pilot,"--output-dir","$Run\sources")
& $Python @SourceArgs
& $Python ".\scripts\ingest_expanded_target_ligands.py" --source-manifest "$Run\sources\source_manifest.json" --registry $Registry --pilot-panel $Pilot --out-dir "$Run\normalized"
~~~

Установить имя выпуска ChEMBL и BindingDB по фактическому официальному файлу, не по дате ответа в чате. Убедиться, что версия и структура SQLite поддерживаются экспортом.

## Два каталога и проверка научных источников

Таблица assay_review_verified.tsv должна содержать реальные источник, релиз, assay ID, метод (binding_direct или biochemical_function), проверяющего, ссылку на первичный источник и approved только после настоящей проверки опыта.

~~~powershell
$IngestId = "ACTUAL_INGEST_BUILD_ID"
$IngestStage = "$Run\normalized\staging\$IngestId"
$Assays = "$Run\assay_review_verified.tsv"
& $Python ".\scripts\build_expanded_target_index.py" --ingest-stage $IngestStage --registry $Registry --reviewed-assays $Assays --output-dir "$Run\index1"
& $Python ".\scripts\build_expanded_target_index.py" --ingest-stage $IngestStage --registry $Registry --reviewed-assays $Assays --output-dir "$Run\index2"

$IndexId = "ACTUAL_INDEX_BUILD_ID"
$Index1 = "$Run\index1\staging\$IndexId"
$Index2 = "$Run\index2\staging\$IndexId"
& $Python ".\scripts\validate_expanded_target_index.py" --stage $Index1
& $Python ".\scripts\validate_expanded_target_index.py" --stage $Index2

& $Python ".\scripts\preflight_tls_v2_real_data.py" --science-root $Science --baseline $Baseline --registry $Registry --registry-manifest "$Run\p02\protein_target_universe_manifest.json" --source-manifest "$Run\sources\source_manifest.json" --assay-review $Assays --out-dir "$Run\preflight" --strict

$ManualReview = "$Run\manual_scientific_review.json"
& $Python ".\scripts\validate_tls_v2_scientific_pilot.py" --preflight "$Run\preflight\preflight_report.json" --pilot-panel $Pilot --ingest-stage $IngestStage --stage $Index1 --second-stage $Index2 --baseline $Baseline --science-root $Science --manual-review $ManualReview --report "$Run\pilot_validation_report.json" --strict
~~~

Все ACTUAL_* ID должны быть заменены выводом программы. При отсутствии файлов выполнение должно блокироваться. Нельзя указывать явно выдуманные одобрения assay и pretend scientific review.

## Контракт ручного научного заключения

Внешний файл manual_scientific_review.json должен содержать:
- build_id из реального build_manifest
- reviewer, review_date, decision=approved, sampling_method
- issues_unresolved_n=0
- reviewed_cases: >=5 реальных независимых строк при достаточном покрытии; иначе все имеющиеся
- Для каждой reviewed_case: source, source_record_id, endpoint, evidence_url первоисточника, decision=approved

Код проверяет ссылки на правильные записи, но не может самостоятельно подтвердить истинность экспертного заключения. Положительный программный статус — PASS_FOR_SCIENTIFIC_REVIEW, не самостоятельная публикация и не клиническая эффективность.

## Запреты и следующий этап

Никакие изменения data/runtime/target_ligand_space/ v1, own_compounds, phenotype, PRISM, DepMap, 65 PYZ-ID и рабочего API на P0.3.6 не допускаются. Результаты остаются в отдельных staging. Имитационный тест отката не подменяет реальный протокол активации/восстановления версии на следующем этапе. Проверку полного архива необходимо проводить отдельно после реального пилота.
