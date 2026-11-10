# MCL Scientific Core v2 — первый безопасный аудит CDK4

**Этап ТЗ 1.1:** P0.0 и диагностическая часть P0.1–P0.3.  
**Статус:** инструмент подготовлен; фактический локальный результат ещё не подтверждён.  
**Назначение:** проверить исходные нормализованные измерения шести лабораторных линий: HCT 116, MCF-7, T-47D, HuH-7, JURKAT и U-937.

## Почему это отдельный аудит

В коде Candidate v2 присутствует ограничение role_rank <= 5, в Own Compounds — model_items[:30], в Uploaded Molecules — selected[:5] по мишеням. Нельзя считать выдачу Candidate v2 полным набором CRISPR и PRISM.

Для CDK4 требуется:
1. сопоставить названия лабораторных культур с DepMap ModelID;
2. получить генетическое измерение ModelID × CDK4 отдельно от препаратов;
3. получить все нормализованные PRISM LFC-наблюдения по каждому известному соединению, связанному с CDK4, **с сохранением observation_id и source_assay_id**;
4. получить строки Candidate v2 с исходной ролью, compound_id и hypothesis_id;
5. найти случаи расхождения ролей разных известных веществ.

## Безопасность

Скрипт scripts/audit_scientific_core_v2_cdk4.py только читает научные Parquet-таблицы из аргумента --data-root и пишет новые отчёты в отдельный каталог --output-root. Он не пересобирает источники, не меняет Git, не перезаписывает Parquet и не удаляет данные. Для записей физической доступности не подставляются выдуманные утверждения о наличии культуры.

Рабочая копия кода: C:\Users\Bogdan\Desktop\mcl-vnext-review.  
Исходные научные данные: C:\Users\Bogdan\Desktop\master-cancer-landscape.

Нельзя выполнять git reset --hard, git clean -fd и массовые перезаписи научных файлов исходного каталога.

## Запуск в новом PowerShell

~~~powershell
$Review = "$env:USERPROFILE\Desktop\mcl-vnext-review"
$Source = "$env:USERPROFILE\Desktop\master-cancer-landscape"
$Python = "$Source\.venv\Scripts\python.exe"

Set-Location $Review
git pull --ff-only
if ($LASTEXITCODE -ne 0) { throw "Git pull failed" }

& $Python -m pytest .\apps\api\tests\test_scientific_core_v2_cdk4_audit.py -q
if ($LASTEXITCODE -ne 0) { throw "Audit unit tests failed" }

& $Python .\scripts\audit_scientific_core_v2_cdk4.py --data-root $Source
if ($LASTEXITCODE -ne 0) { throw "Scientific audit failed" }
~~~

## Выходные файлы

В тестовом каталоге создаётся новая папка outputs/audits/scientific_core_v2/CDK4_YYYYMMDDTHHMMSSZ.

| Файл | Смысл |
|---|---|
| models.tsv | Шесть названий, ModelID, контекст опухоли, CRISPR P(dep) и Gene Effect |
| prism_measurements.tsv | Отдельные измерения ответа на известные препараты с assay ID и observation ID |
| target_annotations.tsv | Исходные аннотации известных препаратов по CDK4 |
| candidate_v2_rows.tsv | Сохранённые (возможно усечённые) роли Candidate v2 |
| candidate_role_trace.tsv | Роль каждого конкретного известного препарата, статус его гипотезы и видимость через старый фильтр |
| legacy_priority_filter_impact.tsv | Для каждой линии: сколько положительных и противоречивых записей скрывает фильтр priority_for_in_vitro |
| target_links.tsv | Связки известных веществ, CDK4 и моделей |
| issues.tsv | Пропуски, дубликаты, неоднозначности |
| input_inventory.tsv | SHA256, размеры и наличие входных файлов |
| manifest.json | Версия аудита, коммит, дата, количества и источники |

**Научное ограничение:** нормализованный PRISM — не оригинальный файл выгрузки. До сверки с сырыми данными источника этот аудит остаётся диагностическим, а не окончательным подтверждением всех измерений.

## Правильная интерпретация

- CRISPR P(dep) > 0.5 означает вероятность генетической зависимости клетки от CDK4 по конкретному релизу DepMap; это не вероятность активности PYZ.
- PRISM LFC относится только к названному известному веществу и клеточной модели.
- Candidate v2 — ограниченный набор производных ролей; отсутствие там модели не означает отсутствие её в DepMap.
- legacy_priority_filter_impact.tsv объясняет, какие положительные строки Candidate v2 скрылись из-за статуса гипотезы; это не оценка активности PYZ и не полный аудит всех PRISM-экспериментов.
- Поле physical_stock_verified=false означает, что отдельная проверка качества и фактического наличия культуры не проводилась в этом инструменте.
- Наличие target-link не доказывает прямое связывание: необходимо проверить тип соответствующих target-evidence.

## Условие завершения

Подтвердить имена/ModelID; получить P(dep) и Gene Effect с релизом; сохранить строки PRISM отдельно по веществам; сравнить с ролями Candidate v2; сформировать список ошибок/пропусков и контрольные суммы. Только после анализа реального отчёта начинать исправление P0.1–P0.3. 

Отчёты содержат исследовательские сведения: не публиковать их без разрешения владельца.
