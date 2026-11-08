# MCL vNext — анализ новой молекулы (MVP)

## Цель

MVP добавляет пользовательский маршрут:

`SMILES → стандартизация → Target Ligand Space → предполагаемые мишени → опухолевые контексты → клеточные модели → экспериментальная панель`.

Главный вопрос системы:

> На каких клеточных моделях наиболее информативно проверить новую молекулу при текущем объёме доказательств?

MCL не выдаёт вероятность клинической эффективности и не утверждает механизм действия до прямой экспериментальной проверки.

## Пользовательский маршрут

Открыть:

`http://127.0.0.1:3000/molecule-analysis`

Поддерживается:

- один SMILES, вставленный вручную;
- один SMILES из первого непустого ряда файла `.smi`, `.smiles` или `.txt`;
- необязательное пользовательское название.

SDF и пакетная загрузка относятся к следующему этапу.

## Стандартизация

Используется RDKit:

`Cleanup → FragmentParent → canonical isomeric SMILES → InChIKey`.

По каноническому SMILES создаётся стабильный идентификатор вида:

`MCL-UPL-XXXXXXXXXXXX`.

Повторная загрузка той же parent-структуры открывает сохранённый анализ вместо создания дубликата.

Результаты сохраняются локально в:

`data/runtime/uploaded_molecules/`.

## Физико-химические характеристики

MVP рассчитывает:

- molecular formula;
- MW и exact MW;
- cLogP;
- TPSA;
- HBD/HBA;
- rotatable bonds;
- rings;
- formal charge;
- QED;
- Lipinski violations;
- Veber;
- synthetic accessibility, если доступен RDKit SA scorer;
- PAINS;
- Brenk.

Это характеристики структуры, а не доказательство эффективности или безопасности.

## Поиск предполагаемых мишеней

Новая структура сравнивается с локальным Target Ligand Space.

Метод:

- Morgan fingerprint;
- radius = 2;
- 2048 bits;
- Tanimoto.

Для каждой текущей мишени сохраняются:

- число экспериментальных лигандов;
- лучший Tanimoto;
- ближайший известный ligand;
- известный Ki/Kd/IC50, если доступен;
- тип endpoint;
- источник;
- число соседей ≥0.35 / ≥0.50 / ≥0.70;
- top-5 структурных соседей.

Классы:

- ≥0.70 — высокая структурная переносимость;
- 0.50–0.70 — умеренная;
- 0.35–0.50 — исследовательская;
- <0.35 — механизм нельзя переносить по 2D-сходству.

Текущий Target Ligand Space ограничен проработанными core-мишенями MCL и не является proteome-wide target prediction.

## Клеточные линии

Для target-гипотез MCL использует существующий Candidate v2 / DepMap / PRISM контекст.

Клеточная модель может получить качественный приоритет:

- высокий;
- средний;
- исследовательский.

Причины показываются явно:

- CRISPR Probability of Dependency;
- Chronos Gene Effect;
- роль Candidate v2;
- наличие фармакологического контекста известных веществ;
- наличие модели в лаборатории.

Клеточная релевантность означает:

> эта модель информативна для проверки предполагаемой мишени.

Она не означает:

> новая молекула уже доказанно активна на этой линии.

## Экспериментальная панель

Target-based панель формируется только если хотя бы одна мишень достигает Tanimoto ≥0.35.

При наличии достаточной опоры система пытается выбрать:

1. до двух положительных опухолевых моделей;
2. отрицательный опухолевый контроль;
3. общий человеческий неопухолевый контроль из лабораторного реестра;
4. первый viability/MTT-пилот.

Если ни одна мишень не достигает 0.35, MCL не выдаёт target-based панель как обоснованную. В этом случае нужен независимый target prediction либо phenotype-first screening.

## Неопухолевый контроль

BJ5ta, если доступна в локальном реестре, используется только как общий человеческий неопухолевый comparator.

Нельзя трактовать BJ5ta как:

- нормальную ткань конкретного органа;
- доказательство безопасности;
- полноценную оценку терапевтического окна.

## ADMET

Универсальная ADMET-модель для произвольной загруженной структуры в MVP не подключена.

ADMET-снимки PYZ не переносятся на новую молекулу.

## Normal Tissue Liability

GTEx / Human Protein Atlas / human genetic constraint / расширенные нормальные клеточные модели в MVP ещё не подключены.

Это следующий крупный биологический слой.

## API

- `POST /api/molecules/analyze`
- `GET /api/molecules/{analysis_id}`
- `GET /api/molecules/{analysis_id}/targets`
- `GET /api/molecules/{analysis_id}/cell-lines`
- `GET /api/molecules/{analysis_id}/experiment`
- `GET /media/molecules/{analysis_id}/structure.svg`

## Научные ограничения

- structural similarity ≠ binding;
- docking ≠ binding evidence;
- CRISPR dependency ≠ pharmacological inhibition;
- PRISM response ≠ target-mediated mechanism;
- known ligand evidence ≠ uploaded molecule evidence;
- cell-line response ≠ in-vivo efficacy;
- in-silico prediction ≠ safety;
- отсутствие similarity ≥0.35 не исключает биологическую активность.

## Следующий этап

После стабилизации MVP:

1. SDF и batch upload;
2. расширение Target Ligand Space;
3. независимый target prediction;
4. Normal Tissue Liability layer;
5. интеграция docking;
6. загрузка собственных MTT/dose-response результатов;
7. IC50, 95% CI и selectivity index внутри MCL.
