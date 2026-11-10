# Scientific Core v2 — первый исправленный маршрут по результатам CDK4 аудита

## Подтверждённые факты на локальном наборе (2026-11-10)

- Для MCF-7 Candidate v2 содержит 4 положительные и 2 противоречивые записи по разным известным веществам; старый фильтр priority_for_in_vitro оставлял 1 противоречивую и исключал 4 положительные.
- Для T-47D: 2 положительные и 3 противоречивые; старый фильтр оставлял 1 противоречивую и исключал 2 положительные.
- Эти исходные роли относятся к реакции на конкретные известные вещества PRISM 24Q2 и не доказывают чувствительность к PYZ.
- Все шесть исследованных линий сопоставлены с DepMap ModelID и генетическими измерениями CDK4.

## Что изменено

1. Модуль full CRISPR context читает полный processed DepMap model × target слой без отсечения по Candidate v2 и без Top-5 / Top-30.
2. CRISPR зависимость, ее отсутствие и пропуск явно разделены; не выводится ложная активность PYZ.
3. Наличие культуры в реестре отдельно от подтверждённого физического состояния, которое не проверено.
4. Старые Candidate v2 ответы не склеиваются для новой PYZ: сохранены идентификаторы референсных веществ, LFC и типы ролей отдельно.
5. Фильтр priority_for_in_vitro больше не удаляет в legacy PYZ-маршруте поддержанные и исследовательские гипотезы известных препаратов.
6. Прежний сокращённый слой остаётся для совместимости, но явно помечен как reference-only и не может выдавать полностью готовую панель.
7. Интерфейс различает генетическую зависимость и ответ на известный препарат.

## Граница этого первого исправления

- Это **не окончательное завершение P0.1–P0.3**. Полное соединение исходных нормализованных PRISM/ChEMBL/BindingDB строк и их первоисточников с независимо измеренными генетическими данными ещё предстоит.
- Научная точность поиска мишеней Tanimoto 0,35 не откалибрована.
- Полный контекст DepMap не является полным массивом всех PRISM-исследований.
- Отсутствие измеренного P(dep) не трактуется как отрицательное измерение.
- Новый полный маршрут не даёт статус лабораторной или фармакологической готовности без проверки.
- Для разработки используются существующие данные DepMap 26Q1 и PRISM 24Q2; новый релиз самопроизвольно не подменяется.

## Безопасный запуск на Windows

**Важно:** виртуальное окружение находится в исходной папке, где может быть установлена старая версия пакета `mcl_api`. Поэтому `PYTHONPATH` необходимо устанавливать **до pytest**, иначе тесты могут запускаться против устаревшего кода и давать ложные регрессии. Рекомендуемый безопасный запуск из тестового каталога после `git pull --ff-only`:

~~~powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\scripts\validate-scientific-core-v2.ps1
~~~

Скрипт перед тестами проверяет фактический путь `mcl_api.__file__`, выполняет тесты и проверяет CDK4/PYZ-035 **без перезаписи карточки молекулы**. Научные Parquet-файлы читаются только из исходной папки.

Старый ручной блок ниже сохранён для диагностики; перед его выполнением `PYTHONPATH` должен быть установлен заранее.


В новой PowerShell-сессии:

~~~powershell
$Review = "$env:USERPROFILE\Desktop\mcl-vnext-review"
$Source = "$env:USERPROFILE\Desktop\master-cancer-landscape"
$Python = "$Source\.venv\Scripts\python.exe"

Set-Location $Review
git pull --ff-only
if ($LASTEXITCODE -ne 0) { throw "Git update failed" }

& $Python -m pytest .\apps\api\tests\test_uploaded_molecules.py .\apps\api\tests\test_scientific_core_v2_cdk4_audit.py -q
if ($LASTEXITCODE -ne 0) { throw "Unit tests failed" }

$env:PYTHONPATH = "$Review\apps\api"
$env:MCL_SCIENTIFIC_DATA_ROOT = $Source

@'
from pathlib import Path
from mcl_api.own_compounds import OwnCompoundStore
from mcl_api.uploaded_molecules import UploadedMoleculeStore

root = Path.cwd()
own = OwnCompoundStore(root)
store = UploadedMoleculeStore(root, own)

full = own.full_target_models("CDK4")
assert full["available"], "Full processed CDK4 context unavailable"
print("=== COMPLETE DEPMAP CDK4 ===")
print("All model rows:", full["total_model_rows"])
print("Measured P(dep):", full["total_measured_models"])
print("Registered lab models:", full["laboratory_registered_models_n"])
for name in ("MCF-7", "T-47D", "HCT 116", "HuH-7", "JURKAT", "U-937"):
    hits = [r for r in full["items"] if r["cell_line_name"].casefold() == name.casefold()]
    for r in hits:
        print(name, r["model_id"], "P(dep)=", r["dependency_probability"], r["genetic_dependency_status"])

print("=== SOURCE-SPECIFIC CANDIDATE V2 CONTEXT ===")
contexts = own._contexts_for_target("CDK4")
for name in ("MCF-7", "T-47D"):
    seen = []
    for context in contexts:
        for model in context.get("models", []):
            if str(model.get("cell_line_name") or "").casefold() == name.casefold():
                for obs in model.get("reference_drug_observations", []):
                    seen.append((obs.get("preferred_name"), obs.get("reference_model_role"), obs.get("prism_lfc")))
    print(name, sorted(set(seen), key=str))

print("=== NEW PYZ-035 ROUTE ===")
registry = own.registry()
row = registry.loc[registry["own_compound_id"].eq("PYZ-035")].iloc[0]
result = store.analyze(str(row["standardized_smiles"]), "PYZ-035", force_recompute=True)
panel = result["experiment_panel"]
print("Source:", result["cell_line_recommendations"][0].get("source_layer")
      if result["cell_line_recommendations"] else "not found")
print("Cell models:", len(result["cell_line_recommendations"]))
print("Panel:", panel.get("status"), panel.get("panel_type"))
print("Genetic candidates:", panel.get("genetic_candidate_models_n"))
print("Genetic comparators:", panel.get("genetic_comparator_models_n"))
print("Laboratory ready:", panel.get("laboratory_ready"))
print("Missing:", panel.get("missing_components"))
'@ | & $Python -

$Exit = $LASTEXITCODE
Remove-Item Env:\PYTHONPATH -ErrorAction SilentlyContinue
Remove-Item Env:\MCL_SCIENTIFIC_DATA_ROOT -ErrorAction SilentlyContinue
if ($Exit -ne 0) { throw "Real data scientific check failed" }
~~~

Файлы в научной исходной папке при этом читаются без записи. В тестовой папке создаётся или заменяется только производная карточка анализа PYZ-035; для чтения исходных структур и данных использовать test worktree как кодовую базу.

## Что требуется для приёмки первого исправления

- В полном контексте должны присутствовать шесть известных ModelID с фактическими P(dep).
- Число полных строк по CDK4 должно быть проверено относительно обработанного исходного слоя; нельзя ожидать, что оно равно старым 212.
- В reference-only слое MCF-7 и T-47D должны показывать запись ribociclib отдельно от ON123300.
- Новая PYZ-035 карточка не должна содержать измеренную чувствительность PYZ в отсутствие такого эксперимента.
- Панель не должна быть помечена laboratory_ready или complete только по аналогии с известным веществом.
- Все автоматические тесты и сборка сайта должны оставаться успешными.
