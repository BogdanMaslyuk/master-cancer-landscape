# MCL Explorer API

Read-only FastAPI layer over the existing Master Cancer Landscape outputs. The API never recomputes DepMap statistics or pathway enrichment.

## Local setup

```powershell
cd C:\Users\Bogdan\Desktop\master-cancer-landscape
.\.venv\Scripts\python.exe -m pip install -e ".\apps\api[dev]"
.\.venv\Scripts\python.exe -m uvicorn mcl_api.main:app --app-dir .\apps\api --reload --port 8000
```

Open `http://127.0.0.1:8000/docs` for interactive API documentation.

Set `MCL_ROOT` only if the API is launched outside the repository layout.

## Cell-model molecular profile index

The Cancer Atlas distinguishes three evidence levels:

1. patient tumour evidence;
2. molecular context used to stratify models;
3. molecular profile of an individual experimental model.

The raw DepMap mutation file is intentionally not read on every API request. Build a compact local index once from the already downloaded `Model.csv` and `OmicsSomaticMutations.csv`:

```powershell
cd C:\Users\Bogdan\Desktop\master-cancer-landscape
.\.venv\Scripts\python.exe .\scripts\build_depmap_model_profiles.py
```

The script infers the DepMap release from `data/processed/depmap_input_inventory.tsv`, keeps only models already present in the MCL context audit, and creates local processed files for model metadata and mutation profiles. Generated indexes are ignored by Git because they are reproducible from the raw DepMap release.

Until this index is present, the Explorer deliberately shows only the variants that were used to assign a model to an MCL molecular context. It does not pretend that those variants are the complete genetics of the cell line.

Patient cohort genomics is a separate future data layer. Statistics across DepMap cell lines must not be interpreted as mutation frequencies in patients.
