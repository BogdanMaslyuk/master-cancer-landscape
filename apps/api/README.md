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
