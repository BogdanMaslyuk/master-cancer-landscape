# MCL local development workflow

## PowerShell roles

Use the same three-window model:

- PowerShell #1: backend / Python / FastAPI
- PowerShell #2: frontend / Next.js
- PowerShell #3: diagnostics / Git / one-off checks

## Update code

PowerShell #3:

```powershell
cd C:\Users\Bogdan\Desktop\master-cancer-landscape
git switch mcl-explorer-v0.1
git pull
```

## Rebuild Explorer indexes

PowerShell #3:

```powershell
.\scripts\build-explorer.ps1
```

Optional:

```powershell
.\scripts\build-explorer.ps1 -RefreshReference
.\scripts\build-explorer.ps1 -AllowPartialMultiomics
```

## Verify

PowerShell #3:

```powershell
.\scripts\verify.ps1
```

Fast verification without the final Next production build:

```powershell
.\scripts\verify.ps1 -SkipFrontendBuild
```

Do not call the branch clean/stable until full verification finishes successfully.

## Start Explorer

PowerShell #1:

```powershell
.\scripts\start-backend.ps1
```

Leave the window open.

PowerShell #2:

```powershell
.\scripts\start-frontend.ps1
```

The frontend script builds automatically if no production build exists. Leave the window open.

Open `http://127.0.0.1:3000`.

## Diagnostics

PowerShell #3:

```powershell
Invoke-RestMethod http://127.0.0.1:8000/health
```

```powershell
Measure-Command { Invoke-RestMethod "http://127.0.0.1:8000/api/genes/search?page_size=3" }
```

```powershell
Measure-Command { Invoke-RestMethod "http://127.0.0.1:8000/api/gene-matrix?limit=3" }
```

## Build-time / runtime rule

If a page request takes tens of seconds or minutes, do not solve it with a larger timeout. Check whether scientific computation accidentally moved into the HTTP request path.

Expensive ontology projection, genome-wide statistics, enrichment and source downloads belong in build scripts. The API serves materialized results.

## Git rule

Development continues on `mcl-explorer-v0.1`. Do not merge into `main` automatically.
