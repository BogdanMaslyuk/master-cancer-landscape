Master Cancer Landscape — Milestone 2A+2B patch

IMPORTANT: apply this patch OVER your existing master-cancer-landscape folder.
Do not delete the existing folder first. Your data/raw, data/processed, outputs/qc and HGNC checkpoint must remain in place.

PowerShell after extracting/overwriting files:

  .\.venv\Scripts\Activate.ps1
  python -m pip install -e ".[dev,analysis]"
  pytest -q
  mcl status
  mcl validate-ot-diseases

Expected tests: 15 passed.
Stop after validate-ot-diseases and review the 4 mappings before running fetch-opentargets.
