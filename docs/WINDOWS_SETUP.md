# Windows / PyCharm setup

Recommended Python: **3.12 or 3.13**. Do not start this project on Python 3.14 until all scientific dependencies are verified.

## PowerShell

```powershell
cd master-cancer-landscape
py -3.12 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
pip install -e ".[dev,analysis]"

mcl status
mcl fetch-hgnc --snapshot-date 2026-09-26
mcl normalize-targets --source-release 2026-09-26
mcl qc-targets
```

## Expected Milestone 1 outputs

```text
data/processed/target_identifiers.parquet
data/processed/target_identifiers.tsv
data/processed/provenance.parquet
data/processed/provenance.tsv
outputs/qc/target_normalization_qc.json
outputs/qc/target_normalization_qc.tsv
outputs/reports/run_manifest.json
```

The raw HGNC files are saved under `data/raw/hgnc/<retrieval-date>/` and must not be edited.
If the same snapshot already exists, the downloader stops instead of overwriting it.

## PyCharm

1. Open the `master-cancer-landscape` folder as a project.
2. Select `.venv/Scripts/python.exe` as the interpreter.
3. Keep `src/` as the source-root through the editable installation (`pip install -e`).
4. Run tests from the project root with `pytest -q`.
