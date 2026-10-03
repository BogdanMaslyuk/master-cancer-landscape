param(
    [string]$SnapshotDate = "2026-09-26"
)
$ErrorActionPreference = "Stop"
if (-not (Test-Path ".\.venv\Scripts\python.exe")) {
    throw "Virtual environment not found. Run scripts/setup_windows.ps1 first."
}
.\.venv\Scripts\Activate.ps1
mcl fetch-hgnc --snapshot-date $SnapshotDate
mcl normalize-targets --source-release $SnapshotDate
mcl qc-targets
