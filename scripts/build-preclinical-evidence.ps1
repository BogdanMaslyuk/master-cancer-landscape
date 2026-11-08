param(
    [string]$Release = "v1.0",
    [switch]$ForceDownload,
    [string]$DownloadUrl = "",
    [string]$ArchivePath = ""
)

$ErrorActionPreference = "Stop"
$Root = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$Python = Join-Path $Root ".venv\Scripts\python.exe"

if (-not (Test-Path $Python)) {
    throw "Python environment not found: $Python"
}

Set-Location $Root

Write-Host "=== Preclinical Database source ===" -ForegroundColor Cyan
$FetchArgs = @(".\scripts\fetch_preclinical_db.py", "--release", $Release)
if ($ForceDownload) { $FetchArgs += "--force" }
if ($DownloadUrl) { $FetchArgs += @("--url", $DownloadUrl) }
if ($ArchivePath) { $FetchArgs += @("--archive", $ArchivePath) }
& $Python @FetchArgs
if ($LASTEXITCODE -ne 0) {
    throw "Preclinical Database source step failed with exit code $LASTEXITCODE"
}

Write-Host "`n=== MCL Preclinical Evidence ===" -ForegroundColor Cyan
& $Python ".\scripts\build_preclinical_oncology_evidence.py" --release $Release
if ($LASTEXITCODE -ne 0) {
    throw "Preclinical evidence build failed with exit code $LASTEXITCODE"
}

Write-Host "`nPreclinical evidence build completed." -ForegroundColor Green
Write-Host "Review data\runtime\preclinical\preclinical_target_summary.tsv and preclinical_compound_summary.tsv."
