# MCL E1-P0-B2: one central project, no manual unpacking.
param(
  [string]$Root = (Join-Path $env:USERPROFILE "MCL_Project"),
  [string]$BindingDBPath = "",
  [string]$NeighborsPath = ""
)
$ErrorActionPreference = "Stop"
$Repo = Join-Path $Root "repositories\master-cancer-landscape"
$Python = Join-Path $Repo ".venv\Scripts\python.exe"
$Audit = Join-Path $PSScriptRoot "audit_pyz_e1_p0b_full_bindingdb.py"
$Reports = Join-Path $Root "reports"
$DataDir = Join-Path $Root "data\raw\bindingdb"
foreach ($Required in @($Root, $Repo, $Python, $Audit, $Reports)) {
    if (-not (Test-Path -LiteralPath $Required)) { throw "Missing required path: $Required" }
}
Write-Host "MCL / E1-P0-B2: centralized BindingDB audit" -ForegroundColor Cyan
Write-Host "Project: $Root"
if (-not $NeighborsPath) {
    $Found = @(Get-ChildItem -LiteralPath $Reports -Recurse -File -Filter "neighbors_audit.tsv" -ErrorAction Stop |
        Sort-Object LastWriteTime -Descending)
    if ($Found.Count -eq 0) { throw "No E1-P0-B neighbors_audit.tsv in MCL_Project\reports. Restore the original report." }
    $NeighborsPath = $Found[0].FullName
}
if (-not (Test-Path -LiteralPath $NeighborsPath -PathType Leaf)) { throw "Missing neighbors: $NeighborsPath" }
if (-not $BindingDBPath) {
    $Folders = @($DataDir, (Join-Path $env:USERPROFILE "Downloads"), (Join-Path $Root "documents")) |
        Where-Object { Test-Path -LiteralPath $_ }
    $Candidates = @()
    foreach ($Folder in $Folders) {
        $Candidates += @(Get-ChildItem -LiteralPath $Folder -Recurse -File -ErrorAction SilentlyContinue |
            Where-Object { $_.Name -match '^BindingDB_All.*\.(zip|tsv)$' })
    }
    $Candidates = @($Candidates | Sort-Object FullName -Unique)
    if ($Candidates.Count -eq 0) {
        throw "No full BindingDB_All*.zip or BindingDB_All*.tsv in Downloads/MCL_Project. Pass -BindingDBPath."
    }
    if ($Candidates.Count -gt 1) {
        $Pinned = @($Candidates | Where-Object { $_.DirectoryName -eq $DataDir })
        if ($Pinned.Count -eq 1) {
            $BindingDBPath = $Pinned[0].FullName
        } else {
            $Choices = ($Candidates | ForEach-Object { $_.FullName }) -join [Environment]::NewLine
            throw ("Several BindingDB versions found. Specify -BindingDBPath: " + [Environment]::NewLine + $Choices)
        }
    } else { $BindingDBPath = $Candidates[0].FullName }
}
if (-not (Test-Path -LiteralPath $BindingDBPath -PathType Leaf)) { throw "Missing BindingDB: $BindingDBPath" }
$BindingDBPath = (Resolve-Path -LiteralPath $BindingDBPath).Path
$SourceHash = (Get-FileHash -LiteralPath $BindingDBPath -Algorithm SHA256).Hash
if (-not (Test-Path -LiteralPath $DataDir)) {
    New-Item -ItemType Directory -Path $DataDir -Force | Out-Null
}
$Canonical = Join-Path $DataDir (Split-Path -Path $BindingDBPath -Leaf)
if ([IO.Path]::GetFullPath($Canonical) -ne [IO.Path]::GetFullPath($BindingDBPath)) {
    if (Test-Path -LiteralPath $Canonical) {
        if ((Get-FileHash -LiteralPath $Canonical -Algorithm SHA256).Hash -ne $SourceHash) {
            throw "Different source with same filename already pinned in MCL_Project. No overwrite."
        }
    } else {
        Write-Host "Copying BindingDB into MCL_Project; the original remains untouched..."
        Copy-Item -LiteralPath $BindingDBPath -Destination $Canonical -ErrorAction Stop
        if ((Get-FileHash -LiteralPath $Canonical -Algorithm SHA256).Hash -ne $SourceHash) {
            throw "Copied BindingDB SHA mismatch; audit blocked."
        }
    }
}
Write-Host "Neighbors: $NeighborsPath"
Write-Host "BindingDB: $Canonical"
Write-Host "BindingDB SHA-256: $SourceHash"
& $Python -c "import rdkit; print('RDKit:', rdkit.__version__)"
if ($LASTEXITCODE -ne 0) { throw "RDKit is missing in MCL Python environment" }
$Stamp = Get-Date -Format "yyyyMMdd_HHmmss_fff"
$Output = Join-Path $Reports ("E1_P0B2_FULL_" + $Stamp)
& $Python $Audit --neighbors $NeighborsPath --bindingdb $Canonical --out-dir $Output --require-ten
if ($LASTEXITCODE -ne 0) { throw "BindingDB scan failed; scientific sources were not changed" }
$ReportPath = Join-Path $Output "bindingdb_full_report.json"
$Report = Get-Content -LiteralPath $ReportPath -Raw -Encoding UTF8 | ConvertFrom-Json
if ($Report.scientific_status -ne "NO_PRIMARY_ASSAYS_APPROVED" -or $Report.neighbors -ne 10) {
    throw "Scientific safeguards failed"
}
$Manifest = [ordered]@{
    stage = "E1_P0B2"
    project_root = $Root
    script_path = $Audit
    git_commit = (& git -C $PSScriptRoot rev-parse HEAD)
    neighbors_path = $NeighborsPath
    neighbors_sha256 = (Get-FileHash -LiteralPath $NeighborsPath -Algorithm SHA256).Hash
    source_path = $BindingDBPath
    canonical_bindingdb_path = $Canonical
    bindingdb_sha256 = $SourceHash
    scientific_status = $Report.scientific_status
}
$ManifestPath = Join-Path $Output "run_manifest.json"
$Manifest | ConvertTo-Json -Depth 5 | Set-Content -LiteralPath $ManifestPath -Encoding UTF8
$ResultZip = Join-Path $Reports ("MCL_E1_P0B2_RESULT_" + $Stamp + ".zip")
Compress-Archive -LiteralPath @(
    $ReportPath,
    (Join-Path $Output "bindingdb_full_matches.tsv"),
    (Join-Path $Output "bindingdb_full_neighbor_coverage.tsv"),
    $ManifestPath
) -DestinationPath $ResultZip -CompressionLevel Fastest
Write-Host "Rows scanned: $($Report.stats.rows_scanned)" -ForegroundColor Green
Write-Host "Exact neighbors: $($Report.neighbors_exact_endpoint_coverage)"
Write-Host "Connectivity-only neighbors: $($Report.neighbors_connectivity_only_coverage)"
Write-Host "Primary assays approved: 0"
if (-not (Test-Path -LiteralPath $ResultZip -PathType Leaf)) {
    throw "ZIP result archive was not created"
}
Write-Host "Result: $ResultZip" -ForegroundColor Green
Write-Host "Archive saved. No Explorer launch is required." -ForegroundColor Green
