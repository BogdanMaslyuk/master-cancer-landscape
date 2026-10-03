param(
    [string]$CancerId = "CANCER-001",
    [ValidateSet("primary", "kras-wt", "other-kras")]
    [string]$Comparison = "kras-wt"
)

$ErrorActionPreference = "Stop"

mcl sync-depmap --offline
mcl validate-depmap-inputs
mcl analyze-depmap-genome-wide --cancer-id $CancerId --comparison $Comparison
mcl qc-depmap-genome-wide --cancer-id $CancerId --comparison $Comparison

$report = Join-Path $PSScriptRoot "..\outputs\reports\depmap_explorer_${CancerId}__${Comparison}.html"
if (Test-Path $report) {
    Write-Host "Explorer: $report"
}
