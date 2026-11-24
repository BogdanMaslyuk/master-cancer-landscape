[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [string]$P09Zip,

    [Parameter(Mandatory = $true)]
    [ValidatePattern('^[0-9a-f]{64}$')]
    [string]$P09ExpectedSHA256,

    [string]$Root = "$env:USERPROFILE\MCL_Project",
    [string]$Python = '',
    [string]$Curation = '',
    [string]$Output = '',
    [string]$OutputZip = ''
)

$ErrorActionPreference = 'Stop'
$Repo = Split-Path -Parent $PSScriptRoot
$Builder = Join-Path $PSScriptRoot 'build_primary_publication_review_v3_p010.py'

if (-not (Test-Path -LiteralPath $Builder -PathType Leaf)) {
    throw "P0.10 builder missing: $Builder"
}
if (-not (Test-Path -LiteralPath $P09Zip -PathType Leaf)) {
    throw "Validated P0.9 ZIP missing: $P09Zip"
}
if ([string]::IsNullOrWhiteSpace($Curation)) {
    $Curation = Join-Path $Repo 'docs\research\evidence\tls_v3_p010_publication_reviews.json'
}
if (-not (Test-Path -LiteralPath $Curation -PathType Leaf)) {
    throw "P0.10 curation journal missing: $Curation"
}
if ([string]::IsNullOrWhiteSpace($Python)) {
    $Python = Join-Path $Root 'repositories\master-cancer-landscape\.venv\Scripts\python.exe'
}
if (-not (Test-Path -LiteralPath $Python -PathType Leaf)) {
    throw "Python interpreter missing: $Python"
}

$ActualP09Sha256 = (Get-FileHash -LiteralPath $P09Zip -Algorithm SHA256).Hash.ToLowerInvariant()
if ($ActualP09Sha256 -ne $P09ExpectedSHA256) {
    throw "P0.9 ZIP SHA-256 mismatch: expected $P09ExpectedSHA256, actual $ActualP09Sha256"
}
$CurationSha256 = (Get-FileHash -LiteralPath $Curation -Algorithm SHA256).Hash.ToLowerInvariant()

$Commit = (& git -C $Repo rev-parse --short=8 HEAD).Trim()
if ($LASTEXITCODE -ne 0 -or [string]::IsNullOrWhiteSpace($Commit)) {
    throw 'Cannot resolve current Git commit'
}
$Stamp = [DateTime]::UtcNow.ToString('yyyyMMddTHHmmssfffZ')
$Reports = Join-Path $Root 'reports'
if ([string]::IsNullOrWhiteSpace($Output)) {
    $Output = Join-Path $Reports ("TLS_V3_P010_{0}_{1}" -f $Stamp, $Commit)
}
if ([string]::IsNullOrWhiteSpace($OutputZip)) {
    $OutputZip = "$Output.zip"
}
if (Test-Path -LiteralPath $Output) { throw "Output already exists: $Output" }
if (Test-Path -LiteralPath $OutputZip) { throw "Output ZIP already exists: $OutputZip" }

Write-Host '=== MCL TLS V3 P0.10: primary-publication review ===' -ForegroundColor Cyan
Write-Host "P0.9 ZIP: $P09Zip"
Write-Host "P0.9 SHA-256: $ActualP09Sha256"
Write-Host "Curation SHA-256: $CurationSha256"
Write-Host "Output: $Output"

& $Python $Builder `
    --p09-zip $P09Zip `
    --expected-p09-sha256 $P09ExpectedSHA256 `
    --curation $Curation `
    --output-dir $Output `
    --output-zip $OutputZip `
    --git-commit $Commit
if ($LASTEXITCODE -ne 0) {
    throw "TLS V3 P0.10 failed (exit code $LASTEXITCODE)"
}

$SummaryPath = Join-Path $Output 'p010_summary.json'
if (-not (Test-Path -LiteralPath $SummaryPath -PathType Leaf)) {
    throw 'P0.10 summary missing despite zero exit code'
}
if (-not (Test-Path -LiteralPath $OutputZip -PathType Leaf)) {
    throw 'P0.10 ZIP missing despite zero exit code'
}
$Summary = Get-Content -LiteralPath $SummaryPath -Raw | ConvertFrom-Json
if ($Summary.scientific_acceptance -ne 'NOT_APPROVED' -or $Summary.merge_permitted -ne $false) {
    throw 'P0.10 acceptance safeguards mismatch'
}
if ($Summary.multi_gigabyte_sources_opened -ne $false -or $Summary.network_used_by_builder -ne $false) {
    throw 'P0.10 source isolation safeguards mismatch'
}
$FinalP09Sha256 = (Get-FileHash -LiteralPath $P09Zip -Algorithm SHA256).Hash.ToLowerInvariant()
$FinalCurationSha256 = (Get-FileHash -LiteralPath $Curation -Algorithm SHA256).Hash.ToLowerInvariant()
if ($FinalP09Sha256 -ne $ActualP09Sha256) { throw 'P0.9 input ZIP changed during P0.10' }
if ($FinalCurationSha256 -ne $CurationSha256) { throw 'P0.10 curation journal changed during build' }

$OutputSha256 = (Get-FileHash -LiteralPath $OutputZip -Algorithm SHA256).Hash.ToLowerInvariant()
Write-Host "P0.10 ZIP: $OutputZip" -ForegroundColor Green
Write-Host "P0.10 ZIP SHA-256: $OutputSha256"
Write-Host 'scientific_acceptance=NOT_APPROVED'
Write-Host 'merge_permitted=false'
