[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)]
    [string]$P08Zip,

    [Parameter(Mandatory = $true)]
    [ValidatePattern('^[0-9a-f]{64}$')]
    [string]$P08ExpectedSHA256,

    [string]$Root = "$env:USERPROFILE\MCL_Project",
    [string]$Python = '',
    [string]$Output = '',
    [string]$OutputZip = ''
)

$ErrorActionPreference = 'Stop'
$Repo = Split-Path -Parent $PSScriptRoot
$Builder = Join-Path $PSScriptRoot 'build_known_ligand_catalog_v3_p09.py'

if (-not (Test-Path -LiteralPath $Builder -PathType Leaf)) {
    throw "P0.9 builder missing: $Builder"
}
if (-not (Test-Path -LiteralPath $P08Zip -PathType Leaf)) {
    throw "Validated P0.8 ZIP missing: $P08Zip"
}
if ([string]::IsNullOrWhiteSpace($Python)) {
    $Python = Join-Path $Root 'repositories\master-cancer-landscape\.venv\Scripts\python.exe'
}
if (-not (Test-Path -LiteralPath $Python -PathType Leaf)) {
    throw "Python interpreter missing: $Python"
}

$ActualP08Sha256 = (Get-FileHash -LiteralPath $P08Zip -Algorithm SHA256).Hash.ToLowerInvariant()
if ($ActualP08Sha256 -ne $P08ExpectedSHA256) {
    throw "P0.8 ZIP SHA-256 mismatch: expected $P08ExpectedSHA256, actual $ActualP08Sha256"
}

$Commit = (& git -C $Repo rev-parse --short=8 HEAD).Trim()
if ($LASTEXITCODE -ne 0 -or [string]::IsNullOrWhiteSpace($Commit)) {
    throw 'Cannot resolve current Git commit'
}
$Stamp = [DateTime]::UtcNow.ToString('yyyyMMddTHHmmssfffZ')
$Reports = Join-Path $Root 'reports'
if ([string]::IsNullOrWhiteSpace($Output)) {
    $Output = Join-Path $Reports ("TLS_V3_P09_{0}_{1}" -f $Stamp, $Commit)
}
if ([string]::IsNullOrWhiteSpace($OutputZip)) {
    $OutputZip = "$Output.zip"
}
if (Test-Path -LiteralPath $Output) { throw "Output already exists: $Output" }
if (Test-Path -LiteralPath $OutputZip) { throw "Output ZIP already exists: $OutputZip" }

Write-Host '=== MCL TLS V3 P0.9: validated P0.8 known-ligand catalog ===' -ForegroundColor Cyan
Write-Host "P0.8 ZIP: $P08Zip"
Write-Host "P0.8 SHA-256: $ActualP08Sha256"
Write-Host "Output: $Output"

& $Python $Builder `
    --p08-zip $P08Zip `
    --expected-p08-sha256 $P08ExpectedSHA256 `
    --output-dir $Output `
    --output-zip $OutputZip `
    --git-commit $Commit
if ($LASTEXITCODE -ne 0) {
    throw "TLS V3 P0.9 failed (exit code $LASTEXITCODE)"
}

$SummaryPath = Join-Path $Output 'p09_summary.json'
if (-not (Test-Path -LiteralPath $SummaryPath -PathType Leaf)) {
    throw 'P0.9 summary missing despite zero exit code'
}
if (-not (Test-Path -LiteralPath $OutputZip -PathType Leaf)) {
    throw 'P0.9 ZIP missing despite zero exit code'
}
$Summary = Get-Content -LiteralPath $SummaryPath -Raw | ConvertFrom-Json
if ($Summary.scientific_acceptance -ne 'NOT_APPROVED' -or $Summary.merge_permitted -ne $false) {
    throw 'P0.9 acceptance safeguards mismatch'
}
$FinalP08Sha256 = (Get-FileHash -LiteralPath $P08Zip -Algorithm SHA256).Hash.ToLowerInvariant()
if ($FinalP08Sha256 -ne $ActualP08Sha256) {
    throw 'P0.8 input ZIP changed during P0.9'
}
$OutputSha256 = (Get-FileHash -LiteralPath $OutputZip -Algorithm SHA256).Hash.ToLowerInvariant()
Write-Host "P0.9 ZIP: $OutputZip" -ForegroundColor Green
Write-Host "P0.9 ZIP SHA-256: $OutputSha256"
Write-Host 'scientific_acceptance=NOT_APPROVED'
Write-Host 'merge_permitted=false'
