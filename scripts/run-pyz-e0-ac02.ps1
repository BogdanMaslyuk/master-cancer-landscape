# MCL E0 / AC-02; Windows PowerShell 5.1 UTF-8 BOM.
param(
    [string]$Root = (Join-Path $env:USERPROFILE 'MCL_Project'),
    [string]$ScienceRoot = '',
    [string]$ComparisonRoot = '',
    [string]$BaselineHashes = ''
)
$ErrorActionPreference = 'Stop'
$Root = [IO.Path]::GetFullPath($Root)
$Reports = Join-Path $Root 'reports'
if (-not $ScienceRoot) { $ScienceRoot = Join-Path $Root 'repositories\mcl-pyz-v2-review' }
$ScienceRoot = [IO.Path]::GetFullPath($ScienceRoot)
if (-not $ComparisonRoot) { $ComparisonRoot = Join-Path $Root 'repositories\master-cancer-landscape' }
$ComparisonRoot = [IO.Path]::GetFullPath($ComparisonRoot)
if (-not (Test-Path -LiteralPath (Join-Path $ComparisonRoot 'data\processed') -PathType Container)) {
    throw "Comparison source folder missing: $ComparisonRoot"
}
$Registry = Join-Path $ScienceRoot 'data\runtime\explorer'
if (-not (Test-Path -LiteralPath $Registry -PathType Container)) {
    throw "Explorer scientific source folder missing: $Registry"
}
if (-not (Test-Path -LiteralPath $Reports -PathType Container)) {
    throw "Missing MCL_Project reports folder: $Reports"
}
$Audit = Join-Path $PSScriptRoot 'audit_pyz_e0_ac02.py'
if (-not (Test-Path -LiteralPath $Audit -PathType Leaf)) { throw 'Missing AC-02 audit script' }
$PythonCandidates = @(
    (Join-Path $Root 'repositories\master-cancer-landscape\.venv\Scripts\python.exe'),
    (Join-Path $ScienceRoot '.venv\Scripts\python.exe'),
    (Join-Path $env:USERPROFILE 'Desktop\master-cancer-landscape\.venv\Scripts\python.exe')
)
$Python = $PythonCandidates | Where-Object { Test-Path -LiteralPath $_ -PathType Leaf } | Select-Object -First 1
if (-not $Python) { throw 'Python venv not found inside existing MCL installation' }
if (-not $BaselineHashes) {
    $Candidates = @(
        (Join-Path $Root 'repositories\mcl-pyz-e0-audit\outputs\audits\pyz_e0'),
        (Join-Path $env:USERPROFILE 'Desktop\mcl-pyz-e0-audit\outputs\audits\pyz_e0')
    )
    foreach ($Dir in $Candidates) {
        if (-not (Test-Path -LiteralPath $Dir -PathType Container)) { continue }
        $Found = Get-ChildItem -LiteralPath $Dir -Directory | Sort-Object Name -Descending |
            ForEach-Object { Join-Path $_.FullName 'file_checksums.sha256' } |
            Where-Object { Test-Path -LiteralPath $_ -PathType Leaf } |
            Select-Object -First 1
        if ($Found) { $BaselineHashes = $Found; break }
    }
}
$Stamp = Get-Date -Format 'yyyyMMdd_HHmmss_fff'
$Out = Join-Path $Reports ('AC02_GENETIC_' + $Stamp)
$Zip = Join-Path $Reports ('MCL_E0_AC02_GENETIC_' + $Stamp + '.zip')
$Arguments = @($Audit, '--root', $ScienceRoot, '--comparison-root', $ComparisonRoot, '--out', $Out)
if ($BaselineHashes) { $Arguments += @('--baseline-hashes', $BaselineHashes) }
else { Write-Warning 'Historical E0 checksum baseline unavailable; comparisons marked NOT_AVAILABLE.' }
Write-Host '=== MCL E0 / AC-02: genetic data read-only audit ===' -ForegroundColor Cyan
Write-Host "Scientific root: $ScienceRoot"
Write-Host "Comparison root (read-only): $ComparisonRoot"
Write-Host "Python: $Python"
Write-Host "Historical checksum record: $BaselineHashes"
Write-Host 'PS1 and PS2 are not required for this file-only audit.'
& $Python @Arguments
$ResultCode = $LASTEXITCODE
if ($ResultCode -notin @(0,2)) { throw "AC-02 failed before report completion (exit $ResultCode)" }
if (-not (Test-Path -LiteralPath (Join-Path $Out 'AC02_REPORT.json') -PathType Leaf)) {
    throw 'AC-02 report was not created'
}
$Report = Get-Content -LiteralPath (Join-Path $Out 'AC02_REPORT.json') -Raw -Encoding UTF8 | ConvertFrom-Json
if ($Report.merge_permitted -ne $false -or $Report.ac02_scientific_acceptance -ne 'NOT_APPROVED') {
    throw 'Refusing to publish an automatically accepted scientific report'
}
Compress-Archive -Path (Join-Path $Out '*') -DestinationPath $Zip -CompressionLevel Optimal
if (-not (Test-Path -LiteralPath $Zip -PathType Leaf)) { throw 'AC-02 ZIP was not created' }
Write-Host "Report ZIP: $Zip" -ForegroundColor Green
Write-Host "Technical status: $($Report.technical_status)" -ForegroundColor Yellow
Write-Host 'Scientific approval: NOT APPROVED. No branch merged or scientific data changed.' -ForegroundColor Yellow
if ($ResultCode -eq 2) { Write-Warning 'Automated blockers detected. Attach ZIP for scientific diagnosis; do NOT rebuild sources automatically.' }