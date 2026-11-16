# MCL E0 AC-05: real browser check in one MCL_Project. Windows PowerShell 5.1 compatible.
param(
    [string]$Root = (Join-Path $env:USERPROFILE 'MCL_Project'),
    [string]$UiBase = 'http://127.0.0.1:3100',
    [string]$ApiBase = 'http://127.0.0.1:8100',
    [string]$Expected = ''
)
$ErrorActionPreference = 'Stop'
$Root = [System.IO.Path]::GetFullPath($Root)
$Reports = Join-Path $Root 'reports'
$Audit = Join-Path $PSScriptRoot 'audit_pyz_e0_ac05_browser.cjs'
$ScienceRoot = Join-Path $Root 'repositories\mcl-pyz-v2-review'
if (-not (Test-Path -LiteralPath $Audit -PathType Leaf)) { throw 'Audit script missing' }
if (-not (Test-Path -LiteralPath $Reports -PathType Container)) { throw 'MCL_Project\reports missing' }
if (-not (Test-Path -LiteralPath (Join-Path $ScienceRoot 'config\own_compounds.tsv') -PathType Leaf)) {
    throw 'The central MCL scientific worktree is missing. Stop to avoid checking a different MCL.'
}
$Node = Get-Command node.exe -ErrorAction SilentlyContinue
if (-not $Node) { throw 'Node.js missing. Install Node.js, then retry.' }
if (-not ($UiBase -eq 'http://127.0.0.1:3100' -or $UiBase -eq 'http://localhost:3100')) {
    throw 'Only the local MCL UI on port 3100 is allowed by this launcher.'
}
if (-not ($ApiBase -eq 'http://127.0.0.1:8100' -or $ApiBase -eq 'http://localhost:8100')) {
    throw 'Only the local MCL API on port 8100 is allowed by this launcher.'
}
# Do not start or stop either service; verify the API serves the expected scientific copy.
try {
    $Health = Invoke-RestMethod -Uri ($ApiBase + '/health') -Method Get -TimeoutSec 15
} catch { throw 'MCL API (PS1, port 8100) is not running. Launch PS1 and PS2 from MCL_Project, then retry.' }
$LiveRoot = [System.IO.Path]::GetFullPath([string]$Health.mcl_root).TrimEnd('\')
$ExpectedRoot = [System.IO.Path]::GetFullPath($ScienceRoot).TrimEnd('\')
if (-not [string]::Equals($LiveRoot, $ExpectedRoot, [StringComparison]::OrdinalIgnoreCase)) {
    throw "API reads a different scientific copy: $LiveRoot ; expected $ExpectedRoot"
}
try {
    $UiCheck = Invoke-WebRequest -Uri ($UiBase + '/pyz') -UseBasicParsing -Method Get -TimeoutSec 45
} catch { throw 'MCL website (PS2, port 3100) is not reachable. Launch PS2 and retry.' }
if ($UiCheck.StatusCode -ne 200) { throw "MCL /pyz returned HTTP $($UiCheck.StatusCode)" }
$ToolsCandidates = @(
    (Join-Path $Root 'reports\_tools\ac05'),
    (Join-Path $Root 'repositories\mcl-pyz-e0-audit\outputs\audits\pyz_e0\_ui_tools'),
    (Join-Path $env:USERPROFILE 'Desktop\mcl-pyz-e0-audit\outputs\audits\pyz_e0\_ui_tools')
)
$Tools = $null
foreach ($Candidate in $ToolsCandidates) {
    if (Test-Path -LiteralPath (Join-Path $Candidate 'node_modules\playwright-core\package.json')) {
        $Tools = $Candidate
        break
    }
}
if (-not $Tools) {
    $Tools = $ToolsCandidates[0]
    $Npm = Get-Command npm.cmd -ErrorAction SilentlyContinue
    if (-not $Npm) { throw 'playwright-core is not installed and npm.cmd is missing.' }
    Write-Host 'Installing playwright-core into reports\_tools\ac05 (no browser download)...' -ForegroundColor Cyan
    & npm.cmd install --prefix $Tools --no-save --ignore-scripts --no-audit --no-fund 'playwright-core@1.56.1'
    if ($LASTEXITCODE -ne 0) { throw 'playwright-core installation failed. Scientific files unchanged.' }
    if (-not (Test-Path -LiteralPath (Join-Path $Tools 'node_modules\playwright-core\package.json'))) {
        throw 'playwright-core was not installed'
    }
}
if (-not $Expected) {
    $Roots = @(
        (Join-Path $Root 'repositories\mcl-pyz-e0-audit\outputs\audits\pyz_e0'),
        (Join-Path $env:USERPROFILE 'Desktop\mcl-pyz-e0-audit\outputs\audits\pyz_e0')
    )
    foreach ($R in $Roots) {
        if (-not (Test-Path -LiteralPath $R -PathType Container)) { continue }
        $Found = Get-ChildItem -LiteralPath $R -Directory | Sort-Object Name -Descending |
            ForEach-Object { Join-Path $_.FullName 'closeout\ui_expected.json' } |
            Where-Object { Test-Path -LiteralPath $_ -PathType Leaf } | Select-Object -First 1
        if ($Found) { $Expected = $Found; break }
    }
}
$Stamp = Get-Date -Format 'yyyyMMdd_HHmmss_fff'
$Out = Join-Path $Reports ('AC05_BROWSER_' + $Stamp)
$Zip = Join-Path $Reports ('MCL_E0_AC05_BROWSER_' + $Stamp + '.zip')
$Args = @($Audit, '--url', $UiBase, '--out', $Out, '--tools', $Tools)
if ($Expected) { $Args += @('--expected', $Expected) }
Write-Host '=== MCL E0 / AC-05: real browser verification ===' -ForegroundColor Cyan
Write-Host "Origin: $UiBase"
Write-Host "Browser tools: $Tools"
Write-Host "Output: $Out"
& $Node.Source @Args
$AuditExit = $LASTEXITCODE
if (Test-Path -LiteralPath $Out -PathType Container) {
    Compress-Archive -Path (Join-Path $Out '*') -DestinationPath $Zip -CompressionLevel Optimal
    if (-not (Test-Path -LiteralPath $Zip -PathType Leaf)) { throw 'Report ZIP was not created' }
    Write-Host "Report ZIP: $Zip" -ForegroundColor Green
}
if ($AuditExit -ne 0) { throw "Browser check was not completed (exit $AuditExit). Attach report ZIP if one was created." }
Write-Host 'AC-05 browser check complete. Manual review and E0 acceptance are still pending.' -ForegroundColor Yellow