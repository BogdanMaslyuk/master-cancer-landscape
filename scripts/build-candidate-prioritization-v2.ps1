param(
    [switch]$SkipAudit,
    [switch]$SkipBenchmark
)

$ErrorActionPreference = "Stop"
$Root = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$Python = Join-Path $Root ".venv\Scripts\python.exe"
$Processed = Join-Path $Root "data\processed"
$RawDepMap = Join-Path $Root "data\raw\depmap"

if (-not (Test-Path $Python)) {
    throw "Python environment not found: $Python"
}

Set-Location $Root

function Run-Step([string]$Label, [string[]]$Arguments) {
    Write-Host "`n=== $Label ===" -ForegroundColor Cyan
    & $Python @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "$Label failed with exit code $LASTEXITCODE"
    }
}

if (-not $SkipAudit) {
    Run-Step "Audit pharmacology layer" @(".\scripts\audit_pharmacology_layer.py")
}

Run-Step "Build pinned DepMap multi-omics indexes" @(".\scripts\build_depmap_multiomics.py", "--allow-partial")

$Dependency = Join-Path $Processed "depmap_model_gene_dependency.parquet"
if (-not (Test-Path $Dependency)) {
    Write-Host "`nCandidate v2 requires CRISPR Probability of Dependency for binary dependent/non-dependent calls." -ForegroundColor Red
    Write-Host "Expected raw source in the pinned DepMap release: CRISPRGeneDependency.csv" -ForegroundColor Yellow
    Write-Host "Current DepMap raw root: $RawDepMap" -ForegroundColor Yellow
    throw "Missing depmap_model_gene_dependency.parquet. Do not fall back to Gene Effect for Candidate v2."
}

Run-Step "Rebuild gene dependency summary with Probability of Dependency" @(".\scripts\build_gene_dependency_summary.py")
Run-Step "Build model × target molecular context" @(".\scripts\build_depmap_target_context.py")

$ContextManifest = Join-Path $Processed "depmap_model_target_context_manifest.json"
if (Test-Path $ContextManifest) {
    $Context = Get-Content $ContextManifest -Raw | ConvertFrom-Json
    Write-Host "`n=== Molecular context availability ===" -ForegroundColor Cyan
    Write-Host ("Dependency probability: " + $Context.available_layers.dependency_probability)
    Write-Host ("RNA expression:        " + $Context.available_layers.expression_log2_tpm1)
    Write-Host ("Relative copy number:  " + $Context.available_layers.copy_number_relative)
    if ($Context.mutation_source) {
        Write-Host ("Mutation source:       " + $Context.mutation_source)
    } else {
        Write-Host "Mutation source:       MISSING" -ForegroundColor Yellow
        Write-Host "Candidate v2 can still build, but mutation context will be absent until a 26Q1 mutation source is added." -ForegroundColor Yellow
    }
}

Run-Step "Build lineage-adjusted pharmacology × CRISPR concordance v2" @(".\scripts\build_pharmacology_target_concordance_v2.py")
Run-Step "Build methodology-corrected candidate hypotheses v2" @(".\scripts\build_candidate_hypotheses_v2.py")

if (-not $SkipBenchmark) {
    Run-Step "Known mechanism benchmark on Candidate v2" @(".\scripts\build_known_mechanism_benchmark.py")
}

Write-Host "`nMCL Candidate Prioritization v2 build completed." -ForegroundColor Green
Write-Host "Restart backend/frontend so /hypotheses switches from v1 to v2."
Write-Host "Review outputs\qc\known_mechanism_benchmark.tsv before using novel candidates as laboratory priorities."
