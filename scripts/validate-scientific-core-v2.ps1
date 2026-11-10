# Scientific Core v2: safe validation using this Git worktree's mcl_api sources.
# Reads the original data checkout. Does not modify Git, Parquet files or cached molecules.
param(
    [string]$SourceRoot = (Join-Path $env:USERPROFILE "Desktop\master-cancer-landscape"),
    [switch]$SkipRealData
)

$ErrorActionPreference = "Stop"
$Review = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$Source = (Resolve-Path $SourceRoot).Path
$ApiRoot = (Resolve-Path (Join-Path $Review "apps\api")).Path
$Python = Join-Path $Source ".venv\Scripts\python.exe"
if (-not (Test-Path $Python)) { throw "Python environment not found: $Python" }

$PreviousPath = (Get-Location).Path
$PreviousPythonPath = [Environment]::GetEnvironmentVariable("PYTHONPATH", "Process")
$PreviousDataRoot = [Environment]::GetEnvironmentVariable("MCL_SCIENTIFIC_DATA_ROOT", "Process")

function Assert-Exit([string]$Label) {
    if ($LASTEXITCODE -ne 0) { throw "$Label failed with exit code $LASTEXITCODE" }
}

try {
    Set-Location $Review
    # The interpreter is from another checkout and may contain an older editable
    # installation. Override PYTHONPATH BEFORE running any test.
    $env:PYTHONPATH = $ApiRoot
    # Tests use temporary isolated fixtures. Never point them at the real data.
    Remove-Item Env:\MCL_SCIENTIFIC_DATA_ROOT -ErrorAction SilentlyContinue

    Write-Host ""
    Write-Host "=== Verify active mcl_api source ===" -ForegroundColor Cyan
    @'
import pathlib
import sys
import mcl_api
from mcl_api.own_compounds import OwnCompoundStore

expected = (pathlib.Path(sys.argv[1]) / "apps" / "api" / "mcl_api").resolve()
actual = pathlib.Path(mcl_api.__file__).resolve().parent
print("Actual:", actual)
print("Expected:", expected)
if actual != expected:
    raise SystemExit("WRONG SOURCE: Python imported an older mcl_api.")
if not hasattr(OwnCompoundStore, "full_target_models"):
    raise SystemExit("WRONG SOURCE: missing full_target_models.")
print("SOURCE CHECK PASS")
'@ | & $Python - $Review
    Assert-Exit "Python source check"

    Write-Host ""
    Write-Host "=== Run Scientific Core v2 tests ===" -ForegroundColor Cyan
    & $Python -m pytest .\apps\api\tests\test_uploaded_molecules.py .\apps\api\tests\test_scientific_core_v2_cdk4_audit.py -q
    Assert-Exit "Scientific Core v2 tests"

    if (-not $SkipRealData) {
        # Only the real-data diagnostic gets the external source location.
        $env:MCL_SCIENTIFIC_DATA_ROOT = $Source
        Write-Host ""
        Write-Host "=== Validate full CDK4 source and PYZ-035 without changing cached results ===" -ForegroundColor Cyan
        @'
from pathlib import Path
from mcl_api.own_compounds import OwnCompoundStore
from mcl_api.uploaded_molecules import UploadedMoleculeStore

root = Path.cwd()
own = OwnCompoundStore(root)
store = UploadedMoleculeStore(root, own)
full = own.full_target_models("CDK4")
if not full["available"]:
    raise SystemExit("Full processed CDK4 source unavailable. Check data root.")

print("FULL CDK4 DATA")
print("Total model rows:", full["total_model_rows"])
print("Measured P(dep):", full["total_measured_models"])
print("Registered laboratory models:", full["laboratory_registered_models_n"])

expected = {"ACH-000971", "ACH-000019", "ACH-000147", "ACH-000480", "ACH-000995", "ACH-000406"}
found = {row["model_id"] for row in full["items"]}
missing = sorted(expected - found)
if missing:
    raise SystemExit("Missing CDK4 reference models: " + ", ".join(missing))
for row in full["items"]:
    if row["model_id"] in expected:
        print(row["cell_line_name"], row["model_id"],
              "P(dep)=", row["dependency_probability"],
              "status=", row["genetic_dependency_status"])

print("\nCANDIDATE V2 (KNOWN DRUGS ONLY)")
contexts = own._contexts_for_target("CDK4")
for name in ("MCF-7", "T-47D"):
    seen = []
    for context in contexts:
        for model in context.get("models", []):
            if str(model.get("cell_line_name") or "").casefold() == name.casefold():
                for obs in model.get("reference_drug_observations", []):
                    seen.append((
                        str(obs.get("preferred_name") or obs.get("compound_id")),
                        str(obs.get("reference_model_role")),
                    ))
    print(name, ":", sorted(set(seen)))
    if not seen:
        raise SystemExit("Missing per-compound reference trace for " + name)

registry = own.registry()
hit = registry.loc[registry["own_compound_id"].astype(str).eq("PYZ-035")]
if hit.empty:
    raise SystemExit("PYZ-035 not found in worktree registry")
mol, canonical, inchikey = store._standardize(str(hit.iloc[0]["standardized_smiles"]))
targets = store._target_search(mol)
match = next((t for t in targets if t["target_gene"] == "CDK4"), None)
if match is None:
    raise SystemExit("CDK4 not in current target-ligand space; check runtime source")
cells = store._cell_lines(targets)
if not cells:
    raise SystemExit("No relevant cells returned; review structural gate and source files")
assert all(c["source_layer"] == "full_processed_depmap_target_context" for c in cells), (
    "The analysis fell back to truncated Candidate v2"
)
assert all(c["direct_uploaded_molecule_response"] == "absent_not_measured" for c in cells)
panel = store._experiment_panel(targets, cells)
assert panel["complete"] is False
assert panel["laboratory_ready"] is False

print("\nPYZ-035: READ-ONLY ANALYSIS")
print("CDK4 best Tanimoto:", match["best_tanimoto"])
print("Cell models:", len(cells))
print("Panel type:", panel.get("panel_type"))
print("Status:", panel.get("status"))
print("Genetically dependent candidate models:", panel.get("genetic_candidate_models_n"))
print("Genetically nondependent comparator models:", panel.get("genetic_comparator_models_n"))
print("Laboratory ready:", panel.get("laboratory_ready"))
print("Missing:", panel.get("missing_components"))
print("\nVALIDATION COMPLETE")
'@ | & $Python -
        Assert-Exit "Full data validation"
    }
    Write-Host ""
    Write-Host "Scientific Core v2 validation: PASS" -ForegroundColor Green
} finally {
    if ($null -eq $PreviousPythonPath) {
        Remove-Item Env:\PYTHONPATH -ErrorAction SilentlyContinue
    } else {
        $env:PYTHONPATH = $PreviousPythonPath
    }
    if ($null -eq $PreviousDataRoot) {
        Remove-Item Env:\MCL_SCIENTIFIC_DATA_ROOT -ErrorAction SilentlyContinue
    } else {
        $env:MCL_SCIENTIFIC_DATA_ROOT = $PreviousDataRoot
    }
    Set-Location $PreviousPath
}
