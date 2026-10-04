# MCL reproducible environment

Architecture v1 uses two independent dependency locks:

- Python 3.13: `constraints/python-3.13.txt`
- Explorer frontend: `apps/explorer/package-lock.json`

The Python `pyproject.toml` files remain the human-readable declaration of supported dependency ranges. The constraints snapshot pins the concrete resolved versions used by CI and local Architecture v1 verification.

## Clean Windows setup

Use PowerShell #3 (diagnostics / Git / builds):

```powershell
cd C:\Users\Bogdan\Desktop\master-cancer-landscape
Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
.\scripts\bootstrap.ps1
```

The bootstrap script:

1. creates `.venv` with Python 3.13 when it is missing;
2. pins pip to the Architecture v1 bootstrap version;
3. installs the scientific core and API against `constraints/python-3.13.txt`;
4. runs `pip check`;
5. verifies that every direct Python dependency declared by both `pyproject.toml` files is represented in the constraints snapshot;
6. installs the frontend with `npm ci` from `package-lock.json`.

To deliberately rebuild the Python environment from zero:

```powershell
.\scripts\bootstrap.ps1 -Recreate
```

This deletes only the local `.venv`; project data and runtime indexes are not removed.

## Verification

After bootstrap:

```powershell
.\scripts\verify.ps1
```

The verification gate checks the Python constraint contract and `pip check` before scientific/API/frontend tests.

## Updating Python dependencies

Dependency updates are explicit Architecture changes:

1. change supported ranges in the relevant `pyproject.toml`;
2. resolve and test a new Python 3.13 environment;
3. update `constraints/python-3.13.txt` with the tested resolved versions;
4. run `scripts/verify_python_constraints.py`;
5. run the full `scripts/verify.ps1` gate;
6. commit the declaration and constraint changes together.

Do not casually run an unconstrained `pip install --upgrade` in CI and call the result reproducible.

## Scope of the Python constraints snapshot

The committed snapshot pins the complete package set observed in the successful Linux CI environment. A small number of OS-specific helper packages may be added by pip on another platform when a dependency declares a platform marker. `uvloop` is explicitly disabled on Windows. The important invariant is that all MCL-declared direct dependencies and the shared scientific/API graph are pinned and continuously checked.
