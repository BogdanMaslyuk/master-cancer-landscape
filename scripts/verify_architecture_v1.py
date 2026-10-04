from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

REQUIRED_PATHS = [
    "docs/ARCHITECTURE.md",
    "docs/DEVELOPMENT.md",
    "docs/REPRODUCIBILITY.md",
    "constraints/python-3.13.txt",
    "apps/explorer/package-lock.json",
    "apps/explorer/lib/generated/api-types.ts",
    "apps/api/mcl_api/routers",
    "apps/api/mcl_api/services",
    "apps/api/mcl_api/repositories",
    "apps/api/mcl_api/schemas",
    "scripts/bootstrap.ps1",
    "scripts/build-explorer.ps1",
    "scripts/verify.ps1",
    "scripts/verify_runtime_indexes.py",
    "scripts/generate_frontend_api_types.py",
]


def main() -> None:
    missing = [relative for relative in REQUIRED_PATHS if not (ROOT / relative).exists()]
    if missing:
        raise SystemExit(
            "Architecture v1 repository contract is incomplete. Missing: " + ", ".join(missing)
        )

    architecture = (ROOT / "docs" / "ARCHITECTURE.md").read_text(encoding="utf-8")
    required_statements = [
        "Status: Architecture v1 complete",
        "data/runtime/explorer/",
        "router -> service -> repository",
        "constraints/python-3.13.txt",
        "mcl-architecture-v1",
    ]
    absent = [statement for statement in required_statements if statement not in architecture]
    if absent:
        raise SystemExit(
            "Architecture v1 documentation contract is incomplete. Missing statements: "
            + ", ".join(absent)
        )

    print("Architecture v1 repository contract: PASS")


if __name__ == "__main__":
    main()
