from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
API_DIR = ROOT / "apps" / "api"
OUTPUT = ROOT / "apps" / "explorer" / "lib" / "generated" / "api-types.ts"
if str(API_DIR) not in sys.path:
    sys.path.insert(0, str(API_DIR))

from mcl_api.main import app  # noqa: E402


EXCLUDED_SCHEMAS = {"HTTPValidationError", "ValidationError"}
HEADER = """// AUTO-GENERATED FILE. DO NOT EDIT MANUALLY.\n// Source of truth: FastAPI OpenAPI schema (`mcl_api.main:app`).\n// Regenerate with: .\\.venv\\Scripts\\python.exe .\\scripts\\generate_frontend_api_types.py\n\n"""


def _ref_name(ref: str) -> str:
    return ref.rsplit("/", 1)[-1]


def _unique(parts: list[str]) -> list[str]:
    seen: set[str] = set()
    output: list[str] = []
    for part in parts:
        if part not in seen:
            seen.add(part)
            output.append(part)
    return output


def _render(schema: dict[str, Any], indent: int = 0) -> str:
    if "$ref" in schema:
        return _ref_name(str(schema["$ref"]))

    if "enum" in schema:
        return " | ".join(json.dumps(value, ensure_ascii=False) for value in schema["enum"])

    for keyword in ("anyOf", "oneOf"):
        if keyword in schema:
            return " | ".join(_unique([_render(item, indent) for item in schema[keyword]]))

    if "allOf" in schema:
        return " & ".join(_unique([_render(item, indent) for item in schema["allOf"]]))

    schema_type = schema.get("type")
    if schema_type == "null":
        return "null"
    if schema_type == "string":
        return "string"
    if schema_type in {"integer", "number"}:
        return "number"
    if schema_type == "boolean":
        return "boolean"
    if schema_type == "array":
        item_type = _render(schema.get("items") or {}, indent)
        if " | " in item_type or " & " in item_type:
            item_type = f"({item_type})"
        return f"{item_type}[]"

    if schema_type == "object" or "properties" in schema or "additionalProperties" in schema:
        properties = schema.get("properties") or {}
        additional = schema.get("additionalProperties")
        if not properties:
            if additional is True or additional is None:
                return "Record<string, unknown>"
            if isinstance(additional, dict):
                return f"Record<string, {_render(additional, indent)}>"
            return "{}"

        required = set(schema.get("required") or [])
        lines = ["{"]
        if additional is True:
            lines.append("  " * (indent + 1) + "[key: string]: unknown;")
        elif isinstance(additional, dict):
            lines.append(
                "  " * (indent + 1)
                + f"[key: string]: {_render(additional, indent + 1)};"
            )

        for name, property_schema in properties.items():
            optional = "" if name in required else "?"
            lines.append(
                "  " * (indent + 1)
                + f"{json.dumps(name, ensure_ascii=False)}{optional}: {_render(property_schema, indent + 1)};"
            )
        lines.append("  " * indent + "}")
        return "\n".join(lines)

    return "unknown"


def generate() -> str:
    openapi = app.openapi()
    schemas = (openapi.get("components") or {}).get("schemas") or {}
    selected = {
        name: schema
        for name, schema in schemas.items()
        if name not in EXCLUDED_SCHEMAS
    }
    chunks = [HEADER]
    for name in sorted(selected):
        chunks.append(f"export type {name} = {_render(selected[name])};\n\n")
    return "".join(chunks)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--check",
        action="store_true",
        help="Fail if the committed generated TypeScript contract is stale.",
    )
    args = parser.parse_args()

    generated = generate()
    if args.check:
        if not OUTPUT.exists():
            raise SystemExit(f"Generated API types are missing: {OUTPUT.relative_to(ROOT)}")
        current = OUTPUT.read_text(encoding="utf-8")
        if current != generated:
            raise SystemExit(
                "Generated frontend API types are stale. Run "
                ".\\.venv\\Scripts\\python.exe .\\scripts\\generate_frontend_api_types.py"
            )
        print("Frontend API type contract: PASS")
        return

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(generated, encoding="utf-8", newline="\n")
    print(f"Wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
