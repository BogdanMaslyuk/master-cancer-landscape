"""Regression guard for Windows PowerShell 5.1 UTF-8 decoding.

PowerShell 5.1 interprets a UTF-8 script lacking BOM through the legacy ANSI
code page, corrupting Cyrillic strings and causing a ParserError.
"""
from pathlib import Path


def test_e0_acceptance_gate_ps51_has_utf8_bom():
    script = (Path(__file__).resolve().parents[1] / "scripts"
              / "run-pyz-e0-acceptance-gate.ps1")
    raw = script.read_bytes()
    assert raw.startswith(b"\xef\xbb\xbf"), "Windows PowerShell 5.1 requires a UTF-8 BOM"
    decoded = raw.decode("utf-8-sig")
    assert not decoded.startswith("\ufeff"), "Expected exactly one UTF-8 BOM"
    assert "\ufffd" not in decoded, "Replacement characters detected"
    assert "Э0" in decoded, "Cyrillic log messages have been corrupted"
    assert "$InputPath" in decoded
