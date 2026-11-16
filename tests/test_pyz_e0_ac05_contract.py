"""AC-05 audit is read-only and compatible with Windows PowerShell 5.1."""
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_windows_powershell_launcher_bom():
    content = (ROOT / "scripts/run-pyz-e0-ac05-browser.ps1").read_bytes()
    assert content.startswith(b"\xef\xbb\xbf")
    text = content.decode("utf-8-sig")
    assert "MCL_Project" in text
    assert "--no-save" in text
    assert "Compress-Archive" in text


def test_script_never_approves_science_or_calls_remote_site():
    text = (ROOT / "scripts/audit_pyz_e0_ac05_browser.cjs").read_text(encoding="utf-8")
    assert "scientific_acceptance: 'NOT_ACCEPTED'" in text
    assert "merge_permitted: false" in text
    assert "localhostUrl" in text
    assert "root_overflow_cases" in text
    assert "isMobile: width === 390" in text
    assert "fs.existsSync(out)" in text