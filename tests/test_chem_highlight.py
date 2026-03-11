# tests/test_chem_highlight.py
import json
import subprocess
import sys
from tests.conftest import ASPIRIN_SMILES, SCRIPTS_DIR, SCRIPT_TIMEOUT

SCRIPT = "chem_highlight.py"
CARBOXYL_SMARTS = "C(=O)O"
AROMATIC_SMARTS = "c1ccccc1"


def run_highlight_to_file(smiles, smarts, out_path, extra_args=None):
    args = [sys.executable, str(SCRIPTS_DIR / SCRIPT),
            "--smiles", smiles, "--out", str(out_path)]
    if smarts:
        args += ["--smarts", smarts]
    if extra_args:
        args += extra_args
    result = subprocess.run(args, capture_output=True, text=True, timeout=SCRIPT_TIMEOUT)
    assert result.returncode == 0, result.stderr
    return json.loads(result.stdout)


def run_highlight_stdout(smiles, smarts=None):
    """Returns raw stdout (SVG string)."""
    args = [sys.executable, str(SCRIPTS_DIR / SCRIPT), "--smiles", smiles]
    if smarts:
        args += ["--smarts", smarts]
    result = subprocess.run(args, capture_output=True, text=True, timeout=SCRIPT_TIMEOUT)
    assert result.returncode == 0, result.stderr
    return result.stdout


def test_svg_file_created(tmp_path):
    out = tmp_path / "mol.svg"
    run_highlight_to_file(ASPIRIN_SMILES, CARBOXYL_SMARTS, out)
    assert out.exists()
    assert out.stat().st_size > 0


def test_svg_file_is_valid_xml(tmp_path):
    out = tmp_path / "mol.svg"
    run_highlight_to_file(ASPIRIN_SMILES, CARBOXYL_SMARTS, out)
    content = out.read_text()
    assert "<svg" in content


def test_json_metadata_fields(tmp_path):
    out = tmp_path / "mol.svg"
    meta = run_highlight_to_file(ASPIRIN_SMILES, CARBOXYL_SMARTS, out)
    for f in ["smiles", "smarts", "format", "output", "width", "height"]:
        assert f in meta


def test_json_metadata_values(tmp_path):
    out = tmp_path / "mol.svg"
    meta = run_highlight_to_file(ASPIRIN_SMILES, CARBOXYL_SMARTS, out)
    assert meta["format"] == "svg"
    assert meta["smarts"] == CARBOXYL_SMARTS
    assert meta["width"] == 400
    assert meta["height"] == 300


def test_stdout_is_svg_without_out():
    svg = run_highlight_stdout(ASPIRIN_SMILES, CARBOXYL_SMARTS)
    assert svg.strip().startswith("<?xml") or "<svg" in svg


def test_no_smarts_still_draws(tmp_path):
    """Drawing without SMARTS should succeed (no highlight, plain molecule)."""
    out = tmp_path / "plain.svg"
    meta = run_highlight_to_file(ASPIRIN_SMILES, smarts=None, out_path=out)
    assert out.exists()
    assert meta["smarts"] is None


def test_custom_dimensions(tmp_path):
    out = tmp_path / "mol.svg"
    meta = run_highlight_to_file(
        ASPIRIN_SMILES, AROMATIC_SMARTS, out,
        extra_args=["--width", "600", "--height", "400"],
    )
    assert meta["width"] == 600
    assert meta["height"] == 400


def test_non_matching_smarts_still_exits_zero(tmp_path):
    """Non-matching SMARTS should warn but not crash."""
    out = tmp_path / "mol.svg"
    result = subprocess.run(
        [sys.executable, str(SCRIPTS_DIR / SCRIPT),
         "--smiles", ASPIRIN_SMILES, "--smarts", "[Si]", "--out", str(out)],
        capture_output=True, text=True, timeout=SCRIPT_TIMEOUT,
    )
    assert result.returncode == 0
    assert out.exists()
