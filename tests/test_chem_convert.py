import subprocess, sys, json
from tests.conftest import run_script, ASPIRIN_SMILES

SCRIPT = "chem_convert.py"

def test_smiles_to_inchi():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN_SMILES, "--to", "inchi"])
    assert "inchi" in result
    assert result["inchi"].startswith("InChI=")

def test_smiles_to_smiles_roundtrip():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN_SMILES, "--to", "smiles"])
    assert "smiles" in result
    assert len(result["smiles"]) > 0

def test_smiles_to_svg():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN_SMILES, "--to", "svg"])
    assert "svg" in result
    assert "<svg" in result["svg"]

def test_invalid_smiles_exits_nonzero():
    from tests.conftest import SCRIPTS_DIR
    proc = subprocess.run(
        [sys.executable, str(SCRIPTS_DIR / "chem_convert.py"),
         "--smiles", "NOT_A_SMILES", "--to", "inchi"],
        capture_output=True, text=True,
    )
    assert proc.returncode != 0
