# tests/test_chem_3d.py
import subprocess, sys
from pathlib import Path
from tests.conftest import run_script, SCRIPTS_DIR, ASPIRIN_SMILES

SCRIPT = "chem_3d.py"

def test_generate_conformers_json():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN_SMILES,
                                 "--conformers", "5"])
    assert result["n_conformers"] == 5
    assert "energies_kcal" in result

def test_generate_conformers_sdf(tmp_path):
    out_file = tmp_path / "confs.sdf"
    proc = subprocess.run(
        [sys.executable, str(SCRIPTS_DIR / SCRIPT),
         "--smiles", ASPIRIN_SMILES, "--conformers", "3",
         "--out", str(out_file)],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0
    assert out_file.exists()
    content = out_file.read_text()
    assert content.count("$$$$") == 3

def test_minimize_mmff():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN_SMILES,
                                 "--conformers", "3", "--minimize", "mmff94"])
    assert result["n_conformers"] == 3
    assert all(e is not None for e in result["energies_kcal"])
