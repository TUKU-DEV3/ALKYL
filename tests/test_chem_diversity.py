# tests/test_chem_diversity.py
import json
import subprocess
import sys
from tests.conftest import ASPIRIN_SMILES, CAFFEINE_SMILES, SCRIPTS_DIR, SCRIPT_TIMEOUT

SCRIPT = "chem_diversity.py"


def make_library(tmp_path, mols):
    """Write a .smi file: 'SMILES name' per line."""
    p = tmp_path / "lib.smi"
    p.write_text("\n".join(f"{s} {n}" for s, n in mols))
    return str(p)


def run_diversity(args):
    result = subprocess.run(
        [sys.executable, str(SCRIPTS_DIR / SCRIPT)] + args,
        capture_output=True, text=True, timeout=SCRIPT_TIMEOUT,
    )
    assert result.returncode == 0, result.stderr
    return json.loads(result.stdout)


def test_output_fields(tmp_path):
    lib = make_library(tmp_path, [
        (ASPIRIN_SMILES, "aspirin"),
        (CAFFEINE_SMILES, "caffeine"),
    ])
    result = run_diversity(["--input", lib, "--n", "2"])
    for f in ["n_requested", "n_selected", "n_library", "selected"]:
        assert f in result


def test_selected_item_fields(tmp_path):
    lib = make_library(tmp_path, [
        (ASPIRIN_SMILES, "aspirin"),
        (CAFFEINE_SMILES, "caffeine"),
    ])
    result = run_diversity(["--input", lib, "--n", "2"])
    for item in result["selected"]:
        for k in ["idx", "name", "smiles"]:
            assert k in item


def test_selects_correct_n(tmp_path):
    lib = make_library(tmp_path, [
        (ASPIRIN_SMILES, "aspirin"),
        (CAFFEINE_SMILES, "caffeine"),
        ("CC(=O)Oc1ccccc1", "phenyl_acetate"),
        ("CC(C)Cc1ccc(C(C)C(=O)O)cc1", "ibuprofen_like"),
    ])
    result = run_diversity(["--input", lib, "--n", "2"])
    assert result["n_selected"] == 2
    assert len(result["selected"]) == 2


def test_no_duplicates(tmp_path):
    lib = make_library(tmp_path, [
        (ASPIRIN_SMILES, "aspirin"),
        (CAFFEINE_SMILES, "caffeine"),
        ("CC", "ethane"),
        ("CCC", "propane"),
    ])
    result = run_diversity(["--input", lib, "--n", "3"])
    indices = [s["idx"] for s in result["selected"]]
    assert len(indices) == len(set(indices))


def test_n_larger_than_library(tmp_path):
    """When n >= library size, return all molecules."""
    lib = make_library(tmp_path, [
        (ASPIRIN_SMILES, "aspirin"),
        (CAFFEINE_SMILES, "caffeine"),
    ])
    result = run_diversity(["--input", lib, "--n", "10"])
    assert result["n_selected"] == 2
    assert result["n_library"] == 2


def test_n_library_count(tmp_path):
    lib = make_library(tmp_path, [
        (ASPIRIN_SMILES, "aspirin"),
        (CAFFEINE_SMILES, "caffeine"),
        ("CC", "ethane"),
    ])
    result = run_diversity(["--input", lib, "--n", "2"])
    assert result["n_library"] == 3
    assert result["n_requested"] == 2
