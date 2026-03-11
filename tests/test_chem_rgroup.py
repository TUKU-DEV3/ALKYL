# tests/test_chem_rgroup.py
import json
import subprocess
import sys
from tests.conftest import SCRIPTS_DIR, SCRIPT_TIMEOUT

SCRIPT = "chem_rgroup.py"

# Phenylacetic acid series: para-substituted, all match c1ccc([*:1])cc1
PARA_SERIES = [
    ("c1ccc(CC(=O)O)cc1", "unsubst"),
    ("c1ccc(F)cc1CC(=O)O", "fluoro"),
    ("c1ccc(Cl)cc1CC(=O)O", "chloro"),
    ("c1ccc(OC)cc1CC(=O)O", "methoxy"),
    ("c1ccc(N)cc1CC(=O)O", "amino"),
]
CORE_PARA = "c1ccc([*:1])cc1"

# One molecule that will NOT match a pyridine core
NON_MATCHING = ("CC(=O)O", "acetic_acid")


def make_smi(tmp_path, mols):
    p = tmp_path / "lib.smi"
    p.write_text("\n".join(f"{s} {n}" for s, n in mols))
    return str(p)


def run_rgroup(args):
    result = subprocess.run(
        [sys.executable, str(SCRIPTS_DIR / SCRIPT)] + args,
        capture_output=True, text=True, timeout=SCRIPT_TIMEOUT,
    )
    assert result.returncode == 0, result.stderr
    return json.loads(result.stdout)


def test_output_top_level_fields(tmp_path):
    lib = make_smi(tmp_path, PARA_SERIES)
    d = run_rgroup(["--input", lib, "--core", CORE_PARA])
    for f in ["core_smarts", "n_input", "n_matched", "n_unmatched",
              "rgroup_labels", "decomposition", "unmatched"]:
        assert f in d


def test_all_matched(tmp_path):
    lib = make_smi(tmp_path, PARA_SERIES)
    d = run_rgroup(["--input", lib, "--core", CORE_PARA])
    assert d["n_matched"] == len(PARA_SERIES)
    assert d["n_unmatched"] == 0


def test_decomposition_entry_fields(tmp_path):
    lib = make_smi(tmp_path, PARA_SERIES)
    d = run_rgroup(["--input", lib, "--core", CORE_PARA])
    for entry in d["decomposition"]:
        assert "name" in entry
        assert "smiles" in entry
        assert "core_match" in entry


def test_rgroup_labels_present(tmp_path):
    lib = make_smi(tmp_path, PARA_SERIES)
    d = run_rgroup(["--input", lib, "--core", CORE_PARA])
    assert len(d["rgroup_labels"]) >= 1
    assert "R1" in d["rgroup_labels"]


def test_unmatched_counted(tmp_path):
    lib = make_smi(tmp_path, PARA_SERIES + [NON_MATCHING])
    d = run_rgroup(["--input", lib, "--core", "n1ccccc1"])  # pyridine core — nothing matches
    assert d["n_unmatched"] == len(PARA_SERIES) + 1
    assert d["n_matched"] == 0


def test_unmatched_entry_fields(tmp_path):
    lib = make_smi(tmp_path, PARA_SERIES + [NON_MATCHING])
    d = run_rgroup(["--input", lib, "--core", "n1ccccc1"])
    for entry in d["unmatched"]:
        assert "name" in entry
        assert "smiles" in entry


def test_core_smarts_in_output(tmp_path):
    lib = make_smi(tmp_path, PARA_SERIES)
    d = run_rgroup(["--input", lib, "--core", CORE_PARA])
    assert d["core_smarts"] == CORE_PARA


def test_output_file(tmp_path):
    lib = make_smi(tmp_path, PARA_SERIES)
    out = tmp_path / "rgroups.json"
    result = subprocess.run(
        [sys.executable, str(SCRIPTS_DIR / SCRIPT),
         "--input", lib, "--core", CORE_PARA, "--out", str(out)],
        capture_output=True, text=True, timeout=SCRIPT_TIMEOUT,
    )
    assert result.returncode == 0
    assert out.exists()
    d = json.loads(out.read_text())
    assert "decomposition" in d
