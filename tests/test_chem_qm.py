# tests/test_chem_qm.py
import json, subprocess, sys
from pathlib import Path
from tests.conftest import run_script, SCRIPTS_DIR, ASPIRIN_SMILES

SCRIPT = "chem_qm.py"

def test_generate_orca_input(tmp_path):
    out_file = tmp_path / "aspirin.inp"
    proc = subprocess.run(
        [sys.executable, str(SCRIPTS_DIR / SCRIPT),
         "--smiles", ASPIRIN_SMILES,
         "--engine", "orca", "--task", "opt",
         "--method", "B3LYP", "--basis", "6-31G*",
         "--charge", "0", "--mult", "1",
         "--out", str(out_file)],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0
    content = out_file.read_text()
    assert "B3LYP" in content
    assert "6-31G*" in content
    assert " Opt" in content

def test_generate_gaussian_input(tmp_path):
    out_file = tmp_path / "aspirin.gjf"
    proc = subprocess.run(
        [sys.executable, str(SCRIPTS_DIR / SCRIPT),
         "--smiles", ASPIRIN_SMILES,
         "--engine", "gaussian", "--task", "sp",
         "--method", "B3LYP", "--basis", "6-31G*",
         "--out", str(out_file)],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0
    content = out_file.read_text()
    assert "#P B3LYP" in content
    assert "6-31G*" in content

def test_stdout_json_when_no_out():
    """When --out is not given, print JSON metadata to stdout."""
    result = run_script(SCRIPT, [
        "--smiles", ASPIRIN_SMILES,
        "--engine", "orca", "--task", "sp",
        "--method", "HF", "--basis", "STO-3G",
    ])
    assert "engine" in result
    assert result["engine"] == "orca"


# ---------------------------------------------------------------------------
# IR parsing tests (ORCA block format)
# ---------------------------------------------------------------------------

ORCA_IR_OUTPUT = """\
$vibrational_frequencies
7
0: 0.0
1: 0.0
2: 0.0
3: 0.0
4: 0.0
5: 0.0
6: 1750.32
$end

$ir_spectrum
7
0: 0.0
1: 0.0
2: 0.0
3: 0.0
4: 0.0
5: 0.0
6: 234.5
$end
"""

ORCA_IR_IMAGINARY = """\
$vibrational_frequencies
7
0: 0.0
1: 0.0
2: 0.0
3: 0.0
4: 0.0
5: 0.0
6: -150.32
$end

$ir_spectrum
7
0: 0.0
1: 0.0
2: 0.0
3: 0.0
4: 0.0
5: 0.0
6: 180.0
$end
"""


def test_parse_ir_output_fields(tmp_path):
    log = tmp_path / "freq.log"
    log.write_text(ORCA_IR_OUTPUT)
    result = run_script(SCRIPT, ["--parse", str(log), "--parse-ir"])
    for f in ["source", "frequencies_cm1", "n_imaginary", "n_modes"]:
        assert f in result


def test_parse_ir_frequencies(tmp_path):
    log = tmp_path / "freq.log"
    log.write_text(ORCA_IR_OUTPUT)
    result = run_script(SCRIPT, ["--parse", str(log), "--parse-ir"])
    freqs = result["frequencies_cm1"]
    assert any(abs(f["freq_cm1"] - 1750.32) < 0.1 for f in freqs)


def test_parse_ir_no_imaginary(tmp_path):
    log = tmp_path / "freq.log"
    log.write_text(ORCA_IR_OUTPUT)
    result = run_script(SCRIPT, ["--parse", str(log), "--parse-ir"])
    assert result["n_imaginary"] == 0


def test_parse_ir_detects_imaginary(tmp_path):
    log = tmp_path / "freq_imag.log"
    log.write_text(ORCA_IR_IMAGINARY)
    result = run_script(SCRIPT, ["--parse", str(log), "--parse-ir"])
    assert result["n_imaginary"] == 1


def test_parse_ir_intensity(tmp_path):
    log = tmp_path / "freq.log"
    log.write_text(ORCA_IR_OUTPUT)
    result = run_script(SCRIPT, ["--parse", str(log), "--parse-ir"])
    freqs = result["frequencies_cm1"]
    # Mode 6 has intensity 234.5
    mode6 = next(f for f in freqs if abs(f["freq_cm1"] - 1750.32) < 0.1)
    assert abs(mode6["intensity"] - 234.5) < 0.01
