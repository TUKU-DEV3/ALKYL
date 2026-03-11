# tests/test_chem_admet.py
import json
import subprocess
import sys
from tests.conftest import ASPIRIN_SMILES, CAFFEINE_SMILES, SCRIPTS_DIR, SCRIPT_TIMEOUT

SCRIPT = "chem_admet.py"

# Lipophilic tertiary amine — should flag hERG risk
VERAPAMIL_SMILES = "COc1ccc(CCN(C)CCCC(C#N)(c2ccc(OC)c(OC)c2)C(C)C)cc1OC"


def run_admet(smiles):
    result = subprocess.run(
        [sys.executable, str(SCRIPTS_DIR / SCRIPT), "--smiles", smiles],
        capture_output=True, text=True, timeout=SCRIPT_TIMEOUT,
    )
    assert result.returncode == 0, result.stderr
    return json.loads(result.stdout)


def test_output_top_level_fields():
    d = run_admet(ASPIRIN_SMILES)
    for f in ["smiles", "solubility", "bbb_penetration", "herg_risk",
              "pgp_substrate", "ppb", "disclaimer"]:
        assert f in d


def test_solubility_fields():
    d = run_admet(ASPIRIN_SMILES)
    s = d["solubility"]
    assert "logS_esol" in s
    assert "class" in s
    assert "method" in s
    assert isinstance(s["logS_esol"], float)


def test_solubility_class_valid():
    d = run_admet(ASPIRIN_SMILES)
    valid_classes = {"Highly soluble", "Soluble", "Moderately soluble",
                     "Poorly soluble", "Insoluble"}
    assert d["solubility"]["class"] in valid_classes


def test_bbb_fields():
    d = run_admet(CAFFEINE_SMILES)
    b = d["bbb_penetration"]
    assert "prediction" in b
    assert "score" in b
    assert "criteria" in b
    assert b["prediction"] in {"likely+", "uncertain", "likely-"}


def test_herg_fields():
    d = run_admet(ASPIRIN_SMILES)
    h = d["herg_risk"]
    assert "risk" in h
    assert "alerts" in h
    assert h["risk"] in {"low", "medium", "high"}
    assert isinstance(h["alerts"], list)


def test_herg_low_risk_for_aspirin():
    d = run_admet(ASPIRIN_SMILES)
    assert d["herg_risk"]["risk"] == "low"


def test_herg_elevated_risk_for_lipophilic_amine():
    d = run_admet(VERAPAMIL_SMILES)
    assert d["herg_risk"]["risk"] in {"medium", "high"}


def test_pgp_fields():
    d = run_admet(ASPIRIN_SMILES)
    p = d["pgp_substrate"]
    assert "substrate_likely" in p
    assert "criteria" in p
    assert isinstance(p["substrate_likely"], bool)


def test_ppb_fields():
    d = run_admet(ASPIRIN_SMILES)
    ppb = d["ppb"]
    assert "high_ppb_likely" in ppb
    assert "predicted_ppb" in ppb
    assert isinstance(ppb["high_ppb_likely"], bool)


def test_output_file(tmp_path):
    out = tmp_path / "admet.json"
    result = subprocess.run(
        [sys.executable, str(SCRIPTS_DIR / SCRIPT),
         "--smiles", ASPIRIN_SMILES, "--out", str(out)],
        capture_output=True, text=True, timeout=SCRIPT_TIMEOUT,
    )
    assert result.returncode == 0
    assert out.exists()
    d = json.loads(out.read_text())
    assert "solubility" in d
