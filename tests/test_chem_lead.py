# tests/test_chem_lead.py
import json
import subprocess
import sys
from tests.conftest import ASPIRIN_SMILES, CAFFEINE_SMILES, SCRIPTS_DIR, SCRIPT_TIMEOUT

SCRIPT = "chem_lead.py"

SERIES = [
    (ASPIRIN_SMILES, "aspirin", 500, 1),
    ("c1ccc(F)cc1CC(=O)O", "fluoro", 200, 1),
    ("c1ccc(Cl)cc1CC(=O)O", "chloro", 150, 1),
    (CAFFEINE_SMILES, "caffeine", 80, 2),
    ("c1ccc(OC)cc1CC(=O)O", "methoxy", 50, 2),
]


def make_csv(tmp_path, rows, extra_cols=None):
    """Write CSV with smiles, name, IC50, round columns."""
    header = "smiles,name,IC50,round"
    lines = [header]
    for smi, name, ic50, rnd in rows:
        lines.append(f"{smi},{name},{ic50},{rnd}")
    p = tmp_path / "leads.csv"
    p.write_text("\n".join(lines))
    return str(p)


def run_lead(args):
    result = subprocess.run(
        [sys.executable, str(SCRIPTS_DIR / SCRIPT)] + args,
        capture_output=True, text=True, timeout=SCRIPT_TIMEOUT,
    )
    assert result.returncode == 0, result.stderr
    return json.loads(result.stdout)


def test_output_top_level_fields(tmp_path):
    csv = make_csv(tmp_path, SERIES)
    d = run_lead(["--csv", csv, "--activity-col", "IC50", "--unit", "nm"])
    for f in ["activity_column", "unit", "n_input", "n_computed", "n_errors",
              "thresholds_reference", "compounds", "round_summary", "errors"]:
        assert f in d


def test_all_compounds_computed(tmp_path):
    csv = make_csv(tmp_path, SERIES)
    d = run_lead(["--csv", csv, "--activity-col", "IC50", "--unit", "nm"])
    assert d["n_computed"] == len(SERIES)
    assert d["n_errors"] == 0


def test_compound_fields(tmp_path):
    csv = make_csv(tmp_path, SERIES)
    d = run_lead(["--csv", csv, "--activity-col", "IC50", "--unit", "nm"])
    for c in d["compounds"]:
        for f in ["name", "smiles", "activity_raw", "activity_unit", "metrics", "thresholds"]:
            assert f in c


def test_metrics_fields(tmp_path):
    csv = make_csv(tmp_path, SERIES)
    d = run_lead(["--csv", csv, "--activity-col", "IC50", "--unit", "nm"])
    for c in d["compounds"]:
        m = c["metrics"]
        for f in ["pIC50", "MW", "LogP", "TPSA", "n_heavy", "LE", "LLE", "BEI", "LELP"]:
            assert f in m


def test_le_positive(tmp_path):
    """LE must be positive for active compounds (IC50 < 1 mM → pIC50 > 0)."""
    csv = make_csv(tmp_path, SERIES)
    d = run_lead(["--csv", csv, "--activity-col", "IC50", "--unit", "nm"])
    for c in d["compounds"]:
        assert c["metrics"]["LE"] > 0


def test_lle_equals_pic50_minus_logp(tmp_path):
    csv = make_csv(tmp_path, SERIES)
    d = run_lead(["--csv", csv, "--activity-col", "IC50", "--unit", "nm"])
    for c in d["compounds"]:
        m = c["metrics"]
        expected_lle = round(m["pIC50"] - m["LogP"], 4)
        assert abs(m["LLE"] - expected_lle) < 1e-3


def test_round_summary_keys(tmp_path):
    csv = make_csv(tmp_path, SERIES)
    d = run_lead(["--csv", csv, "--activity-col", "IC50",
                  "--unit", "nm", "--round-col", "round"])
    assert "1" in d["round_summary"]
    assert "2" in d["round_summary"]


def test_round_summary_fields(tmp_path):
    csv = make_csv(tmp_path, SERIES)
    d = run_lead(["--csv", csv, "--activity-col", "IC50",
                  "--unit", "nm", "--round-col", "round"])
    for rnd_data in d["round_summary"].values():
        assert "mean_pIC50" in rnd_data
        assert "mean_LE" in rnd_data
        assert "n" in rnd_data


def test_thresholds_flags(tmp_path):
    csv = make_csv(tmp_path, SERIES)
    d = run_lead(["--csv", csv, "--activity-col", "IC50", "--unit", "nm"])
    for c in d["compounds"]:
        t = c["thresholds"]
        assert "LE_acceptable" in t
        assert "LLE_acceptable" in t
        assert isinstance(t["LE_acceptable"], bool)


def test_pic50_unit(tmp_path):
    """Compounds provided as pIC50 directly should compute same metrics."""
    lines = ["smiles,name,pIC50"]
    lines.append(f"{ASPIRIN_SMILES},aspirin,6.3")
    p = tmp_path / "pic50.csv"
    p.write_text("\n".join(lines))
    d = run_lead(["--csv", str(p), "--activity-col", "pIC50", "--unit", "pic50"])
    assert d["n_computed"] == 1
    assert abs(d["compounds"][0]["metrics"]["pIC50"] - 6.3) < 1e-3


def test_output_file(tmp_path):
    csv = make_csv(tmp_path, SERIES)
    out = tmp_path / "metrics.json"
    result = subprocess.run(
        [sys.executable, str(SCRIPTS_DIR / SCRIPT),
         "--csv", csv, "--activity-col", "IC50", "--unit", "nm", "--out", str(out)],
        capture_output=True, text=True, timeout=SCRIPT_TIMEOUT,
    )
    assert result.returncode == 0
    assert out.exists()
    d = json.loads(out.read_text())
    assert "compounds" in d
