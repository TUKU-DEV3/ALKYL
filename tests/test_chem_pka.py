# tests/test_chem_pka.py
from tests.conftest import run_script

SCRIPT = "chem_pka.py"
ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"   # carboxylic acid pKa ~3.5 → ionized at pH 7.4
ANILINE = "Nc1ccccc1"                # aniline pKaH ~4.6 → neutral amine at pH 7.4


def test_output_fields():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN, "--ph", "7.4"])
    for f in ["input_smiles", "ph", "dominant_form",
              "ionization_groups", "net_charge_at_ph"]:
        assert f in result, f"Missing field: {f}"


def test_aspirin_ionized_at_physiological_ph():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN, "--ph", "7.4"])
    assert result["net_charge_at_ph"] < 0, (
        f"Expected negative charge at pH 7.4, got {result['net_charge_at_ph']}")
    assert "[O-]" in result["dominant_form"], (
        f"Expected carboxylate [O-] in dominant form: {result['dominant_form']}")


def test_aspirin_neutral_at_low_ph():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN, "--ph", "2.0"])
    assert result["net_charge_at_ph"] == 0, (
        f"Expected neutral charge at pH 2.0, got {result['net_charge_at_ph']}")


def test_aniline_neutral_at_physiological():
    # aniline pKaH ~4.6 → amine fully deprotonated (neutral) at pH 7.4
    result = run_script(SCRIPT, ["--smiles", ANILINE, "--ph", "7.4"])
    assert result["net_charge_at_ph"] == 0, (
        f"Expected aniline neutral at pH 7.4, got {result['net_charge_at_ph']}")


def test_ionization_groups_list():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN, "--ph", "7.4"])
    assert isinstance(result["ionization_groups"], list)
    types = [g["type"] for g in result["ionization_groups"]]
    assert "carboxylic_acid" in types, f"carboxylic_acid not found in types: {types}"


def test_ionization_group_fields():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN, "--ph", "7.4"])
    for group in result["ionization_groups"]:
        for key in ["type", "pka_range", "count", "state_at_ph", "charge_contribution"]:
            assert key in group, f"Missing key '{key}' in group: {group}"


def test_carboxylic_acid_state_at_high_ph():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN, "--ph", "7.4"])
    ca_groups = [g for g in result["ionization_groups"]
                 if g["type"] == "carboxylic_acid"]
    assert len(ca_groups) >= 1
    assert ca_groups[0]["state_at_ph"] == "ionized"


def test_carboxylic_acid_state_at_low_ph():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN, "--ph", "2.0"])
    ca_groups = [g for g in result["ionization_groups"]
                 if g["type"] == "carboxylic_acid"]
    assert len(ca_groups) >= 1
    assert ca_groups[0]["state_at_ph"] == "neutral"


def test_ph_stored_in_output():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN, "--ph", "5.5"])
    assert result["ph"] == 5.5


def test_default_ph_is_7_4():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN])
    assert result["ph"] == 7.4
