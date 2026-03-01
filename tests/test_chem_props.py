# tests/test_chem_props.py
from tests.conftest import run_script, ASPIRIN_SMILES, CAFFEINE_SMILES

SCRIPT = "chem_props.py"

def test_basic_descriptors():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN_SMILES,
                                 "--descriptors", "mw,logp,hbd,hba,tpsa"])
    assert abs(result["mw"] - 180.16) < 0.1
    assert "logp" in result
    assert result["hbd"] == 1
    assert result["hba"] == 3

def test_lipinski_pass():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN_SMILES, "--lipinski"])
    assert result["lipinski"]["pass"] is True
    assert result["lipinski"]["violations"] == []

def test_pains_clean():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN_SMILES, "--pains"])
    assert "pains" in result
    assert isinstance(result["pains"]["alerts"], list)

def test_morgan_fingerprint():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN_SMILES,
                                 "--fingerprint", "morgan", "--radius", "2",
                                 "--nbits", "1024"])
    assert "fingerprint" in result
    assert result["fingerprint"]["type"] == "morgan"
    assert len(result["fingerprint"]["bits"]) > 0

def test_all_descriptors():
    result = run_script(SCRIPT, ["--smiles", CAFFEINE_SMILES,
                                 "--descriptors", "all"])
    for key in ["mw", "logp", "hbd", "hba", "tpsa", "rotbonds", "rings"]:
        assert key in result, f"Missing descriptor: {key}"
