# tests/test_chem_cluster.py
import json
import subprocess
import sys
from tests.conftest import ASPIRIN_SMILES, CAFFEINE_SMILES, SCRIPTS_DIR, SCRIPT_TIMEOUT

SCRIPT = "chem_cluster.py"

PHENYL_LIB = [
    ("c1ccccc1CC(=O)O", "phenylacetic"),
    ("c1ccc(F)cc1CC(=O)O", "fluoro_phenylacetic"),
    ("c1ccc(Cl)cc1CC(=O)O", "chloro_phenylacetic"),
    (ASPIRIN_SMILES, "aspirin"),
    (CAFFEINE_SMILES, "caffeine"),
]


def make_smi(tmp_path, mols):
    p = tmp_path / "lib.smi"
    p.write_text("\n".join(f"{s} {n}" for s, n in mols))
    return str(p)


def run_cluster(args):
    result = subprocess.run(
        [sys.executable, str(SCRIPTS_DIR / SCRIPT)] + args,
        capture_output=True, text=True, timeout=SCRIPT_TIMEOUT,
    )
    assert result.returncode == 0, result.stderr
    return json.loads(result.stdout)


def test_output_top_level_fields(tmp_path):
    lib = make_smi(tmp_path, PHENYL_LIB)
    d = run_cluster(["--input", lib, "--cutoff", "0.5"])
    for f in ["n_molecules", "n_clusters", "n_singletons", "cutoff_distance",
              "similarity_threshold", "fingerprint", "clusters"]:
        assert f in d


def test_cluster_item_fields(tmp_path):
    lib = make_smi(tmp_path, PHENYL_LIB)
    d = run_cluster(["--input", lib, "--cutoff", "0.5"])
    for c in d["clusters"]:
        assert "cluster_id" in c
        assert "size" in c
        assert "centroid" in c
        assert "members" in c
        for k in ["pos", "idx", "name", "smiles"]:
            assert k in c["centroid"]


def test_total_members_equals_n_molecules(tmp_path):
    lib = make_smi(tmp_path, PHENYL_LIB)
    d = run_cluster(["--input", lib, "--cutoff", "0.5"])
    total = sum(c["size"] for c in d["clusters"])
    assert total == d["n_molecules"]


def test_tight_cutoff_produces_more_clusters(tmp_path):
    """Lower distance cutoff (tighter) → more clusters. Higher cutoff (looser) → fewer clusters."""
    lib = make_smi(tmp_path, PHENYL_LIB)
    d_tight = run_cluster(["--input", lib, "--cutoff", "0.2"])  # distance ≤ 0.2, sim ≥ 0.8
    d_loose = run_cluster(["--input", lib, "--cutoff", "0.7"])  # distance ≤ 0.7, sim ≥ 0.3
    assert d_tight["n_clusters"] >= d_loose["n_clusters"]


def test_singletons_count(tmp_path):
    """All-different molecules at tight cutoff → all singletons."""
    lib = make_smi(tmp_path, [
        (ASPIRIN_SMILES, "aspirin"),
        (CAFFEINE_SMILES, "caffeine"),
    ])
    d = run_cluster(["--input", lib, "--cutoff", "0.01"])
    assert d["n_singletons"] == d["n_clusters"]


def test_maccs_fingerprint(tmp_path):
    lib = make_smi(tmp_path, PHENYL_LIB)
    d = run_cluster(["--input", lib, "--cutoff", "0.4", "--fingerprint", "maccs"])
    assert d["fingerprint"] == "maccs"
    assert d["n_clusters"] >= 1


def test_output_file(tmp_path):
    lib = make_smi(tmp_path, PHENYL_LIB)
    out = tmp_path / "clusters.json"
    result = subprocess.run(
        [sys.executable, str(SCRIPTS_DIR / SCRIPT),
         "--input", lib, "--cutoff", "0.5", "--out", str(out)],
        capture_output=True, text=True, timeout=SCRIPT_TIMEOUT,
    )
    assert result.returncode == 0
    assert out.exists()
    d = json.loads(out.read_text())
    assert "clusters" in d
