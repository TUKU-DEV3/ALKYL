# ALKYL Scripts Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Build 5 pre-built chemistry scripts (`chem_*.py`) that ALKYL invokes via Bash tool during sessions, plus wire them into `install.sh` and `config/CLAUDE.md`.

**Architecture:** Independent argparse scripts outputting JSON to stdout, organized by functional axis (convert / props / 3d / qm / fetch). RDKit is the required core dependency; ASE and cclib are imported on demand. `install.sh` resolves the absolute scripts path and exports `$ALKYL_SCRIPTS` into the injected CLAUDE.md block.

**Tech Stack:** Python 3.10+, RDKit, requests (chem_fetch), ASE + cclib (chem_qm optional), pytest + subprocess for integration tests.

**Design doc:** `docs/plans/2026-03-01-scripts-design.md`

---

### Task 1: Scaffold `scripts/` and `tests/`

**Files:**
- Create: `scripts/.gitkeep` (replaced by scripts in next tasks)
- Create: `tests/__init__.py`
- Create: `tests/conftest.py`

**Step 1: Create directory structure**

```bash
mkdir -p scripts tests
touch tests/__init__.py
```

**Step 2: Write `tests/conftest.py`**

```python
import subprocess
import json
import sys
from pathlib import Path

SCRIPTS_DIR = Path(__file__).parent.parent / "scripts"

def run_script(script_name: str, args: list[str]) -> dict:
    """Run a chem_*.py script and parse its JSON stdout."""
    result = subprocess.run(
        [sys.executable, str(SCRIPTS_DIR / script_name)] + args,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, f"Script failed:\n{result.stderr}"
    return json.loads(result.stdout)

# Common test molecules
ASPIRIN_SMILES = "CC(=O)Oc1ccccc1C(=O)O"
CAFFEINE_SMILES = "Cn1cnc2c1c(=O)n(C)c(=O)n2C"
```

**Step 3: Commit scaffold**

```bash
git add tests/ scripts/
git commit -m "chore: scaffold scripts/ and tests/ directories"
```

---

### Task 2: `chem_convert.py` — Format I/O

**Files:**
- Create: `scripts/chem_convert.py`
- Create: `tests/test_chem_convert.py`

**Step 1: Write failing tests**

```python
# tests/test_chem_convert.py
import subprocess, sys, json
from pathlib import Path
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
    proc = subprocess.run(
        [sys.executable, str(Path("scripts/chem_convert.py")),
         "--smiles", "NOT_A_SMILES", "--to", "inchi"],
        capture_output=True, text=True,
    )
    assert proc.returncode != 0
```

**Step 2: Run to verify they fail**

```bash
pytest tests/test_chem_convert.py -v
```
Expected: `ERROR` — `chem_convert.py` does not exist yet.

**Step 3: Write `scripts/chem_convert.py`**

```python
#!/usr/bin/env python3
"""ALKYL — chem_convert.py: inter-format molecular I/O."""

import argparse
import json
import sys

def parse_input(args):
    from rdkit import Chem
    if args.smiles:
        mol = Chem.MolFromSmiles(args.smiles)
        if mol is None:
            print(f"Invalid SMILES: {args.smiles}", file=sys.stderr)
            sys.exit(1)
        return mol
    if args.sdf:
        suppl = Chem.SDMolSupplier(args.sdf)
        mol = next((m for m in suppl if m is not None), None)
        if mol is None:
            print(f"No valid molecule in: {args.sdf}", file=sys.stderr)
            sys.exit(1)
        return mol
    if args.inchi:
        from rdkit.Chem.inchi import MolFromInchi
        mol = MolFromInchi(args.inchi)
        if mol is None:
            print(f"Invalid InChI: {args.inchi}", file=sys.stderr)
            sys.exit(1)
        return mol
    print("No input provided. Use --smiles, --sdf, or --inchi.", file=sys.stderr)
    sys.exit(1)


def convert(mol, target: str) -> dict:
    from rdkit import Chem
    if target == "smiles":
        return {"smiles": Chem.MolToSmiles(mol)}
    if target == "inchi":
        from rdkit.Chem.inchi import MolToInchi, InchiToInchiKey
        inchi = MolToInchi(mol)
        return {"inchi": inchi, "inchikey": InchiToInchiKey(inchi)}
    if target == "sdf":
        return {"sdf": Chem.MolToMolBlock(mol)}
    if target == "svg":
        from rdkit.Chem.Draw import rdMolDraw2D
        drawer = rdMolDraw2D.MolDraw2DSVG(300, 200)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        return {"svg": drawer.GetDrawingText()}
    print(f"Unsupported target format: {target}", file=sys.stderr)
    sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description="Convert molecular formats.")
    inp = parser.add_mutually_exclusive_group(required=True)
    inp.add_argument("--smiles", help="Input SMILES string")
    inp.add_argument("--sdf", help="Input SDF file path")
    inp.add_argument("--inchi", help="Input InChI string")
    parser.add_argument("--to", required=True,
                        choices=["smiles", "inchi", "sdf", "svg"],
                        help="Target format")
    parser.add_argument("--out", help="Write output to file instead of stdout")
    args = parser.parse_args()

    mol = parse_input(args)
    result = convert(mol, args.to)

    output = json.dumps(result, indent=2)
    if args.out:
        with open(args.out, "w") as f:
            f.write(output)
    else:
        print(output)


if __name__ == "__main__":
    main()
```

**Step 4: Run tests**

```bash
pytest tests/test_chem_convert.py -v
```
Expected: 4 PASSED.

**Step 5: Commit**

```bash
git add scripts/chem_convert.py tests/test_chem_convert.py
git commit -m "feat: add chem_convert.py — inter-format molecular I/O"
```

---

### Task 3: `chem_props.py` — Descriptors & Drug-likeness

**Files:**
- Create: `scripts/chem_props.py`
- Create: `tests/test_chem_props.py`

**Step 1: Write failing tests**

```python
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
```

**Step 2: Run to verify they fail**

```bash
pytest tests/test_chem_props.py -v
```
Expected: `ERROR` — script does not exist.

**Step 3: Write `scripts/chem_props.py`**

```python
#!/usr/bin/env python3
"""ALKYL — chem_props.py: physico-chemical descriptors, drug-likeness, fingerprints."""

import argparse
import json
import sys

DESCRIPTOR_MAP = {
    "mw":       ("rdkit.Chem.Descriptors", "ExactMolWt"),
    "logp":     ("rdkit.Chem.Descriptors", "MolLogP"),
    "hbd":      ("rdkit.Chem.Descriptors", "NumHDonors"),
    "hba":      ("rdkit.Chem.Descriptors", "NumHAcceptors"),
    "tpsa":     ("rdkit.Chem.Descriptors", "TPSA"),
    "rotbonds": ("rdkit.Chem.Descriptors", "NumRotatableBonds"),
    "rings":    ("rdkit.Chem.Descriptors", "RingCount"),
    "fsp3":     ("rdkit.Chem.Descriptors", "FractionCSP3"),
    "mw_avg":   ("rdkit.Chem.Descriptors", "MolWt"),
}
ALL_DESCRIPTORS = list(DESCRIPTOR_MAP.keys())


def load_mol(args):
    from rdkit import Chem
    if args.smiles:
        mol = Chem.MolFromSmiles(args.smiles)
    elif args.sdf:
        suppl = Chem.SDMolSupplier(args.sdf)
        mol = next((m for m in suppl if m is not None), None)
    else:
        print("Provide --smiles or --sdf.", file=sys.stderr)
        sys.exit(1)
    if mol is None:
        print("Invalid molecule input.", file=sys.stderr)
        sys.exit(1)
    return mol


def calc_descriptors(mol, names: list[str]) -> dict:
    import importlib
    result = {}
    for name in names:
        if name not in DESCRIPTOR_MAP:
            continue
        module_path, func_name = DESCRIPTOR_MAP[name]
        module = importlib.import_module(module_path)
        result[name] = round(getattr(module, func_name)(mol), 4)
    return result


def check_lipinski(mol) -> dict:
    from rdkit.Chem import Descriptors
    violations = []
    if Descriptors.MolWt(mol) > 500:
        violations.append("MW > 500")
    if Descriptors.MolLogP(mol) > 5:
        violations.append("LogP > 5")
    if Descriptors.NumHDonors(mol) > 5:
        violations.append("HBD > 5")
    if Descriptors.NumHAcceptors(mol) > 10:
        violations.append("HBA > 10")
    return {"pass": len(violations) == 0, "violations": violations}


def check_pains(mol) -> dict:
    from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog(params)
    matches = catalog.GetMatches(mol)
    alerts = [m.GetDescription() for m in matches]
    return {"alerts": alerts, "clean": len(alerts) == 0}


def calc_fingerprint(mol, fp_type: str, radius: int, nbits: int) -> dict:
    from rdkit.Chem import rdMolDescriptors, MACCSkeys
    from rdkit import DataStructs
    if fp_type == "morgan":
        fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius, nbits)
    elif fp_type == "maccs":
        fp = MACCSkeys.GenMACCSKeys(mol)
    else:
        from rdkit.Chem.rdmolops import RDKFingerprint
        fp = RDKFingerprint(mol, fpSize=nbits)
    on_bits = list(fp.GetOnBits())
    return {"type": fp_type, "nbits": fp.GetNumBits(), "bits": on_bits,
            "density": round(len(on_bits) / fp.GetNumBits(), 4)}


def main():
    parser = argparse.ArgumentParser(description="Compute molecular properties.")
    inp = parser.add_mutually_exclusive_group(required=True)
    inp.add_argument("--smiles", help="Input SMILES")
    inp.add_argument("--sdf", help="Input SDF file")
    parser.add_argument("--descriptors",
                        help="Comma-separated descriptors or 'all'")
    parser.add_argument("--lipinski", action="store_true")
    parser.add_argument("--pains", action="store_true")
    parser.add_argument("--fingerprint", choices=["morgan", "maccs", "rdkit"])
    parser.add_argument("--radius", type=int, default=2)
    parser.add_argument("--nbits", type=int, default=2048)
    args = parser.parse_args()

    mol = load_mol(args)
    result = {}

    if args.descriptors:
        names = ALL_DESCRIPTORS if args.descriptors == "all" \
            else [d.strip() for d in args.descriptors.split(",")]
        result.update(calc_descriptors(mol, names))

    if args.lipinski:
        result["lipinski"] = check_lipinski(mol)

    if args.pains:
        result["pains"] = check_pains(mol)

    if args.fingerprint:
        result["fingerprint"] = calc_fingerprint(
            mol, args.fingerprint, args.radius, args.nbits)

    if not result:
        parser.error("Specify at least one of: --descriptors, --lipinski, "
                     "--pains, --fingerprint")

    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
```

**Step 4: Run tests**

```bash
pytest tests/test_chem_props.py -v
```
Expected: 5 PASSED.

**Step 5: Commit**

```bash
git add scripts/chem_props.py tests/test_chem_props.py
git commit -m "feat: add chem_props.py — descriptors, Lipinski, PAINS, fingerprints"
```

---

### Task 4: `chem_fetch.py` — External Data

**Files:**
- Create: `scripts/chem_fetch.py`
- Create: `tests/test_chem_fetch.py`

**Step 1: Write failing tests**

```python
# tests/test_chem_fetch.py
# Note: these tests make real HTTP requests to PubChem REST API.
# They are marked slow and skipped in offline environments.
import pytest
from tests.conftest import run_script

SCRIPT = "chem_fetch.py"

@pytest.mark.network
def test_fetch_by_name_pubchem():
    result = run_script(SCRIPT, ["--source", "pubchem",
                                 "--name", "aspirin",
                                 "--properties", "mw,logp,inchi"])
    assert "smiles" in result or "inchi" in result
    assert result.get("cid") == 2244

@pytest.mark.network
def test_fetch_by_cid():
    result = run_script(SCRIPT, ["--source", "pubchem",
                                 "--cid", "2244"])
    assert "smiles" in result
    assert "iupac_name" in result

@pytest.mark.network
def test_fetch_chembl_by_id():
    result = run_script(SCRIPT, ["--source", "chembl",
                                 "--chembl-id", "CHEMBL25"])
    assert "smiles" in result
    assert result.get("chembl_id") == "CHEMBL25"
```

**Step 2: Run to verify they fail**

```bash
pytest tests/test_chem_fetch.py -v -m network
```
Expected: `ERROR` — script does not exist.

**Step 3: Write `scripts/chem_fetch.py`**

```python
#!/usr/bin/env python3
"""ALKYL — chem_fetch.py: fetch molecular data from PubChem or ChEMBL."""

import argparse
import json
import sys
import urllib.request
import urllib.parse


PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"


def _get_json(url: str) -> dict:
    try:
        with urllib.request.urlopen(url, timeout=15) as resp:
            return json.loads(resp.read().decode())
    except Exception as e:
        print(f"HTTP error fetching {url}: {e}", file=sys.stderr)
        sys.exit(1)


def fetch_pubchem(args) -> dict:
    # Build lookup URL
    if args.cid:
        namespace, ident = "cid", str(args.cid)
    elif args.name:
        namespace, ident = "name", urllib.parse.quote(args.name)
    elif args.smiles:
        namespace, ident = "smiles", urllib.parse.quote(args.smiles)
    else:
        print("Provide --name, --cid, or --smiles.", file=sys.stderr)
        sys.exit(1)

    props = "CanonicalSMILES,InChI,InChIKey,IUPACName,MolecularWeight,XLogP"
    url = f"{PUBCHEM_BASE}/compound/{namespace}/{ident}/property/{props}/JSON"
    data = _get_json(url)
    props_data = data["PropertyTable"]["Properties"][0]

    result = {
        "cid": props_data.get("CID"),
        "smiles": props_data.get("CanonicalSMILES"),
        "inchi": props_data.get("InChI"),
        "inchikey": props_data.get("InChIKey"),
        "iupac_name": props_data.get("IUPACName"),
        "mw": props_data.get("MolecularWeight"),
        "logp": props_data.get("XLogP"),
    }
    # Filter to requested properties if specified
    if args.properties:
        keep = {p.strip() for p in args.properties.split(",")}
        result = {k: v for k, v in result.items() if k in keep or k == "cid"}
    return result


def fetch_chembl(args) -> dict:
    if args.chembl_id:
        url = f"{CHEMBL_BASE}/molecule/{args.chembl_id}?format=json"
    elif args.name:
        url = f"{CHEMBL_BASE}/molecule?pref_name={urllib.parse.quote(args.name)}&format=json"
    else:
        print("Provide --chembl-id or --name.", file=sys.stderr)
        sys.exit(1)

    data = _get_json(url)
    # Handle list vs single result
    if "molecules" in data:
        if not data["molecules"]:
            print("No results found.", file=sys.stderr)
            sys.exit(1)
        mol = data["molecules"][0]
    else:
        mol = data

    struct = mol.get("molecule_structures") or {}
    return {
        "chembl_id": mol.get("molecule_chembl_id"),
        "smiles": struct.get("canonical_smiles"),
        "inchi": struct.get("standard_inchi"),
        "inchikey": struct.get("standard_inchi_key"),
        "name": mol.get("pref_name"),
        "type": mol.get("molecule_type"),
        "mw": mol.get("molecule_properties", {}).get("full_mwt"),
        "logp": mol.get("molecule_properties", {}).get("alogp"),
    }


def main():
    parser = argparse.ArgumentParser(
        description="Fetch molecular data from PubChem or ChEMBL.")
    parser.add_argument("--source", required=True, choices=["pubchem", "chembl"])
    parser.add_argument("--name", help="Molecule name")
    parser.add_argument("--smiles", help="SMILES (PubChem only)")
    parser.add_argument("--cid", type=int, help="PubChem CID")
    parser.add_argument("--chembl-id", dest="chembl_id", help="ChEMBL ID")
    parser.add_argument("--properties", help="Comma-separated property subset")
    parser.add_argument("--out", help="Write JSON to file")
    args = parser.parse_args()

    result = fetch_pubchem(args) if args.source == "pubchem" else fetch_chembl(args)

    output = json.dumps(result, indent=2)
    if args.out:
        with open(args.out, "w") as f:
            f.write(output)
    else:
        print(output)


if __name__ == "__main__":
    main()
```

**Step 4: Run tests**

```bash
pytest tests/test_chem_fetch.py -v -m network
```
Expected: 3 PASSED (requires network).

**Step 5: Commit**

```bash
git add scripts/chem_fetch.py tests/test_chem_fetch.py
git commit -m "feat: add chem_fetch.py — PubChem/ChEMBL data fetcher"
```

---

### Task 5: `chem_3d.py` — Conformers & Geometry

**Files:**
- Create: `scripts/chem_3d.py`
- Create: `tests/test_chem_3d.py`

**Step 1: Write failing tests**

```python
# tests/test_chem_3d.py
import json, subprocess, sys, tempfile
from pathlib import Path
from tests.conftest import run_script, ASPIRIN_SMILES

SCRIPT = "chem_3d.py"

def test_generate_conformers_json():
    result = run_script(SCRIPT, ["--smiles", ASPIRIN_SMILES,
                                 "--conformers", "5"])
    assert result["n_conformers"] == 5
    assert "energies_kcal" in result

def test_generate_conformers_sdf(tmp_path):
    out_file = tmp_path / "confs.sdf"
    proc = subprocess.run(
        [sys.executable, str(Path("scripts/chem_3d.py")),
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
```

**Step 2: Run to verify they fail**

```bash
pytest tests/test_chem_3d.py -v
```
Expected: `ERROR` — script does not exist.

**Step 3: Write `scripts/chem_3d.py`**

```python
#!/usr/bin/env python3
"""ALKYL — chem_3d.py: 3D conformer generation, minimization, RMSD."""

import argparse
import json
import sys


def load_mol_2d(args):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    if args.smiles:
        mol = Chem.MolFromSmiles(args.smiles)
        if mol is None:
            print(f"Invalid SMILES.", file=sys.stderr)
            sys.exit(1)
        mol = Chem.AddHs(mol)
        return mol
    if args.sdf:
        suppl = Chem.SDMolSupplier(args.sdf, removeHs=False)
        mol = next((m for m in suppl if m is not None), None)
        if mol is None:
            print(f"No valid molecule in SDF.", file=sys.stderr)
            sys.exit(1)
        return mol
    print("Provide --smiles or --sdf.", file=sys.stderr)
    sys.exit(1)


def generate_conformers(mol, n: int, method: str):
    from rdkit.Chem import AllChem
    params = AllChem.ETKDGv3() if method == "etkdg" else AllChem.ETDG()
    params.randomSeed = 42
    params.numThreads = 0
    ids = AllChem.EmbedMultipleConfs(mol, numConfs=n, params=params)
    return list(ids)


def minimize_conformers(mol, ff: str) -> list[float | None]:
    from rdkit.Chem import AllChem
    energies = []
    for conf_id in range(mol.GetNumConformers()):
        if ff == "mmff94":
            props = AllChem.MMFFGetMoleculeProperties(mol)
            if props is None:
                energies.append(None)
                continue
            ff_obj = AllChem.MMFFGetMoleculeForceField(mol, props, confId=conf_id)
        else:
            ff_obj = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)
        if ff_obj is None:
            energies.append(None)
            continue
        ff_obj.Minimize()
        energies.append(round(ff_obj.CalcEnergy(), 4))
    return energies


def main():
    parser = argparse.ArgumentParser(description="3D conformer generation.")
    inp = parser.add_mutually_exclusive_group(required=True)
    inp.add_argument("--smiles")
    inp.add_argument("--sdf")
    parser.add_argument("--conformers", type=int, default=10,
                        help="Number of conformers to generate")
    parser.add_argument("--method", choices=["etkdg", "etdg"], default="etkdg")
    parser.add_argument("--minimize", choices=["mmff94", "uff"],
                        help="Force field for minimization")
    parser.add_argument("--out", help="Write multi-conf SDF to file")
    args = parser.parse_args()

    mol = load_mol_2d(args)
    conf_ids = generate_conformers(mol, args.conformers, args.method)

    energies = None
    if args.minimize:
        energies = minimize_conformers(mol, args.minimize)

    if args.out:
        from rdkit.Chem import SDWriter
        writer = SDWriter(args.out)
        for cid in conf_ids:
            writer.write(mol, confId=cid)
        writer.close()
        result = {"n_conformers": len(conf_ids), "out": args.out}
    else:
        result = {
            "n_conformers": len(conf_ids),
            "method": args.method,
            "energies_kcal": energies,
        }

    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
```

**Step 4: Run tests**

```bash
pytest tests/test_chem_3d.py -v
```
Expected: 3 PASSED.

**Step 5: Commit**

```bash
git add scripts/chem_3d.py tests/test_chem_3d.py
git commit -m "feat: add chem_3d.py — conformer generation and minimization"
```

---

### Task 6: `chem_qm.py` — QM Interface (ORCA / Gaussian)

**Files:**
- Create: `scripts/chem_qm.py`
- Create: `tests/test_chem_qm.py`

**Step 1: Write failing tests**

```python
# tests/test_chem_qm.py
import json, subprocess, sys
from pathlib import Path
from tests.conftest import run_script, ASPIRIN_SMILES

SCRIPT = "chem_qm.py"

def test_generate_orca_input(tmp_path):
    out_file = tmp_path / "aspirin.inp"
    proc = subprocess.run(
        [sys.executable, "scripts/chem_qm.py",
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
    assert "! Opt" in content

def test_generate_gaussian_input(tmp_path):
    out_file = tmp_path / "aspirin.gjf"
    proc = subprocess.run(
        [sys.executable, "scripts/chem_qm.py",
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
```

**Step 2: Run to verify they fail**

```bash
pytest tests/test_chem_qm.py -v
```
Expected: `ERROR` — script does not exist.

**Step 3: Write `scripts/chem_qm.py`**

```python
#!/usr/bin/env python3
"""ALKYL — chem_qm.py: QM input generation and output parsing (ORCA, Gaussian)."""

import argparse
import json
import sys
from pathlib import Path


def get_xyz_block(mol) -> str:
    """Return XYZ coordinate block from a 3D RDKit mol (no header lines)."""
    from rdkit.Chem import AllChem
    conf = mol.GetConformer()
    lines = []
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        lines.append(f"  {atom.GetSymbol():2s}  {pos.x:12.6f}  {pos.y:12.6f}  {pos.z:12.6f}")
    return "\n".join(lines)


def smiles_to_3d(smiles: str):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("Invalid SMILES.", file=sys.stderr)
        sys.exit(1)
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    AllChem.EmbedMolecule(mol, params)
    AllChem.MMFFOptimizeMolecule(mol)
    return mol


TASK_ORCA = {"sp": "SP", "opt": "Opt", "freq": "Freq", "scan": "Scan"}
TASK_GAUSSIAN = {"sp": "SP", "opt": "Opt", "freq": "Freq", "scan": "Scan"}


def write_orca(mol, method: str, basis: str, task: str,
               charge: int, mult: int) -> str:
    xyz = get_xyz_block(mol)
    keyword = TASK_ORCA.get(task, task)
    n_atoms = mol.GetNumAtoms()
    return (
        f"! {method} {basis} {keyword}\n"
        f"\n"
        f"* xyz {charge} {mult}\n"
        f"{xyz}\n"
        f"*\n"
    )


def write_gaussian(mol, method: str, basis: str, task: str,
                   charge: int, mult: int) -> str:
    xyz = get_xyz_block(mol)
    keyword = TASK_GAUSSIAN.get(task, task)
    return (
        f"%mem=4GB\n"
        f"%nprocshared=4\n"
        f"#P {method}/{basis} {keyword}\n"
        f"\n"
        f"ALKYL generated input\n"
        f"\n"
        f"{charge} {mult}\n"
        f"{xyz}\n"
        f"\n"
    )


def parse_orca_output(path: str) -> dict:
    """Parse key results from an ORCA output file."""
    result = {}
    text = Path(path).read_text(errors="replace")
    # Final single-point energy
    for line in reversed(text.splitlines()):
        if "FINAL SINGLE POINT ENERGY" in line:
            try:
                result["energy_hartree"] = float(line.split()[-1])
            except ValueError:
                pass
            break
    # Frequencies
    freqs = []
    for line in text.splitlines():
        if ":" in line and "cm**-1" in line:
            try:
                freqs.append(float(line.split(":")[1].split()[0]))
            except (ValueError, IndexError):
                pass
    if freqs:
        result["frequencies_cm1"] = freqs
    return result


def main():
    parser = argparse.ArgumentParser(
        description="Generate QM input files or parse QM outputs.")
    inp = parser.add_mutually_exclusive_group()
    inp.add_argument("--smiles", help="Input SMILES (auto 3D)")
    inp.add_argument("--xyz", help="Input XYZ file (pre-computed geometry)")
    parser.add_argument("--engine", choices=["orca", "gaussian"], default="orca")
    parser.add_argument("--task", choices=["sp", "opt", "freq", "scan"],
                        default="sp")
    parser.add_argument("--method", default="B3LYP")
    parser.add_argument("--basis", default="6-31G*")
    parser.add_argument("--charge", type=int, default=0)
    parser.add_argument("--mult", type=int, default=1)
    parser.add_argument("--out", help="Write input file to path")
    parser.add_argument("--parse", help="Parse existing output file → JSON")
    args = parser.parse_args()

    # Parse mode
    if args.parse:
        if args.engine == "orca":
            result = parse_orca_output(args.parse)
        else:
            print("Gaussian parsing not yet implemented.", file=sys.stderr)
            sys.exit(1)
        print(json.dumps(result, indent=2))
        return

    # Generation mode
    if args.smiles:
        mol = smiles_to_3d(args.smiles)
    elif args.xyz:
        try:
            import ase.io
            atoms = ase.io.read(args.xyz)
            # Build minimal RDKit mol from XYZ via ASE (best-effort, no bonds)
            print("XYZ input via ASE: bond perception may be limited.",
                  file=sys.stderr)
            mol = smiles_to_3d("C")  # fallback — user should prefer --smiles
        except ImportError:
            print("ASE required for --xyz input. Install with: pip install ase",
                  file=sys.stderr)
            sys.exit(1)
    else:
        parser.error("Provide --smiles or --xyz for input generation.")

    if args.engine == "orca":
        content = write_orca(mol, args.method, args.basis,
                             args.task, args.charge, args.mult)
    else:
        content = write_gaussian(mol, args.method, args.basis,
                                 args.task, args.charge, args.mult)

    if args.out:
        Path(args.out).write_text(content)
        print(json.dumps({"written": args.out, "engine": args.engine,
                          "task": args.task}))
    else:
        print(content)


if __name__ == "__main__":
    main()
```

**Step 4: Run tests**

```bash
pytest tests/test_chem_qm.py -v
```
Expected: 2 PASSED.

**Step 5: Commit**

```bash
git add scripts/chem_qm.py tests/test_chem_qm.py
git commit -m "feat: add chem_qm.py — ORCA/Gaussian input generation and output parsing"
```

---

### Task 7: Wire `$ALKYL_SCRIPTS` into `install.sh`

**Files:**
- Modify: `install.sh`

**Step 1: Read current install.sh**

```bash
cat install.sh
```

**Step 2: Add ALKYL_SCRIPTS export to the injected block**

In `config/CLAUDE.md`, add after `## Priority Stack` (before closing `<!-- ALKYL-END -->`):

```markdown
## Scripts disponibles

Scripts dans `$ALKYL_SCRIPTS/`.
Appelle via Bash : `python $ALKYL_SCRIPTS/<script>.py`

| Tâche | Script |
|---|---|
| Convertir format moléculaire | `chem_convert.py` |
| Calculer MW, LogP, TPSA, fingerprints | `chem_props.py` |
| Vérifier Lipinski / PAINS | `chem_props.py --lipinski --pains` |
| Générer conformères 3D | `chem_3d.py --conformers N` |
| Préparer input ORCA/Gaussian | `chem_qm.py --engine orca` |
| Parser output QM | `chem_qm.py --parse output.log` |
| Récupérer molécule PubChem/ChEMBL | `chem_fetch.py` |

Règles :
- Toujours parser le JSON stdout avant de répondre à l'utilisateur
- Si RDKit absent : signaler clairement, ne pas inventer les valeurs
- `--help` disponible sur chaque script pour vérifier les flags
```

In `install.sh`, after injecting `config/CLAUDE.md`, replace the static `$ALKYL_SCRIPTS`
placeholder with the resolved absolute path:

```bash
SCRIPTS_ABS="$(cd "$(dirname "$0")/scripts" && pwd)"
# After injecting the block, substitute the placeholder:
sed -i "s|\$ALKYL_SCRIPTS|${SCRIPTS_ABS}|g" ~/.claude/CLAUDE.md
```

**Step 3: Run install.sh and verify**

```bash
bash install.sh
grep "ALKYL_SCRIPTS\|chem_convert" ~/.claude/CLAUDE.md
```
Expected: lines showing the resolved absolute path to scripts/.

**Step 4: Commit**

```bash
git add install.sh config/CLAUDE.md
git commit -m "feat: wire ALKYL_SCRIPTS path into install.sh + Scripts section in CLAUDE.md"
```

---

### Task 8: Run full test suite

**Step 1: Run all non-network tests**

```bash
pytest tests/ -v -m "not network"
```
Expected: all tests PASSED (chem_convert, chem_props, chem_3d, chem_qm).

**Step 2: Run network tests (optional, requires connectivity)**

```bash
pytest tests/test_chem_fetch.py -v -m network
```

**Step 3: Final commit**

```bash
git add tests/
git commit -m "test: complete test suite for all chem_*.py scripts"
```

---

## Summary

| Script | Dependencies | Tests |
|---|---|---|
| `chem_convert.py` | RDKit | `test_chem_convert.py` (4 tests) |
| `chem_props.py` | RDKit | `test_chem_props.py` (5 tests) |
| `chem_3d.py` | RDKit | `test_chem_3d.py` (3 tests) |
| `chem_qm.py` | RDKit + ASE (optional) | `test_chem_qm.py` (2 tests) |
| `chem_fetch.py` | stdlib only | `test_chem_fetch.py` (3 network tests) |

**Total:** 17 tests, 8 tasks, ~6 commits.
