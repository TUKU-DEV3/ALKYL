# Docking — Protein Preparation

Full workflow: raw PDB → docking-ready PDBQT receptor.

---

## Step 1 — Download and Inspect

```python
import urllib.request
from pathlib import Path

def fetch_pdb(pdb_id: str, out_dir: str = ".") -> Path:
    """Download PDB structure."""
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    out = Path(out_dir) / f"{pdb_id.upper()}.pdb"
    urllib.request.urlretrieve(url, out)
    return out

# Example: CDK2 with inhibitor
pdb_path = fetch_pdb("1HCL")

# Or use ALKYL script to fetch
# python scripts/chem_fetch.py --source pdb --id 1HCL
```

**Inspect first:**
```bash
grep "^REMARK" 1HCL.pdb | head -30    # resolution, R-factor, ligands
grep "^HET " 1HCL.pdb                 # heteroatoms (ligands, cofactors, metals)
grep "^SEQRES" 1HCL.pdb | head -5     # sequence coverage
grep "^MISSING" 1HCL.pdb              # missing residues/atoms
```

**Key decisions before prep:**
- Keep or remove: co-crystallized ligand (keep as reference for box), water molecules (case-by-case), cofactors (always keep if catalytic)
- Note metal ions (Zn, Mg, Ca) — require special treatment

---

## Step 2 — Fix Structure with pdbfixer

```python
from pdbfixer import PDBFixer
from openmm.app import PDBFile

def fix_pdb(input_pdb: str, output_pdb: str,
            ph: float = 7.4,
            keep_water: bool = False) -> None:
    """
    Fix PDB: add missing residues, atoms, hydrogens.
    """
    fixer = PDBFixer(filename=input_pdb)

    # 1. Find and add missing residues (loops)
    fixer.findMissingResidues()
    # Remove terminal missing residues (often disordered tails)
    chains = list(fixer.topology.chains())
    keys = fixer.missingResidues.keys()
    terminal_keys = [k for k in keys
                     if k[1] == 0 or
                     k[1] == len(list(chains[k[0]].residues()))]
    for k in terminal_keys:
        del fixer.missingResidues[k]

    # 2. Find and add missing heavy atoms
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    # 3. Remove water (optional)
    if not keep_water:
        fixer.removeHeterogens(keepWater=False)
    else:
        fixer.removeHeterogens(keepWater=True)

    # 4. Add hydrogens at specified pH
    fixer.addMissingHydrogens(ph)

    # Write
    with open(output_pdb, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)

    print(f"Fixed PDB written to {output_pdb}")

fix_pdb("1HCL.pdb", "1HCL_fixed.pdb", ph=7.4, keep_water=False)
```

**When to keep water:** Crystallographic waters in the binding site that mediate key H-bonds (check publications). Use `keep_water=True` + manually curate after.

---

## Step 3 — Protonation State Assignment (propka)

pdbfixer assigns H at a uniform pH. For accurate protonation of histidine (HID/HIE/HIP), aspartate, glutamate, and lysine:

```bash
# Run propka on the fixed structure
propka3 1HCL_fixed.pdb

# Read 1HCL_fixed.pka — check pKa values for binding site residues
# Key residues with shifted pKa in active sites:
#   His catalytic: often HIP (doubly protonated) or specific tautomer
#   Asp/Glu catalytic: may be protonated (neutral) at pH 7.4 if pKa > 6
#   Lys: typically NH3+ at pH 7.4
```

```python
import subprocess, re

def read_propka_results(pka_file: str) -> dict:
    """Parse propka .pka output → dict {residue: pKa}."""
    results = {}
    with open(pka_file) as f:
        in_summary = False
        for line in f:
            if "SUMMARY OF THIS PREDICTION" in line:
                in_summary = True
                continue
            if in_summary and line.strip():
                parts = line.split()
                if len(parts) >= 3:
                    try:
                        res_id = f"{parts[0]}{parts[1]}{parts[2]}"
                        pka = float(parts[3])
                        results[res_id] = pka
                    except (ValueError, IndexError):
                        pass
    return results

pka_data = read_propka_results("1HCL_fixed.pka")
print("Shifted pKa residues:")
for res, pka in pka_data.items():
    if pka < 6.0 or pka > 9.0:
        print(f"  {res}: pKa = {pka:.1f}")
```

**Histidine tautomers** — must decide manually:
- HID: H on Nδ (more common)
- HIE: H on Nε
- HIP: both N protonated (positive charge)

```bash
# In the fixed PDB, rename HIS → HID/HIE/HIP as needed
sed -i 's/^ATOM.\{9\}HIS A  64/ATOM      X HIE A  64/' 1HCL_fixed.pdb
```

---

## Step 4 — Separate Receptor from Ligand

```python
from rdkit import Chem

def split_receptor_ligand(pdb_file: str,
                           ligand_resname: str,
                           receptor_out: str,
                           ligand_out: str) -> None:
    """Split PDB into receptor-only and ligand-only files."""
    receptor_lines = []
    ligand_lines = []

    with open(pdb_file) as f:
        for line in f:
            record = line[:6].strip()
            if record in ("ATOM", "TER", "END"):
                receptor_lines.append(line)
            elif record == "HETATM":
                resname = line[17:20].strip()
                if resname == ligand_resname:
                    ligand_lines.append(line)
                elif resname not in ("HOH", "WAT"):
                    # Keep cofactors in receptor
                    receptor_lines.append(line)

    with open(receptor_out, "w") as f:
        f.writelines(receptor_lines)
    with open(ligand_out, "w") as f:
        f.writelines(ligand_lines)

split_receptor_ligand("1HCL_fixed.pdb", "LIG",
                      "receptor.pdb", "ligand_ref.pdb")
```

---

## Step 5 — Convert Receptor to PDBQT

```bash
# OpenBabel: PDB → PDBQT (adds partial charges, atom types)
obabel receptor.pdb -O receptor.pdbqt -xr
# -xr : rigid receptor (no rotatable bonds)

# Alternative: AutoDock Tools (pythonsh prepare_receptor4.py)
# More accurate atom typing for Mg2+, Zn2+ etc.
pythonsh prepare_receptor4.py -r receptor.pdb -o receptor.pdbqt -A hydrogens
```

---

## Step 6 — Define Docking Box

**From co-crystallized ligand (best method):**

```python
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

def box_from_ligand(ligand_pdb: str, padding: float = 8.0) -> dict:
    """Compute docking box centered on reference ligand."""
    mol = Chem.MolFromPDBFile(ligand_pdb, removeHs=False)
    if mol is None:
        raise ValueError(f"Could not parse {ligand_pdb}")

    conf = mol.GetConformer()
    coords = np.array([conf.GetAtomPosition(i)
                       for i in range(mol.GetNumAtoms())])

    center = coords.mean(axis=0)
    size = coords.max(axis=0) - coords.min(axis=0) + padding

    return {
        "center_x": round(float(center[0]), 2),
        "center_y": round(float(center[1]), 2),
        "center_z": round(float(center[2]), 2),
        "size_x":   round(float(size[0]), 1),
        "size_y":   round(float(size[1]), 1),
        "size_z":   round(float(size[2]), 1),
    }

box = box_from_ligand("ligand_ref.pdb", padding=10.0)
print(box)
# → {'center_x': 10.5, 'center_y': -2.3, 'center_z': 14.1,
#    'size_x': 22.0, 'size_y': 20.0, 'size_z': 18.0}

# Write Vina config
with open("box.conf", "w") as f:
    for k, v in box.items():
        f.write(f"{k} = {v}\n")
```

**From fpocket (when no co-crystal):**

```bash
# Detect pockets ranked by druggability score
fpocket -f receptor.pdb

# Output: receptor_out/pockets/pocket1_atm.pdb (best pocket)
# Use pocket1 center for docking box
```

```python
def box_from_fpocket(pocket_pdb: str, padding: float = 5.0) -> dict:
    """Derive box from fpocket pocket atoms."""
    coords = []
    with open(pocket_pdb) as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append([x, y, z])

    coords = np.array(coords)
    center = coords.mean(axis=0)
    size = coords.max(axis=0) - coords.min(axis=0) + padding
    return {
        "center_x": round(float(center[0]), 2),
        "center_y": round(float(center[1]), 2),
        "center_z": round(float(center[2]), 2),
        "size_x":   max(round(float(size[0]), 1), 15.0),
        "size_y":   max(round(float(size[1]), 1), 15.0),
        "size_z":   max(round(float(size[2]), 1), 15.0),
    }
```

---

## Common Problems and Fixes

| Problem | Symptom | Fix |
|---------|---------|-----|
| Missing loop near binding site | Hole in surface, wrong contacts | pdbfixer addMissingResidues or model manually |
| His tautomer wrong | Ligand clashes or loses H-bond | Check propka, set HID/HIE manually |
| Metal not typed correctly | PDBQT error or bad poses | Use AutoDockTools prepare_receptor4.py with `-M` flag |
| Co-solvent molecules (DMSO, glycerol) left in pocket | Ligand can't fit | Remove HETATM records with those resnames before prep |
| Disulfide bridge H missing | Steric clashes | Fix with `fixer.addMissingAtoms()` — checks SS bonds |
| PDB has multiple models/chains | Vina only accepts 1 chain | Extract chain A: `grep "^ATOM.*A  " receptor.pdb` |

---

## Validation — Redocking Test

Before screening, validate the prep by redocking the co-crystallized ligand:

```python
def validate_redocking(vina_bin: str, receptor_pdbqt: str,
                       ligand_pdbqt: str, ref_ligand_pdb: str,
                       box: dict) -> float:
    """Redock reference ligand and compute RMSD vs crystal pose."""
    import subprocess, tempfile
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolAlign

    with tempfile.NamedTemporaryFile(suffix=".pdbqt", delete=False) as tmp:
        out_pdbqt = tmp.name

    cmd = [vina_bin,
           "--receptor", receptor_pdbqt,
           "--ligand",   ligand_pdbqt,
           "--out",      out_pdbqt,
           "--exhaustiveness", "32",
           "--num_modes", "1"]
    for k, v in box.items():
        cmd += [f"--{k}", str(v)]

    subprocess.run(cmd, check=True, capture_output=True)

    # Parse best pose RMSD vs crystal
    # Convert out_pdbqt → SDF via obabel, then RDKit RMSD
    subprocess.run(["obabel", out_pdbqt, "-O", "best_pose.sdf"], check=True)

    ref = Chem.MolFromPDBFile(ref_ligand_pdb, removeHs=True)
    docked = Chem.SDMolSupplier("best_pose.sdf", removeHs=True)[0]

    if ref is None or docked is None:
        return float("nan")

    rmsd = rdMolAlign.CalcRMS(docked, ref)
    print(f"Redocking RMSD: {rmsd:.2f} Å  ({'PASS' if rmsd < 2.0 else 'FAIL'})")
    return rmsd

# Acceptable: RMSD < 2.0 Å
```
