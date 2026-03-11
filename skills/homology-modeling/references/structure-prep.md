# Structure Preparation for MD

After building a homology model (or loading a crystal structure), preparation for MD requires: fixing missing atoms/residues, assigning protonation states, capping termini, and handling special residues. Errors here propagate through the entire simulation.

---

## pdbfixer — The Essential Fixer

```bash
pip install pdbfixer
```

```python
from pdbfixer import PDBFixer
from openmm.app import PDBFile

fixer = PDBFixer(filename='model.pdb')

# 1. Find and add missing residues (in gaps or termini)
fixer.findMissingResidues()
# Remove terminal missing residues (often disordered tails — don't add them)
chain_ids = list(fixer.missingResidues.keys())
for key in list(fixer.missingResidues.keys()):
    chain = key[0]
    res_idx = key[1]
    chain_len = sum(1 for r in fixer.topology.residues() if r.chain.id == chain.id)
    # Remove if at N- or C-terminus
    if res_idx == 0 or res_idx >= chain_len:
        del fixer.missingResidues[key]

# 2. Find missing heavy atoms
fixer.findMissingAtoms()
fixer.addMissingAtoms()

# 3. Add hydrogens at target pH
fixer.addMissingHydrogens(pH=7.4)

# 4. Remove water (keep=True if you want crystallographic waters)
fixer.removeHeterogens(keepWater=False)

# 5. Replace nonstandard residues (e.g. MSE → MET)
fixer.findNonstandardResidues()
fixer.replaceNonstandardResidues()

# Save
with open('model_fixed.pdb', 'w') as f:
    PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)
```

---

## propka3 — pKa-Based Protonation

```bash
pip install propka
```

```python
import subprocess
import re

def run_propka(pdb_file, pH=7.4):
    """Run propka and parse predicted pKa values."""
    result = subprocess.run(
        ['propka3', pdb_file],
        capture_output=True, text=True
    )
    # Output: model.pka file
    pka_file = pdb_file.replace('.pdb', '.pka')

    residues = {}
    with open(pka_file) as f:
        for line in f:
            # Parse lines like: ASP  35 A    3.80  ...
            m = re.match(r'(\w{3})\s+(\d+)\s+(\w)\s+([\d.]+)', line)
            if m:
                resname, resid, chain, pka = m.groups()
                residues[(resname, int(resid), chain)] = float(pka)

    return residues

pka_dict = run_propka('model_fixed.pdb', pH=7.4)

# Identify titratable residues that deviate from default
for (resname, resid, chain), pka in sorted(pka_dict.items()):
    if resname == 'HIS':
        if pka > 7.4:
            state = 'HIP (protonated, +1)'
        else:
            state = 'HID or HIE (neutral)'
        print(f"HIS {resid}{chain}: pKa={pka:.1f} → {state}")
    elif resname == 'ASP' and pka > 5.5:
        print(f"ASP {resid}{chain}: pKa={pka:.1f} → ASH (protonated, neutral)")
    elif resname == 'GLU' and pka > 5.5:
        print(f"GLU {resid}{chain}: pKa={pka:.1f} → GLH (protonated, neutral)")
    elif resname == 'LYS' and pka < 9.0:
        print(f"LYS {resid}{chain}: pKa={pka:.1f} → LYN (neutral, deprotonated)")
    elif resname == 'CYS' and pka < 7.0:
        print(f"CYS {resid}{chain}: pKa={pka:.1f} → CYM (deprotonated, -1)")
```

---

## Histidine Tautomer Assignment

HIS has three states; AMBER uses explicit residue names:

| AMBER name | State | ε-N (NE2) | δ-N (ND1) |
|-----------|-------|-----------|-----------|
| HID | neutral, δ-H | — | H |
| HIE | neutral, ε-H | H | — |
| HIP | protonated (+1) | H | H |

```python
from Bio.PDB import PDBParser, PDBIO

def rename_his_tautomers(pdb_file, pka_dict, pH=7.4, output='model_his.pdb'):
    """Rename HIS residues to HID/HIE/HIP based on pKa and environment."""
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure('model', pdb_file)

    for model in struct:
        for chain in model:
            for residue in chain:
                if residue.resname != 'HIS':
                    continue
                resid = residue.id[1]
                chain_id = chain.id
                pka = pka_dict.get(('HIS', resid, chain_id), 6.5)

                if pka > pH:
                    residue.resname = 'HIP'  # protonated
                else:
                    # Check which N is closer to H-bond acceptor
                    # Default: HIE (ε-H exposed to solvent)
                    residue.resname = 'HIE'

    io = PDBIO()
    io.set_structure(struct)
    io.save(output)
    return output
```

---

## Disulfide Bond Detection

```python
from Bio.PDB import PDBParser
import numpy as np

def find_disulfides(pdb_file, dist_threshold=2.2):
    """
    Find CYS–CYS pairs within S–S bond distance.
    Returns list of (chain1, resid1, chain2, resid2, distance) tuples.
    """
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure('model', pdb_file)

    cys_sg = []
    for model in struct:
        for chain in model:
            for res in chain:
                if res.resname in ('CYS', 'CYM') and 'SG' in res:
                    cys_sg.append((chain.id, res.id[1], res['SG'].coord))

    disulfides = []
    for i in range(len(cys_sg)):
        for j in range(i+1, len(cys_sg)):
            c1, r1, sg1 = cys_sg[i]
            c2, r2, sg2 = cys_sg[j]
            dist = np.linalg.norm(sg1 - sg2)
            if dist < dist_threshold:
                disulfides.append((c1, r1, c2, r2, dist))
                print(f"Disulfide: CYS{r1}{c1} — CYS{r2}{c2}  d={dist:.2f} Å")

    return disulfides

ss_bonds = find_disulfides('model_fixed.pdb')
```

### tleap Disulfide Patching (AMBER)
```bash
# In tleap script:
bond complex.A.184.SG complex.A.196.SG
```

### OpenMM Disulfide (automatic)
```python
# pdbfixer handles disulfides automatically when addMissingAtoms() is called
# OpenMM ff14SB recognizes CYX residue name for disulfide-bonded CYS
```

---

## Terminal Capping (ACE/NME)

Protein termini that are internal fragments need capping to avoid spurious charge:

```python
# In tleap (AMBER):
# N-terminus cap (ACE = acetyl):
#   loadoff ace.lib
# C-terminus cap (NME = N-methylamide):
#   loadoff nme.lib
# Then: set complex.A.1 tail complex.A.1.C
#       set ace.1 head ace.1.N

# With pdbfixer: termini are automatically treated as charged (NH3+/COO-)
# For capped fragments, manually add ACE/NME residues before running fixer

# OpenFF approach — specify terminal patches in SystemGenerator
from openmmforcefields.generators import SystemGenerator
generator = SystemGenerator(
    forcefields=['amber/ff14SB.xml', 'amber/tip3p_standard.xml'],
    forcefield_kwargs={'constraints': app.HBonds}
)
# Capping is handled by ff14SB.xml for standard termini
```

---

## Crystallographic Waters

```python
def decide_waters(pdb_file, ligand_resname='LIG', keep_radius=5.0):
    """
    Keep crystallographic waters within keep_radius Å of ligand.
    Remove all others.
    """
    from Bio.PDB import PDBParser, PDBIO, Select
    import numpy as np

    parser = PDBParser(QUIET=True)
    struct = parser.get_structure('s', pdb_file)

    # Get ligand atom positions
    lig_pos = []
    for model in struct:
        for chain in model:
            for res in chain:
                if res.resname == ligand_resname:
                    for atom in res:
                        lig_pos.append(atom.coord)
    lig_pos = np.array(lig_pos)

    class WaterSelector(Select):
        def accept_residue(self, residue):
            if residue.resname not in ('HOH', 'WAT', 'TIP'):
                return True  # keep non-water
            # Keep water if close to ligand
            for atom in residue:
                dists = np.linalg.norm(lig_pos - atom.coord, axis=1)
                if dists.min() <= keep_radius:
                    return True
            return False

    io = PDBIO()
    io.set_structure(struct)
    io.save('model_waters.pdb', WaterSelector())
```

---

## Preparation Checklist

```
Before MD:
  ✓ Remove all HETATM except ligand (cofactors need separate parameterization)
  ✓ Missing internal residues → pdbfixer (loops modeled or flagged)
  ✓ Missing heavy atoms → pdbfixer.addMissingAtoms()
  ✓ Protonation states → propka3 at target pH (7.4 for cytoplasmic)
  ✓ HIS tautomers: check H-bond environment (HID/HIE/HIP)
  ✓ Disulfide bonds: detect SG–SG < 2.2 Å, use CYX in AMBER
  ✓ Net charge: count Asp/Glu/Lys/Arg/His, verify integer
  ✓ Sequence coverage: note any missing loops modeled or skipped
  ✓ Termini: charged (default) or capped (ACE/NME for fragments)
  After solvation:
  ✓ Minimize 500 steps with restraints on heavy atoms
  ✓ Heat gradually: 0 K → 300 K over 100 ps
  ✓ Equilibrate NPT 1 ns before production
```
