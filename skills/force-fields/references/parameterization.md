# Parameterization — GAFF2, CGenFF, Charge Methods

Parameterizing a small molecule = assigning force field types + computing partial charges. Three main workflows:

| Workflow | FF | Charges | Tools |
|---|---|---|---|
| GAFF2 | AMBER ecosystem | AM1-BCC | antechamber + tleap / acpype |
| CGenFF | CHARMM ecosystem | RESP-like | CHARMM-GUI (web) + ffTK |
| OpenFF Sage | SMIRNOFF | AM1-BCC(ELF10) | openff-toolkit (see openff-smirnoff.md) |

---

## GAFF2 with antechamber (AMBER Tools)

GAFF2 = General AMBER Force Field 2. Best for drug-like organics in AMBER/GROMACS simulations.

### Step 1 — Prepare input (need 3D structure)

```bash
# Convert SMILES → 3D SDF → mol2 (via RDKit or OpenBabel)
python -c "
from rdkit import Chem
from rdkit.Chem import AllChem
mol = Chem.MolFromSmiles('c1ccc(cc1)CC(=O)O')
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
AllChem.MMFFOptimizeMolecule(mol)
Chem.MolToMolFile(mol, 'lig.sdf')
"

# Or use obabel
obabel -:'c1ccc(cc1)CC(=O)O' --gen3d -O lig.mol2
```

### Step 2 — Run antechamber (assigns GAFF2 atom types + AM1-BCC charges)

```bash
antechamber -i lig.sdf -fi sdf \
            -o lig.mol2 -fo mol2 \
            -c bcc \       # AM1-BCC charges
            -s 2 \         # verbosity
            -at gaff2 \    # atom type: gaff2
            -nc 0          # net charge (integer)

# Check for missing parameters
parmchk2 -i lig.mol2 -f mol2 -o lig.frcmod -s gaff2
```

### Step 3 — Build AMBER topology (tleap)

```bash
# lig.leap
cat > lig.leap <<EOF
source leaprc.gaff2
loadamberparams lig.frcmod
LIG = loadmol2 lig.mol2
saveamberparm LIG lig.prmtop lig.inpcrd
quit
EOF

tleap -f lig.leap
```

### Step 4 — Full protein-ligand system (tleap)

```bash
cat > complex.leap <<EOF
source leaprc.protein.ff14SB
source leaprc.water.tip3p
source leaprc.gaff2

loadamberparams lig.frcmod
LIG = loadmol2 lig.mol2

complex = loadpdb protein.pdb
complex = combine {complex LIG}

solvatebox complex TIP3PBOX 10.0   # 10 Å padding
addIons complex Na+ 0              # neutralize
addIons complex Na+ 10             # add 150 mM

saveamberparm complex complex.prmtop complex.inpcrd
savepdb complex complex.pdb
quit
EOF

tleap -f complex.leap
```

---

## acpype — GAFF2 to GROMACS/AMBER (Python wrapper)

acpype wraps antechamber and generates GROMACS (`.itp`/`.gro`/`.top`) or AMBER (`.prmtop`/`.inpcrd`) topologies.

```bash
pip install acpype
# or: conda install -c conda-forge acpype

# Generate GROMACS topology
acpype -i lig.sdf -b lig -n 0 -a gaff2 -c bcc

# Output: lig.acpype/
#   lig_GMX.itp     ← GROMACS include topology
#   lig_GMX.gro     ← GROMACS coordinates
#   lig_GMX.top     ← standalone topology
#   lig_AC.frcmod   ← AMBER force field modification
#   lig_AC.prmtop   ← AMBER topology
#   lig_AC.inpcrd   ← AMBER coordinates
```

```python
# Python API
import subprocess

result = subprocess.run([
    'acpype',
    '-i', 'lig.sdf',
    '-b', 'lig',
    '-n', '0',        # net charge
    '-a', 'gaff2',
    '-c', 'bcc',
], capture_output=True, text=True)

print(result.stdout)
```

---

## CGenFF (CHARMM General Force Field)

CHARMM counterpart to GAFF2. Best parameterized via **CHARMM-GUI** (web service) or **CGenFF program** (license required).

### CHARMM-GUI Workflow
1. Go to https://www.charmm-gui.org → Ligand Reader & Modeler
2. Upload SMILES or SDF
3. Download `.str` (stream file containing topology + parameters)
4. Use with `app.CharmmParameterSet` in OpenMM

```python
import openmm.app as app

# Load CHARMM36m + CGenFF ligand parameters
params = app.CharmmParameterSet(
    'charmm36.prm',   # protein parameters
    'lig.str',        # CGenFF ligand stream file
)
psf = app.CharmmPsfFile('complex.psf')
system = psf.createSystem(params, nonbondedMethod=app.PME, ...)
```

### ffTK (Force Field Toolkit) — Penalty Scores
CGenFF assigns a **penalty score** per parameter:
- Score 0 → direct match, high confidence
- Score < 10 → minor analog transfer, generally usable
- Score 10–50 → moderate transfer, test carefully
- Score > 50 → requires QM optimization

```python
# Check penalty in .str file
# Lines look like: ! penalty= 35.000
import re

with open('lig.str') as f:
    content = f.read()

penalties = re.findall(r'penalty=\s*([\d.]+)', content)
max_penalty = max(float(p) for p in penalties)
print(f"Max penalty: {max_penalty}")  # >50 → needs optimization
```

---

## Partial Charge Methods

### AM1-BCC (Recommended for drug-like molecules)
- Semi-empirical AM1 electrostatic potential + BCC (bond charge correction)
- Fast, good accuracy for drug-like organics
- Available in antechamber (free) and OpenFF toolkit

```bash
antechamber -i lig.sdf -fi sdf -o lig.mol2 -fo mol2 -c bcc -nc 0
```

### RESP (Restrained Electrostatic Potential)
- Fit charges to QM (HF/6-31G*) electrostatic potential
- More accurate, especially for polar molecules
- Requires Gaussian or ORCA QM output

```bash
# Step 1: QM ESP calculation (Gaussian)
# In .com file: #HF/6-31G* Pop=MK IOp(6/33=2)

# Step 2: RESP fitting (AMBER Tools)
espgen -i lig.gesp -o lig.esp
respgen -i lig.ac -o lig.respin1 -f resp1
respgen -i lig.ac -o lig.respin2 -f resp2
resp -O -i lig.respin1 -o lig.respout1 -e lig.esp -t qout_stage1
resp -O -i lig.respin2 -o lig.respout2 -e lig.esp -q qout_stage1 -t qout_stage2
```

### RESP2
- Two-phase RESP: gas + implicit solvent (δ=0.6 mixing)
- Better for condensed-phase simulations
- Reference: Schauperl et al. 2020, DOI: 10.1038/s42004-020-0291-4

---

## ParmEd — Topology Manipulation

```python
import parmed as pmd

# Load AMBER topology
struct = pmd.load_file('complex.prmtop', 'complex.inpcrd')

# Inspect
print(struct.atoms[:5])
print(struct.residues[:3])
print(struct.bonds[:5])

# Change charge on atom
struct[0].charge = -0.31

# Save as GROMACS topology
struct.save('output.top')
struct.save('output.gro')

# Save as AMBER
struct.save('output.prmtop')
struct.save('output.inpcrd')

# Strip water (keep protein + ligand)
stripped = struct['!:WAT,Na+,Cl-']
stripped.save('nowater.prmtop')
stripped.save('nowater.inpcrd')

# Merge two structures
complex = protein_struct + ligand_struct

# Fix box vectors
struct.box = [50.0, 50.0, 50.0, 90.0, 90.0, 90.0]  # Å and degrees
```

---

## Validation Checklist

After parameterization, always verify:

```python
# 1. Charges sum to net charge
total_q = sum(a.charge for a in struct.atoms)
print(f"Total charge: {total_q:.4f}")  # should be integer

# 2. No missing parameters (antechamber)
# Check: parmchk2 output has no ATTN warnings

# 3. Quick energy minimization (should complete without NaN)
# Run 100 steps with OpenMM, check energy is finite

# 4. Check for close contacts (bad initial geometry)
state = simulation.context.getState(getEnergy=True)
print(state.getPotentialEnergy())  # should not be >> 1e6 kJ/mol

# 5. Validate atom types assigned correctly
for atom in struct.atoms[:10]:
    print(f"{atom.name:4s} type={atom.type:6s} charge={atom.charge:.4f}")
```
