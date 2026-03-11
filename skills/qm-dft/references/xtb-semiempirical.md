# xTB / Semi-Empirical Methods

xTB (extended tight-binding) by Grimme group. GFN2-xTB is the workhorse: DFT-level geometry and thermochemistry at a fraction of the cost. Free, open-source (xtb 6.7, CREST 3.0).

---

## Method Overview

| Method | Type | Atoms | Notes |
|--------|------|-------|-------|
| GFN-FF | Force field | All | Fastest, pre-screening |
| GFN1-xTB | Semi-empirical | All | First GFN; mostly replaced by GFN2 |
| GFN2-xTB | Semi-empirical | All | Standard; good geometry + thermochem |
| GFN2-xTB//GFN-FF | Dual-level | All | FF geometry + xTB energy (fast) |
| GFN0-xTB | Non-self-consistent | All | Very fast, geometry only |

---

## CLI Usage

### Basic Single Point
```bash
xtb mol.xyz --gfn 2 --sp
# Output: xtb.out (energy, dipole, charges)
```

### Geometry Optimization
```bash
xtb mol.xyz --opt --gfn 2
# Output: xtbopt.xyz (optimized), xtbopt.log

# Optimization level (default: normal)
xtb mol.xyz --opt crude    # fast, loose
xtb mol.xyz --opt normal   # standard
xtb mol.xyz --opt tight    # tighter convergence
xtb mol.xyz --opt verytight
```

### Frequency Calculation
```bash
xtb mol.xyz --hess --gfn 2
# Output: hessian (force constants), vibspectrum (frequencies), g98.out
```

### Solvation
```bash
# ALPB (Analytical Linearized Poisson-Boltzmann) — recommended
xtb mol.xyz --opt --alpb water
xtb mol.xyz --opt --alpb methanol
xtb mol.xyz --opt --alpb dmso
xtb mol.xyz --opt --alpb acetonitrile

# GBSA (Generalized Born)
xtb mol.xyz --opt --gbsa water

# Available solvents: water, methanol, ethanol, acetonitrile, dmso,
# thf, acetone, chloroform, toluene, benzene, hexane, ccl4, cs2
```

### Charge and Multiplicity
```bash
xtb mol.xyz --opt --chrg -1 --uhf 1   # anion (charge=-1, 1 unpaired electron)
xtb mol.xyz --opt --chrg 0 --uhf 0    # neutral singlet (default)
```

### Partial Charges Output
```bash
xtb mol.xyz --sp --gfn 2
# Charges in: charges file (one per atom, in Mulliken or CM5)
# Also in xtb.out under "Mulliken/CM5 charges"
```

---

## Python API (tblite)

```python
from tblite.interface import Calculator
import numpy as np

# Molecule in atomic units (Bohr)
BOHR = 1.8897259886  # Å → Bohr

numbers = np.array([6, 6, 8, 1, 1, 1, 1])   # atomic numbers
positions = np.array([
    [0.0,  0.0,  0.0],
    [1.54, 0.0,  0.0],
    # ... Å, will be converted
]) * BOHR   # must be in Bohr

calc = Calculator(
    method='GFN2-xTB',
    numbers=numbers,
    positions=positions,  # Bohr
    charge=0,
    uhf=0,               # number of unpaired electrons
)

# Set solvation
calc.set('solvent', 'water')   # ALPB

# Compute
result = calc.singlepoint()
energy = result.get('energy')      # Hartree
gradient = result.get('gradient')  # Hartree/Bohr
charges = result.get('charges')    # Mulliken charges

print(f"E = {energy:.8f} Eh = {energy * 627.51:.2f} kcal/mol")
```

### ASE + xTB (via xtb-python or tblite)
```python
# Preferred: use ASE calculator (see ase skill)
from ase.calculators.xtb import XTB  # if xtb-python installed
# or:
from tblite.ase import TBLite

from ase.build import molecule
from ase.optimize import BFGS

atoms = molecule('aspirin')
atoms.calc = TBLite(method='GFN2-xTB', charge=0, uhf=0)

opt = BFGS(atoms, trajectory='xtb_opt.traj')
opt.run(fmax=0.05)  # eV/Å

print(atoms.get_potential_energy())  # eV
```

---

## CREST — Conformer/Ensemble Search

CREST (Conformer-Rotamer Ensemble Sampling Tool). Uses GFN2-xTB + iMTD-GC algorithm. Finds all low-energy conformers within a given energy window.

### Basic Conformer Search
```bash
crest mol.xyz --T 8          # use 8 threads
# Output: crest_conformers.xyz (ensemble), crest_best.xyz (lowest energy)

# With solvation
crest mol.xyz --T 8 --alpb water

# Energy window (default 6 kcal/mol above global minimum)
crest mol.xyz --T 8 --ewin 3.0   # tighter: only within 3 kcal/mol
```

### Parse CREST Output (Python)
```python
from xtb.utils import get_method  # if xtb-python available
# Or parse directly:

def parse_crest_ensemble(xyz_file):
    """Parse multi-conformer XYZ file from CREST."""
    conformers = []
    with open(xyz_file) as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        n = int(lines[i].strip())
        comment = lines[i+1].strip()
        # comment contains energy in Hartree
        try:
            energy = float(comment.split()[0])
        except (ValueError, IndexError):
            energy = None
        coords = []
        for j in range(i+2, i+2+n):
            parts = lines[j].split()
            coords.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
        conformers.append({'energy': energy, 'coords': coords, 'n_atoms': n})
        i += n + 2
    return conformers

conformers = parse_crest_ensemble('crest_conformers.xyz')
print(f"{len(conformers)} conformers found")
for i, c in enumerate(conformers[:5]):
    dE = (c['energy'] - conformers[0]['energy']) * 627.51  # kcal/mol
    print(f"Conf {i+1}: ΔE = {dE:.2f} kcal/mol")
```

### Protonation/Deprotonation Ensemble (pKa-related)
```bash
# Find all reasonable protonation sites
crest mol.xyz --protonate
# Output: protonated.xyz

crest mol.xyz --deprotonate
# Output: deprotonated.xyz
```

### Tautomer Search
```bash
crest mol.xyz --tautomerize
# Output: tautomers.xyz (ranked by energy)
```

---

## pKa Estimation with xTB

GFN2-xTB + ALPB solvation for fast pKa estimation. Not as accurate as DFT-SMD but useful for screening.

```bash
# pKa via protonation free energy (manual workflow)
# 1. Optimize neutral form
xtb mol.xyz --opt --alpb water
mv xtbopt.xyz mol_neutral.xyz

# 2. Generate deprotonated form with CREST
crest mol_neutral.xyz --deprotonate --alpb water

# 3. Optimize deprotonated form
xtb deprotonated.xyz --opt --alpb water
mv xtbopt.xyz mol_deprot.xyz

# 4. Compute free energies
xtb mol_neutral.xyz --hess --alpb water   # G_neutral
xtb mol_deprot.xyz  --hess --alpb water   # G_deprot

# 5. ΔG = G_deprot + G_H+ - G_neutral; pKa = ΔG / (RT·ln10)
# G_H+(aq) ≈ -270.3 kcal/mol (literature value for computational pKa)
```

---

## xTB Output Parsing

```bash
# Key lines in xtb.out:
grep 'TOTAL ENERGY' xtb.out        # E (Hartree)
grep 'TOTAL ENTHALPY' xtb.out      # H (Hartree)
grep 'TOTAL FREE ENERGY' xtb.out   # G (Hartree, includes RRHO)
grep 'DIPOLE' xtb.out              # dipole moment (Debye)
```

```python
def parse_xtb_output(outfile):
    result = {}
    with open(outfile) as f:
        for line in f:
            if '| TOTAL ENERGY' in line:
                result['energy'] = float(line.split()[3])
            elif '| TOTAL ENTHALPY' in line:
                result['enthalpy'] = float(line.split()[3])
            elif '| TOTAL FREE ENERGY' in line:
                result['gibbs'] = float(line.split()[4])
            elif 'molecular dipole:' in line:
                pass  # next lines have components
    return result
```

---

## Limitations of xTB / GFN2

- No excited states (no TD-xTB for UV-Vis)
- Less accurate for transition metals (GFN2 parametrized for Z=1–86)
- Torsion barriers less reliable than DFT (use DFT for rotamer profiles)
- No analytical NMR shifts
- Polarizabilities approximate
- **Always validate important results with at least one DFT single-point**
