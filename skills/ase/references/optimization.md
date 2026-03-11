# ASE — Geometry Optimization

## Overview

```
workflow: set calculator → attach optimizer → run(fmax) → read energy/positions
```

`fmax` = max force on any atom in eV/Å. Typical: 0.01–0.05 eV/Å.

---

## Optimizers

```python
from ase.optimize import BFGS, LBFGS, FIRE, GPMin
from ase.optimize.sciopt import SciPyFminBFGS, SciPyFminCG
```

| Optimizer | Best for | Notes |
|-----------|----------|-------|
| `BFGS` | Most molecular cases | Stable, memory-cheap |
| `LBFGS` | Large systems | Limited-memory BFGS |
| `FIRE` | NEB, periodic cells | Works with constraints |
| `GPMin` | Expensive calculators | Gaussian process surrogate |
| `SciPyFminBFGS` | Small systems | SciPy backend |

---

## Basic Optimization

```python
from ase.build import molecule
from ase.calculators.orca import ORCA
from ase.optimize import BFGS

atoms = molecule('H2O')
atoms.calc = ORCA(charge=0, mult=1, orcasimpleinput='B3LYP def2-SVP')

opt = BFGS(
    atoms,
    trajectory='opt.traj',   # save all steps
    logfile='opt.log',        # '-' for stdout
    maxstep=0.2               # Å, max atomic displacement per step
)
opt.run(fmax=0.01)            # converge when max force < 0.01 eV/Å

print(f"Energy: {atoms.get_potential_energy():.4f} eV")
print(f"Forces:\n{atoms.get_forces()}")
```

---

## Constraints

```python
from ase.constraints import (
    FixAtoms, FixBondLength, FixBondLengths,
    FixCartesian, FixInternals, FixLinearTriatomic
)

# Fix specific atoms (e.g., bottom slab layers)
atoms.set_constraint(FixAtoms(indices=[0, 1, 2]))
atoms.set_constraint(FixAtoms(mask=[True, True, False, False]))

# Fix bond length between atoms 0 and 1
atoms.set_constraint(FixBondLength(0, 1))

# Fix z-coordinate only
atoms.set_constraint(FixCartesian(0, mask=[0, 0, 1]))  # fix z of atom 0

# Multiple constraints
from ase.constraints import FixAtoms, FixBondLength
atoms.set_constraint([FixAtoms([0, 1]), FixBondLength(2, 3)])

# Remove constraints
atoms.set_constraint()
```

---

## Unit Cell Optimization (variable cell)

```python
from ase.filters import ExpCellFilter, FrechetCellFilter

# ExpCellFilter — optimize positions AND cell simultaneously
atoms.calc = calc
ecf = ExpCellFilter(atoms, scalar_pressure=0.0)   # eV/Å³
opt = BFGS(ecf)
opt.run(fmax=0.01)

# FrechetCellFilter — better conditioned, recommended for hard systems
fcf = FrechetCellFilter(atoms)
opt = BFGS(fcf)
opt.run(fmax=0.01)

# Position-only (fixed cell, default BFGS)
opt = BFGS(atoms)
opt.run(fmax=0.01)
```

---

## Convergence and Logging

```python
# Run with step limit
opt.run(fmax=0.05, steps=500)

# Attach observer to run function every N steps
def print_energy():
    print(atoms.get_potential_energy())

opt.attach(print_energy, interval=5)
opt.run(fmax=0.01)

# Check if converged
converged = opt.converged()
nsteps    = opt.nsteps
```

---

## Reading Optimization Trajectory

```python
from ase import io
import matplotlib.pyplot as plt

frames = io.read('opt.traj', index=':')
energies = [f.get_potential_energy() for f in frames]

plt.plot(energies)
plt.xlabel('Step')
plt.ylabel('Energy (eV)')
plt.title('Optimization convergence')
plt.show()
```

---

## Full Example: ORCA + BFGS + xTB pre-optimization

```python
from ase.build import molecule
from ase.optimize import BFGS
from tblite.ase import TBLite
from ase.calculators.orca import ORCA, OrcaProfile

# 1. Fast pre-optimization with GFN2-xTB
atoms = molecule('aspirin')
atoms.calc = TBLite(method='GFN2-xTB', verbosity=0)
BFGS(atoms).run(fmax=0.05)

# 2. Refine with DFT
atoms.calc = ORCA(
    charge=0, mult=1,
    orcasimpleinput='PBE0 def2-TZVP D3BJ TightSCF'
)
opt = BFGS(atoms, trajectory='aspirin_pbe0.traj', logfile='-')
opt.run(fmax=0.01)

energy = atoms.get_potential_energy()
print(f"Final energy: {energy:.6f} eV")
```

---

## Equation of State (EOS)

```python
from ase.eos import EquationOfState

volumes = []
energies = []
for scale in [0.95, 0.97, 0.99, 1.00, 1.01, 1.03, 1.05]:
    a = atoms.copy()
    a.set_cell(atoms.cell * scale, scale_atoms=True)
    a.calc = calc
    energies.append(a.get_potential_energy())
    volumes.append(a.get_volume())

eos = EquationOfState(volumes, energies, eos='birchmurnaghan')
v0, e0, B = eos.fit()
print(f"V0={v0:.2f} Å³, E0={e0:.4f} eV, B={B/units.GPa:.1f} GPa")
eos.plot('eos.png')
```
