# ASE — Calculators

## Calculator Pattern

```python
atoms.calc = SomeCalculator(...)   # attach
e = atoms.get_potential_energy()   # eV
f = atoms.get_forces()             # eV/Å, shape (N, 3)
s = atoms.get_stress()             # Voigt notation, eV/Å³
```

All calculators expose the same interface — swap freely.

---

## EMT — Embedded Atom Method (built-in, no install)

Fast toy potential for FCC metals (Al, Cu, Ag, Au, Ni, Pd, Pt) and their alloys.
**Only for testing workflows — not production quality.**

```python
from ase.calculators.emt import EMT

atoms.calc = EMT()
print(atoms.get_potential_energy())
```

---

## ORCA — DFT / Coupled Cluster (ASE 3.22+)

```python
from ase.calculators.orca import ORCA, OrcaProfile

# Simple setup (ORCA must be in PATH or use OrcaProfile)
calc = ORCA(
    charge=0,
    mult=1,
    orcasimpleinput='B3LYP def2-TZVP D3BJ Opt',
    orcablocks='%pal nprocs 8 end\n%maxcore 2000'
)

# Custom executable path
profile = OrcaProfile(command='/path/to/orca')
calc = ORCA(
    profile=profile,
    charge=0,
    mult=1,
    orcasimpleinput='PBE0 def2-SVP TightSCF',
    orcablocks='%pal nprocs 4 end'
)

atoms.calc = calc
energy = atoms.get_potential_energy()  # triggers ORCA run
```

### Common ORCA inputs

```python
# DFT geometry optimization
orcasimpleinput='B3LYP def2-TZVP D3BJ Opt TightOpt'

# Single-point energy
orcasimpleinput='DLPNO-CCSD(T) def2-TZVPP TightSCF'

# Frequency calculation
orcasimpleinput='B3LYP def2-SVP Freq'

# NMR
orcasimpleinput='B3LYP def2-TZVP NMR'

# Solvation (CPCM)
orcablocks='%cpcm\n  smd True\n  SMDsolvent "water"\nend'
```

> Note: ORCA handles geometry optimization internally via `Opt` keyword. For ASE-driven optimization, omit `Opt` and use an ASE optimizer.

---

## xTB / TBLite — Semi-empirical GFN (fast, no license)

```python
# Via tblite-python (recommended, lightweight)
from tblite.ase import TBLite

calc = TBLite(
    method='GFN2-xTB',   # or 'GFN1-xTB', 'IPEA1-xTB'
    charge=0,
    multiplicity=1,
    verbosity=0
)
atoms.calc = calc

# Alternatively via xtb-python
from xtb.ase.calculator import XTB
calc = XTB(method='GFN2-xTB')
atoms.calc = calc
```

**GFN2-xTB** is the recommended default for organic/drug-like molecules.

---

## GPAW — Plane-wave DFT (open source, periodic systems)

```python
from gpaw import GPAW, PW

calc = GPAW(
    mode=PW(500),           # plane-wave cutoff 500 eV
    xc='PBE',
    kpts=(4, 4, 1),         # k-point sampling
    txt='gpaw.out'
)
atoms.calc = calc
```

---

## Lennard-Jones — Simple pair potential

```python
from ase.calculators.lj import LennardJones

calc = LennardJones(
    sigma=3.4,    # Å
    epsilon=0.01  # eV
)
atoms.calc = calc
```

---

## LAMMPS — Classical MD interoperability

```python
from ase.calculators.lammpslib import LAMMPSlib

cmds = ["pair_style eam/alloy",
        "pair_coeff * * Cu_u3.eam Cu"]

calc = LAMMPSlib(
    lmpcmds=cmds,
    atom_types={'Cu': 1},
    log_file='lammps.log'
)
atoms.calc = calc
```

---

## Quantum ESPRESSO — Plane-wave DFT (periodic)

```python
from ase.calculators.espresso import Espresso, EspressoProfile

profile = EspressoProfile(command='pw.x', pseudo_dir='/path/to/pseudos')

calc = Espresso(
    profile=profile,
    input_data={
        'system': {
            'ecutwfc': 60,      # Ry
            'ecutrho': 480,
            'occupations': 'smearing',
            'smearing': 'mv',
            'degauss': 0.02,
        },
        'control': {'calculation': 'scf'},
    },
    pseudopotentials={'Cu': 'Cu.pbe-dn-kjpaw_psl.1.0.0.UPF'},
    kpts=(4, 4, 4),
)
atoms.calc = calc
```

---

## Wrapper Calculators

### SumCalculator — combine potentials

```python
from ase.calculators.mixing import SumCalculator
from ase.calculators.dftd3 import DFTD3

# DFT + dispersion correction
dft  = GPAW(mode=PW(500), xc='PBE')
disp = DFTD3(damping='bj', xc='pbe')
calc = SumCalculator([dft, disp])
atoms.calc = calc
```

### SocketIOCalculator — persistent server

```python
from ase.calculators.socketio import SocketIOCalculator

# Avoids relaunching DFT code per step
with SocketIOCalculator(calc, unixsocket='ase_socket') as sio_calc:
    atoms.calc = sio_calc
    # Run MD or optimization
```

---

## Calculator Configuration (`~/.config/ase/config.ini`)

```ini
[ORCA]
command = /opt/orca/orca

[espresso]
command = mpirun -np 4 pw.x -in PREFIX.pwi > PREFIX.pwo
pseudo_dir = /home/user/pseudos/
```

```bash
# Check configuration
ase info --calculators
```

---

## Common Properties

```python
atoms.get_potential_energy()      # total energy, eV
atoms.get_forces()                # (N,3) array, eV/Å
atoms.get_stress()                # Voigt stress (6,), eV/Å³
atoms.get_dipole_moment()         # Debye (if supported)
atoms.get_charges()               # partial charges (if supported)
atoms.get_magnetic_moment()       # total magnetic moment
atoms.get_magnetic_moments()      # per-atom moments
```
