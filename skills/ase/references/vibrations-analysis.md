# ASE — Vibrations, Phonons & Thermochemistry

## Overview

- **Vibrations**: molecular normal modes, IR intensities, ZPE (finite differences)
- **Phonons**: periodic solids, band structure, DOS
- **Thermochemistry**: enthalpy, entropy, free energy, heat capacity

---

## Molecular Vibrations

```python
from ase.vibrations import Vibrations

# Atoms must be relaxed (fmax < 0.01 eV/Å) before vibrational analysis
vib = Vibrations(
    atoms,
    name='vib',          # prefix for cache files (vib.*.json)
    delta=0.01,          # finite displacement step, Å
    nfree=2,             # 2 (central diff) or 4 (more accurate)
    indices=None         # None = all atoms, or list of indices to displace
)

vib.run()               # compute force constants (2*3*N DFT calls)
vib.summary()           # print frequencies to stdout

# Access results
freqs = vib.get_frequencies()          # cm⁻¹, complex for imaginary
modes = vib.get_mode(i)                # (N, 3) displacement vector for mode i
zpe   = vib.get_zero_point_energy()    # eV

# Identify imaginary modes (saddle point check)
for i, f in enumerate(freqs):
    if f.imag > 0:
        print(f"Mode {i}: imaginary {f.imag:.1f}i cm⁻¹")
```

### Caching and Restart

```python
# Results are cached in vib.*.json files
# To restart from existing cache:
vib = Vibrations(atoms, name='vib')
# vib.run() will skip already-computed displacements

# Clean cache
vib.clean()
```

---

## Thermochemical Analysis (IdealGas)

```python
from ase.thermochemistry import IdealGasThermo

vib_energies = vib.get_energies()   # eV, real part only

thermo = IdealGasThermo(
    vib_energies=vib_energies,
    potentialenergy=atoms.get_potential_energy(),
    atoms=atoms,
    geometry='nonlinear',    # 'monatomic', 'linear', 'nonlinear'
    symmetrynumber=1,        # molecular symmetry number (Cs=1, C2v=2, C3v=3…)
    spin=0,                  # total spin quantum number S (0=singlet)
)

T = 298.15   # K
P = 101325   # Pa (1 atm)

H  = thermo.get_enthalpy(temperature=T, verbose=True)      # eV
S  = thermo.get_entropy(temperature=T, pressure=P, verbose=True)  # eV/K
G  = thermo.get_gibbs_energy(temperature=T, pressure=P)    # eV
Cv = thermo.get_heat_capacity_V(temperature=T)             # eV/K
ZPE = thermo.get_ZPE_correction()                         # eV
```

### Reaction Thermodynamics

```python
# ΔG = G_products - G_reactants
dG = G_product - G_reactant1 - G_reactant2
print(f"ΔG = {dG:.3f} eV = {dG * 96.485:.1f} kJ/mol")
```

---

## IR Spectrum

```python
from ase.vibrations.infrared import Infrared

# Requires dipole moment from calculator (e.g., ORCA with IR keyword)
ir = Infrared(atoms, name='ir', delta=0.01)
ir.run()
ir.summary()

# Get intensities and frequencies
freqs = ir.get_frequencies()         # cm⁻¹
intens = ir.get_intensities()        # (km/mol) per mode

# Write spectrum
ir.write_spectra('ir_spectrum.dat', type='Gaussian', width=20)

# Plot
import matplotlib.pyplot as plt
spectrum = ir.get_spectrum(start=400, end=4000, width=20, type='Gaussian')
x, y = spectrum
plt.plot(x, y)
plt.xlabel('Wavenumber (cm⁻¹)')
plt.ylabel('Intensity (km/mol)')
plt.show()
```

> For ORCA IR: `orcasimpleinput='B3LYP def2-SVP NumFreq'` or use `Freq` with analytic gradients.

---

## Phonons — Periodic Solids

```python
from ase.phonons import Phonons

# Setup
ph = Phonons(
    atoms,          # bulk/periodic structure
    calc,           # calculator
    supercell=(2, 2, 2),   # supercell for finite displacements
    delta=0.05      # Å displacement
)

# Compute force constants (many DFT calls)
ph.run()

# Read results and impose symmetry
ph.read(acoustic=True)   # acoustic=True: enforce sum rules

# Phonon band structure
from ase.dft.kpoints import bandpath
path = atoms.cell.bandpath('GXMGR', npoints=100)
bs = ph.get_band_structure(path)
bs.plot(emin=-50, emax=200, filename='phonon_bands.png')

# Phonon DOS
dos = ph.get_dos(kpts=(20, 20, 20)).sample_grid(npts=200, width=10)
dos.plot(filename='phonon_dos.png')

# Debye-Waller factors, force constants
ph.get_force_constant()   # (3N, 3N) matrix

# Cleanup
ph.clean()
```

---

## VibrationsData API (modern, ASE 3.20+)

```python
from ase.vibrations import VibrationsData
import numpy as np

# Construct from atoms + Hessian
hessian = ...   # (3N, 3N) array in eV/Å²
vib_data = VibrationsData(atoms, hessian)

# Access
freqs = vib_data.get_frequencies()          # cm⁻¹
modes = vib_data.get_modes(all_atoms=True)  # (3N, N, 3)
zpe   = vib_data.get_zero_point_energy()    # eV

# Thermochemistry via VibrationsData
from ase.thermochemistry import IdealGasThermo
energies = vib_data.get_energies()
```

---

## Quick Reference: Frequency Ranges

| Bond/mode | Typical range (cm⁻¹) |
|-----------|----------------------|
| O-H stretch | 3200–3700 |
| N-H stretch | 3300–3500 |
| C-H stretch | 2800–3100 |
| C≡N stretch | 2100–2260 |
| C=O stretch | 1650–1800 |
| C=C stretch | 1500–1650 |
| C-N stretch | 1000–1350 |
| C-C stretch | 800–1200 |
| Skeletal bending | < 700 |

---

## Workflow: Full Thermochemistry from ORCA

```python
from ase.build import molecule
from ase.calculators.orca import ORCA
from ase.optimize import BFGS
from ase.vibrations import Vibrations
from ase.thermochemistry import IdealGasThermo

mol = molecule('H2O')

# 1. Optimize geometry
mol.calc = ORCA(charge=0, mult=1,
               orcasimpleinput='B3LYP def2-TZVP D3BJ TightSCF TightOpt')
BFGS(mol).run(fmax=0.01)
E0 = mol.get_potential_energy()

# 2. Vibrational frequencies (SP calculator — no Opt keyword)
mol.calc = ORCA(charge=0, mult=1,
               orcasimpleinput='B3LYP def2-TZVP D3BJ TightSCF')
vib = Vibrations(mol, name='h2o_vib', nfree=4)
vib.run()
vib_energies = vib.get_energies()

# 3. Thermochemistry
thermo = IdealGasThermo(
    vib_energies=vib_energies,
    potentialenergy=E0,
    atoms=mol,
    geometry='nonlinear',
    symmetrynumber=2,   # C2v water
    spin=0
)

print(f"ZPE  = {thermo.get_ZPE_correction():.4f} eV")
print(f"H    = {thermo.get_enthalpy(298):.4f} eV")
print(f"G    = {thermo.get_gibbs_energy(298, 101325):.4f} eV")
```
