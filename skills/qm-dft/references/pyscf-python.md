# PySCF — Python-Native Quantum Chemistry

PySCF 2.7 (2024). Pure Python with NumPy/SciPy. No external QM binary required. Supports HF, DFT, MP2, CCSD(T), CASSCF, NMR, and more. Scales to ~500 atoms for DFT.

```bash
pip install pyscf
pip install pyscf-dftd3   # D3 dispersion correction
conda install -c conda-forge pyscf
```

---

## Molecule Setup

```python
from pyscf import gto, scf, dft

# Define molecule
mol = gto.Mole()
mol.atom = """
C   0.000   0.000   0.000
H   0.000   0.000   1.089
H   1.026   0.000  -0.363
H  -0.513   0.888  -0.363
H  -0.513  -0.888  -0.363
"""                              # Angstrom by default
mol.basis = 'def2-SVP'
mol.charge = 0
mol.spin = 0                     # 2S (singlet=0, doublet=1, triplet=2)
mol.symmetry = False
mol.verbose = 3                  # 0=silent, 5=debug
mol.build()

# From SMILES (via RDKit → XYZ)
from rdkit import Chem
from rdkit.Chem import AllChem
rdmol = Chem.AddHs(Chem.MolFromSmiles('CCO'))
AllChem.EmbedMolecule(rdmol, AllChem.ETKDGv3())
AllChem.MMFFOptimizeMolecule(rdmol)

xyz_block = Chem.MolToXYZBlock(rdmol)
mol = gto.Mole()
mol.atom = xyz_block
mol.basis = 'def2-SVP'
mol.build()
```

---

## HF and DFT

```python
from pyscf import scf, dft

# Restricted HF
mf = scf.RHF(mol)
mf.conv_tol = 1e-9
e_hf = mf.kernel()
print(f"HF energy: {e_hf:.8f} Eh")

# Restricted DFT (KS-DFT)
mf = dft.RKS(mol)
mf.xc = 'B3LYP'
mf.grids.level = 4        # integration grid quality (3=medium, 4=fine, 5=very fine)
e_dft = mf.kernel()

# With dispersion (D3BJ)
from pyscf import dftd3
mf = dft.RKS(mol).density_fit()
mf.xc = 'B3LYP'
mf = dftd3.dftd3(mf, xc='b3lyp', version='d3bj')  # attach D3BJ
e_d3 = mf.kernel()

# Density fitting (RI-J approximation, ~10× speedup)
mf = dft.RKS(mol).density_fit()
mf.xc = 'PBE0'
mf.with_df.auxbasis = 'def2-SVP-jkfit'
e = mf.kernel()

# Unrestricted (open-shell)
mf = dft.UKS(mol)
mf.xc = 'B3LYP'
e = mf.kernel()
```

---

## Basis Set Handling

```python
mol.basis = 'def2-SVP'       # Ahlrichs (recommended)
mol.basis = '6-31G*'         # Pople
mol.basis = 'cc-pVTZ'        # Dunning
mol.basis = 'aug-cc-pVDZ'    # with diffuse

# Mixed basis sets
mol.basis = {
    'C': 'def2-TZVP',
    'H': 'def2-SVP',
    'N': 'def2-TZVP',
    'O': 'def2-TZVP',
}
mol.build()
```

---

## Post-HF Correlation Methods

```python
from pyscf import mp, cc

# MP2
mf = scf.RHF(mol).run()
mp2 = mp.MP2(mf)
e_corr, t2 = mp2.kernel()
e_mp2 = mf.e_tot + e_corr
print(f"MP2 energy: {e_mp2:.8f} Eh")

# CCSD
mycc = cc.CCSD(mf)
mycc.kernel()
e_ccsd = mycc.e_tot
print(f"CCSD energy: {e_ccsd:.8f} Eh")

# CCSD(T) perturbative triples
from pyscf.cc import ccsd_t
e_t = ccsd_t.kernel(mycc, mycc.ao2mo())
e_ccsdt = e_ccsd + e_t
print(f"CCSD(T) energy: {e_ccsdt:.8f} Eh")
```

---

## Geometry Optimization (via geomeTRIC)

```python
# PySCF does not include optimizer; use geomeTRIC
pip install geometric

from pyscf import gto, dft
from pyscf.geomopt.geometric_solver import optimize

mol = gto.Mole(atom='...', basis='def2-SVP').build()
mf = dft.RKS(mol)
mf.xc = 'B3LYP'

mol_eq = optimize(mf)           # returns optimized Mole object
print(mol_eq.atom_coords())     # Bohr
print(mol_eq.atom_coords() * 0.529177)  # Å
```

---

## Molecular Properties

### Dipole Moment
```python
mf = dft.RKS(mol).run()
dm = mf.make_rdm1()            # density matrix
dip = mf.dip_moment()          # Debye, array [x, y, z]
print(f"μ = {np.linalg.norm(dip):.3f} D")
```

### Mulliken and Löwdin Population Analysis
```python
from pyscf.lo import pop_analysis

dm = mf.make_rdm1()
# Mulliken
charges_mulliken = mol.atom_charges() - mf.mulliken_pop(dm)[1]

# Löwdin
charges_lowdin = mf.mulliken_pop_meta_lowdin_ao(dm)[1]
```

### ESP Charges (for RESP-like)
```python
from pyscf.tools import cubegen
from pyscf import df

# Write ESP to cube file
mf = dft.RKS(mol).run()
cubegen.mep(mol, 'esp.cube', mf.make_rdm1(), nx=60, ny=60, nz=60)
# Then fit with external RESP code (e.g., resp package or antechamber)
```

---

## NMR Chemical Shifts

```python
from pyscf.prop import nmr

mf = dft.RKS(mol)
mf.xc = 'B3LYP'
mf.basis = 'def2-TZVP'  # larger basis for NMR
mf.kernel()

msc = nmr.RKS(mf)
shielding = msc.kernel()   # isotropic shielding tensor (ppm)

# Convert to chemical shifts (relative to TMS reference)
# δ(H) ≈ σ(TMS,H) - σ(H); σ(TMS,H) ≈ 31.7 ppm at B3LYP/def2-TZVP
TMS_H_ref = 31.7   # calibrated reference (method-dependent)
TMS_C_ref = 183.6  # for ¹³C

for i, atom in enumerate(mol.atom_symbols()):
    sigma_iso = (shielding[i][0,0] + shielding[i][1,1] + shielding[i][2,2]) / 3
    if atom == 'H':
        delta = TMS_H_ref - sigma_iso
        print(f"H{i}: δ = {delta:.2f} ppm")
    elif atom == 'C':
        delta = TMS_C_ref - sigma_iso
        print(f"C{i}: δ = {delta:.2f} ppm")
```

---

## TD-DFT (Excited States)

```python
from pyscf import tdscf

mf = dft.RKS(mol)
mf.xc = 'CAM-B3LYP'    # range-separated for excited states
mf.kernel()

td = tdscf.TDDFT(mf)
td.nstates = 10
td.kernel()

# Excitation energies and oscillator strengths
for i, (e, f) in enumerate(zip(td.e, td.oscillator_strength())):
    e_ev = e * 27.2114   # Hartree → eV
    lam = 1240 / e_ev    # nm
    print(f"State {i+1}: E={e_ev:.3f} eV ({lam:.0f} nm), f={f:.4f}")
```

---

## Solvation (DDCOSMO / PCM)

```python
from pyscf import solvent

# Polarizable Continuum Model
mf = dft.RKS(mol)
mf.xc = 'B3LYP'
mf = mf.ddCOSMO()       # domain-decomposed COSMO
mf.with_solvent.eps = 78.4  # water dielectric constant
mf.kernel()

# PCM
mf = mf.PCM()
mf.with_solvent.method = 'IEF-PCM'   # or 'C-PCM', 'SS(V)PE'
mf.with_solvent.eps = 78.4
```

---

## Key Notes

- PySCF uses **atomic units (Bohr, Hartree)** internally; inputs in Å are auto-converted
- `mol.spin` = **2S** (number of unpaired electrons): singlet=0, doublet=1, triplet=2
- DFT grids: `mf.grids.level = 4` is usually sufficient; level 5 for benchmarks
- NMR reference values (σ_TMS) are **method and basis set dependent** — must recalibrate or use literature values for the exact method used
- No geometry optimizer built-in → use `geomeTRIC` (preferred) or `berny`
- For large molecules (>200 atoms), use density fitting (`mf.density_fit()`) and linear-scaling methods
