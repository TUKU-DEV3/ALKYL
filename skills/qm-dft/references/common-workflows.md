# Common QM Workflows — Recipes

Standard calculation protocols for drug-like molecules. Always cite the method/basis in results.

---

## 1. Geometry Opt → Frequency → Single Point (Standard Protocol)

```bash
# Step 1: Pre-optimize with xTB (fast)
xtb mol.xyz --opt tight --gfn 2 --alpb water
mv xtbopt.xyz mol_preopt.xyz

# Step 2: DFT optimization
cat > opt.inp << 'EOF'
! B3LYP D3BJ def2-SVP Opt TightSCF TightOpt
%pal nprocs 8 end
%maxcore 4000
%geom MaxIter 300 end

* xyzfile 0 1 mol_preopt.xyz
EOF
orca opt.inp > opt.out

# Step 3: Verify no imaginary frequencies
cat > freq.inp << 'EOF'
! B3LYP D3BJ def2-SVP Freq
%pal nprocs 8 end
%maxcore 4000
%freq Temp 298.15 end

* xyzfile 0 1 opt_opt.xyz
EOF
orca freq.inp > freq.out

# Step 4: Single-point energy with larger basis
cat > sp.inp << 'EOF'
! B3LYP D3BJ def2-TZVP TightSCF RIJCOSX def2/J
%pal nprocs 8 end
%maxcore 4000

* xyzfile 0 1 opt_opt.xyz
EOF
orca sp.inp > sp.out
```

```python
# Parse results
e_opt = parse_orca_energy('opt.out')
freqs = parse_orca_frequencies('freq.out')
thermo = parse_orca_thermo('freq.out')
e_sp = parse_orca_energy('sp.out')

n_imag = sum(1 for f in freqs if f < 0)
assert n_imag == 0, f"ERROR: {n_imag} imaginary frequencies — not a true minimum"

# Final Gibbs free energy (single-point energy + thermal correction from freq)
G_thermal_corr = thermo['G'] - e_opt  # correction at opt level
G_final = e_sp + G_thermal_corr       # Hartree
print(f"G(final) = {G_final:.8f} Eh = {G_final * 627.51:.2f} kcal/mol")
```

---

## 2. Reaction Barrier (Transition State)

```bash
# Step 1: Guess TS geometry (manually or from NEB/xTB path)
# xTB pathway scan:
xtb reactant.xyz --path product.xyz --input scan.inp

# Step 2: ORCA TS optimization
cat > ts.inp << 'EOF'
! B3LYP D3BJ def2-SVP OptTS Freq TightSCF
%pal nprocs 8 end
%maxcore 4000
%geom
  Calc_Hess true    # compute initial Hessian
  MaxIter 200
end
%freq Temp 298.15 end

* xyzfile 0 1 ts_guess.xyz
EOF
orca ts.inp > ts.out
```

```python
freqs = parse_orca_frequencies('ts.out')
n_imag = sum(1 for f in freqs if f < 0)
assert n_imag == 1, f"TS has {n_imag} imaginary freqs (expected 1)"
print(f"Imaginary freq: {min(f for f in freqs if f < 0):.1f} cm⁻¹")

# Reaction energetics
e_reactant = parse_orca_energy('reactant_sp.out')
e_ts = parse_orca_energy('ts_sp.out')
e_product = parse_orca_energy('product_sp.out')

Ea = (e_ts - e_reactant) * 627.51   # kcal/mol
dH = (e_product - e_reactant) * 627.51

print(f"Ea = {Ea:.1f} kcal/mol")
print(f"ΔH = {dH:.1f} kcal/mol")
```

### IRC (Intrinsic Reaction Coordinate) Validation
```
! B3LYP D3BJ def2-SVP IRC TightSCF
%irc
  MaxIter 40
  StepSize 0.1
  Direction both      # forward and backward
end

* xyzfile 0 1 ts.xyz
```

---

## 3. RESP Charges (for Force Field Parameterization)

```bash
# ORCA ESP calculation
cat > resp.inp << 'EOF'
! HF 6-31G* TightSCF
%pal nprocs 4 end
%maxcore 4000

%output
  Print[P_Hirshfeld] 1
end

%elprop
  Dipole true
end

* xyzfile 0 1 mol_opt.xyz
EOF
orca resp.inp > resp.out

# Generate CHELPG charges (directly in ORCA)
cat > chelpg.inp << 'EOF'
! HF 6-31G* CHELPG TightSCF
%pal nprocs 4 end
%maxcore 4000

* xyzfile 0 1 mol_opt.xyz
EOF
orca chelpg.inp > chelpg.out
# Look for "CHELPG Charges" in output
```

```python
def parse_chelpg(outfile):
    """Parse CHELPG charges from ORCA output."""
    charges = []
    in_block = False
    with open(outfile) as f:
        for line in f:
            if 'CHELPG Charges' in line:
                in_block = True
                next(f)  # skip dashes
            elif in_block:
                parts = line.strip().split()
                if len(parts) == 3 and parts[0].isdigit():
                    charges.append(float(parts[2]))
                elif line.strip() == '' and charges:
                    break
    return charges

charges = parse_chelpg('chelpg.out')
print(f"Sum of charges: {sum(charges):.4f}")  # should be integer
```

---

## 4. UV-Vis Spectrum (TD-DFT)

```bash
cat > tddft.inp << 'EOF'
! CAM-B3LYP def2-TZVP TightSCF
%pal nprocs 8 end
%maxcore 4000

%tddft
  nroots 20
  maxdim 5
  triplets false
end

* xyzfile 0 1 mol_opt.xyz
EOF
orca tddft.inp > tddft.out
```

```python
import re
import numpy as np
import matplotlib.pyplot as plt

def parse_tddft(outfile):
    """Parse excited states from ORCA TD-DFT output."""
    states = []
    with open(outfile) as f:
        content = f.read()

    # State blocks
    pattern = r'STATE\s+(\d+):\s+E=\s*([\d.]+)\s+au\s+([\d.]+)\s+eV\s+([\d.]+)\s+nm\s+f=([\d.]+)'
    for m in re.finditer(pattern, content):
        states.append({
            'state': int(m.group(1)),
            'E_au': float(m.group(2)),
            'E_eV': float(m.group(3)),
            'wavelength_nm': float(m.group(4)),
            'osc_strength': float(m.group(5)),
        })
    return states

states = parse_tddft('tddft.out')

# Broadened spectrum (Gaussian broadening)
def broaden_spectrum(states, wmin=200, wmax=800, sigma=20):
    """Simulate UV-Vis spectrum with Gaussian line shapes."""
    wavelengths = np.linspace(wmin, wmax, 1000)
    epsilon = np.zeros_like(wavelengths)
    for s in states:
        lam = s['wavelength_nm']
        f = s['osc_strength']
        if wmin <= lam <= wmax:
            # Convert f to ε (approximate): ε_max = 2.174e8 * f / (σ * √(2π))
            epsilon += f * np.exp(-0.5 * ((wavelengths - lam) / sigma)**2)
    return wavelengths, epsilon

lam, eps = broaden_spectrum(states)
plt.plot(lam, eps)
plt.xlabel('Wavelength (nm)')
plt.ylabel('Oscillator strength (broadened)')
plt.title('Simulated UV-Vis spectrum')
plt.savefig('uvvis.png', dpi=150, bbox_inches='tight')
```

---

## 5. NMR Shift Prediction

```bash
# ORCA NMR (GIAO method)
cat > nmr.inp << 'EOF'
! B3LYP def2-TZVP NMR TightSCF
%pal nprocs 8 end
%maxcore 4000

%eprnmr
  Nuclei = all H {shift}
  Nuclei = all C {shift}
  Gauge GIAO
end

* xyzfile 0 1 mol_opt.xyz
EOF
orca nmr.inp > nmr.out
```

```python
def parse_orca_nmr(outfile):
    """Parse NMR shielding tensors from ORCA output."""
    shifts = {}
    in_block = False
    with open(outfile) as f:
        for line in f:
            if 'CHEMICAL SHIELDING SUMMARY' in line:
                in_block = True
            elif in_block and re.match(r'\s+\d+', line):
                parts = line.split()
                if len(parts) >= 5:
                    idx = int(parts[0])
                    atom = parts[1]
                    sigma_iso = float(parts[4])  # isotropic shielding
                    shifts[idx] = (atom, sigma_iso)
            elif in_block and 'CHEMICAL SHIELDING' not in line and line.strip() == '':
                if shifts:
                    break

    # Reference values for B3LYP/def2-TZVP (method-dependent!)
    refs = {'H': 31.7, 'C': 183.6}
    results = []
    for idx, (atom, sigma) in shifts.items():
        if atom in refs:
            delta = refs[atom] - sigma
            results.append((idx, atom, delta))
    return results

nmr_results = parse_orca_nmr('nmr.out')
for idx, atom, delta in nmr_results:
    print(f"Atom {idx} ({atom}): δ = {delta:.2f} ppm")
```

---

## 6. Conformer Ranking with xTB → DFT

```python
import subprocess
from pathlib import Path

def rank_conformers_dft(crest_file, method='B3LYP D3BJ def2-SVP', nprocs=8, top_n=10):
    """Rank CREST conformers with DFT single points."""
    conformers = parse_crest_ensemble(crest_file)
    energies = []

    for i, conf in enumerate(conformers[:top_n]):
        xyz_file = f'conf_{i}.xyz'
        # Write XYZ
        with open(xyz_file, 'w') as f:
            f.write(f"{conf['n_atoms']}\n\n")
            for atom, x, y, z in conf['coords']:
                f.write(f"{atom:3s} {x:12.6f} {y:12.6f} {z:12.6f}\n")

        # ORCA input
        inp = f"""! {method} TightSCF
%pal nprocs {nprocs} end
%maxcore 4000

* xyzfile 0 1 {xyz_file}
"""
        inp_file = f'conf_{i}.inp'
        with open(inp_file, 'w') as f:
            f.write(inp)

        r = subprocess.run(['orca', inp_file], capture_output=True, text=True)
        with open(f'conf_{i}.out', 'w') as f:
            f.write(r.stdout)

        e = parse_orca_energy(f'conf_{i}.out')
        energies.append((i, e))

    # Sort by DFT energy
    energies.sort(key=lambda x: x[1])
    e0 = energies[0][1]
    print("DFT conformer ranking:")
    for rank, (i, e) in enumerate(energies):
        dE = (e - e0) * 627.51
        print(f"Rank {rank+1}: conf_{i} ΔE = {dE:.2f} kcal/mol")

    return energies
```

---

## 7. Reaction Free Energy (ΔG)

```python
HARTREE_TO_KCAL = 627.51

def reaction_dG(reactant_out, product_out):
    """Compute ΔG for A → B using thermochemistry from ORCA Freq jobs."""
    G_r = parse_orca_thermo(reactant_out)['G']
    G_p = parse_orca_thermo(product_out)['G']
    dG = (G_p - G_r) * HARTREE_TO_KCAL
    print(f"ΔG = {dG:.2f} kcal/mol")
    return dG

# For multi-step: ΔG = G(products) - G(reactants)
# Include explicit water/proton if pH matters
# For proton: G(H+) ≈ G(H+ in water) ≈ -6.28 kcal/mol (standard state correction)
```

---

## Best Practices Checklist

```
Before calculation:
  ✓ Pre-optimize geometry with xTB --gfn 2 --opt
  ✓ Check charge and spin multiplicity
  ✓ Choose functional appropriate for task (see dft-theory.md)
  ✓ Include D3BJ dispersion unless justified otherwise

After optimization:
  ✓ Check all vibrational frequencies are real (no imaginary)
  ✓ Verify geometry is chemically sensible (bond lengths, angles)
  ✓ For TS: exactly 1 imaginary frequency + confirm with IRC

After single point:
  ✓ Check SCF convergence (search "SCF CONVERGED" in output)
  ✓ Report level of theory: method/basis//opt_method/opt_basis

Reporting:
  ✓ Energies: relative, not absolute (ΔE or ΔG)
  ✓ Always specify temperature for G (default 298.15 K)
  ✓ Units: kcal/mol or kJ/mol; never Hartree in final results
```
