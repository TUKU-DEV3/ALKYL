# ORCA — Practical Guide

ORCA 6.0 (2024). Free for academic use. Input: plain text `.inp`. Output: `.out` text file. Run via subprocess in Python workflows. See also: `chem_qm.py` (ALKYL script for input gen + output parsing).

---

## Input File Structure

```
! <keywords>          ← method line: functional, basis, job type, misc options
%pal nprocs N end     ← parallelization (N CPUs)
%maxcore M end        ← memory per core in MB (total = N × M)
%<block> ... end      ← optional blocks for fine control

* <coord_type> <charge> <multiplicity>
  <coordinates>
*
```

### Coordinate types
```
* xyz 0 1             ← inline XYZ (charge=0, spin=singlet)
* xyzfile 0 1 file.xyz   ← external XYZ file (recommended)
* int 0 1             ← Z-matrix
```

---

## Common Keywords

### Single Point
```
! B3LYP D3BJ def2-SVP TightSCF
! PBE0 D3BJ def2-TZVP TightSCF RIJCOSX def2/J
! ωB97X-D def2-TZVP TightSCF
! r2SCAN-3c            ← composite, includes D4+gCP+basis
```

### Geometry Optimization
```
! B3LYP D3BJ def2-SVP Opt TightSCF TightOpt
! r2SCAN-3c Opt
! B3LYP D3BJ def2-SVP Opt Freq   ← opt + frequencies in one job
```

### Frequency Calculation (must be at stationary point)
```
! B3LYP D3BJ def2-SVP Freq
%freq
  Temp 298.15     ← temperature for thermochemistry (K)
  Pressure 1.0    ← pressure (atm)
end
```

### Transition State Optimization
```
! B3LYP D3BJ def2-SVP OptTS Freq NumFreq
%geom
  Calc_Hess true   ← compute Hessian at start (recommended for TS)
  MaxIter 200
end
```

### TD-DFT (UV-Vis excited states)
```
! CAM-B3LYP def2-TZVP TightSCF
%tddft
  nroots 10        ← number of excited states
  maxdim 5         ← Davidson expansion space multiplier
  triplets false   ← singlet excitations only
end
```

### NMR Chemical Shifts
```
! B3LYP def2-TZVP NMR
%eprnmr
  Nuclei = all H {shift}     ← ¹H shifts
  Nuclei = all C {shift}     ← ¹³C shifts
  Gauge GIAO                 ← gauge-including AOs (default, most accurate)
end
```

### DLPNO-CCSD(T) Single Point
```
! DLPNO-CCSD(T) def2-TZVPP def2-TZVPP/C TightSCF
%mdci
  TCutPNO 1e-7    ← tighter PNO cutoff for improved accuracy
end
```

---

## Parallelization + Memory

```
%pal nprocs 8 end        ← 8 CPUs (use all available cores)
%maxcore 4000 end        ← 4 GB per core → 32 GB total

# Rule: maxcore × nprocs < 90% of available RAM
# For 64 GB RAM, 8 cores: maxcore 7000
```

---

## Full Example: Optimization → Single Point

```python
import subprocess
import os

def write_orca_input(name, xyz_file, method_opt, method_sp, charge=0, mult=1, nprocs=4):
    """Write opt then SP input files."""
    opt_inp = f"""\
! {method_opt} Opt TightSCF TightOpt
%pal nprocs {nprocs} end
%maxcore 4000
%geom MaxIter 200 end

* xyzfile {charge} {mult} {xyz_file}
"""
    sp_inp = f"""\
! {method_sp} TightSCF
%pal nprocs {nprocs} end
%maxcore 4000

* xyzfile {charge} {mult} {name}_opt.xyz
"""
    with open(f'{name}_opt.inp', 'w') as f: f.write(opt_inp)
    with open(f'{name}_sp.inp', 'w') as f: f.write(sp_inp)

write_orca_input(
    'lig', 'lig.xyz',
    method_opt='B3LYP D3BJ def2-SVP',
    method_sp='B3LYP D3BJ def2-TZVP RIJCOSX def2/J',
    nprocs=8
)

# Run optimization
r = subprocess.run(['orca', 'lig_opt.inp'], capture_output=True, text=True)
with open('lig_opt.out', 'w') as f: f.write(r.stdout)

# Run single point on optimized geometry
r = subprocess.run(['orca', 'lig_sp.inp'], capture_output=True, text=True)
with open('lig_sp.out', 'w') as f: f.write(r.stdout)
```

---

## Output Parsing

```python
def parse_orca_energy(outfile):
    """Extract final single-point energy (Hartree)."""
    energy = None
    with open(outfile) as f:
        for line in f:
            if 'FINAL SINGLE POINT ENERGY' in line:
                energy = float(line.split()[-1])
    return energy

def parse_orca_frequencies(outfile):
    """Extract vibrational frequencies (cm⁻¹)."""
    freqs = []
    in_freq_block = False
    with open(outfile) as f:
        for line in f:
            if 'VIBRATIONAL FREQUENCIES' in line:
                in_freq_block = True
            elif in_freq_block and 'cm**-1' in line:
                parts = line.split()
                if len(parts) >= 3:
                    freqs.append(float(parts[1]))
            elif in_freq_block and line.strip() == '':
                if freqs:
                    in_freq_block = False
    return freqs

def parse_orca_thermo(outfile):
    """Extract thermochemistry: ZPE, H, G (Hartree)."""
    thermo = {}
    with open(outfile) as f:
        for line in f:
            if 'Zero point energy' in line:
                thermo['ZPE'] = float(line.split()[-2])
            elif 'Total enthalpy' in line:
                thermo['H'] = float(line.split()[-2])
            elif 'Final Gibbs free energy' in line:
                thermo['G'] = float(line.split()[-2])
    return thermo

def check_convergence(outfile):
    """Check if optimization converged."""
    with open(outfile) as f:
        content = f.read()
    return 'THE OPTIMIZATION HAS CONVERGED' in content

def check_ts(outfile):
    """Check if TS has exactly one imaginary frequency."""
    freqs = parse_orca_frequencies(outfile)
    n_imag = sum(1 for f in freqs if f < 0)
    return n_imag == 1, n_imag

# Use ALKYL script for comprehensive parsing:
# python chem_qm.py --parse lig_opt.out
# python chem_qm.py --parse lig_opt.out --parse-ir
```

---

## Output Files

| File | Content |
|------|---------|
| `job.out` | Main output (energies, geometry, properties) |
| `job_opt.xyz` | Optimized geometry (XYZ format) |
| `job.hess` | Hessian matrix (for frequency restart) |
| `job.gbw` | Wavefunction/MO coefficients |
| `job.densities` | Electron density files |
| `job.engrad` | Energy + gradient |
| `job.trj` | Optimization trajectory (XYZ) |
| `job_property.txt` | Machine-readable properties |

---

## Solvation

```
! B3LYP D3BJ def2-SVP Opt CPCM(water)

# SMD solvation (for ΔG_solv)
! B3LYP D3BJ def2-TZVP CPCM
%cpcm
  smd true
  solvent "water"
end

# Common solvent strings: "water", "methanol", "dmso", "thf",
#   "acetonitrile", "chloroform", "benzene", "hexane"
```

---

## Constraints

```
%geom
  Constraints
    {B 0 1 1.5 C}      ← fix bond between atoms 0 and 1 at 1.5 Å
    {A 0 1 2 120.0 C}  ← fix angle at 120°
    {D 0 1 2 3 0.0 C}  ← fix dihedral at 0°
    {C 0 C}            ← fix Cartesian position of atom 0
  end
end
```

---

## Common Issues

| Problem | Cause | Fix |
|---------|-------|-----|
| SCF does not converge | Bad geometry / charged system | `SlowConv` keyword; use `! SlowConv` |
| Opt does not converge | Flat PES, wrong initial geom | Pre-optimize with xTB; increase `MaxIter` |
| Imaginary frequency at minimum | Opt stopped early | Restart from `job_opt.xyz` with `TightOpt` |
| Multiple imaginary freqs at TS | Wrong TS guess | Confirm with `IRC` calculation |
| Memory error | maxcore too high | Reduce `%maxcore`; check available RAM |
| `ORCA: command not found` | ORCA not in PATH | `export PATH=/path/to/orca:$PATH` |
| TD-DFT dark state | Only state of different symmetry | Check `triplets true`; try more roots |
