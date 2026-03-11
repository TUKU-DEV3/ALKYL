# Coarse-Grained Theory

## Resolution Levels

| Level | Atoms/bead | Timestep | Timescale | Example FF |
|-------|-----------|----------|-----------|-----------|
| United-atom (UA) | 1 (no H) | 2 fs | ns–µs | GROMOS, OPLS-UA |
| Coarse-grained (CG) | 3–5 heavy | 20–40 fs | µs–ms | MARTINI, SDK |
| Ultra-CG / mesoscale | 10–100+ | 100 fs–ps | ms–s | DPD, one-bead lipid |
| Brownian dynamics | — | µs | s–min | Diffusion-limited |

**MARTINI 3** (2021) is the current standard for biomolecular CG simulations.

## Mapping Schemes

A mapping defines which atoms constitute each CG bead.
Principle: beads should be chemically meaningful (e.g., one bead per sugar ring).

```python
# Example MARTINI 3 mapping for alanine sidechain
# N-Cα-C=O backbone → 1 backbone bead (BB)
# Cβ sidechain → 1 SC1 bead (for Ala: backbone only)

# Tools for mapping:
# - CGSmiles: SMILES-like notation for CG molecules
# - Polyply: automated mapping generator
# - martinize2 --mapping for proteins
```

### Center of Geometry (COG) vs Center of Mass (COM)
```python
import numpy as np

def bead_position_cog(atom_positions):
    """Center of geometry — simplest mapping."""
    return atom_positions.mean(axis=0)

def bead_position_com(atom_positions, atom_masses):
    """Center of mass — physically more rigorous."""
    return np.average(atom_positions, weights=atom_masses, axis=0)
```

## MARTINI 3 Bead Types

MARTINI 3 has a systematic bead type hierarchy based on partitioning free energies.

### Regular beads (R = regular, S = small, T = tiny)
```
R-beads: 4 heavy atoms  (default)
S-beads: 3 heavy atoms  (prefix S, e.g. SP1)
T-beads: 2 heavy atoms  (prefix T, e.g. TC5)
```

### Polarity scale (0=nonpolar → 6=polar, d=donor, a=acceptor)
```
C1, C2, C3, C4, C5, C6   → nonpolar (hydrophobic), increasing polarity
N1–N6                      → intermediate (no H-bond)
P1–P6                      → polar (H-bond capable)
Q1–Q5d/a                   → charged

Examples:
  C1: alkane (dodecane-like)
  C5: aromatic (benzene-like)
  N3a: ester oxygen (non-HBD)
  P2: amide-like
  P4: water (4 water molecules)
  Qd: positively charged (lysine)
  Qa: negatively charged (glutamate)
```

## Parameterization Methods

### 1. MARTINI 3 Official Parameters (preferred)
Use pre-parameterized molecules from the MARTINI library.

```bash
# Check MARTINI 3 small molecule library
# https://cgmartini.nl/index.php/downloads/molecules

# Install: pip install martini3-small-molecule-lib
```

### 2. Boltzmann Inversion (BI)
Derive bonded parameters from all-atom reference distributions:
```
V_bonded(r) = -kT ln P_ref(r)
```
```python
import numpy as np
from scipy.stats import gaussian_kde

def boltzmann_inversion(values, T=300.0, n_bins=100):
    """
    Derive CG potential from AA reference distribution.
    values: 1D array of bond lengths, angles, or dihedrals from AA traj
    Returns: bin_centers, potential (kJ/mol)
    """
    kB = 8.314e-3   # kJ/(mol·K)
    kT = kB * T

    # KDE to get smooth probability density
    kde = gaussian_kde(values, bw_method='silverman')
    x = np.linspace(values.min(), values.max(), n_bins)
    P = kde(x)
    P = np.clip(P, 1e-10, None)   # avoid log(0)

    V = -kT * np.log(P)
    V -= V.min()   # set minimum to 0
    return x, V

def fit_harmonic(x, V):
    """Fit Boltzmann-inverted potential to harmonic V = 0.5*k*(x-x0)^2."""
    from scipy.optimize import curve_fit
    harmonic = lambda x, k, x0: 0.5 * k * (x - x0)**2
    popt, _ = curve_fit(harmonic, x, V, p0=[100, x.mean()])
    k, x0 = popt
    return k, x0   # force constant (kJ/mol/nm² or /rad²), equilibrium value
```

### 3. Force Matching / Multiscale CG (MS-CG)
Minimize difference between CG forces and projected AA forces.
Implemented in VOTCA software:
```bash
# VOTCA-CSG
csg_map --top aa.tpr --trj aa.xtc --cg mapping.xml --out cg.xtc
csg_stat --options ibm.xml    # iterative Boltzmann matching
```

## Time Scaling in MARTINI

MARTINI CG dynamics are faster than real-time due to smoothed energy landscape.
Empirical correction: **CG time ≈ 4× real time** (Marrink 2004).

```python
def effective_time(cg_time_ns, scaling=4.0):
    """Convert CG simulation time to effective real time."""
    return cg_time_ns * scaling

# BUT: this factor is NOT universal
# Membrane lateral diffusion: ×4 correction reasonable
# Protein conformational dynamics: scaling factor uncertain (×2–10)
# Always state both CG time and effective time in publications
```

## MARTINI vs All-Atom: When to Choose

| Use MARTINI when | Use All-atom when |
|-----------------|------------------|
| Lipid bilayer self-assembly | Sub-Å precision needed |
| Protein-membrane insertion | Explicit water H-bonds critical |
| µs+ conformational sampling | QM/MM interface needed |
| Large assemblies (>10k atoms) | Binding free energy (FEP) |
| Membrane protein lateral diffusion | NMR chemical shifts |
| Lipid nanoparticle structure | Exact koff (kinetics) |

## Implicit Solvent CG (DRY MARTINI)

Removes explicit water beads — faster but loses water-mediated interactions.
```
# DRY MARTINI: all beads rescaled to account for implicit water
# Useful for: even larger systems, faster equilibration
# Not recommended for: membrane-water interface analysis, ion permeation
```
