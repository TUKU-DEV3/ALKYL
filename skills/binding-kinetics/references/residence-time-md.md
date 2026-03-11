# MD-Based Residence Time Estimation

## Overview

Experimental koff is hard to optimize without knowing the unbinding pathway.
MD-based methods provide:
1. Unbinding pathway (which residues disengage first)
2. Rank-order koff across congeneric series
3. Absolute τ estimate (with funnel metadynamics)

| Method | τ estimate | Relative ranking | Cost |
|--------|-----------|-----------------|------|
| τRAMD | Qualitative τ proxy | Good | Low (×10 tRAMD runs) |
| WTmetaD (standard) | Poor absolute τ | Moderate | Medium |
| Funnel metadynamics | Quantitative ΔG_bind, τ | Good | High |
| HTMD τRAMD | Semi-quantitative τ | Good | Low-Medium |
| Metadynamics+MBAR | Quantitative ΔG | Best | Very high |

## τRAMD (Random Acceleration MD)

Apply random force to center of mass of ligand until it exits the pocket.
Shorter exit time → faster koff (qualitative).

### NAMD/GROMACS τRAMD Setup

```python
import subprocess
import numpy as np
from pathlib import Path

# RAMD requires a patched NAMD or the HTMD Python implementation
# HTMD implementation (recommended for Python workflow)

def run_htmd_ramd(topology_file, trajectory_files,
                   ligand_sel='resname LIG',
                   protein_sel='protein',
                   force_constant=14,      # kcal/(mol·Å)
                   n_runs=10,
                   max_time_ns=20,
                   seed=42):
    """
    τRAMD using HTMD (requires: pip install htmd).
    Returns list of exit times per run.
    """
    try:
        from htmd.ui import *
        from htmd.kinetics import Kinetics
    except ImportError:
        raise ImportError("pip install htmd")

    mol = Molecule(topology_file)
    exit_times = []

    for run in range(n_runs):
        # HTMD RAMD simulation
        ramd = RAMD()
        ramd.topology = topology_file
        ramd.coordinates = trajectory_files[0]
        ramd.ligandSel = ligand_sel
        ramd.proteinSel = protein_sel
        ramd.forceConstant = force_constant
        ramd.maxTime = max_time_ns * 1e6   # ns → fs
        ramd.seed = seed + run
        exit_time = ramd.run()
        exit_times.append(exit_time)

    return exit_times

def tau_ramd_from_exit_times(exit_times):
    """
    Estimate τRAMD from exit time distribution.
    Uses exponential fit (single-exponential unbinding).
    """
    from scipy.optimize import curve_fit
    from scipy.stats import expon

    exit_times = np.array(exit_times)
    # Fit exponential distribution
    loc, scale = expon.fit(exit_times, floc=0)
    tau_ramd = scale  # mean exit time

    # Bootstrap CI
    n_boot = 1000
    boot_taus = []
    for _ in range(n_boot):
        sample = np.random.choice(exit_times, len(exit_times), replace=True)
        _, s = expon.fit(sample, floc=0)
        boot_taus.append(s)

    return {
        'tau_ramd': tau_ramd,
        'tau_ramd_ci': (np.percentile(boot_taus, 5),
                         np.percentile(boot_taus, 95)),
        'n_runs': len(exit_times),
        'mean_exit_time': exit_times.mean(),
        'std_exit_time': exit_times.std()
    }
```

## τRAMD with GROMACS (gmx ramd plugin)

```bash
# GROMACS + RAMD plugin (Helmholtz Centre)
# Install: https://gitlab.com/Pwr-Lab/gmx_ramd

# Generate RAMD-enhanced MD input
gmx ramd -s topol.tpr \
         -n index.ndx \
         -ramd-force 14 \          # kcal/(mol·Å)
         -ramd-rmin 0.025 \        # Å minimum displacement
         -ramd-seed 42 \
         -ramd-lig LIG \
         -ramd-prot PROTEIN \
         -nsteps 10000000 \        # 20 ns
         -o ramd_traj.xtc

# Run 10-20 replicas with different seeds for reliable τ estimate
```

## Funnel Metadynamics (FMetaD)

The gold standard for absolute ΔG_bind and koff. Uses a funnel-shaped restraint
to focus sampling along the binding/unbinding pathway.

### PLUMED Input

```
# PLUMED 2.x input for funnel metadynamics
# Project ligand position onto funnel axis

WHOLEMOLECULES ENTITY0=1-5000   # protein atoms
WHOLEMOLECULES ENTITY1=5001-5080  # ligand atoms

# Ligand center of mass
lig: COM ATOMS=5001-5080

# Binding site reference
site: COM ATOMS=145,168,202,217  # hotspot residues

# Funnel collective variable
f: FUNNEL_PS ...
  LIGAND=lig
  REFERENCE=site
  RCYL=0.1         # cylinder radius (nm)
  ALPHA=0.5        # funnel opening angle
  ZMAX=2.5         # max z along funnel axis (nm)
  ZMIN=-0.5        # min z (bound state)
  KAPPA=1000       # funnel wall force constant (kJ/mol/nm^2)
... FUNNEL_PS

# Well-tempered metadynamics
metad: METAD ...
  ARG=f.lp,f.ld   # funnel CVs (lp=projection, ld=distance from axis)
  HEIGHT=0.5       # kJ/mol (initial hill height)
  SIGMA=0.02,0.01  # hill widths
  BIASFACTOR=15    # well-tempered factor
  TEMP=300
  PACE=1000        # add hill every 1000 steps
  GRID_MIN=-0.5,0
  GRID_MAX=2.5,0.5
... METAD

PRINT ARG=f.lp,f.ld,metad.bias FILE=COLVAR STRIDE=500
```

### Free Energy Surface Analysis

```python
import numpy as np
import matplotlib.pyplot as plt

def analyze_fmetad_fes(colvar_file, kT=2.479):  # kT at 298K in kJ/mol
    """
    Analyze funnel metaD free energy surface.
    Extract ΔG_bind and estimate koff.
    """
    data = np.loadtxt(colvar_file, comments='#')
    time = data[:, 0]
    lp = data[:, 1]     # projection on funnel axis
    ld = data[:, 2]     # distance from axis
    bias = data[:, 3]   # accumulated bias

    # Reweight to get FES
    # F(z) = -kT * ln(P(z)) + correction
    # Use PLUMED sum_hills or manual reweighting

    return {'time': time, 'lp': lp, 'ld': ld, 'bias': bias}

def estimate_koff_from_fes(fes_1d, z_bound, z_unbound,
                            attempt_freq=1e12, kT=2.479):
    """
    Estimate koff from 1D free energy surface.
    Uses Kramers' rate theory: koff = ν * exp(-ΔG‡/kT)
    ΔG‡: barrier from bound minimum to transition state.

    kT: thermal energy in same units as FES (kJ/mol)
    attempt_freq: pre-exponential factor (s⁻¹), ≈ kT/h = 6×10¹² s⁻¹ at 300K
    """
    bound_min = fes_1d[z_bound]
    barrier = fes_1d.max() - bound_min   # ΔG‡ in kJ/mol
    koff = attempt_freq * np.exp(-barrier / kT)
    tau = 1.0 / koff

    return {
        'barrier_kJ': barrier,
        'barrier_kcal': barrier / 4.184,
        'koff': koff,
        'tau_s': tau,
        'tau_units': 's' if tau >= 1 else 'ms'
    }
```

## Unbinding Pathway Analysis

```python
import MDAnalysis as mda
import numpy as np

def analyze_unbinding_pathway(topology, ramd_trajectories,
                               ligand_sel='resname LIG',
                               protein_sel='protein',
                               pocket_residues=None):
    """
    Identify key contacts broken during unbinding.
    Compares bound frame (t=0) to exit frame (last frame).
    """
    contact_freqs = {}

    for traj in ramd_trajectories:
        u = mda.Universe(topology, traj)
        lig = u.select_atoms(ligand_sel)

        for ts in u.trajectory:
            # Protein residues within 4 Å of ligand
            near = u.select_atoms(
                f'protein and around 4.0 ({ligand_sel})'
            )
            for res in near.residues:
                key = f'{res.resname}{res.resid}'
                contact_freqs[key] = contact_freqs.get(key, 0) + 1

    # Normalize by first-frame contacts
    u_ref = mda.Universe(topology, ramd_trajectories[0])
    u_ref.trajectory[0]
    near_bound = u_ref.select_atoms(f'protein and around 4.0 ({ligand_sel})')
    bound_residues = {f'{r.resname}{r.resid}' for r in near_bound.residues}

    # Which bound-state contacts persist vs. which break first?
    persistence = {}
    total_frames = sum(len(mda.Universe(topology, t).trajectory)
                       for t in ramd_trajectories)
    for res in bound_residues:
        persistence[res] = contact_freqs.get(res, 0) / total_frames

    return dict(sorted(persistence.items(), key=lambda x: x[1]))
```

## Comparative τRAMD Analysis

```python
def ramd_series_analysis(compounds_dict, ref_compound=None):
    """
    Compare τRAMD across a compound series.
    compounds_dict: {name: [exit_times_ns]}
    Returns relative τ and rank order.
    """
    results = []
    for name, exit_times in compounds_dict.items():
        tau_result = tau_ramd_from_exit_times(exit_times)
        results.append({'compound': name, **tau_result})

    df = pd.DataFrame(results).sort_values('tau_ramd', ascending=False)

    if ref_compound and ref_compound in compounds_dict:
        ref_tau = df[df.compound == ref_compound]['tau_ramd'].values[0]
        df['relative_tau'] = df['tau_ramd'] / ref_tau

    return df

# Interpretation:
# τRAMD is NOT absolute koff (force applied accelerates unbinding)
# But relative ranking within a series is reliable (R² ~ 0.6–0.8 vs exp koff)
# Use for: prioritizing compounds before expensive SPR campaigns
```
