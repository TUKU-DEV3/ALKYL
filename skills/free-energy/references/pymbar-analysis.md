# pymbar — MBAR Analysis, Overlap, Convergence

pymbar 4.x (2023). MBAR estimator, overlap matrix, convergence diagnostics. Install: `pip install pymbar`.

---

## Core MBAR Usage

```python
import numpy as np
from pymbar import MBAR

# u_kln[k, l, n] = reduced potential u_l(x_n^k)
# = H_l(x_n^k) / (kB · T)
# k = state sample was drawn from
# l = state at which energy is evaluated
# n = sample index (within state k)

# Shape: (K, K, N_max) where N_max = max(N_k)
K = 12          # number of lambda windows
N_k = np.array([1000] * K)  # samples per window

# Build u_kln from simulation data
# (typically loaded from openmmtools NetCDF or GROMACS xvg)
u_kln = np.zeros((K, K, N_k.max()))

mbar = MBAR(u_kln, N_k, solver_protocol='robust')

# Compute free energy differences
results = mbar.compute_free_energy_differences()
# results['Delta_f'][i, j]  = ΔG(i→j) in kBT
# results['dDelta_f'][i, j] = uncertainty (1σ) in kBT

kBT = 0.5961   # kcal/mol at 298.15 K

dG_total = results['Delta_f'][0, -1] * kBT
ddG_total = results['dDelta_f'][0, -1] * kBT
print(f"ΔG = {dG_total:.3f} ± {ddG_total:.3f} kcal/mol")

# Profile: ΔG at each lambda
dG_profile = results['Delta_f'][0, :] * kBT
print("Lambda profile (kcal/mol):", dG_profile)
```

---

## Loading Data from Files

### From GROMACS .xvg files
```python
import re
import numpy as np

def read_gromacs_xvg(filename):
    """Read GROMACS dH/dλ or Δu output."""
    data = []
    with open(filename) as f:
        for line in f:
            if line.startswith('#') or line.startswith('@'):
                continue
            values = [float(x) for x in line.split()]
            data.append(values)
    return np.array(data)

# Read energy differences for each lambda window
u_kln_list = []
for k in range(K):
    data = read_gromacs_xvg(f'lambda_{k}/dhdl.xvg')
    u_kln_list.append(data[:, 1:])  # columns: time, dH/dl, Δu_0, Δu_1, ...
```

### From OpenMMTools NetCDF
```python
from openmmtools.multistate import MultiStateReporter
import netCDF4

reporter = MultiStateReporter('simulation.nc', open_mode='r')

# Direct energy matrix access
with netCDF4.Dataset('simulation.nc') as nc:
    energies = nc.variables['energies'][:]  # shape: (iter, replica, state)
    states = nc.variables['states'][:]      # replica → state assignment

# Rearrange to u_kln
# (use MultiStateSamplerAnalyzer for convenience)
```

---

## Overlap Matrix

```python
# Compute and visualize phase-space overlap
overlap_matrix = mbar.compute_overlap()['matrix']  # shape (K, K)

import matplotlib.pyplot as plt
import matplotlib.colors as colors

fig, ax = plt.subplots(figsize=(8, 7))
im = ax.imshow(overlap_matrix, cmap='Blues', vmin=0, vmax=1)
plt.colorbar(im, ax=ax, label='Overlap O_ij')

# Annotate each cell
for i in range(K):
    for j in range(K):
        val = overlap_matrix[i, j]
        color = 'white' if val > 0.5 else 'black'
        ax.text(j, i, f'{val:.2f}', ha='center', va='center',
                fontsize=8, color=color)

ax.set_xlabel('State j')
ax.set_ylabel('State i')
ax.set_title('MBAR Overlap Matrix')
plt.tight_layout()
plt.savefig('overlap_matrix.png', dpi=150, bbox_inches='tight')

# Check adjacent overlap (critical diagnostic)
for i in range(K-1):
    ov = overlap_matrix[i, i+1]
    status = 'OK' if ov > 0.03 else 'POOR — add λ window'
    print(f"O({i},{i+1}) = {ov:.3f}  {status}")
```

---

## Convergence Analysis

### Cumulative ΔG vs. Simulation Time

```python
def convergence_analysis(u_kln_full, N_k_full, fractions=None):
    """Compute ΔG using increasing fractions of data."""
    if fractions is None:
        fractions = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

    dG_vs_time = []
    ddG_vs_time = []

    for frac in fractions:
        N_k = (N_k_full * frac).astype(int)
        u_kln = u_kln_full[:, :, :N_k.max()].copy()
        for k in range(len(N_k)):
            u_kln[k, :, N_k[k]:] = 0.0

        mbar = MBAR(u_kln, N_k)
        res = mbar.compute_free_energy_differences()
        dG_vs_time.append(res['Delta_f'][0, -1] * kBT)
        ddG_vs_time.append(res['dDelta_f'][0, -1] * kBT)

    return fractions, dG_vs_time, ddG_vs_time

fracs, dG_t, ddG_t = convergence_analysis(u_kln, N_k)

plt.figure(figsize=(8, 4))
plt.errorbar(fracs, dG_t, yerr=ddG_t, marker='o', capsize=4)
plt.axhline(dG_t[-1], color='r', linestyle='--', label=f'Final: {dG_t[-1]:.2f} kcal/mol')
plt.xlabel('Fraction of simulation time used')
plt.ylabel('ΔG (kcal/mol)')
plt.title('Convergence: ΔG vs. simulation time')
plt.legend()
plt.tight_layout()
plt.savefig('convergence.png', dpi=150)
```

### Forward/Backward Convergence

```python
# Forward (first half) vs. backward (second half) estimates
N_half = N_k // 2

mbar_fwd = MBAR(u_kln[:, :, :N_half.max()], N_half)
mbar_bwd = MBAR(u_kln[:, :, N_half.max():], N_half)

dG_fwd = mbar_fwd.compute_free_energy_differences()['Delta_f'][0, -1] * kBT
dG_bwd = mbar_bwd.compute_free_energy_differences()['Delta_f'][0, -1] * kBT
hysteresis = abs(dG_fwd - dG_bwd)

print(f"ΔG (forward):  {dG_fwd:.3f} kcal/mol")
print(f"ΔG (backward): {dG_bwd:.3f} kcal/mol")
print(f"Hysteresis:    {hysteresis:.3f} kcal/mol")
if hysteresis > 0.5:
    print("WARNING: high hysteresis — simulation not converged")
```

---

## Autocorrelation Time

```python
from pymbar.timeseries import detect_equilibration, statistical_inefficiency

# For each lambda window, compute autocorrelation of ∂H/∂λ
for k in range(K):
    dhdl = u_kln[k, k, :N_k[k]]   # reduced potential at same state

    # Detect equilibration
    t_eq, g, N_eff = detect_equilibration(dhdl)
    # t_eq: equilibration index
    # g: statistical inefficiency (1 + 2*tau_int)
    # N_eff: effective number of uncorrelated samples

    print(f"Window {k}: t_eq={t_eq}, g={g:.1f}, N_eff={N_eff:.0f}")
    if N_eff < 50:
        print(f"  WARNING: only {N_eff:.0f} effective samples — need longer simulation")
```

---

## Expectations of Observables

```python
# MBAR can estimate expectation values of any observable
# Example: average RMSD at each lambda state

rmsd_values = ...  # shape (K, N_max) — RMSD for each frame
weights = mbar.compute_expectations(rmsd_values, N_k)
# Or use full MBAR reweighting:

A_kn = rmsd_values  # observable matrix
results = mbar.compute_expectations(A_kn, N_k, uncertainty_method='bootstrap')
print(f"⟨RMSD⟩_λ=0: {results['mu'][0]:.3f} ± {results['sigma'][0]:.3f} Å")
```

---

## BAR (Bennett Acceptance Ratio) — Pairwise

```python
from pymbar import BAR

# For a single pair of adjacent states (i, i+1)
# w_F: work values for forward perturbation (i→i+1, sampled at i)
# w_R: work values for reverse perturbation (i+1→i, sampled at i+1)

w_F = u_kln[i, i+1, :N_k[i]] - u_kln[i, i, :N_k[i]]   # kBT units
w_R = u_kln[i+1, i, :N_k[i+1]] - u_kln[i+1, i+1, :N_k[i+1]]

dG_bar, ddG_bar = BAR(w_F, w_R)
print(f"BAR ΔG({i}→{i+1}) = {dG_bar * kBT:.3f} ± {ddG_bar * kBT:.3f} kcal/mol")
```

---

## Summary of Diagnostics

| Diagnostic | Threshold | Action if failing |
|-----------|-----------|------------------|
| Adjacent overlap O_ij | > 0.03 | Add λ windows between i and j |
| N_eff per window | > 50 | Run longer simulation |
| Hysteresis fwd/bwd | < 0.5 kcal/mol | More sampling; HREX |
| ΔG plateaus at 50% time | stable | Otherwise: not converged |
| Cycle closure error | < 1 kcal/mol | Check force field; topology errors |
