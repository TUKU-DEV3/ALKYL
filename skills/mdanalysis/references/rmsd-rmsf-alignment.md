# MDAnalysis — RMSD, RMSF & Alignment

## RMSD — Root Mean Square Deviation

Measures structural deviation from a reference frame over time.

```python
import MDAnalysis as mda
from MDAnalysis.analysis import rms

u = mda.Universe('system.prmtop', 'traj.nc')

# Simple RMSD to first frame
rmsd = rms.RMSD(
    u,
    select='backbone',             # selection for fitting + RMSD
    groupselections=[              # extra groups to track (no fitting)
        'protein and name CA',
        'resname LIG'
    ],
    ref_frame=0                    # reference frame index (default: 0)
)
rmsd.run()

# Results: array of shape (n_frames, 3 + len(groupselections))
# Columns: [frame, time(ps), RMSD_select, RMSD_group1, RMSD_group2, ...]
data = rmsd.results.rmsd
frames    = data[:, 0]
times     = data[:, 1]   # ps
rmsd_bb   = data[:, 2]   # backbone RMSD
rmsd_ca   = data[:, 3]   # CA RMSD
rmsd_lig  = data[:, 4]   # ligand RMSD

# Plot
import matplotlib.pyplot as plt
plt.plot(times, rmsd_bb, label='Backbone')
plt.plot(times, rmsd_ca, label='CA')
plt.xlabel('Time (ps)')
plt.ylabel('RMSD (Å)')
plt.legend()
plt.show()
```

### RMSD to external reference

```python
ref = mda.Universe('reference.pdb')  # crystal structure

rmsd = rms.RMSD(
    u,
    reference=ref,
    select='protein and name CA',
    ref_frame=0
)
rmsd.run()
```

---

## RMSF — Root Mean Square Fluctuation

Per-atom (or per-residue) measure of mobility over the trajectory.
**Always align to average structure first.**

```python
from MDAnalysis.analysis import rms, align

u   = mda.Universe('system.prmtop', 'traj.nc')
ref = mda.Universe('system.prmtop', 'traj.nc')  # same universe for average

# Step 1: compute average structure
avg = align.AverageStructure(
    u, ref,
    select='protein and name CA',
    ref_frame=0
).run()
avg_universe = avg.results.universe

# Step 2: align trajectory to average
aligner = align.AlignTraj(
    u,
    avg_universe,
    select='protein and name CA',
    in_memory=True    # load all frames into RAM — use for small traj
).run()

# Step 3: calculate RMSF
ca = u.select_atoms('protein and name CA')
rmsf = rms.RMSF(ca).run()

rmsf_values = rmsf.results.rmsf   # Å per CA atom

# Plot vs residue number
plt.plot(ca.resids, rmsf_values)
plt.xlabel('Residue ID')
plt.ylabel('RMSF (Å)')
plt.show()

# Export as B-factors for visualization
u.add_TopologyAttr('tempfactors')
for res, val in zip(ca.residues, rmsf_values):
    res.atoms.tempfactors = val
u.select_atoms('protein').write('rmsf.pdb')
```

---

## Structural Alignment

### Align single frame

```python
from MDAnalysis.analysis.align import alignto

# Align mobile to reference in-place
rmsd_before, rmsd_after = alignto(
    mobile,             # AtomGroup or Universe
    reference,          # reference Universe
    select='backbone',
    weights='mass'      # or 'None' for unweighted
)
print(f"RMSD before: {rmsd_before:.3f} Å, after: {rmsd_after:.3f} Å")
```

### Align full trajectory

```python
from MDAnalysis.analysis.align import AlignTraj

aligner = AlignTraj(
    mobile_universe,
    reference_universe,
    select='protein and name CA',
    filename='aligned.xtc',    # write to file
    in_memory=False            # stream (large traj)
)
aligner.run()

# Access RMSD per frame
rmsds = aligner.results.rmsd   # shape (n_frames,)
```

### Compute average structure

```python
from MDAnalysis.analysis.align import AverageStructure

avg = AverageStructure(
    u, u,
    select='protein and name CA',
    ref_frame=0
).run()

avg_u = avg.results.universe
avg_u.atoms.write('average.pdb')
```

---

## Radius of Gyration

```python
# Per-frame — manual loop
rgyrs = []
for ts in u.trajectory:
    rgyrs.append(u.select_atoms('protein').radius_of_gyration())

# Or using custom AnalysisBase pattern
times  = u.trajectory.times      # (n_frames,) array of timestamps in ps
```

---

## RMSD Matrix (pairwise)

Compare every pair of frames — useful for conformational clustering:

```python
from MDAnalysis.analysis.rms import RMSD
import numpy as np

u = mda.Universe('system.prmtop', 'traj.nc')
ca = u.select_atoms('protein and name CA')
n  = len(u.trajectory)

matrix = np.zeros((n, n))
ref_pos = {}

for i, ts in enumerate(u.trajectory):
    ref_pos[i] = ca.positions.copy()

from MDAnalysis.analysis.align import rotation_matrix
for i in range(n):
    for j in range(i, n):
        # Use rms.rmsd (single-call, no fitting)
        r = rms.rmsd(ref_pos[i], ref_pos[j], superposition=True)
        matrix[i, j] = matrix[j, i] = r

# Plot
import seaborn as sns
sns.heatmap(matrix, cmap='viridis')
plt.xlabel('Frame'); plt.ylabel('Frame')
plt.show()
```

---

## Quick Reference: RMSD vs RMSF

| | RMSD | RMSF |
|---|---|---|
| What | Deviation from reference **per frame** | Average fluctuation **per atom** |
| Output shape | `(n_frames,)` | `(n_atoms,)` |
| Requires alignment first | Yes (built-in) | Yes (do manually before) |
| Typical use | Convergence check, stability | Flexible regions, B-factor |
| High value means | Structural drift | Mobile/flexible residue |
