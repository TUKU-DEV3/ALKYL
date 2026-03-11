# MDAnalysis — Structure & Dynamics Analysis

## Dihedral Angles

### Ramachandran Plot (φ/ψ)

```python
from MDAnalysis.analysis.dihedrals import Ramachandran

u = mda.Universe('system.prmtop', 'traj.nc')
protein = u.select_atoms('protein')

rama = Ramachandran(protein).run()

# Results: angles[frame, residue, phi/psi]
angles = rama.results.angles   # shape (n_frames, n_res, 2)

# Plot
rama.plot(ax=None, ref=True)   # ref=True draws Ramachandran regions
import matplotlib.pyplot as plt
plt.savefig('ramachandran.png', dpi=150)
```

### Janin Plot (χ1/χ2 side-chain)

```python
from MDAnalysis.analysis.dihedrals import Janin

janin = Janin(protein).run()
janin.plot()
```

### Custom Dihedral (any 4 atoms)

```python
from MDAnalysis.analysis.dihedrals import Dihedral
import numpy as np

# Select 4 atoms defining the dihedral
ags = [u.select_atoms(f'resid {i} and name {n}')
       for i, n in [(10,'N'), (10,'CA'), (10,'CB'), (10,'CG')]]

dihed = Dihedral(ags).run()
chi1 = dihed.results.angles[:, 0]   # degrees, per frame
```

---

## Secondary Structure (DSSP)

Assigns secondary structure codes per residue per frame.

```python
from MDAnalysis.analysis.dssp import DSSP

u = mda.Universe('protein.pdb', 'traj.xtc')
protein = u.select_atoms('protein')

dssp = DSSP(protein).run()

# Results: (n_frames, n_residues) array of single-char codes
# Codes: H=helix, E=strand, C=coil/other
ss = dssp.results.dssp

# Fraction helix over time
helix_fraction = (ss == 'H').mean(axis=1)

# Per-residue persistence of secondary structure
helix_per_res = (ss == 'H').mean(axis=0)

import matplotlib.pyplot as plt
plt.imshow(ss.T == 'H', aspect='auto', cmap='Blues',
           origin='lower', interpolation='none')
plt.xlabel('Frame'); plt.ylabel('Residue')
plt.colorbar(label='Helix')
plt.show()
```

---

## Principal Component Analysis (PCA)

Identifies dominant modes of motion in the trajectory.

```python
from MDAnalysis.analysis.pca import PCA
from MDAnalysis.analysis import align

u = mda.Universe('system.prmtop', 'traj.nc')

# Align first
aligner = align.AlignTraj(u, u, select='backbone', in_memory=True).run()

# PCA on backbone atoms
ca = u.select_atoms('backbone')
pc = PCA(u, select='backbone', align=False, mean=None, n_components=None).run()

# Variance explained
print(pc.results.cumulated_variance[:10])   # cumulative variance, first 10 PCs

# Project trajectory onto PCs
pcs = pc.transform(ca, n_components=3)    # shape (n_frames, 3)
pc1, pc2, pc3 = pcs[:, 0], pcs[:, 1], pcs[:, 2]

# Plot free energy landscape on PC1/PC2
import matplotlib.pyplot as plt
import numpy as np

h, xe, ye = np.histogram2d(pc1, pc2, bins=50)
F = -np.log(h + 1)   # free energy (kT units, avoid log(0))
plt.contourf(xe[:-1], ye[:-1], F.T, levels=20, cmap='viridis')
plt.xlabel('PC1 (Å)'); plt.ylabel('PC2 (Å)')
plt.colorbar(label='-ln(P) [kT]')
plt.show()
```

---

## Radial Distribution Function (RDF)

```python
from MDAnalysis.analysis.rdf import InterRDF

u = mda.Universe('system.gro', 'traj.xtc')

ox_water = u.select_atoms('resname SOL and name OW')
prot_n   = u.select_atoms('protein and name N O')

rdf = InterRDF(
    ox_water,
    prot_n,
    nbins=75,
    range=(0.0, 10.0),   # Å
    norm='rdf'            # 'rdf' (default), 'density', 'none'
)
rdf.run()

plt.plot(rdf.results.bins, rdf.results.rdf)
plt.xlabel('r (Å)'); plt.ylabel('g(r)')
plt.show()
```

---

## Density Analysis

Compute 3D density of water or ligand around protein.

```python
from MDAnalysis.analysis.density import DensityAnalysis

u = mda.Universe('system.tpr', 'traj.xtc')
water_ox = u.select_atoms('resname SOL and name OW', updating=True)

da = DensityAnalysis(
    water_ox,
    delta=1.0,           # grid resolution, Å
    padding=1.5
)
da.run()

# Export density as DX file (open in VMD/PyMOL)
da.results.density.export('water_density.dx')

# Access grid data
density_grid = da.results.density.grid   # 3D array
```

---

## Mean Square Displacement (MSD) & Diffusion

```python
from MDAnalysis.analysis.msd import EinsteinMSD

u = mda.Universe('system.gro', 'traj.xtc')
water = u.select_atoms('resname SOL and name OW')

msd = EinsteinMSD(
    water,
    select='all',      # within the water group
    msd_type='xyz',    # 'xyz', 'x', 'y', 'z', 'xy', ...
    fft=True           # FFT-based (faster for long traj)
)
msd.run()

msd_vals = msd.results.timeseries   # Å²
times_ns  = msd.results.times * 1e-3  # ps → ns

# Diffusion coefficient (linear regime)
import numpy as np
slope, intercept = np.polyfit(times_ns[10:50], msd_vals[10:50], 1)
D = slope / 6   # 3D: MSD = 6Dt → D in Å²/ns
print(f"D = {D:.4f} Å²/ns = {D * 1e-5:.2e} cm²/s")

plt.plot(times_ns, msd_vals)
plt.xlabel('Time (ns)'); plt.ylabel('MSD (Å²)')
plt.show()
```

---

## Linear Density Profile

1D density along an axis — useful for membrane simulations.

```python
from MDAnalysis.analysis.lineardensity import LinearDensity

u = mda.Universe('membrane.tpr', 'traj.xtc')
lipids = u.select_atoms('resname DPPC')

ld = LinearDensity(
    lipids,
    grouping='residues',   # 'atoms', 'residues', 'segments', 'fragments'
    binsize=0.25           # Å
)
ld.run()

# Profiles along x, y, z
profile = ld.results.z.pos   # density along z-axis

plt.plot(profile)
plt.xlabel('z bin'); plt.ylabel('Density (g/cm³)')
plt.show()
```

---

## Ensemble Similarity (ENCORE)

Compare conformational ensembles (e.g., apo vs holo).

```python
from MDAnalysis.analysis.encore import similarity as es

u1 = mda.Universe('apo.prmtop', 'apo.nc')
u2 = mda.Universe('holo.prmtop', 'holo.nc')

# Harmonic ensemble similarity
HES, _ = es.hes([u1, u2], select='backbone')
print(f"HES = {HES[0,1]:.4f}")  # 0 = identical, high = different

# Clustering ensemble similarity
CES, _ = es.ces([u1, u2], select='backbone')
print(f"CES = {CES[0,1]:.4f}")
```
