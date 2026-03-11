# ASE — NEB & Transition State Search

## Overview

NEB (Nudged Elastic Band) finds minimum energy paths (MEP) between two known structures.
Climbing image NEB (CI-NEB) converges the highest image to the saddle point.

> **Module change**: `from ase.mep import NEB` (ASE 3.22+). Old `ase.neb` is deprecated.

---

## Basic NEB Workflow

```python
from ase import io
from ase.mep import NEB
from ase.optimize import BFGS, FIRE
from ase.calculators.orca import ORCA

# 1. Load optimized endpoints
initial = io.read('reactant.traj')   # fully relaxed
final   = io.read('product.traj')    # fully relaxed

# 2. Create image chain (n_images intermediate + 2 endpoints = n+2 total)
n_images = 5
images = [initial]
images += [initial.copy() for _ in range(n_images)]
images += [final]

# 3. Linear interpolation of intermediate positions
neb = NEB(images, method='improvedtangent')
neb.interpolate()

# 4. Attach calculators to intermediate images only
for img in images[1:-1]:
    img.calc = ORCA(charge=0, mult=1,
                    orcasimpleinput='PBE0 def2-SVP D3BJ TightSCF')

# 5. Optimize
opt = BFGS(neb, trajectory='neb.traj', logfile='neb.log')
opt.run(fmax=0.05)   # eV/Å
```

---

## IDPP Interpolation (recommended over linear)

Better initial guess — minimizes strain in bond lengths.

```python
from ase.mep.neb import idpp_interpolate

neb = NEB(images)
neb.interpolate()          # linear first (required)
idpp_interpolate(          # then refine with IDPP
    images,
    traj='idpp.traj',
    log='idpp.log',
    fmax=0.1,
    steps=200
)
# Or: neb.interpolate(method='idpp')
```

---

## Climbing Image NEB (CI-NEB)

Drives the highest-energy image to the exact saddle point.
**Use FIRE, not BFGS**, for climbing image.

```python
neb = NEB(images, climb=True, method='improvedtangent')
neb.interpolate(method='idpp')

for img in images[1:-1]:
    img.calc = calc

# FIRE is recommended for CI-NEB
opt = FIRE(neb, trajectory='cineb.traj')
opt.run(fmax=0.05)
```

---

## Dynamic NEB (DyNEB) — Skip Converged Images

Skips already-converged images to save DFT calls.

```python
from ase.mep import DyNEB

neb = DyNEB(
    images,
    fmax=0.05,               # must match optimizer fmax
    dynamic_relaxation=True,
    scale_fmax=1.0,
    climb=True,
    method='improvedtangent'
)
neb.interpolate(method='idpp')

for img in images[1:-1]:
    img.calc = calc

opt = BFGS(neb, trajectory='dyneb.traj')
opt.run(fmax=0.05)
```

---

## AutoNEB — Automated Image Management

Automatically adds images near the saddle point for better resolution.

```python
from ase.mep.autoneb import AutoNEB

def attach_calculators(images):
    for img in images:
        img.calc = ORCA(charge=0, mult=1,
                        orcasimpleinput='PBE0 def2-SVP TightSCF')

neb = AutoNEB(
    attach_calculators=attach_calculators,
    prefix='autoneb',        # file prefix for images (autoneb000.traj, …)
    n_simul=4,               # parallel image relaxations
    n_max=9,                 # final total images (including endpoints)
    fmax=0.05,
    maxsteps=10000,
    k=0.1,
    climb=True,
    method='improvedtangent',
    space_energy_ratio=0.5   # 0 = equal spacing, 1 = equal energy gaps
)
```

---

## NEB Methods

| Method | Notes |
|--------|-------|
| `'aseneb'` | Classic ASE NEB |
| `'improvedtangent'` | **Default, recommended** — Henkelman improved tangent |
| `'eb'` | Elastic band (no projection) |
| `'spline'` | Spline interpolation, requires `precon` |
| `'string'` | String method, requires `precon` |

---

## Analysis — Barrier Extraction

```python
from ase import io
from ase.mep.neb import NEBTools

# Load final NEB images
images = io.read('neb.traj', index='-7:')  # last n+2 frames

neb_tools = NEBTools(images)

# Energy barrier and reaction energy
barrier, dE = neb_tools.get_barrier(fit=True, raw=False)
print(f"Forward barrier:  {barrier:.3f} eV")
print(f"Reaction energy:  {dE:.3f} eV")

# Max force
fmax = neb_tools.get_fmax()
print(f"Max force: {fmax:.4f} eV/Å")

# Plot energy profile
import matplotlib.pyplot as plt
fig = neb_tools.plot_band()
plt.savefig('neb_profile.png', dpi=150)
plt.show()
```

---

## Reading NEB Trajectory

```python
from ase import io

# All frames for all images
all_frames = io.read('neb.traj', index=':')

# Last snapshot of each image (n+2 images)
n_images = 7   # 5 intermediate + 2 endpoints
last = io.read('neb.traj', index=f'-{n_images}:')

# Individual image energies
for i, img in enumerate(last):
    print(f"Image {i}: {img.get_potential_energy():.4f} eV")
```

---

## MPI Parallelization

Distribute images across MPI ranks (for HPC):

```python
from ase.parallel import world

for i, img in enumerate(images[1:-1]):
    if i % world.size == world.rank:
        img.calc = ORCA(orcasimpleinput='PBE0 def2-SVP')

neb = NEB(images, parallel=True)
opt = BFGS(neb, trajectory='neb.traj')
opt.run(fmax=0.05)
```

```bash
mpirun -np 5 python neb_script.py   # 1 rank per intermediate image
```

---

## Key Parameters Summary

| Parameter | Recommended value | Notes |
|-----------|------------------|-------|
| `k` | 0.1 eV/Å | Spring constant |
| `fmax` | 0.03–0.05 eV/Å | Convergence threshold |
| `n_images` | 5–9 | More images = smoother profile |
| `climb` | `True` after initial convergence | Two-phase: plain NEB first, then CI |
| `method` | `'improvedtangent'` | Best general choice |
| `interpolate` | `'idpp'` | Better than linear |
