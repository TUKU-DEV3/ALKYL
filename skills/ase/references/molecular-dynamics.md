# ASE — Molecular Dynamics

## Ensembles & Imports

```python
from ase.md.verlet import VelocityVerlet       # NVE
from ase.md.langevin import Langevin            # NVT (recommended)
from ase.md.nvtberendsen import NVTBerendsen    # NVT (simple)
from ase.md.nptberendsen import NPTBerendsen    # NPT (Berendsen)
from ase.md.npt import NPT                     # NPT (Parrinello-Rahman)
from ase.md.velocitydistribution import (
    MaxwellBoltzmannDistribution,
    Stationary, ZeroRotation
)
from ase import units
```

---

## Velocity Initialization

Always initialize velocities before MD:

```python
# Maxwell-Boltzmann distribution at 300 K
MaxwellBoltzmannDistribution(atoms, temperature_K=300, rng=None)

# Remove center-of-mass drift
Stationary(atoms)

# Remove angular momentum (for molecules)
ZeroRotation(atoms)
```

---

## NVE — Velocity Verlet (microcanonical)

Conserves: E, V, N. No thermostat — use for equilibrated systems.

```python
dyn = VelocityVerlet(
    atoms,
    timestep=1.0 * units.fs,   # 1 fs recommended for most cases
    trajectory='nve.traj',
    logfile='nve.log',
    loginterval=10              # write every 10 steps
)

dyn.run(1000)   # 1000 steps = 1 ps at 1 fs timestep
```

---

## NVT — Langevin (canonical, recommended)

Correct NVT ensemble. Adds friction + random force to each atom.

```python
dyn = Langevin(
    atoms,
    timestep=2.0 * units.fs,
    temperature_K=300,
    friction=0.01 / units.fs,   # friction coefficient (adjust: 0.001–0.1)
    fixcm=True,                  # fix center of mass
    rng=None,                    # custom RNG for reproducibility
    trajectory='nvt.traj',
    logfile='nvt.log'
)

dyn.run(5000)   # 10 ps at 2 fs
```

**Choosing friction**: 0.01/fs is a good starting point. Higher → faster equilibration, less physical dynamics. Lower → slower equilibration, better dynamics.

---

## NVT — Berendsen (simple velocity rescaling)

Does NOT reproduce true NVT ensemble. Good for rapid equilibration only.

```python
dyn = NVTBerendsen(
    atoms,
    timestep=1.0 * units.fs,
    temperature_K=300,
    taut=500 * units.fs,        # coupling time constant (500 fs typical)
    fixcm=True,
    trajectory='eq.traj'
)
dyn.run(2000)
```

---

## NPT — Berendsen (pressure + temperature coupling)

```python
dyn = NPTBerendsen(
    atoms,
    timestep=1.0 * units.fs,
    temperature_K=300,
    taut=100 * units.fs,
    pressure_au=1.01325 * units.bar,   # 1 atm
    taup=1000 * units.fs,              # pressure coupling time
    compressibility_au=4.57e-5 / units.bar,  # water: 4.57e-5 bar⁻¹
    trajectory='npt.traj',
    logfile='npt.log'
)
dyn.run(5000)
```

---

## Observers — Callbacks During MD

```python
# Simple observer
def print_temp():
    T = atoms.get_temperature()
    e = atoms.get_potential_energy()
    print(f"T={T:.1f} K  E={e:.4f} eV")

dyn.attach(print_temp, interval=100)   # every 100 steps
dyn.run(10000)

# Attach trajectory writer explicitly
from ase.io.trajectory import Trajectory
traj = Trajectory('md.traj', 'w', atoms)
dyn.attach(traj.write, interval=10)
dyn.run(10000)
traj.close()
```

---

## Monitoring Temperature & Energy

```python
# During or after simulation
T  = atoms.get_temperature()         # kinetic temperature, K
Ek = atoms.get_kinetic_energy()      # kinetic energy, eV
Ep = atoms.get_potential_energy()    # potential energy, eV
Et = Ek + Ep                         # total energy

# From trajectory
from ase import io
frames = io.read('md.traj', index=':')
temps  = [f.get_temperature() for f in frames]
energies = [f.get_potential_energy() for f in frames]
```

---

## Complete Example: Drug-like Molecule NVT Simulation

```python
from ase.build import molecule
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary
from ase import units
from tblite.ase import TBLite
import numpy as np

# Setup
atoms = molecule('aspirin')
atoms.calc = TBLite(method='GFN2-xTB', verbosity=0)

# Initialize velocities
MaxwellBoltzmannDistribution(atoms, temperature_K=300)
Stationary(atoms)

# Equilibration — Berendsen thermostat (fast)
from ase.md.nvtberendsen import NVTBerendsen
eq = NVTBerendsen(atoms, 1.0 * units.fs, temperature_K=300,
                  taut=200 * units.fs)
eq.run(500)   # 500 fs equilibration

# Production — Langevin (correct ensemble)
temps = []
def log_T():
    temps.append(atoms.get_temperature())

dyn = Langevin(atoms, 1.0 * units.fs, temperature_K=300,
               friction=0.01 / units.fs,
               trajectory='prod.traj', logfile='prod.log')
dyn.attach(log_T, interval=10)
dyn.run(5000)   # 5 ps production

print(f"Mean T: {np.mean(temps):.1f} K ± {np.std(temps):.1f} K")
```

---

## Analysis of Trajectories

```python
from ase import io
import numpy as np

frames = io.read('prod.traj', index=':')
N = len(frames)

# RMSD from reference
ref_pos = frames[0].get_positions()
rmsds = []
for f in frames:
    diff = f.get_positions() - ref_pos
    rmsds.append(np.sqrt(np.mean(np.sum(diff**2, axis=1))))

# Radial distribution function
from ase.ga.utilities import get_rdf
distances, rdf = get_rdf(frames, rmax=6.0, nbins=100, elements=('O', 'H'))

# Mean square displacement
displacements = [f.get_positions() - frames[0].get_positions() for f in frames]
msd = [np.mean(np.sum(d**2, axis=1)) for d in displacements]
```

---

## Time Units Reference

```python
from ase import units

1 * units.fs   # 1 femtosecond
1 * units.ps   # = 1000 * units.fs (picosecond)

# Recommended timesteps
# H-containing molecules:  0.5–1.0 fs (H stretching ~10 fs period)
# Metals, no H:            2.0–5.0 fs
# Constrained H (SHAKE):   2.0 fs
```
