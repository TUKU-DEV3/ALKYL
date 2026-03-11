# OpenMM — Python MD Engine

OpenMM 8.x. Core objects: `Topology`, `System`, `Integrator`, `Simulation`. Supports CUDA/OpenCL/CPU. Force fields loaded from XML files bundled with `openmm` or `openmmforcefields`.

---

## Minimal Simulation Pipeline

```python
import openmm.app as app
import openmm as mm
import openmm.unit as unit

# 1. Load structure
pdb = app.PDBFile('input.pdb')

# 2. Load force field
ff = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

# 3. Build system
system = ff.createSystem(
    pdb.topology,
    nonbondedMethod=app.PME,
    nonbondedCutoff=1.0 * unit.nanometer,
    constraints=app.HBonds,       # constrain H-bonds → 2 fs timestep
    hydrogenMass=1.5 * unit.amu,  # HMR → 4 fs timestep
)

# 4. Add barostat (NPT)
system.addForce(mm.MonteCarloBarostat(
    1.0 * unit.bar, 300.0 * unit.kelvin, 25
))

# 5. Integrator
integrator = mm.LangevinMiddleIntegrator(
    300.0 * unit.kelvin,   # temperature
    1.0 / unit.picosecond, # friction coefficient
    2.0 * unit.femtoseconds,  # timestep
)

# 6. Simulation
simulation = app.Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)

# 7. Energy minimization
simulation.minimizeEnergy(maxIterations=500, tolerance=10*unit.kilojoule_per_mole/unit.nanometer)

# 8. Reporters
simulation.reporters.append(app.DCDReporter('output.dcd', 1000))  # every 2 ps
simulation.reporters.append(app.StateDataReporter(
    'log.csv', 1000, step=True, potentialEnergy=True,
    temperature=True, density=True, progress=True,
    totalSteps=500000
))

# 9. NVT equilibration (no barostat active if not added)
simulation.step(50000)   # 100 ps
```

---

## Structure Loading

### From PDB
```python
pdb = app.PDBFile('protein.pdb')
# Access: pdb.topology, pdb.positions
```

### From AMBER prmtop/inpcrd
```python
prmtop = app.AmberPrmtopFile('system.prmtop')
inpcrd = app.AmberInpcrdFile('system.inpcrd')
system = prmtop.createSystem(nonbondedMethod=app.PME, ...)
```

### From GROMACS topology
```python
gro = app.GromacsGroFile('system.gro')
top = app.GromacsTopFile('topol.top', periodicBoxVectors=gro.getPeriodicBoxVectors())
system = top.createSystem(nonbondedMethod=app.PME, ...)
```

### From CHARMM PSF
```python
psf = app.CharmmPsfFile('system.psf')
params = app.CharmmParameterSet('charmm36.prm', 'lig.prm')
system = psf.createSystem(params, nonbondedMethod=app.PME, ...)
```

---

## ForceField XML Files (bundled)

```python
# AMBER
'amber14-all.xml'              # ff14SB proteins + GAFF2 small molecules
'amber14/tip3p.xml'            # TIP3P water
'amber14/tip3pfb.xml'          # TIP3P-FB (more accurate)
'amber14/spce.xml'             # SPC/E water

# CHARMM
'charmm36.xml'                 # CHARMM36 proteins
'charmm36/water.xml'           # TIP3P CHARMM-variant

# Additional (requires openmmforcefields)
'amber/ff14SB.xml'
'amber/ff19SB.xml'
'amber/lipid17.xml'
```

---

## Integrators

| Integrator | Use case |
|-----------|----------|
| `VerletIntegrator(dt)` | NVE (constant energy) |
| `LangevinIntegrator(T, γ, dt)` | NVT, older |
| `LangevinMiddleIntegrator(T, γ, dt)` | NVT, recommended (leapfrog-like, 2020+) |
| `NoseHooverIntegrator(T, γ, dt)` | NVT, canonical ensemble |
| `MTSLangevinIntegrator(T, γ, dt, groups)` | Multiple time step |
| `BrownianIntegrator(T, γ, dt)` | Overdamped Langevin (no inertia) |

```python
# Recommended NVT integrator
integrator = mm.LangevinMiddleIntegrator(
    300 * unit.kelvin,
    1.0 / unit.picosecond,
    2.0 * unit.femtoseconds
)
```

---

## Thermostats and Barostats

```python
# NPT: Monte Carlo barostat (add to system, not integrator)
barostat = mm.MonteCarloBarostat(
    1.0 * unit.bar,           # target pressure
    300.0 * unit.kelvin,      # temperature (must match integrator)
    25                        # MC attempts every N steps
)
system.addForce(barostat)

# Anisotropic barostat (membranes)
barostat = mm.MonteCarloAnisotropicBarostat(
    (1.0, 1.0, 1.0) * unit.bar, 300 * unit.kelvin
)

# Membrane barostat (constant surface tension)
barostat = mm.MonteCarloMembraneBarostat(
    1.0 * unit.bar, 0.0 * unit.bar * unit.nanometer,  # P, gamma
    mm.MonteCarloMembraneBarostat.XYIsotropic,
    mm.MonteCarloMembraneBarostat.ZFree,
    300 * unit.kelvin
)
```

---

## Reporters

```python
# DCD trajectory (compact binary, MDAnalysis/CPPTRAJ compatible)
simulation.reporters.append(app.DCDReporter('traj.dcd', 1000))

# XTC trajectory (GROMACS format, via mdtraj or openmm-mdanalysis)
# Install: pip install mdtraj
from mdtraj.reporters import XTCReporter
simulation.reporters.append(XTCReporter('traj.xtc', 1000))

# State data CSV
simulation.reporters.append(app.StateDataReporter(
    'data.csv', 500,
    step=True, time=True, temperature=True,
    potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
    volume=True, density=True,
    progress=True, remainingTime=True, speed=True,
    totalSteps=500000
))

# Checkpoint (restart)
simulation.reporters.append(app.CheckpointReporter('checkpoint.chk', 10000))
```

---

## Restart from Checkpoint

```python
# Save state (full precision, portable)
simulation.saveState('state.xml')

# Load state
simulation2 = app.Simulation(topology, system, integrator)
simulation2.loadState('state.xml')

# Or from checkpoint (faster but not portable across platforms)
simulation2.loadCheckpoint('checkpoint.chk')

# Set velocities to temperature if no checkpoint
simulation.context.setVelocitiesToTemperature(300 * unit.kelvin)
```

---

## Platform Selection

```python
# Auto-detect best (CUDA > OpenCL > CPU)
platform = mm.Platform.getPlatformByName('CUDA')
properties = {'DeviceIndex': '0', 'Precision': 'mixed'}

simulation = app.Simulation(topology, system, integrator, platform, properties)

# CPU (always available, for testing)
platform = mm.Platform.getPlatformByName('CPU')
properties = {'Threads': '4'}
```

---

## System Modeller

```python
from openmm.app import Modeller

# Combine topologies (protein + ligand)
modeller = Modeller(protein_top, protein_pos)
modeller.add(ligand_top, ligand_pos)

# Add solvent (cubic box, 10 Å padding)
modeller.addSolvent(
    ff,
    model='tip3p',
    padding=1.0 * unit.nanometer,
    ionicStrength=0.15 * unit.molar,  # 150 mM NaCl
    positiveIon='Na+',
    negativeIon='Cl-'
)

topology = modeller.topology
positions = modeller.positions
```

---

## State Queries

```python
state = simulation.context.getState(
    getEnergy=True,
    getPositions=True,
    getVelocities=True,
    getForces=True,
    enforcePeriodicBox=True
)

pe = state.getPotentialEnergy().in_units_of(unit.kilocalories_per_mole)
ke = state.getKineticEnergy().in_units_of(unit.kilocalories_per_mole)
pos = state.getPositions(asNumpy=True)  # shape (N, 3), nanometers
forces = state.getForces(asNumpy=True)  # shape (N, 3), kJ/mol/nm
box = state.getPeriodicBoxVectors()
```

---

## Key Pitfalls

- `LangevinMiddleIntegrator` > `LangevinIntegrator` (better sampling, same cost)
- `HydrogenMass=1.5` amu → 4 fs timestep but ONLY with `constraints=app.HBonds`
- Barostat temperature must **exactly match** integrator temperature
- `enforcePeriodicBox=True` in getState → wraps coordinates, needed for analysis
- Remove barostat forces before NVT production: `system.removeForce(idx)`
