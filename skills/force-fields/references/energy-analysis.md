# Energy Analysis — Decomposition, PME, Cutoffs, Water Models

---

## Energy Decomposition by Force Group

OpenMM assigns forces to numbered groups (0–31). Assign groups explicitly to decompose energy contributions.

```python
import openmm as mm
import openmm.app as app
import openmm.unit as unit

# Assign force groups
for i, force in enumerate(system.getForces()):
    name = type(force).__name__
    if 'Bond' in name:
        force.setForceGroup(0)
    elif 'Angle' in name:
        force.setForceGroup(1)
    elif 'Torsion' in name:
        force.setForceGroup(2)
    elif 'NonbondedForce' in name:
        force.setForceGroup(3)
    elif 'CMMotionRemover' in name:
        force.setForceGroup(4)
    else:
        force.setForceGroup(5)

def get_energy_components(context):
    """Return energy contributions by force group (kcal/mol)."""
    labels = ['Bonds', 'Angles', 'Torsions', 'Nonbonded', 'CMMotion', 'Other']
    components = {}
    for i, label in enumerate(labels):
        state = context.getState(getEnergy=True, groups={i})
        e = state.getPotentialEnergy().in_units_of(unit.kilocalories_per_mole)
        components[label] = e._value
    return components

components = get_energy_components(simulation.context)
for name, e in components.items():
    print(f"{name:12s}: {e:12.3f} kcal/mol")
```

### Decompose Nonbonded (LJ vs Coulomb)

```python
# Separate vdW and electrostatics using CustomNonbondedForce
# The standard NonbondedForce mixes both — decomposition requires modification

# Quick approach: zero out charges or LJ params temporarily
from copy import deepcopy
import openmm as mm

# Clone system for LJ-only
system_lj = deepcopy(system)
for force in system_lj.getForces():
    if isinstance(force, mm.NonbondedForce):
        for i in range(force.getNumParticles()):
            q, sigma, eps = force.getParticleParameters(i)
            force.setParticleParameters(i, 0.0, sigma, eps)  # zero charges
        # Also zero 1-4 electrostatics
        for j in range(force.getNumExceptions()):
            p1, p2, chargeProd, sigma, eps = force.getExceptionParameters(j)
            force.setExceptionParameters(j, p1, p2, 0.0, sigma, eps)
```

---

## Nonbonded Methods

```python
from openmm import app

# Infinite cutoff (slow, vacuum/implicit only)
nonbondedMethod=app.NoCutoff

# Hard cutoff (no PME, non-periodic)
nonbondedMethod=app.CutoffNonPeriodic

# Hard cutoff (periodic box, WRONG for charged systems)
nonbondedMethod=app.CutoffPeriodic

# Particle Mesh Ewald (PME) — REQUIRED for periodic explicit solvent
nonbondedMethod=app.PME

# Ewald (PME without mesh, slower)
nonbondedMethod=app.Ewald
```

### PME Parameters

```python
system = ff.createSystem(
    topology,
    nonbondedMethod=app.PME,
    nonbondedCutoff=1.0 * unit.nanometer,  # real-space cutoff (standard: 0.9-1.2 nm)
    ewaldErrorTolerance=0.0005,             # default 0.0005 (accuracy vs speed)
    constraints=app.HBonds,
)
```

**PME vs Cutoff:**
- PME: correct long-range electrostatics, required for polar/charged systems
- Hard cutoff: introduces discontinuities, artifacts for polar systems
- Always use PME for explicit solvent simulations

### Switching Function

```python
# Prevents discontinuity at cutoff — recommended
for force in system.getForces():
    if isinstance(force, mm.NonbondedForce):
        force.setUseSwitchingFunction(True)
        force.setSwitchingDistance(0.9 * unit.nanometer)  # switch starts here
        # Full cutoff at nonbondedCutoff (e.g. 1.0 nm)
```

---

## Water Models

| Model | Geometry | Charges | Notes |
|-------|---------|---------|-------|
| TIP3P | Rigid 3-site | O, 2H | Standard with AMBER; fast |
| TIP3P-FB | Rigid 3-site | Optimized | Better properties, with ff14SB |
| SPC/E | Rigid 3-site | Fixed | GROMOS/CHARMM |
| TIP4P | 4-site | + dummy (M-site) | Better water structure |
| TIP4P-Ew | 4-site | Re-optimized | For PME; recommended with ff19SB |
| TIP4P-2005 | 4-site | Best general | Best pure water properties |
| TIP5P | 5-site | 2 lone pairs | Best for liquid water, expensive |
| OPC | 4-site | Re-optimized | Best with ff19SB; DOI: 10.1021/jz501780a |
| OPC3 | 3-site | Re-optimized | Faster OPC variant |

```python
# Common water model XML files in OpenMM
'amber14/tip3p.xml'          # TIP3P standard
'amber14/tip3pfb.xml'        # TIP3P-FB
'amber14/spce.xml'           # SPC/E
'amber14/tip4pew.xml'        # TIP4P-Ew (use with ff19SB)
'amber14/opc.xml'            # OPC (most accurate, use with ff19SB)
'amber14/opc3.xml'           # OPC3
```

**Matching recommendations:**
- ff14SB → TIP3P (traditional) or TIP3P-FB
- ff19SB → **OPC** (validated together, DOI: 10.1021/acs.jctc.0c00194)
- CHARMM36m → TIP3P (CHARMM variant) or CHARMM-modified TIP3P

---

## Hydrogen Mass Repartitioning (HMR)

Move mass from heavy atoms to H → larger effective H mass → 4 fs timestep (2× speedup).

```python
# HMR in OpenMM
system = ff.createSystem(
    topology,
    constraints=app.HBonds,
    hydrogenMass=1.5 * unit.amu,  # default H mass ~1.008 Da
    nonbondedMethod=app.PME,
)
# Use dt=4 fs (not 2 fs)
integrator = mm.LangevinMiddleIntegrator(
    300 * unit.kelvin, 1 / unit.picosecond, 4.0 * unit.femtoseconds
)
```

---

## Periodic Box Considerations

```python
# Set box from PDB CRYST1 record (automatic if present)
# Or set manually:
simulation.context.setPeriodicBoxVectors(
    mm.Vec3(5.0, 0, 0) * unit.nanometer,
    mm.Vec3(0, 5.0, 0) * unit.nanometer,
    mm.Vec3(0, 0, 5.0) * unit.nanometer,
)

# Get box vectors
state = simulation.context.getState()
box = state.getPeriodicBoxVectors()
print(box)

# Minimum image convention: box must be > 2 * cutoff in each dimension
# Typical: 10 Å padding → box ~6 nm for a 3 nm protein
```

---

## Energy Minimization

```python
# Default L-BFGS minimizer
simulation.minimizeEnergy(
    tolerance=10 * unit.kilojoule_per_mole / unit.nanometer,  # gradient convergence
    maxIterations=1000,  # 0 = unlimited until convergence
)

# Check energy before and after
state = simulation.context.getState(getEnergy=True)
print(f"After minimization: {state.getPotentialEnergy()}")

# Typical signs of bad starting structure (high energy):
# - Clashes from ligand placement
# - Missing H atoms
# - Bad rotamers from homology model
# → If E > 1e6 kJ/mol after 1000 steps, check structure manually
```

---

## Common Energy Issues

| Symptom | Cause | Fix |
|---------|-------|-----|
| NaN energy | Atom overlap / very bad geometry | Add more minimization steps; check structure |
| E >> 10^6 kJ/mol | Clash not resolved | Use staged minimization; restrain heavy atoms |
| Energy drifts in NVE | Timestep too large | Reduce dt; check H constraints |
| Pressure oscillates | Barostat frequency too high | Increase to every 50 steps |
| Temperature off target | Friction too low | Increase γ (1→5 ps⁻¹) |
| Box collapse | Cutoff > box/2 | Use larger box or smaller cutoff |

### Staged Minimization (for clashed systems)

```python
# Phase 1: restrain heavy atoms, minimize H
from openmm import CustomExternalForce

# Add harmonic restraints to heavy atoms
restraint = mm.CustomExternalForce('k*((x-x0)^2+(y-y0)^2+(z-z0)^2)')
restraint.addGlobalParameter('k', 1000 * unit.kilojoule_per_mole / unit.nanometer**2)
restraint.addPerParticleParameter('x0')
restraint.addPerParticleParameter('y0')
restraint.addPerParticleParameter('z0')

positions = simulation.context.getState(getPositions=True).getPositions()
for i, atom in enumerate(topology.atoms()):
    if atom.element.symbol != 'H':
        restraint.addParticle(i, positions[i].value_in_unit(unit.nanometer))

system.addForce(restraint)
simulation.minimizeEnergy(maxIterations=500)

# Phase 2: reduce restraints, minimize again
simulation.context.setParameter('k', 100 * unit.kilojoule_per_mole / unit.nanometer**2)
simulation.minimizeEnergy(maxIterations=500)

# Phase 3: remove restraints
simulation.context.setParameter('k', 0)
simulation.minimizeEnergy(maxIterations=500)
```

---

## Unit Conversions

```python
import openmm.unit as unit

# Common conversions
1 * unit.kilocalories_per_mole  # → 4.184 kJ/mol
1 * unit.nanometer              # → 10 Å
1 * unit.picosecond             # → 1000 fs
1 * unit.bar                    # → 0.1 MPa

# Convert state energy
state = context.getState(getEnergy=True)
e_kj = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
e_kcal = state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
```
