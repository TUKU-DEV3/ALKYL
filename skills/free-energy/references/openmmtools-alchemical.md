# OpenMMTools — Alchemical Transformations

openmmtools provides alchemical factories and thermodynamic state objects that integrate with OpenMM. Key components: `AlchemicalFactory`, `AlchemicalState`, `ThermodynamicState`, `MultiStateSampler`.

```bash
conda install -c conda-forge openmmtools
pip install openmmtools
```

---

## AlchemicalFactory — Creating Hybrid Systems

```python
from openmmtools.alchemy import (
    AlchemicalFactory,
    AlchemicalRegion,
    AlchemicalState,
)
import openmm.app as app
import openmm as mm
import openmm.unit as unit

# Load system (protein-ligand complex, solvated)
pdb = app.PDBFile('complex_solvated.pdb')
ff = app.ForceField('amber14-all.xml', 'amber14/tip3p_standard.xml')
system = ff.createSystem(
    pdb.topology,
    nonbondedMethod=app.PME,
    nonbondedCutoff=1.0 * unit.nanometer,
    constraints=app.HBonds,
)

# Define alchemical region (ligand atom indices)
ligand_atoms = [a.index for a in pdb.topology.atoms()
                if a.residue.name == 'LIG']

alchemical_region = AlchemicalRegion(
    alchemical_atoms=ligand_atoms,
    softcore_alpha=0.5,          # softcore LJ parameter
    softcore_a=1, softcore_b=1, softcore_c=6,  # Beutler exponents
    annihilate_sterics=True,     # True=annihilate (disappear vdW in vacuum too)
    annihilate_electrostatics=True,
)

# Create alchemical system
factory = AlchemicalFactory()
alchemical_system = factory.create_alchemical_system(
    system,
    alchemical_regions=alchemical_region,
)
```

---

## AlchemicalState — Setting Lambda Values

```python
# Alchemical state controls λ parameters
alchemical_state = AlchemicalState.from_system(alchemical_system)

# Individual lambdas
alchemical_state.lambda_sterics = 1.0       # LJ interactions (1=on, 0=off)
alchemical_state.lambda_electrostatics = 1.0  # charges (1=on, 0=off)

# Apply to context
compound_state = openmmtools.states.CompoundThermodynamicState(
    thermodynamic_state=ThermodynamicState(
        alchemical_system,
        temperature=300 * unit.kelvin,
        pressure=1 * unit.bar,
    ),
    composable_states=[alchemical_state],
)
```

---

## ThermodynamicState

```python
from openmmtools.states import ThermodynamicState, SamplerState

# Thermodynamic state: T, P, system
thermo_state = ThermodynamicState(
    system=alchemical_system,
    temperature=300.0 * unit.kelvin,
    pressure=1.0 * unit.bar,   # NPT (set to None for NVT)
)

# Sampler state: positions, velocities, box
sampler_state = SamplerState(
    positions=pdb.positions,
    box_vectors=thermo_state.system.getDefaultPeriodicBoxVectors(),
)
```

---

## MultiStateSampler (HREX / MBAR)

The `MultiStateSampler` runs multiple λ windows simultaneously with optional Hamiltonian Replica Exchange (HREX) — swaps configurations between adjacent λ states to enhance sampling.

```python
from openmmtools.multistate import (
    MultiStateSampler,
    MultiStateReporter,
)
from openmmtools import mcmc

# Define thermodynamic states for each lambda
lambda_schedule = {
    'lambda_sterics':        [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
    'lambda_electrostatics': [0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0, 1.0, 1.0],
}
n_states = len(lambda_schedule['lambda_sterics'])

thermodynamic_states = []
for i in range(n_states):
    alch_state = AlchemicalState.from_system(alchemical_system)
    alch_state.lambda_sterics = lambda_schedule['lambda_sterics'][i]
    alch_state.lambda_electrostatics = lambda_schedule['lambda_electrostatics'][i]

    thermo = ThermodynamicState(
        system=alchemical_system,
        temperature=300 * unit.kelvin,
        pressure=1 * unit.bar,
    )
    compound = openmmtools.states.CompoundThermodynamicState(
        thermo, composable_states=[alch_state]
    )
    thermodynamic_states.append(compound)

# MCMC move: Langevin dynamics + MC barostat
move = mcmc.LangevinDynamicsMove(
    timestep=2.0 * unit.femtoseconds,
    collision_rate=1.0 / unit.picosecond,
    n_steps=500,                # steps per MBAR energy evaluation
    reassign_velocities=False,
)

# Reporter: saves to NetCDF
reporter = MultiStateReporter(
    storage='simulation.nc',
    checkpoint_interval=50,
    analysis_particle_indices=ligand_atoms,
)

# Initialize sampler
sampler = MultiStateSampler(
    mcmc_moves=move,
    number_of_iterations=2000,  # total iterations (= n_steps × n_iter steps per replica)
    online_analysis_interval=100,
    online_analysis_minimum_iterations=200,
)

sampler.create(
    thermodynamic_states=thermodynamic_states,
    sampler_states=[SamplerState(pdb.positions)] * n_states,
    storage=reporter,
)

# Run
sampler.run()
```

---

## HREX (Hamiltonian Replica Exchange)

Enable replica exchange between adjacent λ states — dramatically improves sampling of slow conformational changes.

```python
from openmmtools.multistate import ReplicaExchangeSampler

sampler = ReplicaExchangeSampler(
    mcmc_moves=move,
    number_of_iterations=2000,
    replica_mixing_scheme='swap-neighbors',  # or 'swap-all'
    online_analysis_interval=50,
)

# Acceptance rate target: 20-40% between adjacent λ
# If < 20%: add more λ windows (worse overlap)
# If > 50%: can remove λ windows (excess overlap)
```

---

## Analyzing NetCDF Output

```python
from openmmtools.multistate import MultiStateReporter
from pymbar import MBAR
import numpy as np

# Load reporter
reporter = MultiStateReporter('simulation.nc', open_mode='r')
analyzer = reporter.read_energies()  # u_kln matrix

# Get energies: shape (K, K, N)
u_kln = reporter.read_energies()
n_states = u_kln.shape[0]

# Equilibration detection and subsampling
from openmmtools.multistate.analyzers import MultiStateSamplerAnalyzer

analyzer = MultiStateSamplerAnalyzer(reporter)
dG, ddG = analyzer.estimate_free_energy()

print(f"ΔG = {dG:.3f} ± {ddG:.3f} kBT")
kBT = 0.5961  # kcal/mol at 298K
print(f"ΔG = {dG * kBT:.2f} ± {ddG * kBT:.2f} kcal/mol")
```

---

## Solvation Free Energy (Simple Example)

```python
# Ligand in waterbox — decouple to compute ΔG_hydration
from openmmtools.alchemy import AlchemicalFactory, AlchemicalRegion

# Build waterbox with ligand
modeller = app.Modeller(lig_topology, lig_positions)
modeller.addSolvent(ff, model='tip3p', padding=1.2*unit.nanometer)

system = ff.createSystem(modeller.topology, ...)

# All ligand atoms are alchemical
lig_atoms = [a.index for a in modeller.topology.atoms()
             if a.residue.chain.id == 'A']

region = AlchemicalRegion(
    alchemical_atoms=lig_atoms,
    annihilate_sterics=False,   # False=decouple (keep intra-molecular)
    annihilate_electrostatics=False,
)
alch_system = AlchemicalFactory().create_alchemical_system(system, region)

# Then run MultiStateSampler as above with decoupling λ schedule
```

---

## Lambda Schedule Recommendations

```python
# Conservative (24 windows, ~4% overlap minimum)
lambda_vdw  = [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35,
               0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75,
               0.8, 0.85, 0.9, 0.95, 1.0]
lambda_elec = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.2,
               0.4, 0.6, 0.8, 1.0, 1.0, 1.0, 1.0, 1.0,
               1.0, 1.0, 1.0, 1.0, 1.0]

# Minimal (12 windows, faster but may have poor overlap at endpoints)
lambda_vdw  = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
lambda_elec = [0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0, 1.0, 1.0]
```

---

## Key Notes

- `annihilate=True` — atom fully vanishes (intra-molecular LJ set to 0); use for **binding** free energies
- `annihilate=False` (decouple) — only inter-molecular interactions turned off; use for **solvation**
- Charges must reach 0 before LJ (or vice versa) to avoid singularities in alchemical endpoint
- `n_steps=500` per iteration × 2 fs = 1 ps per λ evaluation; 2000 iterations = 2 ns per window
