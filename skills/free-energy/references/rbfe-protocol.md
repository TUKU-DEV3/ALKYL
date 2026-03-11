# RBFE Protocol — Relative Binding Free Energy

Relative ΔΔG_bind between two ligands via alchemical mutation. Most efficient for congeneric series (shared scaffold). Typical accuracy: 0.5-1.5 kcal/mol RMSE vs experiment.

---

## Overview

```
Thermodynamic cycle:

  Ligand A (complex)  →[ΔG_alch,prot]→  Ligand B (complex)
       ↑                                       ↑
  ΔG_bind(A)                             ΔG_bind(B)
       ↑                                       ↑
  Ligand A (solvent) →[ΔG_alch,solv]→  Ligand B (solvent)

ΔΔG_bind(A→B) = ΔG_bind(B) - ΔG_bind(A)
               = ΔG_alch,prot - ΔG_alch,solv
```

Run **two** alchemical simulations per edge:
1. Protein complex (ligand A → B inside binding pocket)
2. Solvent (ligand A → B in waterbox)

---

## Ligand Mapping (LOMAP)

Before running FEP, choose which ligand pairs (edges) to compute. A sparse network is more efficient than all-pairs.

```python
# lomap2: compute 2D MCS mapping between ligand pairs
from lomap import DBMolecules, lomap

# Load ligands
db = DBMolecules(directory='ligands/', output=False, time=20, ecr=True)
strict, loose = db.build_matrices()

# Generate optimal network (minimizes edges while maintaining connectivity)
nx_graph = db.build_graph()

# Plot the network
import matplotlib.pyplot as plt
import networkx as nx

pos = nx.spring_layout(nx_graph, seed=42)
nx.draw(nx_graph, pos, with_labels=True, node_color='lightblue',
        edge_color='gray', font_size=8)
labels = nx.get_edge_attributes(nx_graph, 'weight')
nx.draw_networkx_edge_labels(nx_graph, pos, edge_labels=labels)
plt.title('RBFE Ligand Network')
plt.savefig('rbfe_network.png', dpi=150, bbox_inches='tight')
```

---

## Perses — RBFE Pipeline (OpenMM-based)

perses automates protein-ligand RBFE with a single API call.

```python
from perses.app.setup_relative_calculation import RelativeFEPSetup

# Define calculation
setup = RelativeFEPSetup(
    ligand_input='ligands.sdf',          # SDF with all ligands
    old_ligand_index=0,                  # ligand A index in SDF
    new_ligand_index=1,                  # ligand B index in SDF
    forcefield_files=['amber14-all.xml', 'amber14/tip3p_standard.xml'],
    small_molecule_forcefield='openff-2.2.0',
    protein_pdb_filename='protein.pdb',
    phases=['complex', 'solvent'],        # run both legs
    pressure=1.0,                         # atm
    temperature=300.0,                    # K
    n_equilibration_iterations=1000,
    n_production_iterations=5000,
    n_lambda=12,                          # lambda windows per phase
)

# Run
ne_fep = setup.create_and_run()

# Results
dG_complex = ne_fep['complex']['ΔG']   # kcal/mol
dG_solvent = ne_fep['solvent']['ΔG']
ddGbind = dG_complex - dG_solvent
print(f"ΔΔG_bind = {ddGbind:.2f} kcal/mol")
```

---

## Manual RBFE with OpenMMTools

```python
from openmmtools.alchemy import AlchemicalFactory, AlchemicalRegion, AlchemicalState
from openmmtools.multistate import ReplicaExchangeSampler, MultiStateReporter
from openmmtools.states import ThermodynamicState, SamplerState, CompoundThermodynamicState
import openmmtools
import openmm.app as app
import openmm.unit as unit
import numpy as np

def run_rbfe_leg(complex_system, topology, positions, ligand_atoms, phase_name,
                 n_lambda=12, n_iter=2000, output_dir='.'):
    """Run one leg (complex or solvent) of an RBFE calculation."""

    # Lambda schedule (charges off first, on last)
    lambda_sterics =       np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5,
                                     0.6, 0.7, 0.8, 0.9, 1.0, 1.0])
    lambda_electrostatics = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                      0.25, 0.5, 0.75, 1.0, 1.0, 1.0])
    n_lambda = len(lambda_sterics)

    # Alchemical system
    region = AlchemicalRegion(alchemical_atoms=ligand_atoms,
                              annihilate_sterics=True,
                              annihilate_electrostatics=True)
    alch_system = AlchemicalFactory().create_alchemical_system(complex_system, region)

    # Build thermodynamic states
    states = []
    for i in range(n_lambda):
        alch_state = AlchemicalState.from_system(alch_system)
        alch_state.lambda_sterics = lambda_sterics[i]
        alch_state.lambda_electrostatics = lambda_electrostatics[i]
        thermo = ThermodynamicState(alch_system, 300*unit.kelvin, 1*unit.bar)
        compound = CompoundThermodynamicState(thermo, [alch_state])
        states.append(compound)

    # MCMC move
    from openmmtools import mcmc
    move = mcmc.LangevinDynamicsMove(
        timestep=2*unit.femtoseconds,
        collision_rate=1/unit.picosecond,
        n_steps=500,
    )

    # Reporter
    import os
    storage_path = os.path.join(output_dir, f'{phase_name}.nc')
    reporter = MultiStateReporter(storage_path, checkpoint_interval=50)

    # Sampler (HREX)
    sampler = ReplicaExchangeSampler(
        mcmc_moves=move,
        number_of_iterations=n_iter,
        online_analysis_interval=50,
    )
    sampler_states = [SamplerState(positions)] * n_lambda
    sampler.create(states, sampler_states, reporter)
    sampler.run()

    return storage_path

# Run both legs
complex_nc = run_rbfe_leg(complex_system, complex_topology, complex_positions,
                           ligand_atoms, 'complex', output_dir='rbfe_run/')
solvent_nc = run_rbfe_leg(solvent_system, solvent_topology, solvent_positions,
                           lig_atoms_solvent, 'solvent', output_dir='rbfe_run/')
```

---

## Results Analysis

```python
from openmmtools.multistate import MultiStateReporter
from openmmtools.multistate.analyzers import MultiStateSamplerAnalyzer

def get_dG_from_nc(nc_file, kBT=0.5961):
    reporter = MultiStateReporter(nc_file, open_mode='r')
    analyzer = MultiStateSamplerAnalyzer(reporter)
    dG, ddG = analyzer.estimate_free_energy()
    return dG * kBT, ddG * kBT

dG_complex, ddG_complex = get_dG_from_nc('rbfe_run/complex.nc')
dG_solvent, ddG_solvent = get_dG_from_nc('rbfe_run/solvent.nc')

ddGbind = dG_complex - dG_solvent
err_ddGbind = np.sqrt(ddG_complex**2 + ddG_solvent**2)  # propagate errors

print(f"ΔG_complex = {dG_complex:.2f} ± {ddG_complex:.2f} kcal/mol")
print(f"ΔG_solvent = {dG_solvent:.2f} ± {ddG_solvent:.2f} kcal/mol")
print(f"ΔΔG_bind   = {ddGbind:.2f} ± {err_ddGbind:.2f} kcal/mol")
```

---

## Handling Charge Changes

Changing net charge between ligands A and B requires special treatment (alchemical ion):

```python
# Add a neutralizing ion that simultaneously disappears as the charge changes
# Standard approach: co-mutate a Na+ (or Cl-) with the ligand
# Reference: Rocklin et al. 2013, DOI: 10.1063/1.4792208

# In perses: handled automatically when charge_changes=True
setup = RelativeFEPSetup(..., alchemical_charge_change=True)
```

---

## Network Analysis and MLE

For a network of N edges, use Maximum Likelihood Estimation to get consistent ΔG values:

```python
import networkx as nx
import numpy as np

# Build network from computed edges
G = nx.DiGraph()
edges = [
    ('LIG001', 'LIG002', -0.8, 0.2),   # (A, B, ΔΔG, uncertainty)
    ('LIG002', 'LIG003', +1.2, 0.3),
    ('LIG001', 'LIG003', +0.3, 0.25),
]

for a, b, ddg, err in edges:
    G.add_edge(a, b, ddg=ddg, err=err)

# Check cycle closure (should be ~0 for consistent network)
cycles = list(nx.simple_cycles(G))
for cycle in cycles:
    cycle_sum = sum(G[cycle[i]][cycle[(i+1)%len(cycle)]]['ddg']
                    for i in range(len(cycle)))
    print(f"Cycle {' → '.join(cycle)}: closure = {cycle_sum:.3f} kcal/mol")

# For MLE: use cinnabar or freenrgworkflows
# pip install cinnabar
from cinnabar import FEMap
fe_map = FEMap()
for a, b, ddg, err in edges:
    fe_map.add_relative_calculation(a, b, ddg, err)
fe_map.generate_absolute_values()
```

---

## Best Practices

```
Preparation:
  ✓ Dock all ligands into same protein conformation (consistent starting pose)
  ✓ Use LOMAP to plan edge network (≥ 2 paths between each ligand pair)
  ✓ Verify atom mapping covers all heavy atoms correctly
  ✓ Pre-equilibrate each complex (≥ 2 ns NVT/NPT before FEP)

Running:
  ✓ 12-24 λ windows (add more for large perturbations / charge changes)
  ✓ ≥ 5 ns production per window for drug-like molecules
  ✓ Use HREX for protein-ligand complex leg
  ✓ Run complex + solvent legs in parallel

Validation:
  ✓ Check overlap matrix O_ij > 0.03 for all adjacent pairs
  ✓ Check convergence (ΔG plateau in last 50% of data)
  ✓ Verify cycle closure errors < 1 kcal/mol
  ✓ Compare to experimental IC50/Kd where available (target: RMSE ≤ 1.5 kcal/mol)
```
