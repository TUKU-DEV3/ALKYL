# ABFE Protocol — Absolute Binding Free Energy

Compute absolute ΔG_bind for a single ligand without a reference compound. More expensive and noisier than RBFE but independent — no reference ligand needed.

---

## Double Decoupling Scheme

```
ΔG_bind = ΔG_decouple(complex) - ΔG_decouple(solvent) + ΔG_standard_state

Leg 1 (complex):  Ligand (bound) → Ligand (gas phase, restrained in site)
Leg 2 (solvent):  Ligand (free)  → Ligand (gas phase in waterbox)
ΔG_bind = −ΔG_leg1 + ΔG_leg2 + ΔG_std
```

### Standard State Correction
Converts from simulation concentration (1 ligand in box) to standard 1 M:

```python
import numpy as np

def standard_state_correction(box_volume_nm3, temperature=298.15):
    """
    ΔG_std = -kBT ln(V_box / V_std)
    V_std = 1.66054 nm³ (volume per molecule at 1 M = 1/(NA × 1 mol/L))
    """
    kBT = 8.314e-3 * temperature   # kJ/mol
    V_std = 1.66054  # nm³ (1/(6.022e23 × 1000 L/m³ × 1e-27 nm³/m³))
    dG_std = -kBT * np.log(box_volume_nm3 / V_std)
    return dG_std * 0.239006  # kJ/mol → kcal/mol

# Example: 10 nm³ box
box_vol = 10.0**3  # 10 nm side
print(f"ΔG_std = {standard_state_correction(box_vol):.2f} kcal/mol")
```

---

## Restraints

Restraints prevent the ligand from diffusing away when decoupled. Must be removed (or corrected for) analytically.

### Boresch Restraints (most common)
One distance + two angles + three dihedrals → 6 degrees of freedom constrained.

```python
def boresch_restraint_correction(r0, theta_A0, theta_B0,
                                  K_r, K_theta, K_phi,
                                  temperature=298.15):
    """
    Analytical free energy for releasing Boresch restraints.
    Eq. 32 from Boresch et al. J. Phys. Chem. B 107, 9535 (2003).

    r0 (Å), theta_A0/theta_B0 (rad), K_r (kcal/mol/Å²),
    K_theta/K_phi (kcal/mol/rad²)
    """
    import numpy as np
    kBT = 0.001987 * temperature  # kcal/mol

    # Standard state volume (1M): V0 = 1660.54 Å³
    V0 = 1660.54

    numerator = (8 * np.pi**2 * V0 *
                 (K_r * K_theta**2 * K_phi**3)**0.5 *
                 r0**2 *
                 np.sin(theta_A0) * np.sin(theta_B0))
    denominator = (2 * np.pi * kBT)**3

    dG_rest = -kBT * np.log(numerator / denominator)
    return dG_rest

# Example: restraint with moderate force constants
dG_rest = boresch_restraint_correction(
    r0=5.0,       # Å, protein-ligand anchor distance
    theta_A0=np.pi/2, theta_B0=np.pi/2,  # ~90° angles
    K_r=10.0,     # kcal/mol/Å²
    K_theta=10.0, # kcal/mol/rad²
    K_phi=10.0,   # kcal/mol/rad²
)
print(f"Boresch correction: {dG_rest:.2f} kcal/mol")
```

### ABFE with OpenMMTools + Restraints

```python
from openmmtools.alchemy import AlchemicalFactory, AlchemicalRegion
import openmm as mm
import openmm.unit as unit

def add_boresch_restraint(system, protein_atom, ligand_atom,
                           r0=0.5*unit.nanometer,
                           K_r=4184*unit.kilojoule_per_mole/unit.nanometer**2):
    """Add harmonic distance restraint between protein and ligand atoms."""
    restraint = mm.CustomBondForce('0.5*K_r*(r - r0)^2')
    restraint.addGlobalParameter('K_r', K_r)
    restraint.addGlobalParameter('r0', r0)
    restraint.addBond(protein_atom, ligand_atom)
    restraint.setForceGroup(5)  # separate group for energy decomposition
    system.addForce(restraint)
    return system
```

---

## Lambda Schedule for ABFE

ABFE requires more λ windows because the ligand must fully decouple.

```python
# Recommended 20-window schedule
# Phase 1: turn off charges (λ_elec: 1→0)
# Phase 2: turn off LJ (λ_vdw: 1→0)
# Phase 3: release restraints (λ_rest: 1→0)

lambda_schedule = {
    'lambda_electrostatics': [1.0, 0.75, 0.5, 0.25, 0.0,
                               0.0, 0.0, 0.0, 0.0, 0.0,
                               0.0, 0.0, 0.0, 0.0, 0.0,
                               0.0, 0.0, 0.0, 0.0, 0.0],
    'lambda_sterics':        [1.0, 1.0, 1.0, 1.0, 1.0,
                               0.9, 0.8, 0.7, 0.6, 0.5,
                               0.4, 0.3, 0.2, 0.1, 0.0,
                               0.0, 0.0, 0.0, 0.0, 0.0],
    'lambda_restraints':     [1.0, 1.0, 1.0, 1.0, 1.0,
                               1.0, 1.0, 1.0, 1.0, 1.0,
                               1.0, 1.0, 1.0, 1.0, 1.0,
                               0.75, 0.5, 0.25, 0.1, 0.0],
}
```

---

## Full ABFE Formula

```python
def compute_abfe(dG_complex_leg, dG_solvent_leg, dG_restraint_correction,
                 dG_standard_state):
    """
    ΔG_bind = ΔG_decouple(complex) - ΔG_decouple(solvent)
              + ΔG_restraint_release + ΔG_standard_state

    Sign convention (Gilson-Irikura):
    dG_complex_leg > 0 (costs energy to remove ligand from pocket)
    dG_solvent_leg > 0 (costs energy to remove from solvent)
    binding is favorable if dG_complex_leg > dG_solvent_leg
    """
    dG_bind = (-dG_complex_leg + dG_solvent_leg
               + dG_restraint_correction + dG_standard_state)
    return dG_bind

# Example values (kcal/mol)
dG_bind = compute_abfe(
    dG_complex_leg=20.5,       # costs to remove from pocket
    dG_solvent_leg=15.2,       # costs to remove from water
    dG_restraint_correction=3.1,  # Boresch correction
    dG_standard_state=2.8,    # standard state (1M)
)
print(f"ΔG_bind = {dG_bind:.2f} kcal/mol")
# Convert to Kd: Kd = exp(ΔG_bind / kBT) × 1M
import numpy as np
kBT = 0.5961  # kcal/mol at 298K
Kd_M = np.exp(dG_bind / kBT)  # Molar
print(f"Kd = {Kd_M*1e6:.2f} μM")
```

---

## Practical Considerations

### Starting Poses
- Use docking pose (validated by co-crystal if available)
- Pre-equilibrate complex with unrestrained MD (≥ 5 ns) before ABFE
- Cluster trajectory and pick representative poses if multiple binding modes

### Convergence
- ABFE requires more sampling than RBFE due to larger perturbation
- Typical: 10-20 ns per window, 20 windows = 200-400 ns total per leg
- Two legs (complex + solvent) × 200-400 ns = 400-800 ns total

### When to Use ABFE vs. RBFE

| Scenario | Preferred |
|----------|-----------|
| Series of similar compounds | RBFE (faster, accurate) |
| Structurally diverse compounds | ABFE (no reference needed) |
| First-in-class, no reference | ABFE |
| Validate force field or pose | ABFE |
| Prioritize from large set | RBFE with scaffold network |

---

## PMX — Alchemical FEP Toolkit (alternative)

pmx is a Python library for alchemical free energy with GROMACS.

```python
# pip install pmx  (or from github: pmx-gromacs)
from pmx import Model
from pmx.alchemy import AlchemicalModel

# Mutate ligand A to B (single topology)
model = Model('ligA.pdb')
alch = AlchemicalModel(model)
alch.mutate_residue(resname='LIG', to='LIG2', ff='amber99sb-ildn')
alch.write('hybrid.pdb', 'hybrid.itp')

# Then run GROMACS FEP with hybrid topology
# grompp -f fep.mdp -c hybrid.gro -p hybrid.top -o fep.tpr
# mdrun -v -deffnm fep
```

---

## openfe — Campaign Management

openfe provides high-level RBFE/ABFE campaign management with atom mapping and network planning.

```python
# pip install openfe
from openfe import SmallMoleculeComponent, ProteinComponent
from openfe.protocols import openmm_rbfe
from openfe import LigandAtomMapping, LigandNetwork
import openfe.setup

# Load molecules
ligA = SmallMoleculeComponent.from_sdf_file('ligA.sdf')
ligB = SmallMoleculeComponent.from_sdf_file('ligB.sdf')
protein = ProteinComponent.from_pdb_file('protein.pdb')

# Atom mapping (LOMAP-based)
mapper = openfe.setup.LomapAtomMapper(max3d=1.0, element_change=False)
mapping = next(mapper.suggest_mappings(ligA, ligB))

# RBFE protocol
protocol = openmm_rbfe.RelativeHybridTopologyProtocol(
    settings=openmm_rbfe.RelativeHybridTopologyProtocol.default_settings()
)

# Create transformation
transformation = openfe.Transformation(
    stateA=openfe.ChemicalSystem({'ligand': ligA, 'protein': protein}),
    stateB=openfe.ChemicalSystem({'ligand': ligB, 'protein': protein}),
    mapping={'ligand': mapping},
    protocol=protocol,
)

# Run
dag = transformation.create()
results = openfe.execute_dag(dag)
dG = results.get_estimate()
print(f"ΔΔG_bind = {dG}")
```
