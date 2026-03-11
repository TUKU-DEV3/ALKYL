# OpenFF — SMIRNOFF Force Field Toolkit

Open Force Field Initiative. **SMIRNOFF** (SMIRKS Native Open Force Field) uses SMARTS-based parameter assignment — no atom types, parameters attached directly to chemical patterns. Current release: **Sage 2.2.0** (2023).

## Install

```bash
pip install openff-toolkit openmmforcefields
# or
conda install -c conda-forge openff-toolkit openmmforcefields
```

---

## Core Concepts

| AMBER/CHARMM | SMIRNOFF |
|---|---|
| Atom types (CA, CT, NH1…) | SMARTS patterns |
| Hand-crafted type tables | Type assigned by first matching SMARTS |
| Needs antechamber to assign types | Molecule → parameters directly |
| Binary parameter files | `.offxml` (human-readable XML/JSON) |

---

## Molecule Object

```python
from openff.toolkit import Molecule

# From SMILES
mol = Molecule.from_smiles('c1ccc(cc1)CCN')

# From RDKit
from rdkit import Chem
rdmol = Chem.MolFromSmiles('CCO')
mol = Molecule.from_rdkit(rdmol)

# From file (SDF, MOL2, PDB with connectivity)
mol = Molecule.from_file('ligand.sdf')

# Generate 3D conformer (required for parameterization)
mol.generate_conformers(n_conformers=1)

# Assign AM1-BCC charges (default, recommended for drug-like molecules)
mol.assign_partial_charges('am1bcc')
# or: 'am1bccelf10' (extended conformer sampling, more accurate but slower)
# or: 'gasteiger' (fast, less accurate)
# or: 'mmff94'

print(mol.partial_charges)  # unit array
print(mol.to_smiles())      # canonical SMILES
print(mol.n_atoms, mol.n_bonds)
```

---

## Loading Force Fields

```python
from openff.toolkit import ForceField

# Sage 2.2.0 — current recommended
ff = ForceField('openff-2.2.0.offxml')

# Parsley 1.3.1 — previous generation
ff = ForceField('openff-1.3.1.offxml')

# List available
import openff.toolkit
import pkg_resources
# OFFXML files ship with openff-forcefields package:
# conda install -c conda-forge openff-forcefields
```

---

## Parameterize a Small Molecule

```python
from openff.toolkit import Molecule, ForceField
from openff.interchange import Interchange

mol = Molecule.from_smiles('CC(=O)Oc1ccccc1C(=O)O')  # aspirin
mol.generate_conformers(n_conformers=1)

ff = ForceField('openff-2.2.0.offxml')

# Create Interchange (unified topology + parameters)
interchange = Interchange.from_smirnoff(
    force_field=ff,
    topology=mol.to_topology(),
)

# Export to OpenMM
omm_system = interchange.to_openmm()
omm_topology = interchange.to_openmm_topology()
positions = interchange.positions.to_openmm()

# Export to GROMACS
interchange.to_gromacs('ligand')  # generates ligand.top, ligand.gro

# Export to AMBER
interchange.to_amber('ligand')    # ligand.prmtop, ligand.inpcrd
```

---

## Protein-Ligand System (openmmforcefields)

```python
from openmmforcefields.generators import SystemGenerator
from openff.toolkit import Molecule
import openmm.app as app
import openmm.unit as unit

# Ligand(s) to parameterize
ligand = Molecule.from_smiles('c1ccc2c(c1)CC(N2)C(=O)O')
ligand.generate_conformers(n_conformers=1)

# SystemGenerator bridges protein (AMBER XML) + ligand (SMIRNOFF)
generator = SystemGenerator(
    forcefields=[
        'amber/ff14SB.xml',
        'amber/tip3p_standard.xml',
    ],
    small_molecule_forcefield='openff-2.2.0',  # or 'gaff-2.11'
    molecules=[ligand],
    forcefield_kwargs={
        'nonbondedMethod': app.PME,
        'nonbondedCutoff': 1.0 * unit.nanometer,
        'constraints': app.HBonds,
    },
    cache='ff_cache.json',  # caches parameterization
)

# Load pre-solvated complex PDB (protein + ligand)
pdb = app.PDBFile('complex_solvated.pdb')
system = generator.create_system(pdb.topology, molecules=[ligand])
```

---

## Charge Methods Comparison

| Method | Speed | Quality | Use Case |
|--------|-------|---------|----------|
| `am1bcc` | Fast (~1s) | Good | Standard drug-like molecules |
| `am1bccelf10` | Slow (~30s) | Better | Flexible/charged molecules |
| `gasteiger` | Instant | Poor | Screening only |
| `mmff94` | Fast | Moderate | Non-polar molecules |
| RESP | Very slow | Best | High-accuracy, QM-based |
| RESP2 | Slow | Best+ | Condensed-phase corrected |

```python
# am1bcc uses OpenEye or RDKit-based toolkit
# Default toolkit priority: OpenEye > RDKit + AmberTools > RDKit
from openff.toolkit.utils.toolkits import ToolkitRegistry, RDKitToolkitWrapper, AmberToolsToolkitWrapper

registry = ToolkitRegistry([RDKitToolkitWrapper(), AmberToolsToolkitWrapper()])
mol.assign_partial_charges('am1bcc', toolkit_registry=registry)
```

---

## Topology from Multiple Molecules

```python
from openff.toolkit import Topology, Molecule

protein_mol = Molecule.from_pdb('protein.pdb')  # requires connectivity
ligand = Molecule.from_smiles('CCO')
water = Molecule.from_smiles('O')

topology = Topology.from_molecules([protein_mol, ligand] + [water] * 1000)
```

---

## Inspecting Parameters

```python
ff = ForceField('openff-2.2.0.offxml')

# What bond parameters match a SMILES?
mol = Molecule.from_smiles('c1ccccc1')
labels = ff.label_molecules(mol.to_topology())[0]

for bond, params in labels['Bonds'].items():
    print(f"Bond {bond}: k={params.k}, length={params.length}")

for angle, params in labels['Angles'].items():
    print(f"Angle {angle}: k={params.k}, angle={params.angle}")

for torsion, params in labels['ProperTorsions'].items():
    print(f"Torsion {torsion}: k={[p.k for p in params]}")
```

---

## Common Errors

| Error | Cause | Fix |
|-------|-------|-----|
| `ParameterLookupError` | Molecule has unusual chemistry not covered by Sage | Try older Sage or CGenFF/GAFF2 |
| `ChargeCalculationError` | AM1-BCC failed (unusual atoms) | Use `gasteiger` or specify charges manually |
| Partial charges don't sum to net charge | Numerical drift | Check `mol.total_charge` vs assigned |
| Inconsistent stereo | SMILES stereo not propagated to 3D | Use `mol.generate_conformers()` after stereo assignment |

---

## Version History

| Version | Name | Notes |
|---------|------|-------|
| openff-2.2.0 | Sage 2.2 | Current recommended (2023) |
| openff-2.1.0 | Sage 2.1 | Improved torsions |
| openff-2.0.0 | Sage | First Sage release |
| openff-1.3.1 | Parsley | Previous generation |
| openff-1.0.0 | Parsley | First release (2019) |
