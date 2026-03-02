# RDKit — 3D Conformer Generation

## Overview

```
Workflow: SMILES → AddHs → EmbedMolecule (ETKDG) → Optimize (MMFF/UFF) → RemoveHs
```

---

## Single Conformer — ETKDG (recommended)

```python
from rdkit import Chem
from rdkit.Chem import AllChem

mol = Chem.MolFromSmiles('CC(N)C(=O)O')   # alanine
mol = Chem.AddHs(mol)

# ETKDGv3 — best current method (RDKit >= 2020.09)
params = AllChem.ETKDGv3()
params.randomSeed = 42           # reproducibility
result = AllChem.EmbedMolecule(mol, params)

if result == -1:
    raise RuntimeError("3D embedding failed — try different seed or ETKDGv3()")

# Get coordinates
conf = mol.GetConformer()
pos  = conf.GetAtomPosition(0)   # Point3D object
print(pos.x, pos.y, pos.z)

# Remove Hs for clean output
mol_no_h = Chem.RemoveHs(mol)
```

## Multiple Conformers

```python
mol = Chem.MolFromSmiles('CCCCCCCO')
mol = Chem.AddHs(mol)

params = AllChem.ETKDGv3()
params.randomSeed       = 42
params.numThreads       = 0      # 0 = use all available CPU cores
params.useSmallRingTorsions = True   # improved for small rings
params.useMacrocycleTorsions = True  # for macrocycles (>12 atoms)

cids = AllChem.EmbedMultipleConfs(mol, numConfs=50, params=params)
print(f"Generated {len(cids)} conformers")

# Access individual conformers
for cid in cids:
    conf = mol.GetConformer(cid)
    # ...

# Number of conformers
mol.GetNumConformers()   # 50
```

---

## Force Field Optimization

### MMFF94 (recommended)

```python
# Single conformer
AllChem.MMFFOptimizeMolecule(mol)
AllChem.MMFFOptimizeMolecule(mol, mmffVariant='MMFF94s')  # 'MMFF94' or 'MMFF94s'

# Multiple conformers (returns list of (not_converged, energy) tuples)
results = AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=0)
for cid, (not_conv, energy) in enumerate(results):
    print(f"Conformer {cid}: E={energy:.3f} kcal/mol, converged={not not_conv}")

# Get MMFF force field object (for custom calculations)
mp = AllChem.MMFFGetMoleculeProperties(mol)
ff = AllChem.MMFFGetMoleculeForceField(mol, mp, confId=0)
ff.Minimize(maxIts=1000)
e = ff.CalcEnergy()
```

### UFF (Universal Force Field — when MMFF unavailable)

```python
AllChem.UFFOptimizeMolecule(mol)
results = AllChem.UFFOptimizeMoleculeConfs(mol, numThreads=0)
```

---

## ETKDGv3 Parameters Reference

```python
params = AllChem.ETKDGv3()

params.randomSeed            = 42         # int, 0xf00d also common
params.numThreads            = 0          # 0 = max threads
params.maxIterations         = 1000       # embedding optimization steps
params.useSmallRingTorsions  = True       # better small ring sampling
params.useMacrocycleTorsions = True       # enable for macrocycles
params.useBasicKnowledge     = True       # use bond length/angle knowledge
params.ETversion             = 2          # 1=ETKDG, 2=ETKDGv3
params.pruneRmsThresh        = 0.5        # prune conformers within RMSD threshold
                                          # (used with EmbedMultipleConfs)
params.embedFragmentsSeparately = False   # separate fragments
```

---

## Conformer Alignment and RMSD

```python
# Align all conformers to the first one
rmslist = []
AllChem.AlignMolConformers(mol, RMSlist=rmslist)
print(f"Max RMSD from ref: {max(rmslist):.3f} A")

# RMSD between two specific conformers
rmsd = AllChem.GetConformerRMS(mol, confId1=0, confId2=1, prealigned=False)

# Align mol2 onto mol1 (return RMSD)
rmsd = AllChem.AlignMol(mol2, mol1)  # mol2 is moved in place
```

---

## 2D → 3D (common workflow)

```python
def smi_to_3d(smiles: str, n_confs: int = 1, seed: int = 42) -> Chem.Mol:
    """SMILES → minimized 3D conformer(s)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")

    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.numThreads = 0

    if n_confs == 1:
        ok = AllChem.EmbedMolecule(mol, params)
        if ok == -1:
            raise RuntimeError(f"Embedding failed for: {smiles}")
        AllChem.MMFFOptimizeMolecule(mol)
    else:
        AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, params=params)
        AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=0)

    return mol

mol_3d = smi_to_3d('CC1CCC(CC1)C(C)(C)O')
```

---

## Writing 3D Molecules

```python
# SDF with 3D coordinates
with Chem.SDWriter('output_3d.sdf') as w:
    for cid in range(mol.GetNumConformers()):
        w.write(mol, confId=cid)

# XYZ format
def mol_to_xyz(mol, conf_id=0) -> str:
    conf  = mol.GetConformer(conf_id)
    lines = [str(mol.GetNumAtoms()), '']
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        lines.append(f"{atom.GetSymbol():2s}  {pos.x:10.4f}  {pos.y:10.4f}  {pos.z:10.4f}")
    return '\n'.join(lines)

xyz_string = mol_to_xyz(mol, conf_id=0)
```

---

## 3D Descriptors (require conformer)

```python
from rdkit.Chem import rdMolDescriptors

# Requires AddHs + EmbedMolecule + MMFFOptimize
mol_h = Chem.AddHs(Chem.MolFromSmiles('CC(N)C(=O)O'))
AllChem.EmbedMolecule(mol_h, AllChem.ETKDGv3())
AllChem.MMFFOptimizeMolecule(mol_h)

# Shape / PMI descriptors
rdMolDescriptors.CalcAsphericity(mol_h)         # 0=sphere, 1=rod
rdMolDescriptors.CalcEccentricity(mol_h)
rdMolDescriptors.CalcSpherocityIndex(mol_h)     # 0=rod, 1=sphere
rdMolDescriptors.CalcInertialShapeFactor(mol_h)
rdMolDescriptors.CalcRadiusOfGyration(mol_h)    # Angstroms
rdMolDescriptors.CalcPMI1(mol_h)                # smallest principal moment
rdMolDescriptors.CalcPMI2(mol_h)
rdMolDescriptors.CalcPMI3(mol_h)                # largest principal moment

# Descriptor vectors (for 3D-QSAR)
rdf     = rdMolDescriptors.CalcRDF(mol_h)       # 210 values
morse   = rdMolDescriptors.CalcMORSE(mol_h)     # 224 values
whim    = rdMolDescriptors.CalcWHIM(mol_h)      # 114 values
getaway = rdMolDescriptors.CalcGETAWAY(mol_h)   # 273 values
autocorr3d = rdMolDescriptors.CalcAUTOCORR3D(mol_h)  # 80 values
```

---

## Working with Conformers

```python
# Get conformer by ID
conf = mol.GetConformer(confId=0)

# Get all coordinates as numpy array
import numpy as np
coords = conf.GetPositions()   # (n_atoms, 3) ndarray

# Set coordinates manually
from rdkit.Geometry import Point3D
conf.SetAtomPosition(idx, Point3D(x, y, z))

# Remove a specific conformer
mol.RemoveConformer(confId=0)

# Remove all conformers
mol.RemoveAllConformers()

# Copy a conformer
new_mol = Chem.RWMol(mol)
new_mol.AddConformer(mol.GetConformer(0), assignId=True)
```

---

## Constrained Embedding (fix part of structure)

```python
from rdkit.Chem import AllChem

ref_mol = Chem.MolFromSmiles('c1ccccc1C(=O)O')
AllChem.Compute2DCoords(ref_mol)   # or use 3D ref

probe = Chem.MolFromSmiles('c1ccccc1C(=O)NC')
AllChem.GenerateDepictionMatching2DStructure(probe, ref_mol)

# For 3D constrained embedding
match = probe.GetSubstructMatch(ref_mol)
AllChem.EmbedMolecule(
    probe,
    coordMap={i: ref_mol.GetConformer().GetAtomPosition(j)
              for i, j in enumerate(match)},
    params=AllChem.ETKDGv3()
)
```
