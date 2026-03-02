# RDKit — Molecular I/O

## SMILES

```python
from rdkit import Chem

# Parse SMILES
mol = Chem.MolFromSmiles('CC(=O)Oc1ccccc1C(=O)O')
if mol is None:
    raise ValueError("Invalid SMILES")

# Generate canonical SMILES
smi = Chem.MolToSmiles(mol)                       # canonical
smi = Chem.MolToSmiles(mol, isomericSmiles=True)  # with stereo (default True)
smi = Chem.MolToSmiles(mol, isomericSmiles=False) # strip stereo
smi = Chem.MolToSmiles(mol, kekuleSmiles=True)    # Kekule form

# Random SMILES (data augmentation / ML)
rand_smi  = Chem.MolToSmiles(mol, doRandom=True, canonical=False)
rand_list = Chem.MolToRandomSmilesVect(mol, 20, randomSeed=42)

# SMARTS
query = Chem.MolFromSmarts('[C:1](=O)[OH]')       # carboxylic acid query
```

## SDF / Molfile

```python
# Read SDF — sequential access (large files)
with Chem.SDMolSupplier('library.sdf') as suppl:
    mols = [m for m in suppl if m is not None]

# Random access by index
suppl = Chem.SDMolSupplier('library.sdf')
first_mol = suppl[0]
n_mols    = len(suppl)

# Read compressed SDF
import gzip
with gzip.open('library.sdf.gz') as f:
    with Chem.ForwardSDMolSupplier(f) as suppl:
        mols = [m for m in suppl if m is not None]

# Read from string buffer
suppl = Chem.SDMolSupplier()
suppl.SetData(open('mol.sdf').read())
mols = [m for m in suppl if m is not None]

# Write SDF
with Chem.SDWriter('output.sdf') as w:
    for mol in mols:
        w.write(mol)

# Write with specific conformer (3D)
w = Chem.SDWriter('output_3d.sdf')
for mol in mols:
    w.write(mol, confId=0)
w.close()

# Write to string buffer (no file)
from io import StringIO
buf = StringIO()
with Chem.SDWriter(buf) as w:
    for mol in mols:
        w.write(mol)
sdf_string = buf.getvalue()
```

## MOL Block

```python
# Mol → MOL block string
mol_block = Chem.MolToMolBlock(mol)

# MOL block → Mol
mol2 = Chem.MolFromMolBlock(mol_block)

# With 2D coordinates
from rdkit.Chem import AllChem
AllChem.Compute2DCoords(mol)
mol_block_2d = Chem.MolToMolBlock(mol)

# V3000 format (large molecules, enhanced stereochemistry)
mol_block_v3k = Chem.MolToV3KMolBlock(mol)
mol_from_v3k  = Chem.MolFromMolBlock(mol_block_v3k)
```

## SMILES File Supplier

```python
# File: one SMILES per line, optional name column
# 'CC(=O)O acetic_acid'
with Chem.SmilesMolSupplier('input.smi', titleLine=False) as suppl:
    mols = [m for m in suppl if m is not None]

with Chem.SmilesWriter('output.smi') as w:
    for mol in mols:
        w.write(mol)
```

## PDB

```python
# Read PDB
mol = Chem.MolFromPDBFile('protein.pdb', removeHs=False)
mol = Chem.MolFromPDBBlock(pdb_string)

# Write PDB (molecule must have 3D coordinates)
Chem.MolToPDBFile(mol, 'output.pdb')
pdb_block = Chem.MolToPDBBlock(mol)
```

## Molecule Properties (SD tags)

```python
# Set properties
mol.SetProp('Name', 'aspirin')
mol.SetDoubleProp('MW', 180.16)
mol.SetIntProp('NumAtoms', mol.GetNumAtoms())

# Get properties
mol.GetProp('Name')       # 'aspirin'
mol.GetDoubleProp('MW')   # 180.16
mol.HasProp('MW')         # True
mol.GetPropNames()        # ['Name', 'MW', 'NumAtoms', ...]

# Iterate all properties
for key in mol.GetPropNames():
    print(key, mol.GetProp(key))

mol.ClearProp('MW')       # delete a property
```

## Binary Serialization

```python
# Binary (compact, fast — for trusted internal data only)
# Note: Only use with molecules from trusted sources
bin_data = mol.ToBinary()          # → bytes
mol2     = Chem.Mol(bin_data)      # reconstruct

# Store in dict / database column
store = {Chem.MolToSmiles(m): m.ToBinary() for m in mols}
mol_back = Chem.Mol(store[smiles_key])
```

## Sanitization and Validation

```python
# Full sanitization (done automatically by MolFromSmiles)
Chem.SanitizeMol(mol)

# Partial sanitization (skip specific flags)
Chem.SanitizeMol(mol,
    Chem.SanitizeFlags.SANITIZE_ALL ^
    Chem.SanitizeFlags.SANITIZE_PROPERTIES)

# Detect chemistry problems without raising
for prob in Chem.DetectChemistryProblems(mol):
    print(prob.GetType(), prob.Message())

# Safe parsing helper
def safe_mol(smi: str):
    mol = Chem.MolFromSmiles(smi, sanitize=False)
    if mol is None:
        return None
    try:
        Chem.SanitizeMol(mol)
        return mol
    except Exception:
        return None
```

## Hydrogens

```python
mol_h    = Chem.AddHs(mol)              # add explicit H (required for 3D)
mol_no_h = Chem.RemoveHs(mol_h)         # back to heavy-atom representation

mol.GetNumAtoms()                        # heavy atoms only
mol.GetNumAtoms(onlyExplicit=False)      # include implicit H
```

## Batch Processing Pattern

```python
from rdkit import Chem
from rdkit.Chem import Descriptors
import pandas as pd

def process_sdf(path: str) -> pd.DataFrame:
    records = []
    with Chem.SDMolSupplier(path) as suppl:
        for i, mol in enumerate(suppl):
            if mol is None:
                continue
            records.append({
                'idx':    i,
                'smiles': Chem.MolToSmiles(mol),
                'name':   mol.GetProp('_Name') if mol.HasProp('_Name') else '',
                'mw':     round(Descriptors.MolWt(mol), 3),
                'logp':   round(Descriptors.MolLogP(mol), 3),
            })
    return pd.DataFrame(records)
```
