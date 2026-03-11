# OpenBabel — pybel Python API

```python
from openbabel import pybel
from openbabel import openbabel as ob   # low-level OBMol API
```

---

## Reading Molecules

```python
# From string
mol = pybel.readstring('smi', 'CCO')               # ethanol
mol = pybel.readstring('smi', 'c1ccccc1')          # benzene
mol = pybel.readstring('inchi', 'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3')
mol = pybel.readstring('mol', open('aspirin.mol').read())

# From file (iterator — lazy, memory-efficient)
for mol in pybel.readfile('sdf', 'library.sdf'):
    print(mol.title, mol.molwt)

for mol in pybel.readfile('smi', 'smiles.smi'):
    pass

for mol in pybel.readfile('pdb', 'protein.pdb'):
    pass

# Collect all
mols = list(pybel.readfile('sdf', 'small_lib.sdf'))
```

---

## Writing Molecules

```python
# To string
smi  = mol.write('can')       # canonical SMILES (includes newline)
smi  = mol.write('smi').strip()
sdf  = mol.write('sdf')
inchi = mol.write('inchi').strip()
ikey  = mol.write('inchikey').strip()
pdb   = mol.write('pdb')

# To file (single molecule)
mol.write('sdf', 'output.sdf', overwrite=True)
mol.write('mol2', 'output.mol2', overwrite=True)
mol.write('pdb', 'output.pdb', overwrite=True)

# To file (multiple molecules)
outfile = pybel.Outputfile('sdf', 'library_out.sdf', overwrite=True)
for mol in mols:
    outfile.write(mol)
outfile.close()

# Context manager (preferred)
with pybel.Outputfile('sdf', 'out.sdf', overwrite=True) as out:
    for mol in mols:
        out.write(mol)
```

---

## Molecule Properties

```python
mol.molwt           # molecular weight (float)
mol.exactmass       # exact monoisotopic mass
mol.formula         # molecular formula string, e.g. 'C9H8O4'
mol.charge          # formal charge (int)
mol.spin            # spin multiplicity
mol.energy          # current force field energy (eV or kcal/mol depending on FF)
mol.title           # molecule name/title
mol.dim             # dimensionality: 0 (no coords), 2, or 3

mol.sssr            # smallest set of smallest rings (list of rings)

# Atom and bond access
mol.atoms           # list of pybel.Atom objects
len(mol.atoms)      # number of atoms
```

---

## Hydrogen Management

```python
mol.addh()         # add all hydrogens (explicit)
mol.removeh()      # remove hydrogens (back to implicit)

# Low-level: add polar H only
mol.OBMol.AddPolarHydrogens()

# Protonation state at given pH (via pybel)
# → Use obabel CLI with -p flag (see obabel-cli.md)
# Python equivalent:
conv = ob.OBConversion()
conv.SetInAndOutFormats('sdf', 'sdf')
mol.OBMol.AddHydrogens(False, True, 7.4)   # (polarOnly, correctForPH, pH)
```

---

## 3D Coordinate Generation

```python
# make3D: adds H, embeds, force field minimization
mol = pybel.readstring('smi', 'CC(=O)Nc1ccc(O)cc1')  # paracetamol
mol.make3D(forcefield='mmff94', steps=500)             # mmff94, uff, gaff, ghemical
print(mol.dim)  # 3

# Local optimization only (molecule must already have 3D coords)
mol.localopt(forcefield='mmff94', steps=1000)

# Check energy
print(mol.energy)
```

**Force field choice:**
- `'mmff94'` — best for drug-like organic molecules
- `'uff'` — universal, covers unusual atoms (organometallics)
- `'gaff'` — AMBER GAFF, for AMBER-compatible systems
- `'ghemical'` — older, less reliable

---

## Descriptors

```python
# Available descriptors
print(pybel.descs)
# ['HBA1', 'HBA2', 'HBD', 'logP', 'MR', 'MW', 'TPSA', 'abonds', 'atoms',
#  'bonds', 'cansmi', 'cansmiNS', 'dbonds', 'HBA2', ...]

# Calculate all descriptors
desc = mol.calcdesc()
print(desc['logP'])
print(desc['TPSA'])
print(desc['MR'])     # molar refractivity
print(desc['MW'])

# Calculate specific descriptors only
desc = mol.calcdesc(['logP', 'TPSA', 'MW'])
```

---

## Fingerprints and Similarity

```python
# Available types: ['FP2', 'FP3', 'FP4', 'MACCS']
fp1 = mol1.calcfp('FP2')    # default: FP2 (path-based, 1021 bits)
fp2 = mol2.calcfp('FP2')

# Tanimoto similarity via | operator
tanimoto = fp1 | fp2        # float in [0, 1]
print(f"Tanimoto: {tanimoto:.3f}")

# Pairwise similarity matrix
mols = list(pybel.readfile('sdf', 'library.sdf'))
fps  = [mol.calcfp('FP2') for mol in mols]

import numpy as np
n = len(fps)
sim_matrix = np.zeros((n, n))
for i in range(n):
    for j in range(i, n):
        sim_matrix[i, j] = sim_matrix[j, i] = fps[i] | fps[j]

# Similarity search: find molecules similar to query
query_fp = query_mol.calcfp('FP2')
similarities = [(i, query_fp | fp) for i, fp in enumerate(fps)]
similarities.sort(key=lambda x: -x[1])
print("Top-5 similar:", similarities[:5])
```

**Fingerprint types:**
| Type | Bits | Description |
|------|------|-------------|
| `FP2` | 1021 | Path-based (linear fragments up to 7 atoms) — most used |
| `FP3` | 55 | SMARTS-based query patterns |
| `FP4` | 307 | More SMARTS patterns, designed for substructure queries |
| `MACCS` | 166 | MDL MACCS keys — functional group based |

---

## SMARTS Substructure Matching

```python
from openbabel import pybel

# Create SMARTS pattern
smarts = pybel.Smarts('[NH2]c1ccccc1')   # aromatic amine

# Search molecule
matches = smarts.findall(mol)
# Returns list of tuples: [(atom_idx1, atom_idx2, ...), ...]
print(f"Matches: {len(matches)}")

# Filter library
amine_mols = []
smarts_amine = pybel.Smarts('[NH2,NH1]')
for mol in pybel.readfile('sdf', 'library.sdf'):
    if smarts_amine.findall(mol):
        amine_mols.append(mol)

# PAINS-like filter (example)
pains_smarts = [
    pybel.Smarts('[CH]=O'),             # aldehyde
    pybel.Smarts('c1ccc(cc1)-c1ccncc1'), # biphenyl-pyridine
]
def is_pains(mol):
    return any(sm.findall(mol) for sm in pains_smarts)
```

---

## Molecule Data (SDF Properties)

```python
# Access SDF data fields
data = mol.data
print(dict(data))           # all key-value pairs from SDF $$$$

# Read specific field
ic50 = data.get('IC50')

# Set field
data['source'] = 'chembl'
data['IC50'] = '0.42'

# Delete field
del data['unwanted_field']
```

---

## Low-level OBMol Access

```python
# Access C++ OBMol for operations not in pybel
obmol = mol.OBMol

# Atom iteration
for atom in ob.OBMolAtomIter(obmol):
    print(atom.GetAtomicNum(), atom.GetX(), atom.GetY(), atom.GetZ())

# Bond iteration
for bond in ob.OBMolBondIter(obmol):
    print(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond.GetBondOrder())

# Aromaticity perception
obmol.SetAromaticPerceived(True)
for atom in ob.OBMolAtomIter(obmol):
    print(atom.IsAromatic())

# Ring membership
for ring in obmol.GetSSSR():
    print(ring.Size(), ring.IsAromatic())
```

---

## Batch Processing Pattern

```python
from openbabel import pybel

input_file  = 'raw_library.sdf'
output_file = 'filtered.sdf'

smarts_filter = pybel.Smarts('[#7,#8]')   # must contain N or O

with pybel.Outputfile('sdf', output_file, overwrite=True) as out:
    for mol in pybel.readfile('sdf', input_file):
        desc = mol.calcdesc(['MW', 'logP', 'TPSA'])
        if (200 <= desc['MW'] <= 500 and
            -0.5 <= desc['logP'] <= 5.0 and
            desc['TPSA'] <= 140 and
            smarts_filter.findall(mol)):
            out.write(mol)
```
