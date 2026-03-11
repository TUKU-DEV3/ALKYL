# OpenBabel — Formats, Fingerprints & Descriptors

## Key Format Codes

### Cheminformatics

| Code | Format | R/W | Notes |
|------|--------|-----|-------|
| `smi` | SMILES | R/W | Daylight SMILES |
| `can` | Canonical SMILES | R/W | Unique, round-trippable |
| `smy` | SMILES (Smiley parser) | R | Strict parser |
| `sdf` / `sd` / `mol` | MDL SDF/MOL | R/W | Standard structure + data |
| `mol2` / `ml2` | Tripos MOL2 | R/W | Includes atom types |
| `inchi` | InChI | R/W | IUPAC identifier |
| `inchikey` | InChIKey | W | Hashed InChI |
| `pdb` / `ent` | Protein Data Bank | R/W | For proteins |
| `xyz` | XYZ Cartesian | R/W | Simple 3D |
| `cif` | Crystallographic CIF | R/W | Crystal structures |
| `cml` | Chemical Markup Language | R/W | XML-based |

### Quantum Chemistry / MD

| Code | Format | Notes |
|------|--------|-------|
| `g09` / `g03` | Gaussian output | Read geometry |
| `orca` | ORCA output | Read geometry |
| `vasp` | VASP POSCAR | Crystal for DFT |
| `gjf` | Gaussian input | Write input file |
| `tmol` | Turbomole | Coord files |
| `nw` | NWChem input | W |
| `mop` | MOPAC input/output | R/W |
| `lmpdat` | LAMMPS data | R/W |
| `gro` | GROMACS GRO | R/W |
| `dcd` | CHARMM DCD traj | R |
| `nc` | AMBER NetCDF | R |

### Spectra & Properties

| Code | Format | Notes |
|------|--------|-------|
| `jdx` | JCAMP-DX | IR/NMR spectra |
| `cdx` / `cdxml` | ChemDraw | R/W |
| `rxn` | MDL Reaction | R/W |
| `rsmi` | Reaction SMILES | R/W |

### Enumerate All at Runtime

```python
from openbabel import pybel

# Input formats
for code, desc in sorted(pybel.informats.items()):
    print(f"{code:12s} {desc}")

# Output formats
for code, desc in sorted(pybel.outformats.items()):
    print(f"{code:12s} {desc}")
```

---

## Fingerprint Types

| Type | Bits | Description | Best for |
|------|------|-------------|----------|
| `FP2` | 1021 | Path-based linear fragments (2–7 atoms) | General similarity |
| `FP3` | 55 | SMARTS query patterns | Substructure similarity |
| `FP4` | 307 | More SMARTS patterns | Functional group similarity |
| `MACCS` | 166 | MDL MACCS structural keys | Drug-like feature comparison |

```python
from openbabel import pybel

mol = pybel.readstring('smi', 'CCO')
fp2  = mol.calcfp('FP2')
fp3  = mol.calcfp('FP3')
maccs = mol.calcfp('MACCS')

# Bit positions
print(fp2.bits)    # list of set bit indices
print(len(fp2.bits))

# Tanimoto (same type only!)
mol2 = pybel.readstring('smi', 'CCCO')
fp2b = mol2.calcfp('FP2')
print(fp2 | fp2b)   # Tanimoto coefficient
```

### FP2 vs RDKit Morgan/ECFP

| | OpenBabel FP2 | RDKit ECFP4 |
|---|---|---|
| Method | Path-based | Circular (Morgan) |
| Bits | 1021 | 2048 (typical) |
| Speed | Fast | Fast |
| Drug-like | Good | Excellent |
| Unusual atoms | Good | Limited |

---

## Molecular Descriptors

```python
from openbabel import pybel

mol = pybel.readstring('smi', 'CC(=O)Oc1ccccc1C(=O)O')   # aspirin

# All descriptors
desc = mol.calcdesc()

# Key descriptors
desc['MW']      # molecular weight
desc['logP']    # Wildman-Crippen LogP
desc['TPSA']    # topological polar surface area (Å²)
desc['MR']      # molar refractivity
desc['HBD']     # H-bond donors (O-H, N-H count)
desc['HBA1']    # H-bond acceptors (N + O count)
desc['HBA2']    # H-bond acceptors (Egel def.)
desc['abonds']  # aromatic bonds
desc['atoms']   # atom count
desc['bonds']   # bond count
desc['dbonds']  # double bonds
desc['sbonds']  # single bonds
desc['tbonds']  # triple bonds
desc['cansmi']  # canonical SMILES

# Specific subset only
desc = mol.calcdesc(['MW', 'logP', 'TPSA', 'HBD', 'HBA1'])
```

### Lipinski Ro5 Filter Example

```python
def lipinski_ro5(mol):
    d = mol.calcdesc(['MW', 'logP', 'HBD', 'HBA1'])
    return (d['MW'] <= 500 and
            d['logP'] <= 5 and
            d['HBD'] <= 5 and
            d['HBA1'] <= 10)

for mol in pybel.readfile('sdf', 'library.sdf'):
    if lipinski_ro5(mol):
        print(f"PASS: {mol.title} MW={mol.calcdesc(['MW'])['MW']:.1f}")
```

---

## Format Interoperability with RDKit

```python
# OpenBabel → RDKit (via SDF string)
from openbabel import pybel
from rdkit import Chem

def pybel_to_rdkit(pybel_mol):
    sdf_str = pybel_mol.write('sdf')
    return Chem.MolFromMolBlock(sdf_str, removeHs=False)

# RDKit → OpenBabel (via SDF string)
from rdkit.Chem import AllChem

def rdkit_to_pybel(rdkit_mol):
    sdf_str = Chem.MolToMolBlock(rdkit_mol)
    return pybel.readstring('sdf', sdf_str)

# Workflow: unusual format → OpenBabel → RDKit
def load_any_to_rdkit(filename, fmt=None):
    """Load CIF/XYZ/etc via OpenBabel, return RDKit mol"""
    if fmt is None:
        fmt = filename.rsplit('.', 1)[-1].lower()
    mols = list(pybel.readfile(fmt, filename))
    if not mols:
        raise ValueError(f"No molecules in {filename}")
    return pybel_to_rdkit(mols[0])
```

---

## Similarity Screening (Library vs Query)

```python
from openbabel import pybel
import heapq

def screen_library(query_smiles, library_sdf, top_n=10, fp_type='FP2'):
    query = pybel.readstring('smi', query_smiles)
    q_fp  = query.calcfp(fp_type)

    results = []
    for mol in pybel.readfile('sdf', library_sdf):
        mol_fp = mol.calcfp(fp_type)
        tc = q_fp | mol_fp
        # Keep top_n by Tanimoto
        if len(results) < top_n:
            heapq.heappush(results, (tc, mol.title))
        elif tc > results[0][0]:
            heapq.heapreplace(results, (tc, mol.title))

    return sorted(results, reverse=True)

hits = screen_library('CC(=O)Nc1ccc(O)cc1', 'chembl.sdf', top_n=20)
for tc, title in hits:
    print(f"{title:20s}  Tc={tc:.3f}")
```
