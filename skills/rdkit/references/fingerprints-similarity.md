# RDKit — Fingerprints & Similarity

## Fingerprint Types at a Glance

| Type | Generator | Bits | Captures | Best for |
|------|-----------|------|----------|----------|
| Morgan (ECFP) | `GetMorganGenerator(radius=2)` | 2048 | Circular atom environments | General similarity, VS |
| Morgan count | `GetSparseCountFingerprint` | variable | Count of environments | ML features |
| RDKit topological | `GetRDKitFPGenerator()` | 2048 | Bond paths | Substructure similarity |
| Atom pair | `GetAtomPairGenerator()` | 2^32 | Atom-pair distances | Scaffold hopping |
| Topological torsion | `GetTopologicalTorsionGenerator()` | 2^32 | 4-atom torsion chains | 3D proxy |
| MACCS | `MACCSkeys.GenMACCSKeys()` | 167 | Structural keys | Quick drug-likeness |
| FCFP | Morgan + feature invars | 2048 | Feature environments | Pharmacophore sim. |

---

## Generator API (recommended, RDKit ≥ 2020)

```python
from rdkit import Chem
from rdkit.Chem import AllChem, MACCSkeys
from rdkit import DataStructs

mols = [Chem.MolFromSmiles(s) for s in [
    'CC(=O)Oc1ccccc1C(=O)O',   # aspirin
    'CC(=O)Nc1ccc(O)cc1',       # paracetamol
    'OC(=O)c1ccccc1O',          # salicylic acid
]]
```

### Morgan Fingerprints (ECFP-like)

```python
# ECFP4 equivalent (radius=2, 2048 bits)
fpgen = AllChem.GetMorganGenerator(radius=2)
fps = [fpgen.GetFingerprint(m) for m in mols]

# ECFP6 equivalent
fpgen6 = AllChem.GetMorganGenerator(radius=3)

# Custom bit size
fpgen = AllChem.GetMorganGenerator(radius=2, fpSize=4096)

# Sparse count fingerprint (dense features for ML)
fpgen = AllChem.GetMorganGenerator(radius=2)
count_fp = fpgen.GetSparseCountFingerprint(mols[0])
count_fp.GetNonzeroElements()   # {bit_id: count, ...}

# Bit fingerprint (2048 bits)
bit_fp = fpgen.GetFingerprint(mols[0])
len(bit_fp)    # 2048
bit_fp.GetNumOnBits()
list(bit_fp.GetOnBits())[:5]
```

### FCFP (Feature-Based Morgan — pharmacophore-aware)

```python
# Feature-based atom invariants (Donor, Acceptor, Aromatic, etc.)
invgen = AllChem.GetMorganFeatureAtomInvGen()
fcfp = AllChem.GetMorganGenerator(radius=2, atomInvariantsGenerator=invgen)

fps_fcfp = [fcfp.GetFingerprint(m) for m in mols]
```

### RDKit Topological Fingerprints

```python
fpgen = AllChem.GetRDKitFPGenerator()
fps_rdk = [fpgen.GetFingerprint(m) for m in mols]

# Custom path length
fpgen = AllChem.GetRDKitFPGenerator(minPath=1, maxPath=7, fpSize=2048)
```

### Atom Pair Fingerprints

```python
fpgen = AllChem.GetAtomPairGenerator()
fps_ap = [fpgen.GetSparseCountFingerprint(m) for m in mols]
fps_ap_bit = [fpgen.GetFingerprint(m) for m in mols]

# Without count simulation
fpgen_nc = AllChem.GetAtomPairGenerator(countSimulation=False)
```

### Topological Torsion Fingerprints

```python
fpgen = AllChem.GetTopologicalTorsionGenerator()
fps_tt = [fpgen.GetSparseCountFingerprint(m) for m in mols]
```

### MACCS Keys (167-bit structural keys)

```python
fps_maccs = [MACCSkeys.GenMACCSKeys(m) for m in mols]
DataStructs.TanimotoSimilarity(fps_maccs[0], fps_maccs[1])
```

---

## Similarity Metrics

```python
from rdkit import DataStructs

fp1, fp2 = fps[0], fps[1]

# Common metrics
DataStructs.TanimotoSimilarity(fp1, fp2)      # Jaccard — default for drug discovery
DataStructs.DiceSimilarity(fp1, fp2)          # 2*intersection/(|A|+|B|)
DataStructs.CosineSimilarity(fp1, fp2)
DataStructs.SokalSimilarity(fp1, fp2)
DataStructs.KulczynskiSimilarity(fp1, fp2)
DataStructs.McConnaugheySimilarity(fp1, fp2)

# Bulk similarity (query vs. library) — fast
query_fp = fpgen.GetFingerprint(mols[0])
library_fps = [fpgen.GetFingerprint(m) for m in mols[1:]]

similarities = DataStructs.BulkTanimotoSimilarity(query_fp, library_fps)
# Returns list of floats, one per molecule
```

### Find Most Similar

```python
def find_similar(query_smi: str, library: list, n=10, radius=2) -> list:
    """Return top-n most similar molecules from library."""
    fpgen = AllChem.GetMorganGenerator(radius=radius)
    qfp   = fpgen.GetFingerprint(Chem.MolFromSmiles(query_smi))
    lfps  = [fpgen.GetFingerprint(m) for m in library]

    sims = DataStructs.BulkTanimotoSimilarity(qfp, lfps)
    ranked = sorted(enumerate(sims), key=lambda x: x[1], reverse=True)
    return [(library[i], sim) for i, sim in ranked[:n]]
```

---

## Diversity Selection (MaxMin Picker)

```python
from rdkit.SimDivFilters.rdSimDivPickers import MaxMinPicker
from rdkit.Chem import AllChem
from rdkit import DataStructs

fpgen = AllChem.GetMorganGenerator(radius=2)
fps   = [fpgen.GetFingerprint(m) for m in mols]

picker = MaxMinPicker()

# Pick 10 maximally diverse molecules
n_pick = min(10, len(fps))
indices = picker.LazyBitVectorPick(fps, len(fps), n_pick, seed=42)
diverse = [mols[i] for i in indices]
```

---

## Bit Explanation & Visualization

```python
from rdkit.Chem import AllChem, Draw

fpgen = AllChem.GetMorganGenerator(radius=2)

# Collect bit info during fingerprint generation
ao = AllChem.AdditionalOutput()
ao.CollectBitInfoMap()
fp = fpgen.GetFingerprint(mols[0], additionalOutput=ao)

bit_info = ao.GetBitInfoMap()
# {bit_id: [(atom_idx, radius), ...], ...}

# Visualize a specific bit
bit_id = list(bit_info.keys())[0]
svg = Draw.DrawMorganBit(mols[0], bit_id, bit_info, useSVG=True)

# Draw multiple bits in a grid
on_bits = list(fp.GetOnBits())[:4]
img = Draw.DrawMorganBits(
    [(mols[0], bit, bit_info) for bit in on_bits],
    molsPerRow=4,
    subImgSize=(200, 200)
)
```

---

## Similarity Maps (atom-level contribution)

```python
from rdkit.Chem.Draw import SimilarityMaps
from rdkit.Chem import Draw

ref_mol = Chem.MolFromSmiles('CCCN(CCCCN1CCN(c2ccccc2OC)CC1)Cc1ccc2ccccc2c1')
mol     = Chem.MolFromSmiles('COc1cccc2cc(C(=O)NCCCCN3CCN(c4cccc5nccnc54)CC3)oc21')

d2d = Draw.MolDraw2DCairo(400, 400)
_, maxweight = SimilarityMaps.GetSimilarityMapForFingerprint(
    ref_mol, mol,
    SimilarityMaps.GetMorganFingerprint,   # fingerprint function
    d2d
)
d2d.FinishDrawing()
png_data = d2d.GetDrawingText()
```

---

## Legacy API (deprecated but still common in older code)

```python
from rdkit.Chem import AllChem

# Old Morgan FP API — use Generator API instead
fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
fp = AllChem.GetMorganFingerprint(mol, radius=2)   # count-based

# Old RDKit FP
from rdkit.Chem import RDKFingerprint
fp = Chem.RDKFingerprint(mol)
```

---

## Fingerprint Conversion

```python
# Bit vector → numpy array
import numpy as np

fp = fpgen.GetFingerprint(mol)
arr = np.zeros(len(fp), dtype=np.uint8)
DataStructs.ConvertToNumpyArray(fp, arr)

# Sparse count → dense array
count_fp = fpgen.GetSparseCountFingerprint(mol)
dense    = count_fp.GetNonzeroElements()   # dict {bit: count}

# Stack many fingerprints as matrix (for sklearn)
def fps_to_matrix(mols, fpgen) -> np.ndarray:
    fps = [fpgen.GetFingerprint(m) for m in mols]
    X   = np.zeros((len(fps), len(fps[0])), dtype=np.uint8)
    for i, fp in enumerate(fps):
        DataStructs.ConvertToNumpyArray(fp, X[i])
    return X
```
