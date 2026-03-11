# RDKit Pharmacophore — Pharm2D and Pharm3D

## ChemicalFeatures — Feature Extraction

```python
import os
from rdkit import Chem, RDConfig
from rdkit.Chem import rdMolChemicalFeatures

# Feature factory (bundled with RDKit)
fdefName = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
factory = rdMolChemicalFeatures.BuildFeatureFactory(fdefName)

mol = Chem.MolFromSmiles('c1ccc(cc1)NC(=O)c2ccc(cc2)O')  # benzamide derivative
feats = factory.GetFeaturesForMol(mol)

for feat in feats:
    family = feat.GetFamily()   # e.g. 'HBDonor', 'Aromatic'
    ftype = feat.GetType()      # e.g. 'SingleAtomDonor'
    atom_ids = list(feat.GetAtomIds())
    print(f"{family:15s} {ftype:25s} atoms={atom_ids}")

# With 3D positions (requires conformer)
from rdkit.Chem import AllChem
mol3d = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol3d, AllChem.ETKDGv3())
AllChem.MMFFOptimizeMolecule(mol3d)

feats3d = factory.GetFeaturesForMol(mol3d)
for feat in feats3d:
    pos = feat.GetPos()   # 3D centroid (requires conformer)
    print(f"{feat.GetFamily():15s} pos=({pos.x:.2f}, {pos.y:.2f}, {pos.z:.2f})")
```

---

## Pharm2D Fingerprints (Gobbi)

2D pharmacophore fingerprints encode pairs/triplets of features + topological distance (bond count). Fast, no conformer needed.

```python
from rdkit.Chem.Pharm2D import Gobbi_Pharm2D, Generate
from rdkit import DataStructs

factory = Gobbi_Pharm2D.factory  # pre-built Gobbi factory

# Generate fingerprint
mol = Chem.MolFromSmiles('c1ccc(cc1)NC(=O)c2ccc(cc2)O')
fp = Generate.Gen2DFingerprint(mol, factory)
# Returns: SparseBitVector (bit-encoded feature pairs + distances)

# Convert to numpy
import numpy as np
arr = np.zeros(fp.GetNumBits(), dtype=np.uint8)
DataStructs.ConvertToNumpyArray(fp, arr)

# Tanimoto similarity between two molecules
mol2 = Chem.MolFromSmiles('c1ccc(cc1)C(=O)Nc2ccc(cc2)N')
fp2 = Generate.Gen2DFingerprint(mol2, factory)
sim = DataStructs.TanimotoSimilarity(fp, fp2)
print(f"Pharm2D Tanimoto: {sim:.3f}")

# Batch similarity matrix
mols = [mol, mol2, Chem.MolFromSmiles('CC(=O)Oc1ccccc1C(=O)O')]
fps = [Generate.Gen2DFingerprint(m, factory) for m in mols]

n = len(fps)
sim_matrix = np.zeros((n, n))
for i in range(n):
    for j in range(n):
        sim_matrix[i, j] = DataStructs.TanimotoSimilarity(fps[i], fps[j])
print(sim_matrix)
```

---

## Pharm3D — 3D Pharmacophore Matching

Pharm3D embeds a pharmacophore query into a molecule's conformer ensemble and returns matching conformers.

### Define a 3D Pharmacophore

```python
from rdkit.Chem.Pharm3D import Pharmacophore, EmbedLib
from rdkit.Geometry import rdGeometry
from rdkit.Chem import rdDistGeom

# Create pharmacophore from feature list
# Each feature: (family, pointList, radius)
pharmacophore = Pharmacophore.Pharmacophore(feats)

# Or build manually:
from rdkit.Chem.Pharm3D import Pharmacophore as Ph3D
from rdkit.Geometry.rdGeometry import Point3D

# Define feature centroids manually (from known co-crystal)
feat_positions = {
    'HBDonor':   Point3D(2.1, 3.4, 0.5),
    'HBAcceptor': Point3D(-1.2, 1.8, 2.1),
    'Aromatic':  Point3D(0.0, 0.0, 0.0),
}

# Create from MolFeature objects
from rdkit.Chem import MolFromSmiles
from rdkit.Chem import rdMolChemicalFeatures as rcf

# Build pharmacophore from extracted features
ref_mol = Chem.MolFromPDBFile('ligand.pdb')
AllChem.EmbedMolecule(ref_mol)
ref_feats = factory.GetFeaturesForMol(ref_mol)
pcophore = Pharmacophore.Pharmacophore(ref_feats)
```

### Set Distance Bounds

```python
# Distance bounds matrix between features
n_feats = len(ref_feats)
bounds = rdDistGeom.GetMoleculeBoundsMatrix(ref_mol)

# Alternatively: compute pairwise distances from 3D positions
import numpy as np
positions = np.array([f.GetPos() for f in ref_feats])

def build_bounds_matrix(positions, tolerance=1.5):
    n = len(positions)
    lower = np.zeros((n, n))
    upper = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            d = np.linalg.norm(positions[i] - positions[j])
            lower[i, j] = lower[j, i] = max(0, d - tolerance)
            upper[i, j] = upper[j, i] = d + tolerance
    return lower, upper

lower, upper = build_bounds_matrix(
    [ref_feats[i].GetPos() for i in range(n_feats)]
)
```

### Match Query Against a Molecule

```python
from rdkit.Chem.Pharm3D import EmbedLib

def match_pharmacophore(query_feats, mol, factory, n_confs=50,
                         rmsTol=2.0, max_matches=10):
    """
    Try to embed pharmacophore into molecule's conformer ensemble.
    Returns list of (match_score, conformer_id, match_indices).
    """
    # Add Hs and generate conformers
    mol_h = Chem.AddHs(mol)
    AllChem.EmbedMultipleConfs(mol_h, numConfs=n_confs,
                               params=AllChem.ETKDGv3())
    AllChem.MMFFOptimizeMoleculeConfs(mol_h)

    mol_feats = factory.GetFeaturesForMol(mol_h)

    # Try to embed
    canMatch, allMatches = EmbedLib.MatchPharmacophoreToMol(
        mol_h, factory, query_feats
    )

    if not canMatch:
        return []

    matches = EmbedLib.MatchMolToPharmacophore(
        mol_h, factory, query_feats,
        maxMatches=max_matches,
        useRandomCoords=False,
    )
    return matches

results = match_pharmacophore(ref_feats, test_mol, factory)
print(f"Matches found: {len(results)}")
```

---

## Full 3D Pharmacophore Screening (EmbedLib)

```python
from rdkit.Chem.Pharm3D import EmbedLib

def screen_pharmacophore(query_feats, library, factory, bounds_matrix,
                          score_threshold=0.5):
    """
    Screen a list of molecules against a 3D pharmacophore.
    Returns hits as (smiles, score).
    """
    hits = []
    for mol in library:
        if mol is None:
            continue
        try:
            mol_h = Chem.AddHs(mol)
            AllChem.EmbedMultipleConfs(mol_h, numConfs=50,
                                       params=AllChem.ETKDGv3())
            AllChem.MMFFOptimizeMoleculeConfs(mol_h)

            canMatch, allMatches = EmbedLib.MatchPharmacophoreToMol(
                mol_h, factory, query_feats
            )
            if not canMatch:
                continue

            failed, matches, _, _ = EmbedLib.EmbedPharmacophore(
                mol_h, query_feats, pcophore,
                count=10, randomSeed=42
            )
            if matches:
                # Use best RMS as score
                score = 1.0 / (1.0 + min(m[0] for m in matches))
                if score >= score_threshold:
                    hits.append((Chem.MolToSmiles(mol), score))
        except Exception:
            continue

    hits.sort(key=lambda x: -x[1])
    return hits
```

---

## Pharmacophore Fingerprint Comparison

```python
from rdkit.Chem.Pharm2D import Gobbi_Pharm2D, Generate
from rdkit import DataStructs
import pandas as pd

def pharm2d_similarity_matrix(smiles_list):
    factory = Gobbi_Pharm2D.factory
    mols = [Chem.MolFromSmiles(s) for s in smiles_list]
    fps = [Generate.Gen2DFingerprint(m, factory) for m in mols if m]

    n = len(fps)
    mat = pd.DataFrame(index=smiles_list[:n], columns=smiles_list[:n])
    for i in range(n):
        for j in range(n):
            mat.iloc[i, j] = DataStructs.TanimotoSimilarity(fps[i], fps[j])
    return mat.astype(float)

# Cluster by pharmacophore similarity
from sklearn.cluster import AgglomerativeClustering
import numpy as np

sim_mat = pharm2d_similarity_matrix(smiles_list)
dist_mat = 1 - sim_mat.values

clustering = AgglomerativeClustering(
    n_clusters=None,
    distance_threshold=0.4,    # = Tanimoto similarity > 0.6
    metric='precomputed',
    linkage='complete'
)
labels = clustering.fit_predict(dist_mat)
print(f"Clusters: {max(labels)+1}")
```

---

## Visualizing Pharmacophore Features (2D)

```python
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG, display

def draw_pharmacophore_2d(mol, factory, size=(400, 300)):
    """Highlight atoms by pharmacophore feature type."""
    feats = factory.GetFeaturesForMol(mol)

    colors = {
        'HBDonor':    (0.0, 0.5, 1.0),    # blue
        'HBAcceptor': (1.0, 0.0, 0.0),    # red
        'Aromatic':   (0.8, 0.5, 0.0),    # orange
        'Hydrophobe': (0.5, 0.8, 0.0),    # green
        'PosIonizable': (0.5, 0.0, 0.8),  # purple
        'NegIonizable': (1.0, 0.6, 0.0),  # amber
    }

    atom_colors = {}
    highlight_atoms = []
    for feat in feats:
        color = colors.get(feat.GetFamily(), (0.7, 0.7, 0.7))
        for idx in feat.GetAtomIds():
            atom_colors[idx] = color
            if idx not in highlight_atoms:
                highlight_atoms.append(idx)

    drawer = rdMolDraw2D.MolDraw2DSVG(*size)
    drawer.drawOptions().addStereoAnnotation = True
    drawer.DrawMolecule(mol, highlightAtoms=highlight_atoms,
                        highlightAtomColors=atom_colors)
    drawer.FinishDrawing()
    return SVG(drawer.GetDrawingText())

display(draw_pharmacophore_2d(mol, factory))
```
