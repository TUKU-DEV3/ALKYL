# Ligand-Based Pharmacophore

Derive a pharmacophore from a set of active ligands when no protein structure is available. Requires 3D alignment of conformers to find the common feature arrangement.

---

## Workflow

```
1. Input: set of actives (SMILES/SDF) + optionally inactives
2. Generate conformer ensembles for each active
3. Align conformers (shape or feature-based)
4. Find common pharmacophore features across all actives
5. Optional: use inactives to add exclusion volumes
6. Validate: score actives vs. inactives (ROC/EF)
```

---

## Step 1 — Conformer Generation

```python
from rdkit import Chem
from rdkit.Chem import AllChem

def generate_conformers(smiles_list, n_confs=100, energy_window=10.0):
    """Generate and prune conformer ensembles."""
    mols = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        mol = Chem.AddHs(mol)

        # ETKDGv3 + MMFF
        params = AllChem.ETKDGv3()
        params.numThreads = 0        # use all CPU cores
        params.pruneRmsThresh = 0.5  # remove similar conformers
        AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, params=params)

        # MMFF minimization + energy window filter
        res = AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=500)
        energies = [e for ok, e in res if ok == 0]
        if not energies:
            mols.append(mol)
            continue

        e_min = min(energies)
        ids_to_keep = [i for i, (ok, e) in enumerate(res)
                       if ok == 0 and (e - e_min) < energy_window]

        # Prune conformers outside energy window
        for cid in reversed(range(mol.GetNumConformers())):
            if cid not in ids_to_keep:
                mol.RemoveConformer(cid)

        mols.append(mol)
        print(f"{Chem.MolToSmiles(mol)}: {mol.GetNumConformers()} confs")

    return mols

actives = generate_conformers(active_smiles, n_confs=100)
```

---

## Step 2 — Extract Features from All Actives

```python
import os
import numpy as np
from rdkit import RDConfig
from rdkit.Chem import rdMolChemicalFeatures

fdefName = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
factory = rdMolChemicalFeatures.BuildFeatureFactory(fdefName)

def get_all_features(mol):
    """Get features for all conformers of a molecule."""
    all_feats = []
    for cid in range(mol.GetNumConformers()):
        feats = factory.GetFeaturesForMol(mol, confId=cid)
        for feat in feats:
            pos = feat.GetPos(cid)
            all_feats.append({
                'family': feat.GetFamily(),
                'pos': np.array([pos.x, pos.y, pos.z]),
                'conf_id': cid,
            })
    return all_feats

# Get feature sets for each active
active_feat_sets = [get_all_features(mol) for mol in actives]
```

---

## Step 3 — Common Feature Pharmacophore

```python
from sklearn.cluster import DBSCAN
from collections import Counter

def find_common_features(feat_sets, family, min_actives_frac=0.7,
                          cluster_eps=2.0, min_samples=1):
    """
    Find feature positions common to at least min_actives_frac of actives.
    Uses DBSCAN to cluster positions across all actives.
    """
    # Collect all feature positions for this family
    positions = []
    mol_origins = []
    for mol_idx, feats in enumerate(feat_sets):
        for feat in feats:
            if feat['family'] == family:
                positions.append(feat['pos'])
                mol_origins.append(mol_idx)

    if not positions:
        return []

    positions = np.array(positions)
    mol_origins = np.array(mol_origins)

    # Cluster nearby positions
    db = DBSCAN(eps=cluster_eps, min_samples=min_samples, metric='euclidean')
    labels = db.fit_predict(positions)

    # Filter clusters that span enough actives
    n_actives = len(feat_sets)
    min_count = int(n_actives * min_actives_frac)
    common = []

    for cluster_id in set(labels):
        if cluster_id == -1:  # noise
            continue
        mask = labels == cluster_id
        mols_in_cluster = set(mol_origins[mask])
        if len(mols_in_cluster) >= min_count:
            centroid = positions[mask].mean(axis=0)
            common.append({
                'family': family,
                'centroid': centroid,
                'n_actives': len(mols_in_cluster),
                'coverage': len(mols_in_cluster) / n_actives,
            })

    return sorted(common, key=lambda x: -x['coverage'])

# Find common features for each family
families = ['HBDonor', 'HBAcceptor', 'Aromatic', 'Hydrophobe',
            'PosIonizable', 'NegIonizable']

pharmacophore = []
for family in families:
    common = find_common_features(active_feat_sets, family,
                                   min_actives_frac=0.8)
    pharmacophore.extend(common)
    for c in common:
        print(f"{family:15s} coverage={c['coverage']:.0%} "
              f"pos=({c['centroid'][0]:.1f}, {c['centroid'][1]:.1f}, {c['centroid'][2]:.1f})")
```

---

## Step 4 — Alignment (Shape-Based)

For proper ligand-based pharmacophore, actives must be aligned in 3D space first.

```python
from rdkit.Chem import rdMolAlign, AllChem

def align_to_reference(ref_mol, probe_mols, ref_conf_id=0):
    """Align each probe to best-matching conformer of reference."""
    aligned = []
    for mol in probe_mols:
        best_rms = float('inf')
        best_cid = -1
        for cid in range(mol.GetNumConformers()):
            try:
                rms = rdMolAlign.CalcRMS(ref_mol, mol,
                                          prbCid=cid, refCid=ref_conf_id)
                if rms < best_rms:
                    best_rms = rms
                    best_cid = cid
            except Exception:
                continue
        if best_cid >= 0:
            aligned.append((mol, best_cid, best_rms))

    return sorted(aligned, key=lambda x: x[2])

# O3A (Open3DAlign) — shape + feature overlay
def o3a_align(ref_mol, probe_mol, ref_conf_id=0, probe_conf_id=0):
    """RDKit Open3DAlign: shape + atom feature overlay."""
    pyO3A = AllChem.GetO3A(probe_mol, ref_mol,
                            prbCid=probe_conf_id, refCid=ref_conf_id)
    score = pyO3A.Score()
    rms = pyO3A.Align()  # modifies probe in place
    return score, rms

# MMFF-based alignment (crippen)
def crippen_align(ref_mol, probe_mol, ref_cid=0, probe_cid=0):
    AllChem.MMFFGetMoleculeForceField(probe_mol, confId=probe_cid)
    rms, _ = rdMolAlign.GetCrippenO3A(probe_mol, ref_mol,
                                        prbCid=probe_cid, refCid=ref_cid).Align()
    return rms
```

---

## Step 5 — Validate with Inactives (ROC/EF)

```python
import numpy as np
from sklearn.metrics import roc_auc_score

def score_mol_against_pharmacophore(mol, pharmacophore, factory,
                                     n_confs=50, tolerance=1.5):
    """Score molecule: fraction of pharmacophore features matched."""
    mol_h = Chem.AddHs(mol)
    AllChem.EmbedMultipleConfs(mol_h, numConfs=n_confs, params=AllChem.ETKDGv3())
    AllChem.MMFFOptimizeMoleculeConfs(mol_h)

    best_score = 0.0
    for cid in range(mol_h.GetNumConformers()):
        feats = factory.GetFeaturesForMol(mol_h, confId=cid)
        feat_pos = {f['family']: [] for f in pharmacophore}
        for feat in feats:
            fam = feat.GetFamily()
            if fam in feat_pos:
                pos = feat.GetPos(cid)
                feat_pos[fam].append(np.array([pos.x, pos.y, pos.z]))

        # Count matched features
        matched = 0
        for pharm_feat in pharmacophore:
            fam = pharm_feat['family']
            centroid = pharm_feat['centroid']
            positions = feat_pos.get(fam, [])
            for pos in positions:
                if np.linalg.norm(pos - centroid) <= tolerance:
                    matched += 1
                    break

        score = matched / len(pharmacophore)
        if score > best_score:
            best_score = score

    return best_score

# Compute scores for actives and inactives
active_scores = [score_mol_against_pharmacophore(m, pharmacophore, factory)
                 for m in actives]
inactive_scores = [score_mol_against_pharmacophore(m, pharmacophore, factory)
                   for m in inactives]

labels = [1] * len(active_scores) + [0] * len(inactive_scores)
scores = active_scores + inactive_scores
auc = roc_auc_score(labels, scores)
print(f"ROC AUC = {auc:.3f}")

# Enrichment factor at 1%
def enrichment_factor(labels, scores, frac=0.01):
    n = len(labels)
    top_n = max(1, int(n * frac))
    sorted_idx = np.argsort(scores)[::-1]
    top_labels = [labels[i] for i in sorted_idx[:top_n]]
    n_actives = sum(labels)
    ef = (sum(top_labels) / top_n) / (n_actives / n)
    return ef

ef1 = enrichment_factor(labels, scores, frac=0.01)
print(f"EF(1%) = {ef1:.1f}")
```

---

## Step 6 — Export and Reuse

```python
import json

def save_pharmacophore(points, outfile='ligand_based.json'):
    """Save ligand-based pharmacophore for later use."""
    data = []
    for p in points:
        data.append({
            'family': p['family'],
            'x': float(p['centroid'][0]),
            'y': float(p['centroid'][1]),
            'z': float(p['centroid'][2]),
            'radius': 1.5,
            'coverage': float(p['coverage']),
        })
    with open(outfile, 'w') as f:
        json.dump(data, f, indent=2)

def load_pharmacophore(infile):
    with open(infile) as f:
        return json.load(f)
```
