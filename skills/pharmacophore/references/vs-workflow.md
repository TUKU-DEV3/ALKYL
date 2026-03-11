# Pharmacophore-Based Virtual Screening Workflow

End-to-end pipeline: pharmacophore query → conformer generation → scoring → hit list.

---

## Full Pipeline

```
Library (SMILES/SDF)
    │
    ▼
Pre-filter (MW, logP, rotatable bonds)        → reject PAINS / Ro5 failures
    │
    ▼
3D Conformer generation (ETKDGv3 + MMFF)      → ~20-50 confs per molecule
    │
    ▼
Pharmacophore matching (RDKit Pharm3D)         → binary pass/fail
    │
    ▼
Pharmacophore scoring (feature overlap)        → ranked score 0-1
    │
    ▼
Exclusion volume filter                        → remove pocket clashes
    │
    ▼
Docking (Vina / Gnina) — optional refinement  → ΔG score
    │
    ▼
Diversity selection (MaxMin)                   → final hit list
```

---

## Step 1 — Pre-filter Library

```python
from rdkit import Chem
from rdkit.Chem import Descriptors, FilterCatalog

def prefilter(mol):
    """Apply Lipinski + PAINS filter."""
    if mol is None:
        return False, "invalid"
    mw = Descriptors.ExactMolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)

    if mw > 600 or logp > 6 or hbd > 5 or hba > 10:
        return False, "Lipinski"

    # PAINS filter
    params = FilterCatalog.FilterCatalogParams()
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog.FilterCatalog(params)
    if catalog.HasMatch(mol):
        return False, "PAINS"

    return True, "OK"

def load_and_filter(sdf_file):
    supplier = Chem.SDMolSupplier(sdf_file, removeHs=False)
    passed = []
    rejected = Counter()
    for mol in supplier:
        ok, reason = prefilter(mol)
        if ok:
            passed.append(mol)
        else:
            rejected[reason] += 1
    print(f"Passed: {len(passed)}, Rejected: {dict(rejected)}")
    return passed
```

---

## Step 2 — Parallel Conformer Generation

```python
from concurrent.futures import ProcessPoolExecutor
from rdkit.Chem import AllChem
from functools import partial

def gen_confs_single(mol, n_confs=50, energy_window=8.0):
    """Generate conformers for one molecule."""
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.numThreads = 1
    params.pruneRmsThresh = 0.5
    AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, params=params)
    if mol.GetNumConformers() == 0:
        return None
    res = AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=200)
    energies = [e for ok, e in res if ok == 0]
    if not energies:
        return mol
    e_min = min(energies)
    keep = [i for i, (ok, e) in enumerate(res)
            if ok == 0 and (e - e_min) < energy_window]
    for cid in reversed(range(mol.GetNumConformers())):
        if cid not in keep:
            mol.RemoveConformer(cid)
    return mol

def parallel_conformers(mols, n_workers=8, **kwargs):
    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        results = list(ex.map(partial(gen_confs_single, **kwargs), mols))
    return [m for m in results if m is not None]

library_3d = parallel_conformers(filtered_mols, n_workers=8)
print(f"3D conformers generated for {len(library_3d)} molecules")
```

---

## Step 3 — Pharmacophore Screening

```python
import os
import numpy as np
from rdkit import RDConfig
from rdkit.Chem import rdMolChemicalFeatures
from tqdm import tqdm

fdefName = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
factory = rdMolChemicalFeatures.BuildFeatureFactory(fdefName)

def score_pharmacophore_match(mol, pharm_points, factory,
                               tolerance=1.5, min_features=None):
    """
    Score molecule against pharmacophore.
    Returns (score, best_conf_id) or (0, -1) if no match.
    score = fraction of required features matched.
    """
    if min_features is None:
        min_features = len(pharm_points)  # all must match

    best_score = 0.0
    best_cid = -1

    for cid in range(mol.GetNumConformers()):
        feats = factory.GetFeaturesForMol(mol, confId=cid)

        # Group feature positions by family
        feat_by_family = {}
        for feat in feats:
            fam = feat.GetFamily()
            pos = feat.GetPos(cid)
            feat_by_family.setdefault(fam, []).append(
                np.array([pos.x, pos.y, pos.z])
            )

        # Count matched pharmacophore points
        matched = 0
        for pp in pharm_points:
            fam = pp['family']
            centroid = np.array([pp['x'], pp['y'], pp['z']])
            radius = pp.get('radius', tolerance)
            positions = feat_by_family.get(fam, [])
            if any(np.linalg.norm(p - centroid) <= radius for p in positions):
                matched += 1

        score = matched / len(pharm_points)
        if matched >= min_features and score > best_score:
            best_score = score
            best_cid = cid

    return best_score, best_cid

def screen_library(library_3d, pharm_points, factory,
                    score_cutoff=0.8, tolerance=1.5):
    """Screen 3D conformer library against pharmacophore."""
    hits = []
    for mol in tqdm(library_3d, desc="Screening"):
        score, best_cid = score_pharmacophore_match(
            mol, pharm_points, factory, tolerance=tolerance
        )
        if score >= score_cutoff and best_cid >= 0:
            hits.append({
                'mol': mol,
                'score': score,
                'conf_id': best_cid,
                'smiles': Chem.MolToSmiles(Chem.RemoveHs(mol)),
            })

    hits.sort(key=lambda x: -x['score'])
    print(f"Hits: {len(hits)} / {len(library_3d)} screened")
    return hits
```

---

## Step 4 — Exclusion Volume Filter

```python
def passes_exclusion_volumes(mol, conf_id, ev_spheres, overlap_tol=0.5):
    """Return True if molecule does not clash with exclusion volumes."""
    conf = mol.GetConformer(conf_id)
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'H':
            continue
        atom_pos = np.array(conf.GetAtomPosition(atom.GetIdx()))
        for ev in ev_spheres:
            dist = np.linalg.norm(atom_pos - np.array(ev['pos']))
            if dist < (ev['radius'] - overlap_tol):
                return False
    return True

# Filter hits
clean_hits = [
    h for h in hits
    if passes_exclusion_volumes(h['mol'], h['conf_id'], ev_spheres)
]
print(f"After EV filter: {len(clean_hits)} hits")
```

---

## Step 5 — Diversity Selection

```python
from rdkit.Chem import AllChem, DataStructs
import numpy as np

def maxmin_diversity(hits, n_select=50):
    """Select diverse subset using MaxMin on Morgan FP."""
    smiles_list = [h['smiles'] for h in hits]
    mols = [Chem.MolFromSmiles(s) for s in smiles_list]
    fps = [AllChem.GetMorganFingerprintAsBitVect(m, 2, 1024) for m in mols if m]

    selected = [0]  # start with best-scoring hit
    while len(selected) < min(n_select, len(hits)):
        max_min_dist = -1
        best_idx = -1
        for i in range(len(fps)):
            if i in selected:
                continue
            min_sim = min(DataStructs.TanimotoSimilarity(fps[i], fps[j])
                          for j in selected)
            if min_sim > max_min_dist:
                max_min_dist = min_sim
                best_idx = i
        selected.append(best_idx)

    return [hits[i] for i in selected]

final_hits = maxmin_diversity(clean_hits, n_select=50)
```

---

## Step 6 — Export and Reporting

```python
from rdkit.Chem import PandasTools
import pandas as pd

def export_hits(hits, output_sdf='hits.sdf', output_csv='hits.csv'):
    """Export pharmacophore hits."""
    writer = Chem.SDWriter(output_sdf)
    records = []

    for h in hits:
        mol = Chem.RemoveHs(h['mol'])
        mol.SetProp('_PharmScore', f"{h['score']:.3f}")
        mol.SetProp('_SMILES', h['smiles'])
        writer.write(mol, confId=h['conf_id'])
        records.append({
            'SMILES': h['smiles'],
            'PharmScore': round(h['score'], 3),
        })

    writer.close()
    df = pd.DataFrame(records)
    df.to_csv(output_csv, index=False)
    print(f"Saved {len(hits)} hits to {output_sdf} and {output_csv}")
    return df

df_hits = export_hits(final_hits)
print(df_hits.head(10))
```

---

## Pharmit — Web-Based Pharmacophore Search

For large-scale VS (>1M compounds) use Pharmit (pharmit.csb.pitt.edu) with ZINC/ChEMBL databases:

```python
import requests
import json

def search_pharmit(pharm_json_file, receptor_pdb=None, max_hits=1000):
    """
    Submit pharmacophore search to Pharmit web server.
    Returns list of SMILES hits with scores.
    """
    with open(pharm_json_file) as f:
        query = json.load(f)

    # Pharmit API endpoint
    url = 'https://pharmit.csb.pitt.edu/pharmit/search'

    payload = {
        'points': query['points'],
        'max-hits': max_hits,
        'min-weight': 150,
        'max-weight': 500,
        'database': 'chembl',  # or 'zinc'
    }

    resp = requests.post(url, json=payload, timeout=300)
    if resp.status_code == 200:
        results = resp.json()
        return results.get('hits', [])
    else:
        raise RuntimeError(f"Pharmit error: {resp.status_code}")

# hits = search_pharmit('pharmacophore.json', max_hits=500)
```

---

## Enrichment Metrics

```python
from sklearn.metrics import roc_auc_score, roc_curve
import matplotlib.pyplot as plt

def compute_enrichment(actives_scores, inactives_scores):
    """Compute ROC AUC and EF at 1%, 5%, 10%."""
    labels = [1]*len(actives_scores) + [0]*len(inactives_scores)
    scores = actives_scores + inactives_scores

    auc = roc_auc_score(labels, scores)

    # Enrichment factors
    n = len(labels)
    n_act = len(actives_scores)
    order = sorted(range(n), key=lambda i: -scores[i])
    sorted_labels = [labels[i] for i in order]

    def ef(frac):
        top = max(1, int(n * frac))
        return (sum(sorted_labels[:top]) / top) / (n_act / n)

    print(f"ROC AUC: {auc:.3f}")
    print(f"EF(1%):  {ef(0.01):.1f}x")
    print(f"EF(5%):  {ef(0.05):.1f}x")
    print(f"EF(10%): {ef(0.10):.1f}x")

    # ROC curve
    fpr, tpr, _ = roc_curve(labels, scores)
    plt.figure(figsize=(6, 5))
    plt.plot(fpr, tpr, lw=2, label=f'AUC = {auc:.3f}')
    plt.plot([0,1],[0,1],'k--', lw=1)
    plt.xlabel('FPR'); plt.ylabel('TPR')
    plt.title('Pharmacophore VS ROC')
    plt.legend(); plt.tight_layout()
    plt.savefig('roc.png', dpi=150)
```
