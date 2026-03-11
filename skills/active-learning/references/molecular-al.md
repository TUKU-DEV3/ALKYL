# Molecular Active Learning

## Molecular Representations for AL

### Morgan Fingerprints (ECFP4) — default
```python
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

def mol_to_ecfp4(smiles, nbits=2048):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=nbits)
    return np.array(fp)

X = np.array([mol_to_ecfp4(s) for s in smiles_list])
```

### Tanimoto Kernel for GP
```python
def tanimoto_kernel(X1, X2):
    """Tanimoto similarity kernel for binary fingerprints."""
    # |X1 ∩ X2| / |X1 ∪ X2|
    dot = X1 @ X2.T
    norm1 = X1.sum(axis=1, keepdims=True)
    norm2 = X2.sum(axis=1, keepdims=True)
    return dot / (norm1 + norm2.T - dot + 1e-9)
```
Positive semi-definite (Haussler 1999) → valid GP kernel.

### Continuous descriptors (for RF/XGBoost)
```python
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
from sklearn.preprocessing import StandardScaler

desc_names = [n for n, _ in Descriptors._descList]
calc = MoleculeDescriptors.MolecularDescriptorCalculator(desc_names)
X_desc = np.array([calc.CalcDescriptors(Chem.MolFromSmiles(s)) for s in smiles])
X_desc = StandardScaler().fit_transform(X_desc)
```

### Graph neural embeddings
```python
# Use pre-trained Chemprop embeddings as features
import chemprop
model = chemprop.models.MoleculeModel.load_from_file('model.pt')
X_gnn = model.fingerprint(smiles_list)   # shape (N, hidden_size)
```

## Batch Active Learning for Molecules

### Complete Pipeline
```python
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import StandardScaler

class MolecularAL:
    def __init__(self, smiles_pool, oracle_fn, k_batch=50, n_rounds=10):
        self.smiles_pool = list(smiles_pool)
        self.oracle = oracle_fn
        self.k = k_batch
        self.n_rounds = n_rounds
        self.labeled_smiles = []
        self.labeled_y = []

    def featurize(self, smiles_list):
        fps = []
        for s in smiles_list:
            mol = Chem.MolFromSmiles(s)
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048)
            fps.append(np.array(fp))
        return np.array(fps)

    def seed(self, seed_smiles):
        """Initialize with a diverse seed set."""
        self.labeled_smiles = list(seed_smiles)
        self.labeled_y = self.oracle(seed_smiles)
        for s in seed_smiles:
            self.smiles_pool.remove(s)

    def predict_with_uncertainty(self, X_pool):
        """Ensemble variance from RF tree predictions."""
        tree_preds = np.array([t.predict(X_pool)
                               for t in self.model.estimators_])
        mu = tree_preds.mean(axis=0)
        sigma = tree_preds.std(axis=0)
        return mu, sigma

    def select_batch(self, smiles_pool, strategy='ucb', kappa=1.0):
        X_pool = self.featurize(smiles_pool)
        mu, sigma = self.predict_with_uncertainty(X_pool)

        if strategy == 'ucb':
            scores = mu + kappa * sigma
        elif strategy == 'uncertainty':
            scores = sigma
        elif strategy == 'greedy':
            scores = mu
        else:
            raise ValueError(f"Unknown strategy: {strategy}")

        # Cluster-then-rank for diversity
        from sklearn.cluster import KMeans
        n_clusters = min(self.k, len(smiles_pool))
        km = KMeans(n_clusters=n_clusters, n_init=3, random_state=42)
        km.fit(X_pool)
        batch = []
        for c in range(n_clusters):
            mask = km.labels_ == c
            if mask.any():
                idx_in_cluster = np.where(mask)[0]
                best = idx_in_cluster[np.argmax(scores[mask])]
                batch.append(smiles_pool[best])
        return batch[:self.k]

    def run(self, seed_smiles, strategy='ucb'):
        self.seed(seed_smiles)
        history = []
        for r in range(self.n_rounds):
            X = self.featurize(self.labeled_smiles)
            y = np.array(self.labeled_y)
            self.model = RandomForestRegressor(n_estimators=100,
                                              random_state=42)
            self.model.fit(X, y)
            batch = self.select_batch(self.smiles_pool, strategy)
            labels = self.oracle(batch)
            self.labeled_smiles += batch
            self.labeled_y += labels
            for s in batch:
                self.smiles_pool.remove(s)
            hit_rate = np.mean(np.array(labels) > threshold)
            history.append({'round': r, 'hit_rate': hit_rate,
                            'n_labeled': len(self.labeled_smiles)})
            print(f"Round {r}: hit_rate={hit_rate:.2f}, "
                  f"best={max(self.labeled_y):.2f}")
        return history
```

## Diversity in AL Batches

### MaxMin Diversity Selection (ALKYL script)
```python
# Use existing chem_diversity.py
import subprocess
result = subprocess.run(
    ['python', 'scripts/chem_diversity.py',
     'candidates.smi', '--n', '50', '--fingerprint', 'morgan'],
    capture_output=True, text=True
)
```

### DPP-based selection
```python
def dpp_sample(quality_scores, similarity_matrix, k):
    """
    Sample k molecules via DPP.
    L_ij = q_i * K_ij * q_j where K is similarity (Tanimoto), q is quality score.
    """
    q = np.sqrt(np.clip(quality_scores, 0, None))
    L = np.outer(q, q) * similarity_matrix
    # Eigen decomposition
    eigvals, eigvecs = np.linalg.eigh(L)
    # Approximate: greedy MAP
    selected = []
    remaining = list(range(len(q)))
    for _ in range(k):
        gains = []
        for i in remaining:
            S = selected + [i]
            sub_L = L[np.ix_(S, S)]
            gains.append(np.linalg.det(sub_L + 1e-9 * np.eye(len(S))))
        best = remaining[np.argmax(gains)]
        selected.append(best)
        remaining.remove(best)
    return selected
```

### Greedy cluster-based (recommended for speed)
```python
from sklearn.cluster import MiniBatchKMeans

def cluster_then_rank(X_pool, scores, k):
    """Cluster pool, pick best-score molecule per cluster."""
    km = MiniBatchKMeans(n_clusters=k, n_init=3, random_state=42)
    km.fit(X_pool)
    selected = []
    for c in range(k):
        mask = km.labels_ == c
        if not mask.any():
            continue
        idx = np.where(mask)[0][np.argmax(scores[mask])]
        selected.append(idx)
    return selected
```

## Handling Class Imbalance (Active vs Inactive)

Drug discovery datasets are typically 1-5% active. Strategies:

### Stratified seed
```python
# Ensure seed has known actives (from literature/prior screen)
known_actives = df[df.active == 1].sample(10).index.tolist()
known_inactives = df[df.active == 0].sample(40).index.tolist()
seed = known_actives + known_inactives
```

### Balanced resampling for model fitting
```python
from sklearn.utils import resample
n_active = (y == 1).sum()
idx_active = np.where(y == 1)[0]
idx_inactive = np.where(y == 0)[0]
# Oversample actives 5x
idx_inactive_down = resample(idx_inactive, n_samples=5*n_active,
                              replace=False, random_state=42)
idx_bal = np.concatenate([idx_active, idx_inactive_down])
model.fit(X[idx_bal], y[idx_bal])
```

### Cost-sensitive loss
```python
from sklearn.ensemble import RandomForestClassifier
model = RandomForestClassifier(class_weight={0: 1, 1: 10}, n_estimators=200)
```

## Tracking AL Progress

```python
import pandas as pd
import matplotlib.pyplot as plt

def plot_al_curve(history, metric='hit_rate'):
    df = pd.DataFrame(history)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    # Learning curve
    axes[0].plot(df['n_labeled'], df[metric], 'b-o')
    axes[0].set_xlabel('Number of labeled compounds')
    axes[0].set_ylabel(metric)
    axes[0].set_title('Active Learning Curve')

    # Cumulative hits
    if 'n_hits' in df:
        axes[1].plot(df['n_labeled'], df['n_hits'], 'g-o')
        axes[1].set_xlabel('Experiments run')
        axes[1].set_ylabel('Cumulative actives found')
        axes[1].set_title('Hit Discovery Curve')

    plt.tight_layout()
    return fig

# Compute enrichment factor at k%
def enrichment_factor(oracle_values, selected_idx, top_pct=0.01):
    threshold = np.percentile(oracle_values, (1-top_pct)*100)
    n_total_hits = (oracle_values >= threshold).sum()
    k = int(len(oracle_values) * top_pct)
    n_selected_hits = (oracle_values[selected_idx] >= threshold).sum()
    random_ef = k / len(oracle_values)
    return (n_selected_hits / k) / (n_total_hits / len(oracle_values))
```

## Practical Tips

| Issue | Fix |
|-------|-----|
| Queries cluster in one region | Add diversity term (cluster-then-rank or DPP) |
| Model slow to retrain | Use incremental learner (SGDRegressor) or warm-start RF |
| Pool too large (>1M) | Pre-filter with fast property filter; subsample for acquisition scoring |
| Imbalanced actives | Stratified seed + cost-sensitive loss + classification AL |
| Early rounds noisy | Use wide uncertainty + high kappa early; reduce kappa over rounds |
| Scaffold hopping needed | Add scaffold diversity filter to batch selection |
