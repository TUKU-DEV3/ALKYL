# Active Learning for Virtual Screening & Docking

## Motivation
Docking 10M compounds is slow (days). AL surrogate model can identify top-1% actives
by docking only 2–5% of the library. The oracle is the docking score (Vina/Gnina).

## Two-Stage Architecture

```
Stage 1: Pre-screen with fast property filters (chem_filter.py)
Stage 2: AL loop with docking oracle on filtered subset
         → model trained on docking scores, queries most informative compounds
```

## Setting Up the Docking Oracle

```python
import subprocess
import pandas as pd
from pathlib import Path

class VinaDockingOracle:
    """Docking oracle wrapping Vina Python API."""
    def __init__(self, receptor_pdbqt, center, box_size,
                 exhaustiveness=8, n_poses=1):
        from vina import Vina
        self.v = Vina(sf_name='vina', cpu=4, verbosity=0)
        self.v.set_receptor(receptor_pdbqt)
        self.v.compute_vina_maps(center=center, box_size=box_size)
        self.exhaustiveness = exhaustiveness

    def __call__(self, smiles_list):
        """Dock list of SMILES, return Vina scores (negative = better)."""
        from meeko import MoleculePreparation
        from rdkit import Chem
        from rdkit.Chem import AllChem

        scores = []
        for smi in smiles_list:
            try:
                mol = Chem.MolFromSmiles(smi)
                mol = Chem.AddHs(mol)
                AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
                AllChem.MMFFOptimizeMolecule(mol)

                prep = MoleculePreparation()
                mol_setups = prep.prepare(mol)
                pdbqt_string = mol_setups[0].write_pdbqt_string()

                self.v.set_ligand_from_string(pdbqt_string)
                self.v.dock(exhaustiveness=self.exhaustiveness)
                score = self.v.score()[0]   # best pose
                scores.append(score)
            except Exception:
                scores.append(0.0)   # failed docking → 0 (neutral)
        return scores
```

## Gnina CNN Oracle (higher accuracy, slower)
```python
import subprocess
import tempfile
from pathlib import Path

class GninaDockingOracle:
    def __init__(self, receptor_pdb, center, box_size, gnina_path='gnina'):
        self.receptor = receptor_pdb
        self.center = center
        self.box = box_size
        self.gnina = gnina_path

    def __call__(self, smiles_list, n_poses=3):
        scores = []
        for smi in smiles_list:
            try:
                # Write SMILES to file
                with tempfile.NamedTemporaryFile(suffix='.smi', mode='w',
                                                 delete=False) as f:
                    f.write(smi)
                    smi_path = f.name

                out_path = smi_path.replace('.smi', '_out.sdf')
                cmd = [
                    self.gnina,
                    '--receptor', self.receptor,
                    '--ligand', smi_path,
                    '--center_x', str(self.center[0]),
                    '--center_y', str(self.center[1]),
                    '--center_z', str(self.center[2]),
                    '--size_x', str(self.box[0]),
                    '--size_y', str(self.box[1]),
                    '--size_z', str(self.box[2]),
                    '--num_modes', str(n_poses),
                    '--out', out_path,
                    '--no_refine'  # faster
                ]
                result = subprocess.run(cmd, capture_output=True,
                                        timeout=120)
                # Parse CNNscore from output
                for line in result.stdout.decode().split('\n'):
                    if 'CNNscore' in line:
                        score = -float(line.split()[-1])  # negate (higher=better)
                        scores.append(score)
                        break
                else:
                    scores.append(0.0)
            except Exception:
                scores.append(0.0)
        return scores
```

## AL Virtual Screening Pipeline

```python
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.ensemble import RandomForestRegressor

def run_al_vs(library_smiles, oracle, seed_size=200, batch_size=200,
              n_rounds=10, strategy='ucb', kappa=2.0):
    """
    Active learning virtual screening.

    Returns:
        DataFrame with all docked compounds + scores, sorted by score.
    """
    import pandas as pd
    from scripts.chem_diversity import maxmin_select  # ALKYL script

    # 1. Featurize pool
    def fps(smiles):
        fps = []
        for s in smiles:
            mol = Chem.MolFromSmiles(s)
            if mol:
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048)
                fps.append(np.array(fp))
            else:
                fps.append(np.zeros(2048))
        return np.array(fps)

    pool = list(library_smiles)

    # 2. Diverse seed selection using MaxMin
    seed_idx = maxmin_select(pool, seed_size)
    seed_smiles = [pool[i] for i in seed_idx]
    labeled_smiles = list(seed_smiles)
    labeled_scores = oracle(seed_smiles)

    for s in seed_smiles:
        pool.remove(s)

    results = pd.DataFrame({
        'smiles': labeled_smiles,
        'score': labeled_scores,
        'round': 0
    })

    # 3. AL loop
    for r in range(1, n_rounds + 1):
        X_lab = fps(labeled_smiles)
        y_lab = np.array(labeled_scores)

        # Train RF ensemble
        rf = RandomForestRegressor(n_estimators=200, min_samples_leaf=2,
                                    random_state=42)
        rf.fit(X_lab, y_lab)

        # Predict pool with uncertainty
        X_pool = fps(pool)
        tree_preds = np.array([t.predict(X_pool) for t in rf.estimators_])
        mu = tree_preds.mean(0)
        sigma = tree_preds.std(0)

        # UCB acquisition
        if strategy == 'ucb':
            acq_scores = mu + kappa * sigma
        elif strategy == 'greedy':
            acq_scores = mu
        else:
            acq_scores = sigma

        # Cluster-then-rank for diversity
        from sklearn.cluster import MiniBatchKMeans
        k = min(batch_size, len(pool))
        km = MiniBatchKMeans(n_clusters=k, n_init=3, random_state=42).fit(X_pool)
        batch_idx = [np.where(km.labels_ == c)[0][
                         np.argmax(acq_scores[km.labels_ == c])]
                     for c in range(k) if (km.labels_ == c).any()]
        batch = [pool[i] for i in batch_idx[:batch_size]]

        # Oracle queries
        batch_scores = oracle(batch)
        labeled_smiles += batch
        labeled_scores += batch_scores
        for s in batch:
            pool.remove(s)

        round_df = pd.DataFrame({
            'smiles': batch, 'score': batch_scores, 'round': r
        })
        results = pd.concat([results, round_df], ignore_index=True)

        hit_pct = np.percentile(labeled_scores, 5)
        n_hits = (np.array(labeled_scores) <= hit_pct).sum()
        print(f"Round {r}: docked {len(labeled_smiles)}, "
              f"best={min(labeled_scores):.2f}, "
              f"top-5% threshold={hit_pct:.2f}")

    return results.sort_values('score')
```

## Metrics for AL Virtual Screening

### BEDROC (Boltzmann-Enhanced Discrimination of ROC)
```python
def bedroc(scores, actives_mask, alpha=20.0):
    """
    scores: array of predicted scores (lower = better for docking)
    actives_mask: boolean array where True = active (known from reference set)
    alpha: controls early enrichment weight (20 = top 8% emphasized)
    """
    n = len(scores)
    n_actives = actives_mask.sum()
    order = np.argsort(scores)   # ascending (best docking first)
    ranks = np.where(actives_mask[order])[0] + 1  # 1-indexed

    ri_sum = np.exp(-alpha * ranks / n).sum()
    ri_random = n_actives / n * (1 - np.exp(-alpha)) / (np.exp(alpha/n) - 1)
    ra = n_actives / n

    bedroc_score = (ri_sum / ri_random * ra * np.sinh(alpha/2) /
                    (np.cosh(alpha/2) - np.cosh(alpha/2 - alpha * ra)) + 1 /
                    (1 - np.exp(alpha * (1 - ra))))
    return bedroc_score
```

### Enrichment Factor
```python
def enrichment_factor(scores, actives_mask, top_pct=0.01):
    """EF at top x%."""
    n = len(scores)
    n_actives = actives_mask.sum()
    k = max(1, int(n * top_pct))
    top_k_idx = np.argsort(scores)[:k]
    hits_in_top_k = actives_mask[top_k_idx].sum()
    ef = (hits_in_top_k / k) / (n_actives / n)
    return ef
```

### Hit Rate per Round (AL efficiency)
```python
def compute_al_efficiency(results_df, score_threshold):
    """Fraction of each round's queries that hit the threshold."""
    by_round = results_df.groupby('round').apply(
        lambda g: (g['score'] <= score_threshold).mean()
    )
    return by_round
```

## Gnina CNN Rescoring After AL

Top-N from AL surrogate → high-accuracy CNN rescore:
```python
# Select top 500 from AL model
al_top = results_df.nsmallest(500, 'score')

# Rescore with Gnina CNN
gnina_oracle = GninaDockingOracle(receptor_pdb, center, box_size)
gnina_scores = gnina_oracle(al_top['smiles'].tolist())
al_top['gnina_score'] = gnina_scores

# Final ranking by CNN score
final_hits = al_top.nsmallest(50, 'gnina_score')
```

## Pre-filtering Before AL (Reduces Pool Size)

```python
import subprocess

# Use ALKYL chem_filter.py to remove Lipinski failures and PAINS
result = subprocess.run(
    ['python', 'scripts/chem_filter.py', 'library.smi',
     '--filters', 'lipinski', 'pains',
     '--out', 'library_filtered.smi'],
    capture_output=True, text=True
)
print(result.stdout)  # summary statistics

# Library size reduction: typically 60-80% pass Lipinski+PAINS
```

## Benchmarking Parameters

| Parameter | Default | Notes |
|-----------|---------|-------|
| Seed size | 200 | 50–500; larger = better initialization |
| Batch size | 200 | 50–500; larger batches → fewer rounds |
| Rounds | 10 | Stop when hit rate plateaus |
| kappa (UCB) | 2.0 | 0.5 exploit → 3.0 explore |
| RF trees | 200 | 100 min for reliable σ |
| Pre-filter | Lipinski+PAINS | Run first, reduces oracle calls 40–80% |

## Computational Cost (rough estimates)

| Component | Per molecule | 1000 compounds |
|-----------|-------------|----------------|
| Fingerprint | 0.1 ms | 0.1 s |
| RF predict | 0.5 ms | 0.5 s |
| Vina docking | 5–30 s | 1–8 h (serial) |
| Gnina docking | 15–60 s | 4–16 h (serial) |
| Vina (4-core) | 2–8 s | 30 min–2 h |

→ AL can reduce required dockings from 1M to 20k–50k (50–200× speedup).
