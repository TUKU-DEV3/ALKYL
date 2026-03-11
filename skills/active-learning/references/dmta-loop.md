# DMTA Loop: Design-Make-Test-Analyze

## Overview

The DMTA cycle is the physical manifestation of active learning in drug discovery:

```
Design ──→ Make ──→ Test ──→ Analyze
  ↑                              │
  └──────────── (retrain) ───────┘
```

| Phase | Tool | Timeline | Cost/cpd |
|-------|------|----------|----------|
| Design | ALKYL generative/AL | 1–3 days | ~$0 |
| Make | CRO synthesis | 2–6 weeks | $500–5k |
| Test | Biochemical assay | 1–2 weeks | $50–500 |
| Analyze | QSAR model retrain | 1–2 days | ~$0 |

**Implication**: batch size must be large enough (50–200 cpds/round) to justify cycle cost.

## Round 0: Hit Identification Seed

```python
# 1. Start from known actives (literature/patent/HTS)
known_actives = load_actives('known_actives.smi')

# 2. Virtual library enumeration (R-group substitution)
from rdkit import Chem
from rdkit.Chem import AllChem

scaffold_smarts = '[*:1]c1ccc([*:2])cc1'  # example scaffold
r_groups_1 = ['NC(=O)', 'CF3', 'OMe', 'CN']
r_groups_2 = ['F', 'Cl', 'Me', 'CF3']

rxn = AllChem.ReactionFromSmarts(scaffold_smarts)
virtual_library = enumerate_library(scaffold_smarts, r_groups_1, r_groups_2)

# 3. Filter (Lipinski + PAINS)
# → chem_filter.py

# 4. Diverse seed: MaxMin on filtered library
# → chem_diversity.py --n 100

# 5. Synthesize seed set, test, feed into AL
```

## Experiment Tracking

```python
import pandas as pd
from pathlib import Path
from datetime import datetime

class DMTAExperimentLog:
    """Tracks all compounds across DMTA rounds."""

    def __init__(self, path='dmta_log.csv'):
        self.path = Path(path)
        if self.path.exists():
            self.df = pd.read_csv(path)
        else:
            self.df = pd.DataFrame(columns=[
                'smiles', 'inchikey', 'round', 'source',
                'synthesis_date', 'assay_date',
                'pIC50', 'pIC50_error', 'selectivity',
                'synthesized', 'tested', 'notes'
            ])

    def add_design(self, smiles_list, round_idx, source='al_ucb'):
        from rdkit import Chem
        from rdkit.Chem.inchi import MolToInchiKey
        rows = []
        for s in smiles_list:
            mol = Chem.MolFromSmiles(s)
            ik = MolToInchiKey(mol) if mol else ''
            rows.append({
                'smiles': s, 'inchikey': ik,
                'round': round_idx, 'source': source,
                'synthesized': False, 'tested': False
            })
        self.df = pd.concat([self.df, pd.DataFrame(rows)],
                            ignore_index=True)
        self._save()

    def record_assay_results(self, inchikey_to_pic50):
        for ik, (pic50, err) in inchikey_to_pic50.items():
            mask = self.df.inchikey == ik
            self.df.loc[mask, 'pIC50'] = pic50
            self.df.loc[mask, 'pIC50_error'] = err
            self.df.loc[mask, 'tested'] = True
            self.df.loc[mask, 'assay_date'] = datetime.today().isoformat()
        self._save()

    def get_labeled(self):
        return self.df[self.df.tested == True].dropna(subset=['pIC50'])

    def _save(self):
        self.df.to_csv(self.path, index=False)
```

## Batch Design: Balancing Exploitation, Exploration, and Diversity

### Recommended batch composition (100 compounds)
```python
def design_dmta_batch(model, pool, labeled_df, batch_size=100):
    """
    Balanced batch for DMTA round.
    40% UCB (model-guided exploitation+exploration)
    20% scaffold hops (diversity)
    20% analogs of best actives (SAR expansion)
    20% random (safety net against model bias)
    """
    X_pool = featurize(pool)
    mu, sigma = model.predict_with_uncertainty(X_pool)
    ucb = mu + 1.5 * sigma

    n_ucb = int(0.4 * batch_size)
    n_hop = int(0.2 * batch_size)
    n_sar = int(0.2 * batch_size)
    n_random = batch_size - n_ucb - n_hop - n_sar

    # UCB top-k (cluster-diversified)
    ucb_batch = cluster_then_rank(X_pool, ucb, n_ucb)

    # Scaffold-diverse exploration (Core-Set)
    hop_batch = coreset_greedy(featurize(labeled_df.smiles),
                               X_pool, n_hop)

    # SAR: analogs of top-N actives via MMPA
    top_actives = labeled_df.nlargest(5, 'pIC50').smiles.tolist()
    sar_batch = generate_mmpa_analogs(top_actives, pool, n_sar)

    # Random
    remaining = [i for i in range(len(pool))
                 if i not in set(ucb_batch + hop_batch + sar_batch)]
    rand_batch = np.random.choice(remaining, n_random, replace=False).tolist()

    all_idx = ucb_batch + hop_batch + sar_batch + rand_batch
    return [pool[i] for i in all_idx[:batch_size]]
```

## Model Update After Assay Results

```python
def update_model_with_new_data(log: DMTAExperimentLog,
                                method='rf_ensemble'):
    """Retrain model incorporating all new experimental data."""
    labeled = log.get_labeled()
    X = featurize(labeled.smiles.tolist())
    y = labeled.pIC50.values

    if method == 'rf_ensemble':
        from sklearn.ensemble import RandomForestRegressor
        model = RandomForestRegressor(n_estimators=200, min_samples_leaf=2,
                                      random_state=42)
        model.fit(X, y)

    elif method == 'gp':
        from sklearn.gaussian_process import GaussianProcessRegressor
        from sklearn.gaussian_process.kernels import PairwiseKernel
        kernel = PairwiseKernel(metric=tanimoto_kernel)
        model = GaussianProcessRegressor(kernel=kernel, alpha=0.3**2,
                                          normalize_y=True)
        model.fit(X, y)

    # Cross-validate on all data
    from sklearn.model_selection import cross_val_score
    cv_r2 = cross_val_score(model, X, y, cv=5, scoring='r2')
    print(f"CV R² = {cv_r2.mean():.3f} ± {cv_r2.std():.3f} "
          f"(n={len(y)} compounds)")

    return model
```

## Stopping Criteria

### 1. Potency target reached
```python
if labeled.pIC50.max() >= 8.0:   # pIC50 ≥ 8 → IC50 ≤ 10 nM
    print("Target potency reached!")
    stop = True
```

### 2. Hit rate plateau
```python
# Fraction of each round with pIC50 > 7 (IC50 < 100 nM)
hit_rates = [labeled[labeled.round == r].pIC50.gt(7).mean()
             for r in range(current_round)]
if len(hit_rates) >= 3 and max(hit_rates[-3:]) - min(hit_rates[-3:]) < 0.05:
    print(f"Hit rate plateau: {hit_rates}")
    stop = True
```

### 3. Model uncertainty collapsed
```python
# When most of pool has low uncertainty, we've learned the landscape
sigma_mean = sigma.mean()
if sigma_mean < sigma_initial * 0.1:
    print("Pool uncertainty collapsed — model saturated")
    stop = True
```

### 4. Budget exhausted
```python
if total_compounds_synthesized >= synthesis_budget:
    stop = True
```

## SAR Visualization After Each Round

```python
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
from rdkit.Chem import Draw

def dmta_sar_dashboard(log: DMTAExperimentLog):
    labeled = log.get_labeled()
    X = featurize(labeled.smiles.tolist())
    y = labeled.pIC50.values
    rounds = labeled.round.values

    # t-SNE embedding
    tsne = TSNE(n_components=2, random_state=42, perplexity=min(30, len(y)-1))
    coords = tsne.fit_transform(X)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Color by pIC50
    sc = axes[0].scatter(coords[:, 0], coords[:, 1], c=y,
                          cmap='RdYlGn', s=40, alpha=0.8)
    plt.colorbar(sc, ax=axes[0], label='pIC50')
    axes[0].set_title(f'Chemical Space (n={len(y)})')

    # Color by round
    sc2 = axes[1].scatter(coords[:, 0], coords[:, 1], c=rounds,
                           cmap='viridis', s=40, alpha=0.8)
    plt.colorbar(sc2, ax=axes[1], label='Round')
    axes[1].set_title('AL Exploration Progress')

    plt.tight_layout()
    return fig

def plot_potency_distribution(log: DMTAExperimentLog):
    """Violin plot of pIC50 per DMTA round."""
    import seaborn as sns
    labeled = log.get_labeled()
    fig, ax = plt.subplots(figsize=(10, 5))
    sns.violinplot(data=labeled, x='round', y='pIC50',
                   cut=0, ax=ax)
    ax.axhline(y=8.0, color='r', linestyle='--', label='Target (pIC50=8)')
    ax.set_title('Potency Evolution Across DMTA Rounds')
    ax.legend()
    return fig
```

## Integration with REINVENT (Generative AL)

Close the loop between DMTA results and molecular generation:

```python
# After round r, use top actives to bias REINVENT prior
# REINVENT TOML: add activity score component from updated QSAR model

toml_template = """
[parameters]
prior = "reinvent4_prior.prior"
agent = "agent.prior"

[[stage]]
type = "reinforcement-learning"
max_steps = 200

[stage.scoring]
type = "custom_product"

[[stage.scoring.component]]
[stage.scoring.component.custom_alerts]
smarts = {pains_smarts}

[[stage.scoring.component]]
[stage.scoring.component.qsar_oracle]
# Updated RF model from latest DMTA round
model_path = "qsar_round_{round}.pkl"
weight = 1.0
"""
# → See generative-design skill for full REINVENT integration
```

## Case Study: Kinase Inhibitor Campaign

```
Round 0 (seed):  100 cpds (literature + MaxMin), best IC50 = 850 nM
Round 1 (DMTA):  80 AL-designed cpds, best = 120 nM (7× improvement)
Round 2 (DMTA):  60 cpds (UCB + scaffold hops), best = 28 nM
Round 3 (DMTA):  40 cpds (SAR expansion around best scaffolds), best = 4.2 nM
Total:           280 compounds to reach 4.2 nM — vs. ~5000 for random screening
```

## Common Pitfalls

| Problem | Cause | Fix |
|---------|-------|-----|
| Model overfits after round 1 | Too few labeled data | Use GP or RF (not deep NN) early |
| AL keeps querying same scaffold | Acquisition ignores diversity | Add cluster-then-rank or DPP |
| Potency does not improve | Model noise > signal | Lower batch size, improve assay QC |
| All novel designs fail synthesis | No synthetical accessibility filter | Apply SA score (SA ≤ 4) before design |
| Selectivity ignored | Single-objective optimization | Multi-objective: pIC50 + selectivity UCB |

## Deliverables Per Round

```
Round {r} Report:
  - Compounds tested: N
  - Hit rate (IC50 < 1 µM): X%
  - Best pIC50: Y
  - Model R² (CV): Z ± σ
  - Selected queries: batch.sdf (with AL scores)
  - SAR insights: key_sar.md
  - t-SNE dashboard: sar_round_{r}.png
  - Stopping criterion: [not reached / potency / plateau / budget]
```
