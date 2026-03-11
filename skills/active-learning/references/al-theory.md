# Active Learning Theory

## Problem Formulation
- **Pool-based AL**: large unlabeled pool U, query sequentially or in batches
- **Stream-based AL**: molecules arrive one-by-one; decide query/discard online
- **Query-by-synthesis**: generate new molecules to query (generative AL)
- Goal: minimize labels needed to reach target performance / find top-k actives

## Query Strategies

### Uncertainty Sampling (exploitation of model doubt)
```python
# Highest variance / entropy
scores = model.predict_uncertainty(X_pool)   # σ, entropy, or interval width
batch_idx = np.argsort(scores)[-k:]
```
- **Regression**: query highest σ (GP posterior std, ensemble std dev)
- **Classification**: query highest entropy H = -Σ p log p, or margin = p1 - p2
- **Conformal**: query widest prediction interval (MAPIE)

Pros: simple, fast. Cons: clusters in high-uncertainty regions (redundant queries).

### Greedy / Expected Improvement (exploitation of predicted value)
```python
# EI: expected gain over best observed y*
from scipy.stats import norm
def expected_improvement(mu, sigma, y_best, xi=0.01):
    z = (mu - y_best - xi) / (sigma + 1e-9)
    return sigma * (z * norm.cdf(z) + norm.pdf(z))
```
Best for optimization (maximize pIC50). Equivalent to BO acquisition function.

### Upper Confidence Bound (UCB)
```python
scores = mu + kappa * sigma    # kappa controls exploration/exploitation tradeoff
# kappa=1.0 → balanced; kappa→0 greedy; kappa→∞ pure exploration
```
- Most popular for drug discovery AL
- Monotone in both predicted value and uncertainty

### Query by Committee (QBC)
```python
# Train ensemble of N models, query where models disagree most
predictions = np.array([m.predict(X_pool) for m in committee])
disagreement = predictions.std(axis=0)   # vote entropy for classification
batch_idx = np.argsort(disagreement)[-k:]
```
Pros: model-agnostic, diversity implicit. Cons: N × train cost.

### Core-Set / Greedy Diversity
```python
# Farthest-first traversal in feature space
from sklearn.metrics import pairwise_distances
def coreset_greedy(X_labeled, X_pool, k):
    selected = []
    dists = pairwise_distances(X_pool, X_labeled).min(axis=1)
    for _ in range(k):
        i = np.argmax(dists)
        selected.append(i)
        d_new = pairwise_distances(X_pool, X_pool[[i]]).squeeze()
        dists = np.minimum(dists, d_new)
    return selected
```
Pure exploration. Use as diversity component in hybrid strategies.

## Acquisition Function Comparison

| Strategy | Explores | Exploits | Requires | Best For |
|----------|----------|----------|----------|---------|
| Uncertainty sampling | ✓ | – | uncertainty model | Coverage, AD extension |
| Greedy EI | – | ✓ | GP/uncertainty | Single-objective optimization |
| UCB | ✓ | ✓ | uncertainty model | General optimization |
| QBC | ✓ | partial | ensemble | Model-agnostic |
| Core-set | ✓✓ | – | features only | Seed set, exploration |
| Thompson sampling | ✓ | ✓ | GP posterior | Exploration-exploitation balance |

## Hybrid Strategies

### ε-greedy
```python
if np.random.rand() < epsilon:
    batch = random_sample(pool, k)    # explore
else:
    batch = top_k_acquisition(pool, k)  # exploit
epsilon = max(epsilon_min, epsilon * decay)
```

### Weighted sum
```python
alpha = 0.5
scores = alpha * normalize(uncertainty) + (1 - alpha) * normalize(predicted_value)
```

### BALD (Bayesian Active Learning by Disagreement)
Maximizes mutual information I(y; θ | x, D):
```python
# For deep ensembles: BALD ≈ H[y|x,D] - E_θ[H[y|x,θ]]
# = predictive entropy - mean entropy of individual models
pred_entropy = entropy(mean_probs)
mean_model_entropy = mean([entropy(p) for p in ensemble_probs])
bald_score = pred_entropy - mean_model_entropy
```

## Batch Mode Selection

**Problem**: top-k greedy from acquisition yields redundant, clustered queries.

### Determinantal Point Process (DPP)
```python
# DPP selects diverse + high-quality subset
# Kernel: L_ij = quality(i) * similarity(i,j) * quality(j)
from dppy.finite_dpps import FiniteDPP
L = quality_diag @ similarity_matrix @ quality_diag
dpp = FiniteDPP('likelihood', L=L)
dpp.sample_exact_k_dpp(size=k)
```

### Greedy Submodular (practical alternative)
```python
# Maximize f(S) = Σ_i max_{j∈S} sim(i,j) + λ Σ_j score(j)
# Greedy achieves 1-1/e approximation
selected = []
for _ in range(k):
    gains = [submodular_gain(selected + [i], pool) for i in candidates]
    selected.append(candidates[np.argmax(gains)])
```

### Cluster-then-rank
```python
from sklearn.cluster import KMeans
clusters = KMeans(n_clusters=k).fit(X_pool)
# Select top-scoring molecule per cluster
batch = [pool[clusters.labels_==c][np.argmax(scores[clusters.labels_==c])]
         for c in range(k)]
```
Simple and effective. Use with Morgan/ECFP4 features.

## Convergence & Stopping

### Performance plateau
```python
if len(rounds) >= 3:
    recent_improvement = metric[-1] - metric[-3]
    if recent_improvement < delta_min:
        stop()
```

### Hit rate saturation
```python
# Stop when hit rate in queries drops below threshold
hit_rate = labeled_df[labeled_df.round == r]['active'].mean()
if hit_rate < 0.05:  # less than 5% hits in last batch → stop
    stop()
```

### Budget constraint
```python
max_experiments = 500   # fixed budget
if len(labeled_pool) >= max_experiments:
    stop()
```

## Evaluation Metrics for AL

| Metric | Formula | Purpose |
|--------|---------|---------|
| Area under learning curve (AULC) | ∫ performance(n_labels) dn | Overall AL efficiency |
| Enrichment Factor (EF) | (hits_queried/k) / (total_hits/N) | How many actives found |
| BEDROC | Boltzmann-enhanced ranking | Early enrichment |
| Time to top-k | Labels to find k best cpds | Practical efficiency |
| Hit rate per round | Actives / batch_size | Convergence signal |

## Benchmarking AL

Use **Guacamol benchmark tasks** or **TDC ADMET benchmarks** with held-out sets:
```python
# Simulate oracle: reveal labels on demand
class SimulatedOracle:
    def __init__(self, df):
        self.df = df.set_index('smiles')
    def query(self, smiles_list):
        return [self.df.loc[s, 'pIC50'] for s in smiles_list]
```

Key: always report AULC ± SD over multiple random seeds (3–5 runs minimum).
