# Uncertainty Integration with Active Learning

## Uncertainty Sources (from uncertainty-qsar skill)
- **GP posterior std** (σ): exact Bayesian uncertainty
- **RF/ensemble std dev**: aleatoric + epistemic mixed
- **Conformal prediction interval width** (MAPIE): coverage-guaranteed
- **MC-Dropout variance**: approximate Bayesian deep learning
- **Deep ensembles**: best empirical uncertainty for NNs

## GP-Based Acquisition (Bayesian Optimization)

### Tanimoto GP + UCB
```python
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import PairwiseKernel
from scipy.stats import norm
import numpy as np

def tanimoto_kernel(X1, X2):
    dot = X1 @ X2.T
    n1 = X1.sum(axis=1, keepdims=True)
    n2 = X2.sum(axis=1, keepdims=True)
    return dot / (n1 + n2.T - dot + 1e-9)

kernel = PairwiseKernel(metric=tanimoto_kernel)
gp = GaussianProcessRegressor(kernel=kernel, alpha=1e-3,
                               normalize_y=True, random_state=42)

# Fit
gp.fit(X_labeled, y_labeled)

# Acquisition
mu, sigma = gp.predict(X_pool, return_std=True)

# UCB
kappa = 2.0   # β^0.5 from GP-UCB theory; set higher early (exploration)
ucb_scores = mu + kappa * sigma

# EI (for maximization)
y_best = y_labeled.max()
xi = 0.01
z = (mu - y_best - xi) / (sigma + 1e-9)
ei_scores = sigma * (z * norm.cdf(z) + norm.pdf(z))

# PI (Probability of Improvement)
pi_scores = norm.cdf((mu - y_best - xi) / (sigma + 1e-9))
```

### GP Hyperparameters Matter
```python
# Noise level α corresponds to experimental measurement error
# For pIC50 assays: typically σ_noise = 0.3–0.5 log units
gp = GaussianProcessRegressor(kernel=kernel, alpha=0.3**2,
                               n_restarts_optimizer=5,
                               normalize_y=True)
# normalize_y=True critical for EI/PI (shifts mean to 0)
```

## Ensemble Uncertainty for AL

### Random Forest Committee
```python
from sklearn.ensemble import RandomForestRegressor

rf = RandomForestRegressor(n_estimators=200, min_samples_leaf=3,
                           random_state=42)
rf.fit(X_labeled, y_labeled)

# Tree-level predictions → uncertainty estimate
tree_preds = np.array([t.predict(X_pool) for t in rf.estimators_])
mu = tree_preds.mean(axis=0)
sigma = tree_preds.std(axis=0)

# UCB acquisition
scores = mu + 1.5 * sigma
```

### Deep Ensemble (for GNN/Chemprop)
```python
# Train N=5 models with different random seeds
ensemble_preds = []
for seed in range(5):
    model = train_model(X_labeled, y_labeled, seed=seed)
    preds = model.predict(X_pool)
    ensemble_preds.append(preds)

ensemble_preds = np.array(ensemble_preds)  # (5, N_pool)
mu = ensemble_preds.mean(axis=0)
sigma = ensemble_preds.std(axis=0)
```

### MC-Dropout (single model, faster)
```python
import torch
import torch.nn as nn

class MCDropoutModel(nn.Module):
    def __init__(self, in_features, hidden=256, p=0.2):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(in_features, hidden), nn.ReLU(), nn.Dropout(p),
            nn.Linear(hidden, hidden), nn.ReLU(), nn.Dropout(p),
            nn.Linear(hidden, 1)
        )

    def forward(self, x):
        return self.net(x).squeeze(-1)

    def predict_mc(self, x, n_samples=50):
        """Keep dropout active during inference."""
        self.train()  # enables dropout
        with torch.no_grad():
            preds = torch.stack([self(x) for _ in range(n_samples)])
        return preds.mean(0), preds.std(0)
```

## Conformal Prediction Integration

### MAPIE for Prediction Intervals
```python
from mapie.regression import MapieRegressor
from sklearn.ensemble import GradientBoostingRegressor

base = GradientBoostingRegressor(n_estimators=200)
mapie = MapieRegressor(base, method='plus', cv=5)
mapie.fit(X_labeled, y_labeled)

# Predict with intervals
mu, intervals = mapie.predict(X_pool, alpha=0.1)  # 90% coverage
lower = intervals[:, 0, 0]
upper = intervals[:, 1, 0]
width = upper - lower   # ← use as uncertainty proxy

# AL acquisition: wide intervals = uncertain
scores_uncertainty = width
scores_ucb = mu + 0.5 * width   # pseudo-UCB using interval width
```

### Mondrian Conformal (activity class)
```python
from mapie.classification import MapieClassifier
from sklearn.svm import SVC

svc = SVC(probability=True)
mapie_clf = MapieClassifier(svc, method='score', cv=5)
mapie_clf.fit(X_labeled, y_labeled_binary)

# Predict sets: compound is uncertain if {0,1} both in set
_, set_preds = mapie_clf.predict(X_pool, alpha=0.1,
                                  include_last_label=True)
uncertain_mask = (set_preds[:, 0, 0] == 1) & (set_preds[:, 1, 0] == 1)
# Query uncertain compounds first
uncertain_idx = np.where(uncertain_mask)[0]
```

## Calibration Before Querying

**Uncalibrated uncertainty → wrong acquisition scores.**

```python
from sklearn.calibration import calibration_curve
from sklearn.isotonic import IsotonicRegression

def calibrate_uncertainty(sigma_raw, X_cal, y_cal, mu_cal):
    """Isotonic regression to calibrate σ against actual errors."""
    errors = np.abs(y_cal - mu_cal)
    iso = IsotonicRegression(increasing=True, out_of_bounds='clip')
    iso.fit(sigma_raw[:len(errors)], errors)
    return iso

# Calibrate on held-out validation set
calibrator = calibrate_uncertainty(sigma_val, X_val, y_val, mu_val)
sigma_calibrated = calibrator.predict(sigma_pool)
```

### Expected Calibration Error (ECE)
```python
def compute_ece(sigma, errors, n_bins=10):
    """Check if 90% PI contains 90% of true values."""
    bins = np.linspace(0, sigma.max(), n_bins+1)
    ece = 0
    for lo, hi in zip(bins[:-1], bins[1:]):
        mask = (sigma >= lo) & (sigma < hi)
        if mask.sum() == 0:
            continue
        avg_error = errors[mask].mean()
        avg_sigma = sigma[mask].mean()
        ece += mask.mean() * abs(avg_error - avg_sigma)
    return ece
```

## Adaptive κ Schedule (UCB)

Start with high exploration (κ=3), anneal to exploitation (κ=0.5):
```python
def get_kappa(round_idx, n_rounds, kappa_start=3.0, kappa_end=0.5):
    """Cosine annealing of kappa."""
    progress = round_idx / max(n_rounds - 1, 1)
    return kappa_end + 0.5 * (kappa_start - kappa_end) * (1 + np.cos(np.pi * progress))
```

## GP vs Conformal vs Ensemble: When to Use

| Method | Use when | Pros | Cons |
|--------|----------|------|------|
| Tanimoto GP | <5k labeled, fingerprint features | True Bayesian, EI exact | O(n³) scaling |
| RF ensemble | General use, tabular | Fast, calibratable | σ not theoretically grounded |
| MAPIE conformal | Coverage guarantee needed | Distribution-free | Wider intervals than GP |
| MC-Dropout | GNN already trained | Cheap inference | Requires dropout architecture |
| Deep ensemble | Best uncertainty for NNs | Empirically best | 5× training cost |

## Full Integration Example

```python
class BayesianMolecularAL:
    def __init__(self, uncertainty_method='rf_ensemble', kappa_start=2.5):
        self.method = uncertainty_method
        self.kappa_start = kappa_start

    def fit_predict(self, X_labeled, y_labeled, X_pool, round_idx, n_rounds):
        kappa = get_kappa(round_idx, n_rounds,
                          kappa_start=self.kappa_start, kappa_end=0.5)

        if self.method == 'gp':
            gp.fit(X_labeled, y_labeled)
            mu, sigma = gp.predict(X_pool, return_std=True)
        elif self.method == 'rf_ensemble':
            rf.fit(X_labeled, y_labeled)
            tree_preds = np.array([t.predict(X_pool) for t in rf.estimators_])
            mu, sigma = tree_preds.mean(0), tree_preds.std(0)
        elif self.method == 'mapie':
            mapie.fit(X_labeled, y_labeled)
            mu, intervals = mapie.predict(X_pool, alpha=0.1)
            sigma = (intervals[:, 1, 0] - intervals[:, 0, 0]) / 2

        return mu + kappa * sigma  # UCB scores
```
