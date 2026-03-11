# Conformal Prediction with MAPIE

## Installation

```bash
pip install mapie
# Optional: pip install uncertainty-toolbox  # calibration plots
```

## Regression with Split Conformal (MAPIE)

```python
from mapie.regression import MapieRegressor
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.pipeline import Pipeline
import numpy as np

# Any sklearn-compatible base model
base = RandomForestRegressor(n_estimators=200, random_state=42)

# Methods:
#   "naive"     — simple split conformal (faster, less efficient)
#   "plus"      — cross-conformal with leave-one-out (CV+, better coverage)
#   "minmax"    — conservative; guaranteed marginal coverage
#   "cv"        — jackknife-based (cross-conformal)
mapie_reg = MapieRegressor(
    estimator=base,
    method="plus",  # recommended: best coverage-efficiency tradeoff
    cv=5,           # k-fold for CV+; or "prefit" if model already fitted
    random_state=42,
)
mapie_reg.fit(X_train, y_train)

# Predict with intervals at multiple confidence levels
alpha = [0.05, 0.10, 0.20]  # 95%, 90%, 80% CI
y_pred, y_pi = mapie_reg.predict(X_test, alpha=alpha)
# y_pred:  (n_samples,)
# y_pi:    (n_samples, 2, n_alpha)  — [lower, upper] for each alpha

for i, a in enumerate(alpha):
    lower = y_pi[:, 0, i]
    upper = y_pi[:, 1, i]
    coverage = np.mean((y_test >= lower) & (y_test <= upper))
    width = np.mean(upper - lower)
    print(f"α={a:.2f}: coverage={coverage:.3f} (target={(1-a):.2f}), "
          f"mean_width={width:.3f}")
```

## Using a Pre-fitted Model ("prefit")

```python
# When base model is already trained on a separate train set
# Provide a calibration set explicitly

mapie_prefit = MapieRegressor(estimator=fitted_model, cv="prefit")
mapie_prefit.fit(X_cal, y_cal)   # calibration only (no retraining)
y_pred, y_pi = mapie_prefit.predict(X_test, alpha=0.10)
```

## Classification with Conformal Sets

```python
from mapie.classification import MapieClassifier
from sklearn.ensemble import RandomForestClassifier

# For binary or multiclass classification
base_clf = RandomForestClassifier(n_estimators=100, random_state=42)

mapie_clf = MapieClassifier(
    estimator=base_clf,
    method="score",        # "score" (softmax-based) or "cumulated_score" (RAPS)
    cv=5,
    random_state=42,
)
mapie_clf.fit(X_train, y_train)

# Predict: returns prediction set (may contain 0, 1, or 2+ classes)
alpha = 0.10
y_pred_class, y_pred_set = mapie_clf.predict(X_test, alpha=alpha, include_last_label=True)
# y_pred_set: (n_samples, n_classes) boolean array
# Row = True for all classes in the conformal prediction set

# Coverage
n_correct = sum(y_pred_set[i, y_test[i]] for i in range(len(y_test)))
coverage = n_correct / len(y_test)
print(f"Classification coverage: {coverage:.3f} (target={(1-alpha):.2f})")

# Average prediction set size (smaller = more informative)
avg_set_size = y_pred_set.sum(axis=1).mean()
print(f"Avg set size: {avg_set_size:.2f}")
```

## Molecular QSAR Example: pIC50 Prediction with Conformal CI

```python
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import numpy as np
import pandas as pd
from mapie.regression import MapieRegressor
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline

def smiles_to_morgan(smiles_list, radius=2, n_bits=2048) -> np.ndarray:
    fps = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius, n_bits)
            fps.append(list(fp))
        else:
            fps.append([0] * n_bits)
    return np.array(fps)

# Prepare features
X_train = smiles_to_morgan(train_smiles)
X_test  = smiles_to_morgan(test_smiles)
y_train = np.array(train_pIC50)
y_test  = np.array(test_pIC50)

# Conformal QSAR model
base_gbm = GradientBoostingRegressor(
    n_estimators=200, max_depth=4, learning_rate=0.05, random_state=42
)
mapie = MapieRegressor(base_gbm, method="plus", cv=5)
mapie.fit(X_train, y_train)

# Predict
alpha = 0.10
y_pred, y_pi = mapie.predict(X_test, alpha=alpha)
lower, upper = y_pi[:, 0, 0], y_pi[:, 1, 0]

# Results table
results = pd.DataFrame({
    "smiles": test_smiles,
    "pIC50_true": y_test,
    "pIC50_pred": y_pred,
    "ci_lower": lower,
    "ci_upper": upper,
    "ci_width": upper - lower,
    "in_ci": (y_test >= lower) & (y_test <= upper),
})
print(f"Coverage: {results['in_ci'].mean():.3f}")
print(f"Mean CI width: {results['ci_width'].mean():.3f}")
```

## Conditional Coverage (Local Conformalization)

Standard conformal gives marginal coverage (average over all test points), but not per-molecule. For local/conditional coverage:

```python
from mapie.regression import MapieRegressor
from mapie.conformity_scores import GammaConformityScore, AbsoluteConformityScore

# GammaConformityScore: multiplicative residuals (better for heteroscedastic data)
# AbsoluteConformityScore: standard |y - ŷ| (default)
mapie_gamma = MapieRegressor(
    estimator=base,
    conformity_score=GammaConformityScore(),  # relative intervals
    method="plus",
    cv=5,
)
mapie_gamma.fit(X_train, y_train)

# Mondrian conformal: stratify calibration by chemical class
# (ensures coverage per scaffold or activity bin)
def mondrian_conformal(X_train, y_train, X_test, y_test,
                        train_groups, test_groups, alpha=0.10):
    """Separate calibration quantile per group."""
    y_preds = {}
    base = GradientBoostingRegressor(n_estimators=100).fit(X_train, y_train)

    # Compute residuals on training set (leave-one-out or held-out cal)
    # Group-specific quantiles
    from sklearn.model_selection import cross_val_predict
    y_train_pred = cross_val_predict(base, X_train, y_train, cv=5)
    residuals = np.abs(y_train - y_train_pred)

    results = []
    for group in np.unique(test_groups):
        # Calibration residuals for this group only
        mask_cal = train_groups == group
        if mask_cal.sum() < 10:
            q = np.quantile(residuals, 1 - alpha)  # fallback: global quantile
        else:
            q = np.quantile(residuals[mask_cal], 1 - alpha)

        mask_test = test_groups == group
        y_g = base.predict(X_test[mask_test])
        results.append((y_g - q, y_g + q, y_test[mask_test]))

    return results
```

## Calibration Diagnostics

```python
import matplotlib.pyplot as plt

def calibration_plot(mapie_model, X_cal, y_cal, n_alpha=20):
    """Reliability diagram: stated confidence vs empirical coverage."""
    alphas = np.linspace(0.02, 0.50, n_alpha)
    coverages = []

    for a in alphas:
        _, y_pi = mapie_model.predict(X_cal, alpha=float(a))
        lower = y_pi[:, 0, 0]
        upper = y_pi[:, 1, 0]
        cov = np.mean((y_cal >= lower) & (y_cal <= upper))
        coverages.append(cov)

    stated = 1 - alphas
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.plot([0, 1], [0, 1], "k--", label="Perfect calibration")
    ax.plot(stated, coverages, "o-", color="steelblue", label="Model")
    ax.fill_between(stated, stated - 0.05, stated + 0.05, alpha=0.2, color="gray",
                    label="±5% tolerance")
    ax.set_xlabel("Stated confidence (1-α)")
    ax.set_ylabel("Empirical coverage")
    ax.set_title("Conformal Calibration Plot")
    ax.legend()
    return fig, dict(zip(stated, coverages))


def compute_ece(coverages: dict) -> float:
    """Expected Calibration Error."""
    stated = np.array(list(coverages.keys()))
    actual = np.array(list(coverages.values()))
    return float(np.mean(np.abs(stated - actual)))
```

## Uncertainty-Toolbox Integration

```python
# pip install uncertainty-toolbox
import uncertainty_toolbox as uct

# uct requires: y_mean (predictions), y_std (std devs), y_true
# Derive std from conformal intervals (approximate):
y_std_approx = (upper - lower) / (2 * 1.645)  # assumes ~90% CI ≈ ±1.645σ

metrics = uct.metrics.get_all_metrics(y_pred, y_std_approx, y_test)
print(f"ECE: {metrics['cal']['ece']:.4f}")
print(f"MACE: {metrics['cal']['mace']:.4f}")
print(f"Sharpness: {metrics['sharp']['sharp']:.4f}")

# Plots
fig, axes = uct.viz.plot_calibration(y_pred, y_std_approx, y_test)
fig, ax  = uct.viz.plot_intervals(y_pred, y_std_approx, y_test, n_subset=100)
```

## Conformal in Active Learning Loop

```python
def active_learning_round(
    labeled_X, labeled_y,
    unlabeled_X, unlabeled_smiles,
    alpha=0.10,
    n_acquire=10,
    strategy="width",  # or "outside_ci"
):
    """
    Select next molecules to test using conformal uncertainty.
    strategy='width': acquire widest CI (most uncertain)
    strategy='outside_ci': acquire when predicted value exceeds threshold but outside CI
    """
    mapie = MapieRegressor(
        GradientBoostingRegressor(n_estimators=100),
        method="plus", cv=5
    )
    mapie.fit(labeled_X, labeled_y)

    y_pred, y_pi = mapie.predict(unlabeled_X, alpha=alpha)
    lower, upper = y_pi[:, 0, 0], y_pi[:, 1, 0]
    widths = upper - lower

    if strategy == "width":
        acq_scores = widths
    else:  # highest predicted value outside current CI
        acq_scores = np.where(y_pred > upper, y_pred, 0.0)

    top_idx = np.argsort(-acq_scores)[:n_acquire]
    return top_idx, unlabeled_smiles[top_idx], y_pred[top_idx], widths[top_idx]
```

## Key Pitfalls

- **"cv=5" refits the model 5×**: use `cv="prefit"` when base model is expensive to train (GNN, large GBM)
- **`method="plus"` (CV+) has marginal coverage guarantee; `method="naive"` does not** — always prefer `"plus"` or `"minmax"`
- **Do NOT calibrate on test data**: the calibration set must be held out before any model tuning
- **Distribution shift invalidates coverage**: conformal relies on exchangeability — temporal/scaffold splits violate this; always check empirical coverage on shifted splits
- **MAPIE version matters**: API changed between v0.6 and v0.8; check `mapie.__version__`
- **Interval width ≠ prediction error**: a wide interval means high uncertainty, not necessarily large error; both informative metrics for drug discovery
