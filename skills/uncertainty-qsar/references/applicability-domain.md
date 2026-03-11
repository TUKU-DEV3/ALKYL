# Applicability Domain (AD) Assessment

## What is the Applicability Domain?

The AD defines the chemical space where a QSAR model's predictions are considered reliable. A molecule outside the AD may still get a prediction, but that prediction should be flagged as potentially unreliable.

**Principle**: predictions are reliable when test molecules are structurally similar to training molecules (interpolation), not when they are distant from the training set (extrapolation).

## Distance-Based AD (Tanimoto/Euclidean)

```python
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, DataStructs
from sklearn.neighbors import NearestNeighbors

def compute_tanimoto_ad(
    train_smiles: list,
    test_smiles: list,
    radius: int = 2,
    n_bits: int = 2048,
    threshold: float = None,  # None = auto-set at 95th percentile of train distances
    k: int = 1,               # k nearest neighbors
) -> dict:
    """
    Applicability domain based on Tanimoto distance to nearest training neighbor.
    Inside AD if distance to nearest train mol < threshold.
    """
    def to_fps(smiles_list):
        fps = []
        for smi in smiles_list:
            mol = Chem.MolFromSmiles(smi)
            if mol:
                fps.append(rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius, n_bits))
            else:
                fps.append(None)
        return fps

    train_fps = to_fps(train_smiles)
    test_fps  = to_fps(test_smiles)

    # For each test FP, compute nearest-neighbor Tanimoto distance in training set
    def nn_tanimoto_distance(query_fp, pool_fps):
        if query_fp is None: return 1.0
        sims = DataStructs.BulkTanimotoSimilarity(query_fp, [fp for fp in pool_fps if fp])
        return 1.0 - max(sims) if sims else 1.0

    # Auto-set threshold: 5th percentile of training set self-distances (LOO)
    if threshold is None:
        train_dists = []
        for i, fp_i in enumerate(train_fps):
            if fp_i is None: continue
            others = [fp for j, fp in enumerate(train_fps) if j != i and fp is not None]
            d = nn_tanimoto_distance(fp_i, others)
            train_dists.append(d)
        threshold = float(np.percentile(train_dists, 95))
        print(f"Auto AD threshold (95th pctile of train LOO distances): {threshold:.4f}")

    # Compute test distances
    test_distances = [nn_tanimoto_distance(fp, [fp for fp in train_fps if fp])
                      for fp in test_fps]
    in_ad = [d <= threshold for d in test_distances]

    return {
        "distances": np.array(test_distances),
        "in_ad": np.array(in_ad),
        "threshold": threshold,
        "n_in_ad": sum(in_ad),
        "fraction_in_ad": sum(in_ad) / len(in_ad),
    }
```

## k-NN Average Distance AD

More robust than single nearest neighbor:

```python
from sklearn.neighbors import NearestNeighbors
import numpy as np

def knn_ad(X_train: np.ndarray, X_test: np.ndarray, k: int = 5,
           threshold: float = None) -> dict:
    """
    AD using mean distance to k nearest training neighbors.
    Works with any feature vector (Morgan FP bits, descriptors, etc.).
    """
    nbrs = NearestNeighbors(n_neighbors=k + 1, metric="jaccard")  # Jaccard = Tanimoto for binary
    nbrs.fit(X_train)

    # Training distances (for threshold setting)
    train_dists, _ = nbrs.kneighbors(X_train)
    train_knn_dists = train_dists[:, 1:].mean(axis=1)  # skip self (dist=0)

    if threshold is None:
        threshold = float(np.percentile(train_knn_dists, 95))
        print(f"kNN AD threshold (95th pctile): {threshold:.4f}")

    # Test distances
    test_dists, _ = nbrs.kneighbors(X_test)
    test_knn_dists = test_dists[:, :k].mean(axis=1)
    in_ad = test_knn_dists <= threshold

    return {
        "knn_distances": test_knn_dists,
        "in_ad": in_ad,
        "threshold": threshold,
        "fraction_in_ad": in_ad.mean(),
    }
```

## Williams Plot (Leverage-Based AD)

The Williams plot uses the hat matrix (statistical leverage) to identify outliers and AD boundaries. Standard QSAR regulatory practice (OECD guidelines).

```python
import numpy as np
import matplotlib.pyplot as plt

def williams_plot(X_train: np.ndarray, X_test: np.ndarray,
                  y_train: np.ndarray, y_pred_train: np.ndarray,
                  y_pred_test: np.ndarray, y_test: np.ndarray,
                  p: int = None) -> tuple:
    """
    Williams plot: leverage h (x-axis) vs standardized residual (y-axis).
    Identifies:
      - Outliers (high residual): outside ±3σ horizontal lines
      - High leverage (influential): beyond h* = 3p/n vertical line
      - Both = "true outliers" — may indicate AD boundary violations

    Args:
        X_train: (n_train, n_features) feature matrix
        p: number of model parameters (defaults to n_features + 1)
    """
    n = X_train.shape[0]
    p = p or (X_train.shape[1] + 1)

    # Hat matrix leverage: h = x (XᵀX)⁻¹ xᵀ
    XtX = X_train.T @ X_train + 1e-8 * np.eye(X_train.shape[1])  # ridge for stability
    XtX_inv = np.linalg.pinv(XtX)

    def compute_leverage(X):
        return np.diag(X @ XtX_inv @ X.T)

    h_train = compute_leverage(X_train)
    h_test  = compute_leverage(X_test)

    # Standardized residuals (using train residuals to estimate σ)
    res_train = y_train - y_pred_train
    sigma = np.std(res_train)
    std_res_train = res_train / sigma
    std_res_test  = (y_test - y_pred_test) / sigma

    # AD threshold: h* = 3p/n
    h_star = 3 * p / n

    # Williams plot
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(h_train, std_res_train, alpha=0.5, color="steelblue", s=20, label="Train")
    ax.scatter(h_test, std_res_test, alpha=0.7, color="orange", s=30, label="Test")

    ax.axhline(3,  color="red", ls="--", alpha=0.7, label="±3σ")
    ax.axhline(-3, color="red", ls="--", alpha=0.7)
    ax.axvline(h_star, color="purple", ls="--", label=f"h*={h_star:.3f} (3p/n)")

    ax.set_xlabel("Leverage (h)")
    ax.set_ylabel("Standardized Residual")
    ax.set_title("Williams Plot (Applicability Domain)")
    ax.legend()

    # Identify test molecules outside AD
    outside_ad = (h_test > h_star) | (np.abs(std_res_test) > 3)
    n_outside = outside_ad.sum()
    ax.set_title(f"Williams Plot — {n_outside}/{len(h_test)} test outside AD")

    return fig, {
        "h_train": h_train, "h_test": h_test,
        "std_res_train": std_res_train, "std_res_test": std_res_test,
        "h_star": h_star, "sigma": sigma,
        "test_outside_ad": outside_ad,
        "fraction_outside": outside_ad.mean(),
    }
```

## GP Posterior Std as AD Indicator

When using a GP model, posterior std directly encodes AD:

```python
def gp_ad_threshold(model, likelihood, X_train, percentile=95):
    """
    Set AD threshold as percentile of GP posterior std on training data.
    """
    result = predict_gp(model, likelihood, X_train)
    train_stds = result["std"]
    threshold = float(np.percentile(train_stds, percentile))
    print(f"GP AD threshold ({percentile}th pctile of train σ): {threshold:.4f}")
    return threshold


def gp_ad_flag(model, likelihood, X_test, threshold):
    """Flag test predictions outside GP-based AD."""
    result = predict_gp(model, likelihood, X_test)
    in_ad = result["std"] <= threshold
    return in_ad, result["std"]
```

## Mahalanobis Distance AD

For descriptor-based features (handles correlations):

```python
from scipy.spatial.distance import mahalanobis

def mahalanobis_ad(X_train: np.ndarray, X_test: np.ndarray,
                   threshold: float = None) -> dict:
    """Mahalanobis distance from training centroid."""
    mu = X_train.mean(axis=0)
    cov = np.cov(X_train.T) + 1e-6 * np.eye(X_train.shape[1])  # regularized
    cov_inv = np.linalg.pinv(cov)

    def dm(x):
        return mahalanobis(x, mu, cov_inv)

    train_dists = np.array([dm(x) for x in X_train])
    test_dists  = np.array([dm(x) for x in X_test])

    if threshold is None:
        threshold = float(np.percentile(train_dists, 95))

    return {
        "distances": test_dists,
        "in_ad": test_dists <= threshold,
        "threshold": threshold,
        "fraction_in_ad": (test_dists <= threshold).mean(),
    }
```

## AD Visualization

```python
import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import PCA

def visualize_ad_pca(X_train, X_test, ad_flags, title="Applicability Domain"):
    """
    PCA plot: training set + test set colored by AD status.
    Green = in AD; Red = outside AD.
    """
    pca = PCA(n_components=2, random_state=42)
    X_all = np.vstack([X_train, X_test])
    Z = pca.fit_transform(X_all)

    Z_train = Z[:len(X_train)]
    Z_test  = Z[len(X_train):]

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(Z_train[:, 0], Z_train[:, 1],
               color="steelblue", alpha=0.4, s=15, label="Training")
    in_ad  = Z_test[ad_flags]
    out_ad = Z_test[~ad_flags]
    ax.scatter(in_ad[:, 0],  in_ad[:, 1],
               color="green", alpha=0.7, s=30, label="Test (in AD)")
    ax.scatter(out_ad[:, 0], out_ad[:, 1],
               color="red",   alpha=0.7, s=30, marker="x", label="Test (outside AD)")

    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.1%})")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.1%})")
    ax.set_title(title)
    ax.legend()
    return fig


def ad_prediction_scatter(y_pred, y_true, in_ad, prop_name="pIC50"):
    """Scatter plot: predicted vs actual, colored by AD status."""
    fig, ax = plt.subplots(figsize=(6, 6))
    in_mask  = np.array(in_ad)
    out_mask = ~in_mask

    ax.scatter(y_true[in_mask], y_pred[in_mask],
               color="steelblue", alpha=0.7, s=30, label="In AD")
    ax.scatter(y_true[out_mask], y_pred[out_mask],
               color="red", alpha=0.5, s=30, marker="x", label="Outside AD")

    lo = min(y_true.min(), y_pred.min()) - 0.5
    hi = max(y_true.max(), y_pred.max()) + 0.5
    ax.plot([lo, hi], [lo, hi], "k--", alpha=0.5)
    ax.set_xlabel(f"Observed {prop_name}")
    ax.set_ylabel(f"Predicted {prop_name}")
    ax.set_title("Predicted vs Actual (AD coloring)")
    ax.legend()

    from sklearn.metrics import mean_absolute_error
    mae_in  = mean_absolute_error(y_true[in_mask],  y_pred[in_mask])  if in_mask.any()  else None
    mae_out = mean_absolute_error(y_true[out_mask], y_pred[out_mask]) if out_mask.any() else None
    print(f"MAE inside AD:  {mae_in:.3f}" if mae_in else "No in-AD samples")
    print(f"MAE outside AD: {mae_out:.3f}" if mae_out else "No out-AD samples")

    return fig
```

## AD Method Comparison

| Method | Input | Cost | Captures | OECD-compliant |
|--------|-------|------|----------|----------------|
| Tanimoto NN distance | Binary FP | O(n) | Structural novelty | Yes |
| k-NN average distance | Any features | O(kn) | Structural similarity | Yes |
| Williams plot (leverage) | Descriptor matrix | O(p²) | Statistical outliers | Yes |
| Mahalanobis distance | Descriptors | O(p²) | Correlated feature space | Yes |
| GP posterior std | Any kernel | O(n³) | Epistemic uncertainty | Best |
| DModX (PCA residuals) | Descriptors | O(np) | Systematic variation | Yes |

## OECD Principle 3 (for Regulatory Submissions)

QSAR models submitted to regulatory agencies (REACH, EPA, etc.) must define an AD according to OECD Principle 3:

```
Required elements:
1. Descriptor range (min/max of training set descriptors)
2. Structural fragment coverage (training set scaffold coverage)
3. Response range (training set activity range)
4. Distance metric definition (which method + threshold)
5. Visual AD representation (Williams plot or PCA)

Typical threshold selection:
  - 95th percentile of training set self-distances (LOO)
  - Williams: h* = 3p/n; residuals |Δ| > 3σ
```

## Key Pitfalls

- **AD threshold is task-specific**: no universal threshold; validate empirically on temporal/scaffold split
- **Tanimoto AD can miss 3D/conformational novelty**: fingerprint similarity ≠ binding mode similarity
- **Williams plot requires linear model assumption**: use with caution for non-linear models (RF, GNN)
- **Large p/n ratio**: when n_features >> n_train, leverage h → 1 for all molecules (uninformative); use PCA-reduced features first
- **Mahalanobis on correlated features**: singular covariance matrix → use regularized inverse (`+ 1e-6 * I`) or PCA first
- **AD ≠ prediction reliability**: a molecule inside AD can still have large error if training data is noisy or if it represents a rare scaffold within the AD
