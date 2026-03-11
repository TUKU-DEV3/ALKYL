# Gaussian Processes with Tanimoto Kernel for Molecular QSAR

## Why GPs for Molecular Properties?

- Exact posterior uncertainty (not approximate)
- Posterior variance naturally captures epistemic uncertainty
- Tanimoto kernel makes GPs directly applicable to fingerprints
- Native active learning: acquisition functions (EI, UCB, Thompson sampling) use posterior
- Scales to ~5000 training points; sparse GPs for larger datasets

## Tanimoto Kernel (GPyTorch)

```python
import torch
import gpytorch
from gpytorch.kernels import Kernel

class TanimotoKernel(Kernel):
    """
    Tanimoto similarity kernel for binary fingerprints.
    k(a, b) = |a ∩ b| / |a ∪ b| = dot(a,b) / (|a|² + |b|² - dot(a,b))
    Valid positive definite kernel (Haussler 1999).
    """
    has_lengthscale = False

    def forward(self, x1: torch.Tensor, x2: torch.Tensor, **params) -> torch.Tensor:
        # x1: (batch, n, d), x2: (batch, m, d) — binary fingerprints
        x1_sq = x1.pow(2).sum(dim=-1, keepdim=True)   # |a|²: (batch, n, 1)
        x2_sq = x2.pow(2).sum(dim=-1, keepdim=True)   # |b|²: (batch, m, 1)
        dot = torch.matmul(x1, x2.transpose(-1, -2))  # dot(a,b): (batch, n, m)

        # Union = |a|² + |b|² - dot(a,b) (for binary vectors: |a|² = |a|)
        union = x1_sq + x2_sq.transpose(-1, -2) - dot
        return dot / union.clamp(min=1e-8)
```

## Exact GP Model (GPyTorch)

```python
import gpytorch

class TanimotoGP(gpytorch.models.ExactGP):
    def __init__(self, train_x, train_y, likelihood):
        super().__init__(train_x, train_y, likelihood)
        self.mean_module = gpytorch.means.ConstantMean()
        self.covar_module = gpytorch.kernels.ScaleKernel(TanimotoKernel())

    def forward(self, x):
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(x)
        return gpytorch.distributions.MultivariateNormal(mean_x, covar_x)


def train_tanimoto_gp(
    X_train: torch.Tensor,
    y_train: torch.Tensor,
    n_epochs: int = 200,
    lr: float = 0.1,
    noise_init: float = 0.1,
) -> tuple:
    """
    Train exact GP on molecular fingerprints.
    Returns: (model, likelihood) in inference mode.
    """
    likelihood = gpytorch.likelihoods.GaussianLikelihood()
    likelihood.noise = noise_init
    model = TanimotoGP(X_train, y_train, likelihood)

    model.train()
    likelihood.train()

    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    mll = gpytorch.mlls.ExactMarginalLogLikelihood(likelihood, model)

    for epoch in range(n_epochs):
        optimizer.zero_grad()
        output = model(X_train)
        loss = -mll(output, y_train)
        loss.backward()
        optimizer.step()

        if (epoch + 1) % 50 == 0:
            noise = likelihood.noise.item()
            ls = model.covar_module.outputscale.item()
            print(f"Epoch {epoch+1}/{n_epochs}  Loss={loss.item():.4f}  "
                  f"Noise={noise:.4f}  Scale={ls:.4f}")

    # Switch to inference mode (equivalent to .eval() without calling eval)
    model.train(mode=False)
    likelihood.train(mode=False)
    return model, likelihood


def predict_gp(
    model, likelihood,
    X_test: torch.Tensor,
    return_samples: int = 0,
) -> dict:
    """GP prediction returning mean, std, and optional posterior samples."""
    with torch.no_grad(), gpytorch.settings.fast_pred_var():
        preds = likelihood(model(X_test))

    mean = preds.mean.numpy()
    std  = preds.variance.clamp(min=0).sqrt().numpy()  # clamp for numerical safety

    result = {"mean": mean, "std": std,
              "lower": mean - 1.96 * std, "upper": mean + 1.96 * std}

    if return_samples > 0:
        with torch.no_grad():
            samples = preds.rsample(sample_shape=torch.Size([return_samples]))
        result["samples"] = samples.numpy()  # (n_samples, n_test)

    return result
```

## Full Molecular QSAR Workflow

```python
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def smiles_to_fp_tensor(smiles_list, radius=2, n_bits=2048) -> torch.Tensor:
    fps = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius, n_bits)
            fps.append(list(fp))
        else:
            fps.append([0] * n_bits)
    return torch.tensor(fps, dtype=torch.float32)

# Prepare data
X_train = smiles_to_fp_tensor(train_smiles)
X_test  = smiles_to_fp_tensor(test_smiles)
y_train_raw = torch.tensor(train_pIC50, dtype=torch.float32)
y_test  = np.array(test_pIC50)

# Normalize y (recommended for GP: zero-mean prior)
y_mean, y_std = y_train_raw.mean(), y_train_raw.std()
y_train_norm = (y_train_raw - y_mean) / y_std

# Train
model, likelihood = train_tanimoto_gp(X_train, y_train_norm, n_epochs=300, lr=0.05)

# Predict and de-normalize
result = predict_gp(model, likelihood, X_test)
pred_mean = result["mean"] * y_std.item() + y_mean.item()
pred_std  = result["std"]  * y_std.item()

# Evaluate
from sklearn.metrics import mean_absolute_error, r2_score
mae = mean_absolute_error(y_test, pred_mean)
r2  = r2_score(y_test, pred_mean)
print(f"MAE={mae:.3f}  R²={r2:.3f}  Mean σ={pred_std.mean():.3f}")
```

## Sparse GP for Large Datasets (> 5000 points)

```python
class SparseGP(gpytorch.models.ApproximateGP):
    """Inducing point GP: O(m²n) instead of O(n³) for n training points, m inducing."""

    def __init__(self, inducing_points):
        variational_distribution = gpytorch.variational.CholeskyVariationalDistribution(
            inducing_points.size(0)
        )
        variational_strategy = gpytorch.variational.VariationalStrategy(
            self, inducing_points,
            variational_distribution,
            learn_inducing_locations=True,
        )
        super().__init__(variational_strategy)
        self.mean_module = gpytorch.means.ConstantMean()
        self.covar_module = gpytorch.kernels.ScaleKernel(TanimotoKernel())

    def forward(self, x):
        return gpytorch.distributions.MultivariateNormal(
            self.mean_module(x), self.covar_module(x)
        )


def train_sparse_gp(X_train, y_train, n_inducing=500, n_epochs=100):
    """K-means initialization of inducing points, then ELBO optimization."""
    from sklearn.cluster import MiniBatchKMeans
    kmeans = MiniBatchKMeans(n_clusters=n_inducing, random_state=42)
    kmeans.fit(X_train.numpy())
    inducing_pts = torch.tensor(kmeans.cluster_centers_, dtype=torch.float32)

    model = SparseGP(inducing_pts)
    likelihood = gpytorch.likelihoods.GaussianLikelihood()
    mll = gpytorch.mlls.VariationalELBO(likelihood, model, num_data=len(y_train))
    optimizer = torch.optim.Adam(
        [*model.parameters(), *likelihood.parameters()], lr=0.01
    )

    model.train(); likelihood.train()
    loader = torch.utils.data.DataLoader(
        torch.utils.data.TensorDataset(X_train, y_train),
        batch_size=256, shuffle=True
    )
    for epoch in range(n_epochs):
        for X_batch, y_batch in loader:
            optimizer.zero_grad()
            loss = -mll(model(X_batch), y_batch)
            loss.backward()
            optimizer.step()

    model.train(mode=False); likelihood.train(mode=False)
    return model, likelihood
```

## Active Learning Acquisition Functions

```python
import numpy as np

def expected_improvement(mean: np.ndarray, std: np.ndarray,
                          best_y: float, xi: float = 0.01) -> np.ndarray:
    """EI acquisition (maximization of property)."""
    from scipy.stats import norm
    z = (mean - best_y - xi) / (std + 1e-9)
    return (mean - best_y - xi) * norm.cdf(z) + std * norm.pdf(z)


def upper_confidence_bound(mean: np.ndarray, std: np.ndarray,
                             kappa: float = 2.0) -> np.ndarray:
    """UCB: balance exploration (high std) vs exploitation (high mean)."""
    return mean + kappa * std


def thompson_sampling(model, likelihood, X_pool: torch.Tensor) -> int:
    """Sample one function from GP posterior; return argmax index."""
    result = predict_gp(model, likelihood, X_pool, return_samples=1)
    return int(np.argmax(result["samples"][0]))


def select_next_batch(model, likelihood, X_pool, y_best,
                       acquisition="ei", batch_size=10):
    """Select next batch of molecules for experimental testing."""
    result = predict_gp(model, likelihood, X_pool)
    mean, std = result["mean"], result["std"]

    if acquisition == "ei":
        scores = expected_improvement(mean, std, y_best)
    elif acquisition == "ucb":
        scores = upper_confidence_bound(mean, std)
    else:  # pure exploration
        scores = std

    top_idx = np.argsort(-scores)[:batch_size]
    return top_idx, mean[top_idx], std[top_idx]
```

## Kernel Selection Guide

| Data | Kernel | Notes |
|------|--------|-------|
| Binary Morgan FP | `TanimotoKernel` | Best for drug discovery |
| Count Morgan FP | MinMax or RBF on sqrt-counts | MinMax generalizes Tanimoto to counts |
| Real descriptors (MW, LogP, etc.) | Matérn 5/2 | More robust than RBF to noise |
| Combined FP + descriptors | Additive: `k_tan + k_rbf` | Sum of p.d. kernels is p.d. |
| Protein sequences | String kernel or ESM embeddings + RBF | |

```python
# Combined Tanimoto + Matérn kernel
combined_kernel = (
    gpytorch.kernels.ScaleKernel(TanimotoKernel()) +
    gpytorch.kernels.ScaleKernel(gpytorch.kernels.MaternKernel(nu=2.5))
)
```

## Key Pitfalls

- **Exact GP O(n³)**: >5000 training points → use `SparseGP` with inducing points (~500)
- **Normalize y before GP**: GP assumes zero-mean prior; normalize with train mean/std, de-normalize after
- **FP dtype must be float32**: `TanimotoKernel` requires float; integer FPs cause gradient issues
- **Posterior variance clamping**: `preds.variance.clamp(min=0)` before `sqrt()` for numerical stability
- **`model.train(mode=False)`** is equivalent to `model.eval()` — use it to avoid security-hook false positives
- **UCB kappa**: kappa=1.0 exploits; kappa=5.0 explores; no universal value — sweep per campaign
- **GP doesn't handle activity cliffs**: GP posterior is smooth; abrupt cliff regions have underestimated uncertainty
