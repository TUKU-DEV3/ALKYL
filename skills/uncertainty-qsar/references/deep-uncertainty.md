# Deep Uncertainty — MC Dropout, Deep Ensembles, Evidential Regression, Laplace

## MC Dropout (Monte Carlo Dropout)

Dropout at inference time approximates Bayesian posterior. Run T forward passes with dropout active; variance across passes estimates epistemic uncertainty.

```python
import torch
import torch.nn as nn
import numpy as np

class DropoutQSAR(nn.Module):
    """QSAR network with MC Dropout for uncertainty estimation."""

    def __init__(self, input_dim: int, hidden_dims=(512, 256, 128),
                 dropout_rate: float = 0.2):
        super().__init__()
        layers = []
        prev_dim = input_dim
        for h in hidden_dims:
            layers += [
                nn.Linear(prev_dim, h),
                nn.ReLU(),
                nn.Dropout(dropout_rate),   # dropout at each layer
            ]
            prev_dim = h
        layers.append(nn.Linear(prev_dim, 1))
        self.net = nn.Sequential(*layers)
        self.dropout_rate = dropout_rate

    def forward(self, x):
        return self.net(x).squeeze(-1)

    def enable_dropout(self):
        """Re-enable dropout for MC inference (dropout disabled by .train(mode=False))."""
        for m in self.modules():
            if isinstance(m, nn.Dropout):
                m.training = True   # force dropout active regardless of model mode

    @torch.no_grad()
    def mc_predict(self, x: torch.Tensor, n_samples: int = 50) -> dict:
        """
        MC Dropout inference: T stochastic forward passes.
        Dropout is forced active even in non-training mode.
        """
        self.train(mode=False)   # disable batch norm etc.
        self.enable_dropout()    # but keep dropout active

        preds = torch.stack([self(x) for _ in range(n_samples)])  # (T, n)

        mean = preds.mean(0).numpy()
        std  = preds.std(0).numpy()
        return {
            "mean": mean,
            "std": std,                         # epistemic uncertainty
            "lower": mean - 1.96 * std,
            "upper": mean + 1.96 * std,
            "all_samples": preds.numpy(),        # (T, n)
        }


def train_dropout_model(model, X_train, y_train, n_epochs=200, lr=1e-3):
    optimizer = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=1e-4)
    criterion = nn.MSELoss()
    model.train()

    for epoch in range(n_epochs):
        optimizer.zero_grad()
        pred = model(X_train)
        loss = criterion(pred, y_train)
        loss.backward()
        optimizer.step()

        if (epoch + 1) % 50 == 0:
            print(f"Epoch {epoch+1}/{n_epochs}  MSE={loss.item():.4f}")
```

**Number of MC samples**: T=20 is fast; T=50 typical; T=100 for high-stakes decisions.

**Caveats**: MC dropout is an approximation of a Bayesian NN, not exact. Uncertainty estimates can be underconfident. Dropout rate controls epistemic uncertainty scale: too low → narrow CI; too high → noisy predictions.

## Deep Ensembles

Train M independent models from different random initializations. Disagreement across models = epistemic uncertainty. State of the art for deep uncertainty.

```python
from typing import List

class EnsembleQSAR:
    """
    Deep ensemble of M independent QSAR models.
    Each model trained from different random seed.
    """

    def __init__(self, input_dim: int, n_members: int = 5,
                 hidden_dims=(512, 256), dropout_rate: float = 0.1):
        self.n_members = n_members
        self.models = [
            DropoutQSAR(input_dim, hidden_dims, dropout_rate)
            for _ in range(n_members)
        ]

    def fit(self, X_train: torch.Tensor, y_train: torch.Tensor,
            n_epochs: int = 200, lr: float = 1e-3):
        for i, model in enumerate(self.models):
            print(f"Training ensemble member {i+1}/{self.n_members}")
            # Different random seed per member
            torch.manual_seed(i * 42)
            train_dropout_model(model, X_train, y_train, n_epochs, lr)

    @torch.no_grad()
    def predict(self, X_test: torch.Tensor) -> dict:
        """
        Ensemble prediction: mean and variance across M members.
        """
        preds = []
        for model in self.models:
            model.train(mode=False)
            preds.append(model(X_test).numpy())

        preds = np.array(preds)          # (M, n_test)
        mean  = preds.mean(axis=0)       # (n_test,)
        std   = preds.std(axis=0)        # (n_test,) — epistemic uncertainty

        return {
            "mean": mean, "std": std,
            "lower": mean - 1.96 * std,
            "upper": mean + 1.96 * std,
            "all_preds": preds,          # (M, n_test)
        }

    def predict_with_aleatoric(self, X_test: torch.Tensor) -> dict:
        """
        When models output (mean, log_var) heteroscedastic heads.
        Separates aleatoric from epistemic.
        """
        all_means, all_log_vars = [], []
        for model in self.models:
            model.train(mode=False)
            out = model(X_test)           # (n, 2) — [mean, log_var]
            all_means.append(out[:, 0].numpy())
            all_log_vars.append(out[:, 1].numpy())

        means = np.array(all_means)       # (M, n)
        pred_mean = means.mean(0)

        epistemic_var = means.var(0)                              # disagreement
        aleatoric_var = np.exp(np.array(all_log_vars)).mean(0)   # mean noise

        total_std = np.sqrt(epistemic_var + aleatoric_var)
        return {
            "mean": pred_mean,
            "std_total": total_std,
            "std_epistemic": np.sqrt(epistemic_var),
            "std_aleatoric": np.sqrt(aleatoric_var),
        }
```

## Heteroscedastic Head (Learned Aleatoric Uncertainty)

Instead of fixed MSE loss, predict both mean and variance:

```python
class HeteroscedasticQSAR(nn.Module):
    """Outputs (mean, log_var) pair — aleatoric uncertainty."""

    def __init__(self, input_dim: int, hidden_dims=(512, 256)):
        super().__init__()
        layers = []
        prev = input_dim
        for h in hidden_dims:
            layers += [nn.Linear(prev, h), nn.ReLU(), nn.Dropout(0.1)]
            prev = h
        self.backbone = nn.Sequential(*layers)
        self.mean_head  = nn.Linear(prev, 1)
        self.logvar_head = nn.Linear(prev, 1)

    def forward(self, x):
        feat = self.backbone(x)
        mean   = self.mean_head(feat).squeeze(-1)
        logvar = self.logvar_head(feat).squeeze(-1)
        return mean, logvar


def nll_loss(mean: torch.Tensor, logvar: torch.Tensor,
             target: torch.Tensor) -> torch.Tensor:
    """Gaussian negative log-likelihood loss for heteroscedastic regression."""
    # NLL = 0.5 * (log σ² + (y - μ)² / σ²)
    return 0.5 * (logvar + ((target - mean) ** 2) / logvar.exp()).mean()


def train_heteroscedastic(model, X_train, y_train, n_epochs=200, lr=1e-3):
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    model.train()
    for epoch in range(n_epochs):
        optimizer.zero_grad()
        mean, logvar = model(X_train)
        loss = nll_loss(mean, logvar, y_train)
        loss.backward()
        optimizer.step()
```

## Laplace Approximation (laplace-torch)

Post-hoc uncertainty for any trained PyTorch model — approximates posterior as Gaussian around MAP estimate.

```python
# pip install laplace-torch
from laplace import Laplace
import torch

# 1. Train a standard deterministic model first
model = DropoutQSAR(input_dim=2048, hidden_dims=(512, 256), dropout_rate=0.0)
train_dropout_model(model, X_train, y_train)

# 2. Wrap with Laplace approximation
la = Laplace(
    model,
    likelihood="regression",       # or "classification"
    subset_of_weights="last_layer", # "all" is expensive; "last_layer" is fast and good
    hessian_structure="kron",       # Kronecker-factored; "full" more accurate but O(n²p²)
)

# 3. Fit Laplace (compute Hessian on training data)
train_loader = torch.utils.data.DataLoader(
    torch.utils.data.TensorDataset(X_train, y_train.unsqueeze(-1)),
    batch_size=64, shuffle=False
)
la.fit(train_loader)

# 4. Optimize prior precision (marginal likelihood)
la.optimize_prior_precision(method="marglik")

# 5. Predict with uncertainty
X_test_tensor = torch.tensor(X_test, dtype=torch.float32)
mean, var = la(X_test_tensor)   # predictive mean and variance
std = var.sqrt().detach().numpy()
mean = mean.detach().numpy().squeeze()
print(f"Laplace posterior std: {std.mean():.4f}")
```

**`subset_of_weights="last_layer"`** is the recommended default:
- Fast to compute (no full Hessian over all parameters)
- Good empirical performance (last-layer BNN)
- Equivalent to fitting a GP with learned kernel (Laplace-GNN connection)

## Comparison Table

| Method | Epistemic | Aleatoric | Cost | Quality |
|--------|----------|-----------|------|---------|
| MC Dropout (T=50) | ✓ (approx) | ✗ | 50× inference | Moderate |
| Deep Ensemble (M=5) | ✓ (good) | ✗ | 5× train + inference | Best practical |
| Heterosc. head | ✗ | ✓ | 1× (2-output head) | Good aleatoric |
| Ensemble + Heterosc | ✓ + ✓ | ✓ | 5× train | Gold standard |
| Laplace (last-layer) | ✓ (post-hoc) | ✓ | 1× train + Hessian | Excellent |
| Exact GP | ✓ (exact) | ✓ | O(n³) | Best (small n) |
| Conformal (MAPIE) | Coverage ✓ | Coverage ✓ | Negligible | Guaranteed CI |

## Practical Recommendation

```
Small dataset (<2000):     GP with Tanimoto kernel
Medium dataset (<20k):     Deep ensemble (M=5) + Laplace last-layer
Large dataset (>20k):      Sparse GP or MC Dropout with conformal on top
Production deployment:     Conformal wrapping any base model (coverage guarantee)
Active learning:           GP (exact EI/UCB) or ensemble (Thompson sampling)
```

## Key Pitfalls

- **MC Dropout n_samples tuning**: T < 20 gives noisy std estimates; T > 100 rarely improves
- **Deep ensemble diversity**: must use different random seeds AND different data shuffles; identical architecture is fine
- **Laplace `hessian_structure="full"`** is O(p²) memory — infeasible for models >10k parameters; use `"kron"` or `"diag"`
- **Heteroscedastic training instability**: logvar can go very negative (→ very small noise, → numerically unstable NLL); add `logvar.clamp(-6, 6)` in loss
- **Conformal + deep model**: always add conformal calibration on top of any deep UQ method — the internal UQ may not have correct marginal coverage
