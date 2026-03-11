# Uncertainty Theory — Epistemic vs Aleatoric, GP Kernels, Conformal Coverage, AD

## Epistemic vs Aleatoric Uncertainty

```
Total uncertainty = Epistemic + Aleatoric

Epistemic (model uncertainty):
  - Caused by limited training data
  - Reducible: more data → less epistemic uncertainty
  - High far from training data (interpolation vs extrapolation)
  - Captured by: GP posterior variance, ensemble disagreement, MC dropout std

Aleatoric (data uncertainty):
  - Intrinsic noise in the measurement (assay variability, experimental error)
  - Irreducible: more data does NOT reduce it
  - Constant or input-dependent (heteroscedastic)
  - Captured by: learned noise σ in GP, NLL-minimizing heteroscedastic heads

For drug discovery:
  - Typical assay σ: IC50 ≈ 0.3-0.5 log units (replicate measurements)
  - Epistemic dominates for novel scaffolds
  - Aleatoric dominates for well-characterized chemical series
```

## Total Predictive Uncertainty Decomposition

For a Bayesian model with parameters θ:

```
E[y|x] = ∫ f(x;θ) p(θ|data) dθ           (predictive mean)

Var[y|x] = E_θ[Var[y|x,θ]] + Var_θ[E[y|x,θ]]
           ↑ aleatoric             ↑ epistemic

For ensemble of M models:
  Epistemic σ² = (1/M) Σ (f_m(x) - ȳ)²
  Aleatoric σ² = (1/M) Σ σ²_noise,m(x)      (if models output noise σ)
```

## Gaussian Process Foundations

A GP defines a distribution over functions: f ~ GP(μ, k), where k is a kernel function.

**Posterior prediction** at test point x*:

```
Posterior mean:     μ*(x*) = k(x*, X) [K(X,X) + σ²I]⁻¹ y
Posterior variance: σ²*(x*) = k(x*,x*) - k(x*,X) [K(X,X) + σ²I]⁻¹ k(X,x*)

Where:
  K(X,X)_ij = k(x_i, x_j)  [n×n kernel matrix]
  k(x*, X) = [k(x*,x_1), ..., k(x*,x_n)]  [1×n cross-kernel]
  σ² = observation noise variance (learned or fixed)
```

**Key property**: σ²*(x*) increases as x* moves away from training data → naturally captures epistemic uncertainty.

## Kernel Validity for Molecular Fingerprints

A kernel k(x,y) must be **positive definite** (all Gram matrices have non-negative eigenvalues).

| Kernel | Valid for fingerprints? | Notes |
|--------|------------------------|-------|
| **Tanimoto** `T(a,b) = |a∩b|/|a∪b|` | Yes (binary FP) | Haussler 1999; valid p.d. kernel |
| **MinMax** `min(a,b)/max(a,b)` | Yes (count FP) | Generalization of Tanimoto |
| **RBF on Euclidean** | Yes (any real FP) | Works on Morgan counts, not bits |
| Cosine similarity | NOT p.d. in general | Can cause negative definiteness |
| Dice similarity | Yes (binary FP) | `2|a∩b|/(|a|+|b|)` |

```python
# Verify positive definiteness:
import numpy as np
from scipy.linalg import eigvalsh

def is_positive_definite(K):
    eigenvalues = eigvalsh(K)
    return bool(np.all(eigenvalues >= -1e-8))
```

## Conformal Prediction — Coverage Guarantee

**Split conformal prediction** (Papadopoulos et al. 2002):

```
Algorithm:
1. Split data: train set D_train, calibration set D_cal (typically 20-30%)
2. Train model f on D_train
3. Compute nonconformity scores on calibration set:
   s_i = |y_i - f(x_i)|   [for regression]
4. Compute quantile: q̂ = quantile(s_1,...,s_n, (1-α)(1 + 1/n))
5. Prediction interval for new x*: Ĉ(x*) = [f(x*) - q̂, f(x*) + q̂]

Coverage guarantee (finite-sample, distribution-free):
  P(y ∈ Ĉ(x)) ≥ 1 - α

Required assumption: exchangeability of (x_i, y_i) — approximately holds
  if calibration and test molecules are from the same distribution.
```

**Key advantage**: no distributional assumptions on the model or data. Works with ANY base model (RF, GNN, GP, linear).

**Width tradeoff**: smaller α (higher confidence) → wider intervals. Interval width reflects model uncertainty.

## Conformal vs Bayesian Credible Interval

| Property | Conformal | Bayesian |
|----------|-----------|---------|
| Coverage guarantee | ✓ Frequentist, finite-sample | ✗ Asymptotic only |
| Model assumptions | None | Requires correct prior/likelihood |
| Computational cost | Low (quantile on calibration) | High (MCMC or variational) |
| Conditional coverage | Not guaranteed per-molecule | Not guaranteed |
| Integration with any model | Yes | Requires Bayesian model |
| Uncertainty reflects data structure | Via nonconformity score | Via posterior |

**Best practice**: use conformal for guaranteed coverage + Bayesian/GP for uncertainty structure (which molecules are uncertain and why).

## Applicability Domain (AD) Theory

The AD defines the chemical space where the model's predictions are reliable.

**Formal definition**: a molecule x* is within AD if its uncertainty (or distance to training) is below a threshold τ:

```
AD(x*) = True  iff  d(x*, X_train) ≤ τ

Where d can be:
- Tanimoto distance: d = 1 - max_{x_i ∈ train} Tanimoto(x*, x_i)
- kNN distance: d = mean of k nearest neighbors' Tanimoto distances
- Leverage (hat value): h = x*(XᵀX)⁻¹x*ᵀ
- Mahalanobis distance: D_M = sqrt((x* - μ)ᵀ Σ⁻¹ (x* - μ))
- GP posterior std: σ*(x*) directly
```

**Threshold setting**: set τ so that X% of training set is within AD (e.g., 95th percentile of training distances).

## Why AD Fails Silently

A molecule can be:
1. **Within AD** (small distance to train) but **wrong prediction** → model memorized noisy data
2. **Outside AD** (large distance) but **correct prediction** → model generalizes well to that scaffold

AD is necessary but not sufficient for reliability. Always combine AD with:
- Calibration check (is stated CI actually achieving correct coverage?)
- Temporal/scaffold split validation (does coverage hold on new chemotypes?)

## Statistical Calibration

**Calibration**: is the model's stated confidence correct?

```
For α ∈ {0.05, 0.10, ..., 0.95}:
  Compute empirical coverage(α) = mean(y ∈ CI_α)

Perfect calibration: empirical coverage ≈ 1 - α for all α

Calibration plot (reliability diagram):
  x-axis: stated confidence (1-α)
  y-axis: empirical coverage
  Ideal: diagonal line
  Overconfident: below diagonal (actual coverage < stated)
  Underconfident: above diagonal (intervals too wide)
```

**Expected Calibration Error (ECE)**:
```
ECE = (1/M) Σ |coverage(α_m) - (1 - α_m)|
```
