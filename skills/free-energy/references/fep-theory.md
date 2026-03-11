# FEP Theory — Thermodynamic Cycles, Alchemical Paths, Estimators

## Free Energy Basics

Free energy differences are **state functions** — path-independent. Only endpoints matter for ΔG. However, the *estimator* efficiency depends on how you sample the path.

```
ΔG(A→B) = -kBT ln ⟨exp(-(H_B - H_A)/kBT)⟩_A    [FEP / Zwanzig equation]
```

For large ΔG: direct FEP has catastrophic variance (rare events dominate). Solution: intermediate λ states.

---

## Thermodynamic Cycles

### Relative Binding Free Energy (RBFE)

```
Ligand A  + Protein  →  Complex A
    |                       |
  ΔG_alch(solv)          ΔG_alch(prot)
    |                       |
Ligand B  + Protein  →  Complex B

ΔΔG_bind(A→B) = ΔG_alch(prot) - ΔG_alch(solv)
```
- Mutate A→B *in situ*, in both the protein complex and in solvent (waterbox)
- Much cheaper than computing two absolute ΔG_bind values
- ΔΔG_bind = ΔG_bind(B) - ΔG_bind(A)

### Absolute Binding Free Energy (ABFE)

```
Ligand (bound) → Ligand (unbound, restrained) → Ligand (free in solution)

ΔG_bind = ΔG_restrain + ΔG_decouple(complex) - ΔG_decouple(solvent) + ΔG_standard_state
```

### Solvation Free Energy

```
Ligand (gas phase) → Ligand (water)

ΔG_solv = ΔG_decouple(vacuum) - ΔG_decouple(water)
```

---

## Alchemical Transformation

A **hybrid Hamiltonian** interpolates between state A and state B:

```
H(λ) = (1 - λ) · H_A + λ · H_B       [linear mixing; naive]
```

For nonbonded terms, softcore potentials prevent singularities when atoms appear/disappear:

### Softcore Lennard-Jones (recommended)
```
V_sc(r, λ) = 4ε·λⁿ · [1/(α(1-λ)ᵐ + (r/σ)⁶)² - 1/(α(1-λ)ᵐ + (r/σ)⁶)]
```
- α = 0.5, n = 1, m = 1 (Beutler parameters — standard)
- Avoids r→0 singularity when λ→0 (atom vanishes)

### Lambda Schedule Design
```
# Electrostatics turned off first (before LJ), turned on last
# Avoids charge singularity under LJ softcore

λ_vdw:   [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
λ_elec:  [0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0, 1.0, 1.0]

# Rule: annihilate charges before decoupling LJ
#        restore LJ before restoring charges
```

---

## Free Energy Estimators

### FEP (Zwanzig, 1954)
```
ΔG(A→B) = -kBT ln ⟨exp(-ΔH/kBT)⟩_A
```
- Forward FEP: sample from A, evaluate B
- Reverse FEP: sample from B, evaluate A
- Problem: **exponential averaging** — dominated by rare low-energy configurations
- Use: only for small perturbations (ΔG ≤ 2 kBT)

### TI (Thermodynamic Integration)
```
ΔG(A→B) = ∫₀¹ ⟨∂H/∂λ⟩_λ dλ
```
- Numerically integrate the ensemble average of ∂H/∂λ at each λ point
- Integrand ⟨∂H/∂λ⟩ must be smooth → need well-designed λ schedule
- Numerical integration: Gaussian quadrature (optimal nodes) or trapezoidal rule
- Use: robust but requires smooth integrand; slower convergence than MBAR

### BAR (Bennett Acceptance Ratio, 1976)
```
ΔG(λ_i → λ_j) from pair of simulations using overlap samples
Solved self-consistently: ΔG = argmin_C [⟨f(H_j - H_i - C)⟩_i / ⟨f(H_i - H_j + C)⟩_j]
where f(x) = 1/(1 + exp(x))   [Fermi function]
```
- Optimal for a **pair** of adjacent states
- Requires bidirectional sampling (simulate at λ_i, evaluate at λ_j and vice versa)

### MBAR (Multistate BAR, Shirts & Chodera 2008)
```
ΔG from K states simultaneously:
-ln Z_k = ln Σ_j Σ_n exp(u_k(x_n^j) - f_j) - f_k
```
- **Global** estimator: uses all data from all K states
- Minimum variance for given simulation data
- Provides free energies between **any** pair of states
- Also estimates expectations of observables and their uncertainties
- **Recommended default estimator**

### Overlap Matrix
```
O_ij = ∫ p_i(x) p_j(x) dx / [∫ p_i(x) dx · ∫ p_j(x) dx]
```
- Quantifies phase-space overlap between adjacent λ windows
- O_ij > 0.03 (3%) required for reliable MBAR
- O_ij ≈ 0.2-0.4 is ideal — add more λ windows if lower

---

## Variance and Efficiency

### Exponential Variance (FEP failure mode)
```
Var[ΔG] ∝ exp(σ²_ΔH / (kBT)²)
```
High variance if ΔH distribution is wide → exponential blow-up → need intermediate λ.

### Statistical Error in MBAR
From bootstrap or analytical MBAR covariance:
```
σ(ΔG) ≈ 1/√N_eff   where N_eff = N / (1 + 2·τ_int)
```
τ_int = integrated autocorrelation time (decorrelation time in MD steps).

### Optimal λ Spacing
- Space λ to equalize variance ⟨δH²⟩ across windows
- In practice: denser λ near endpoints (0 and 1) where structure changes most
- Rule of thumb: ⟨∂H/∂λ⟩ should vary smoothly; |ΔH| between adjacent windows ≈ 2-3 kBT

---

## Convergence Criteria

```
1. Free energy plateau: ΔG(t) flat for last 50% of simulation
2. Forward/backward agreement: |ΔG_fwd - ΔG_bwd| < 0.5 kcal/mol
3. Overlap: O_ij > 0.03 for all adjacent pairs
4. Autocorrelation: effective N_eff > 50 uncorrelated samples per window
5. Hysteresis: ΔG(A→B) + ΔG(B→A) ≈ 0 (cycle closure)
```

---

## Common Problems

| Problem | Symptom | Fix |
|---------|---------|-----|
| Poor overlap | O_ij < 0.03 | Add λ windows between offending pair |
| End-point catastrophe | Large spike in ⟨∂H/∂λ⟩ at λ=0 or 1 | Use softcore vdW + decouple charges first |
| Slow convergence | Large τ_int | Longer simulation; replica exchange |
| Charge change | RBFE with net charge change | Use charge-neutralizing co-transformation; alchemical ion |
| Ring breaking/forming | Non-congeneric series | Use core-hopping protocols; single topology careful design |
| Protein conformational change | Hysteresis between FE legs | Enhanced sampling; multiple starting poses |
