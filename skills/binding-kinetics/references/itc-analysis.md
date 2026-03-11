# ITC Thermogram Analysis

## ITC Basics

Isothermal Titration Calorimetry measures heat released/absorbed per injection.
Directly gives: KD, ΔH, n (stoichiometry).
ΔG and ΔS derived: ΔG = RT ln(KD), ΔS = (ΔH - ΔG) / T.

**Instruments**: MicroCal PEAQ-ITC, TA Instruments Nano ITC, MicroCal iTC200.

**Typical setup**:
- Cell: protein (5–50 µM)
- Syringe: ligand (50–500 µM, 10× cell concentration)
- 19–25 injections, 2–4 µL each
- Temperature: 25°C (standard), but measure at multiple T for ΔCp

## Loading ITC Data

### MicroCal Origin format (.itc file)
```python
import numpy as np
import pandas as pd
from pathlib import Path

def load_microcal_itc(filepath):
    """
    Parse MicroCal .itc raw data file.
    Returns: injections (µL), heats (µcal), time (s), thermogram.
    """
    lines = Path(filepath).read_text().split('\n')

    # Header parsing (MicroCal format)
    params = {}
    data_start = 0
    for i, line in enumerate(lines):
        if line.startswith('@'):
            key, val = line[1:].split('=', 1)
            params[key.strip()] = val.strip()
        elif line.strip() == '$':
            data_start = i + 1
            break

    # Parse data block (time, power µcal/s)
    data = []
    for line in lines[data_start:]:
        parts = line.strip().split()
        if len(parts) >= 2:
            try:
                data.append([float(p) for p in parts[:2]])
            except ValueError:
                continue

    df = pd.DataFrame(data, columns=['time_s', 'power_ucal_s'])
    return df, params

def load_nitpic_integrated(csv_path):
    """
    Load NITPIC-integrated heats (recommended over manual Origin integration).
    Columns: injection_vol, heat_ucal, [concentration_ratio]
    """
    return pd.read_csv(csv_path)
```

## Fitting: One-Site Binding Model

```python
from scipy.optimize import curve_fit
import numpy as np

def one_site_itc_model(molar_ratio, n, KA, dH, dH_dilution,
                        M_total, L_total):
    """
    One-site ITC binding model.
    molar_ratio: [L]/[M] cumulative molar ratio at each injection
    n: stoichiometry (binding sites per macromolecule)
    KA: association constant (M⁻¹)
    dH: enthalpy of binding (kcal/mol)
    dH_dilution: heat of dilution (constant offset)
    M_total: macromolecule concentration in cell (M)
    L_total: cumulative ligand concentration at each injection (M)
    """
    # Fraction of sites occupied (Wiseman isotherm)
    # Θ = fraction bound, derived from quadratic binding equation
    c = n * KA * M_total  # Wiseman c-value (should be 1–1000)
    x = molar_ratio / n   # normalized ligand ratio

    # Quadratic solution for fraction bound
    term = 1 + x + 1/(n * KA * M_total)
    theta = (term - np.sqrt(term**2 - 4*x)) / 2

    dQ_dtheta = n * M_total * dH   # heat per injection (differential)
    # Numerical derivative of theta w.r.t injection
    dQ = np.gradient(theta) * n * M_total * dH

    return dQ + dH_dilution

def fit_itc_wiseman(molar_ratios, heats_per_injection,
                     M_conc, L_conc_per_injection,
                     p0=None):
    """
    Fit integrated ITC heats using Wiseman isotherm.

    Parameters:
        molar_ratios: cumulative [L]/[M] at each injection midpoint
        heats_per_injection: µcal per injection (from integration)
        M_conc: macromolecule concentration [M]
        p0: initial guess [n, KD_M, dH_kcal_mol]
    """
    R_gas = 1.987e-3  # kcal/(mol·K)
    T = 298.15

    def model(x, n, KD, dH, offset):
        KA = 1.0 / KD
        # Wiseman isotherm: cumulative Q at molar ratio x
        r = x / n
        c_param = n * KA * M_conc
        term = 1 + r + 1/c_param
        theta = (term - np.sqrt(np.maximum(term**2 - 4*r, 0))) / 2
        Q_total = n * M_conc * dH * theta   # kcal/L * vol_cell
        return np.gradient(Q_total, x) + offset

    if p0 is None:
        p0 = [1.0, 1e-7, -10.0, 0.0]  # n, KD, dH, offset

    popt, pcov = curve_fit(
        model, molar_ratios, heats_per_injection,
        p0=p0,
        bounds=([0.1, 1e-12, -100, -10], [4.0, 1e-3, 100, 10]),
        maxfev=50000
    )
    perr = np.sqrt(np.diag(pcov))

    n, KD, dH, offset = popt
    dG = R_gas * T * np.log(KD)   # kcal/mol (negative = favorable)
    dS = (dH - dG) / T            # kcal/(mol·K)

    return {
        'n': n, 'n_err': perr[0],
        'KD': KD, 'KD_err': perr[1],
        'KD_nM': KD * 1e9,
        'pKD': -np.log10(KD),
        'dH': dH, 'dH_err': perr[2],       # kcal/mol
        'dG': dG,                           # kcal/mol
        'TdS': T * dS,                      # kcal/mol
        'dS': dS,                           # kcal/(mol·K)
        'offset': offset,
        'c_wiseman': n / KD / M_conc,       # should be 1–1000
    }
```

## Wiseman c-Value

```
c = n × KA × [M] = n / (KD × [M_cell])
```
- c < 1: too tight → sigmoidal curve too sharp, KD poorly defined
- 1 ≤ c ≤ 1000: optimal, good sigmoidal
- c > 1000: too tight → need displacement ITC or lower [M]

```python
def wiseman_c(n, KD, M_conc):
    return n / (KD * M_conc)

# To achieve c=10–100:
# KD = 1 nM → [M] = n / (KD × c) = 1 / (1e-9 × 100) = 10 µM
# KD = 1 µM → [M] = 10 nM (very dilute — noisy)
```

## Thermodynamic Analysis

```python
def thermodynamic_profile(dH, dG, T=298.15):
    """Decompose ΔG into ΔH and -TΔS components."""
    TdS = dH - dG
    dS = TdS / T

    profile = {
        'dG': dG, 'dH': dH, 'dS': dS, 'TdS': TdS, 'T': T,
        'frac_enthalpic': dH / dG if dG != 0 else None,
        'signature': (
            'enthalpy-driven' if dH < 0 and abs(dH) > abs(TdS)
            else 'entropy-driven' if TdS > 0 and abs(TdS) > abs(dH)
            else 'both-favorable' if dH < 0 and TdS < 0
            else 'entropically_compensated'
        )
    }
    return profile

def van_t_hoff_analysis(temps_K, KD_values):
    """
    van't Hoff: ln(KA) vs 1/T → slope = -ΔH/R, intercept = ΔS/R
    temps_K: temperatures in Kelvin
    KD_values: KD at each temperature [M]
    """
    R = 1.987e-3  # kcal/(mol·K)
    KA_values = 1.0 / np.array(KD_values)
    inv_T = 1.0 / np.array(temps_K)
    ln_KA = np.log(KA_values)

    slope, intercept = np.polyfit(inv_T, ln_KA, 1)
    dH = -slope * R
    dS = intercept * R
    dCp = None  # requires multiple temperature points for curvature

    return {'dH_vH': dH, 'dS_vH': dS, 'slope': slope, 'intercept': intercept}
```

## Heat Capacity Change ΔCp

```python
def delta_Cp(temps_K, dH_values):
    """
    ΔCp = dΔH/dT (slope of ΔH vs T).
    Negative ΔCp: hydrophobic burial (typical for protein-ligand).
    Positive ΔCp: exposure of nonpolar surface.
    """
    dCp, intercept = np.polyfit(temps_K, dH_values, 1)
    return dCp   # kcal/(mol·K)

# Typical values:
# Hydrophobic burial: ΔCp = -0.2 to -1.0 kcal/(mol·K)
# Polar contacts only: ΔCp ≈ 0
```

## ITC Visualization

```python
import matplotlib.pyplot as plt

def plot_itc(molar_ratios, heats_per_injection, fit_result=None):
    """
    Standard ITC binding isotherm plot.
    Top panel: raw thermogram (optional).
    Bottom panel: integrated heats + fit.
    """
    fig, ax = plt.subplots(figsize=(8, 6))

    ax.scatter(molar_ratios, heats_per_injection,
               color='black', s=50, zorder=3, label='Data')

    if fit_result:
        x_fit = np.linspace(0, max(molar_ratios), 200)
        # Reconstruct model curve
        n, KD = fit_result['n'], fit_result['KD']
        dH = fit_result['dH']
        KA = 1.0 / KD
        # Wiseman: cumulative Q
        M_conc = 1e-5  # placeholder — recompute with actual [M]
        ax.plot(x_fit, x_fit * 0 + dH * 0.5, '--', alpha=0.3)  # approximate

        textstr = (f"n = {n:.2f} ± {fit_result.get('n_err', 0):.2f}\n"
                   f"KD = {fit_result['KD_nM']:.1f} ± "
                   f"{fit_result.get('KD_err', 0)*1e9:.1f} nM\n"
                   f"ΔH = {dH:.1f} kcal/mol\n"
                   f"ΔG = {fit_result['dG']:.1f} kcal/mol\n"
                   f"-TΔS = {-fit_result['TdS']:.1f} kcal/mol")
        ax.text(0.65, 0.95, textstr, transform=ax.transAxes, fontsize=9,
                va='top', bbox=dict(boxstyle='round', alpha=0.8))

    ax.axhline(y=0, color='grey', linestyle='-', alpha=0.3)
    ax.set_xlabel('Molar Ratio [L]/[M]', fontsize=12)
    ax.set_ylabel('Heat per injection (µcal)', fontsize=12)
    ax.set_title('ITC Binding Isotherm', fontsize=13)
    return fig

def plot_thermodynamic_bars(compounds, dG_list, dH_list, TdS_list):
    """Stacked bar chart: ΔG decomposition into ΔH and -TΔS."""
    fig, ax = plt.subplots(figsize=(10, 5))
    x = np.arange(len(compounds))
    ax.bar(x, dH_list, label='ΔH', color='#4488cc', alpha=0.8)
    ax.bar(x, [-v for v in TdS_list], bottom=dH_list,
           label='-TΔS', color='#cc4444', alpha=0.8)
    ax.plot(x, dG_list, 'ko-', label='ΔG', zorder=3)
    ax.axhline(y=0, color='black', linewidth=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(compounds, rotation=45, ha='right')
    ax.set_ylabel('kcal/mol', fontsize=11)
    ax.set_title('Thermodynamic Signature', fontsize=12)
    ax.legend()
    return fig
```

## Common ITC Problems

| Problem | Likely cause | Fix |
|---------|-------------|-----|
| Flat isotherm | c < 1 (too weak) or c > 1000 (too tight) | Adjust [M], use displacement ITC |
| Noisy baseline | Protein aggregation, air bubbles | Spin protein, degas both solutions |
| Mismatch n ≠ 1 | Incomplete active fraction | Correct [M] by active fraction |
| Endothermic ΔH ≈ 0 | Entropy-driven binding | Normal; consider DSF/SPR to confirm KD |
| Biphasic curve | Two binding sites | Fit two-site model |
