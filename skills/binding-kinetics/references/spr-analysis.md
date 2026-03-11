# SPR Sensorgram Analysis

## SPR Basics

Surface Plasmon Resonance measures mass change on a sensor chip surface.
Signal unit: RU (Response Units) — 1 RU ≈ 1 pg/mm²

**Typical experiment**:
1. Immobilize protein on chip (amine coupling or capture)
2. Flow analyte (drug) at multiple concentrations
3. Record: association (drug on) → dissociation (buffer only) phases
4. Fit to kinetic model → kon, koff, KD

## Theoretical Rmax

```python
def Rmax_theoretical(MW_analyte, MW_ligand, Rl, valency=1):
    """
    Expected maximum response.
    Rl: immobilized ligand response (RU)
    valency: binding sites per ligand molecule (usually 1)
    """
    return (MW_analyte / MW_ligand) * Rl * valency

# If measured Rmax << Rtheor: partial activity (surface heterogeneity)
# If measured Rmax >> Rtheor: non-specific binding
# Acceptable: 50–80% of Rtheor
```

## Loading SPR Data (Biacore T200 / BIAevaluation format)

```python
import pandas as pd
import numpy as np

def load_biacore_csv(filepath):
    """
    Load Biacore exported sensorgram data.
    Format: Time(s), Response(RU) columns per concentration.
    """
    df = pd.read_csv(filepath, skiprows=1)  # skip header line
    # Typical columns: Time, Conc1_RU, Conc2_RU, ...
    # Normalize: subtract reference channel (Fc1)
    return df

def parse_multichannel_spr(filepath, ref_channel='Fc1'):
    """Parse multi-channel SPR with reference subtraction."""
    df = pd.read_csv(filepath, header=0)
    time_col = df.columns[0]
    ref = df[ref_channel] if ref_channel in df.columns else 0
    channels = [c for c in df.columns if c != time_col and c != ref_channel]
    result = pd.DataFrame({'time': df[time_col]})
    for ch in channels:
        result[ch] = df[ch] - ref   # double-reference subtraction
    return result
```

## 1:1 Langmuir Kinetic Model

```
d[RU]/dt = kon * C * (Rmax - RU) - koff * RU
```

Association phase (t₀ to t_off, constant analyte concentration C):
```
RU(t) = Req * (1 - exp(-(kon*C + koff) * t))
Req = Rmax * C / (C + KD)
ks = kon*C + koff    (observed rate)
```

Dissociation phase (buffer, C=0):
```
RU(t) = RU(t_off) * exp(-koff * (t - t_off))
```

```python
from scipy.optimize import curve_fit
import numpy as np

def langmuir_association(t, kon, koff, Rmax, C):
    """1:1 Langmuir association phase."""
    ks = kon * C + koff
    Req = Rmax * C / (C + koff/kon)
    return Req * (1 - np.exp(-ks * t))

def langmuir_dissociation(t, koff, R0, t_off):
    """1:1 Langmuir dissociation phase (starts at t_off with response R0)."""
    return R0 * np.exp(-koff * (t - t_off))

def fit_spr_dissociation(time, response, t_off, p0=None):
    """
    Fit dissociation phase to extract koff.
    More reliable than full kinetic fit for slow binders.
    """
    mask = time >= t_off
    t_diss = time[mask]
    r_diss = response[mask]
    R0_init = r_diss[0]

    if p0 is None:
        p0 = [1e-3, R0_init]  # [koff, R0]

    popt, pcov = curve_fit(
        lambda t, koff, R0: langmuir_dissociation(t, koff, R0, t_off),
        t_diss, r_diss,
        p0=p0,
        bounds=([1e-7, 0], [1.0, R0_init * 2]),
        maxfev=10000
    )
    perr = np.sqrt(np.diag(pcov))
    koff, R0 = popt
    return {
        'koff': koff,
        'koff_err': perr[0],
        'R0': R0,
        't_half': np.log(2) / koff,
        'RT': 1.0 / koff
    }

def fit_spr_full_kinetics(time, response, concentrations,
                           t_assoc, t_dissoc, Rmax_init=None):
    """
    Global fit across multiple concentrations.
    time: 1D array (common time axis)
    response: 2D array (n_conc × n_timepoints)
    concentrations: list of analyte concentrations [M]
    """
    from scipy.optimize import minimize

    def model(params, t, C, t_off):
        kon, koff, Rmax = params
        result = np.zeros_like(t)
        assoc_mask = t <= t_off
        diss_mask = t > t_off

        # Association
        ks = kon * C + koff
        Req = Rmax * C / (C + koff/kon)
        result[assoc_mask] = Req * (1 - np.exp(-ks * t[assoc_mask]))

        # Dissociation
        R_at_toff = Req * (1 - np.exp(-ks * t_off))
        result[diss_mask] = R_at_toff * np.exp(-koff * (t[diss_mask] - t_off))
        return result

    def objective(params):
        if any(p <= 0 for p in params):
            return 1e10
        total_residuals = 0
        for i, C in enumerate(concentrations):
            pred = model(params, time, C, t_dissoc)
            total_residuals += np.sum((response[i] - pred)**2)
        return total_residuals

    Rmax_guess = Rmax_init or response.max()
    result = minimize(
        objective,
        x0=[1e5, 1e-3, Rmax_guess],
        method='Nelder-Mead',
        options={'maxiter': 10000, 'xatol': 1e-8, 'fatol': 1e-8}
    )

    kon, koff, Rmax = result.x
    return {
        'kon': kon, 'koff': koff,
        'KD': koff / kon,
        'pKD': -np.log10(koff / kon),
        'RT': 1.0 / koff,
        't_half': np.log(2) / koff,
        'Rmax': Rmax,
        'success': result.success,
        'residuals': result.fun
    }
```

## Equilibrium Analysis (Steady-State)

For fast kinetics (equilibrium reached during injection):

```python
def fit_equilibrium_spr(concentrations, responses_eq):
    """
    Fit steady-state Req vs C to extract KD.
    Req = Rmax * C / (KD + C)
    """
    def langmuir_eq(C, Rmax, KD):
        return Rmax * C / (KD + C)

    popt, pcov = curve_fit(
        langmuir_eq,
        concentrations, responses_eq,
        p0=[max(responses_eq), np.median(concentrations)],
        bounds=([0, 1e-12], [max(responses_eq)*2, 1e-3])
    )
    Rmax, KD = popt
    perr = np.sqrt(np.diag(pcov))
    return {
        'KD': KD, 'KD_err': perr[1],
        'pKD': -np.log10(KD),
        'Rmax': Rmax,
        'KD_nM': KD * 1e9
    }
```

## Two-State Binding Model

For induced-fit or conformational change:
```
D + T ⇌ DT* ⇌ DT**
    k1/k-1    k2/k-2
```
Observed as biphasic association and slow dissociation.

```python
def two_state_model(t, k1, k_minus1, k2, k_minus2, Rmax, C, t_off):
    """
    Two-state (induced fit) SPR model.
    More complex fit — use only when 1:1 residuals are systematic.
    """
    # Rate matrix eigenvalue solution
    a = k1 * C + k_minus1 + k2 + k_minus2
    b = k_minus1 * k_minus2 + k1 * C * k_minus2 + k1 * C * k2

    discriminant = a**2 - 4*b
    if discriminant < 0:
        return np.zeros_like(t)

    s1 = (-a + np.sqrt(discriminant)) / 2
    s2 = (-a - np.sqrt(discriminant)) / 2

    # Association phase
    assoc = t <= t_off
    A = Rmax * k1 * C / (s1 - s2)
    R_assoc = A * ((s1 + k_minus1 + k2) * (np.exp(s1*t[assoc]) - 1) / s1 -
                   (s2 + k_minus1 + k2) * (np.exp(s2*t[assoc]) - 1) / s2)

    # Apparent rates (observable)
    KD_app = (k_minus1 * k_minus2) / (k1 * (k2 + k_minus2))
    return R_assoc   # simplified; full model requires numerical ODE
```

## Quality Checks

```python
def spr_quality_checks(fit_result, concentrations, responses_eq=None):
    """Check common SPR data quality issues."""
    issues = []

    kon, koff = fit_result['kon'], fit_result['koff']
    KD = fit_result['KD']

    # Mass transport limitation: kon > 10^6 for slow dissociation check
    if kon > 5e6:
        issues.append("WARNING: High kon — check for mass transport limitation "
                      "(reduce flow rate or immobilization level)")

    # Biphasic dissociation: check with residuals
    if fit_result.get('residuals', 0) > 1.0:
        issues.append("Poor fit quality — consider two-state model or "
                      "rebinding artifact")

    # Rmax consistency
    if responses_eq is not None:
        Req_max = max(responses_eq)
        if fit_result['Rmax'] / Req_max < 0.5:
            issues.append("Rmax << max(Req) — possible surface heterogeneity")

    if not issues:
        issues.append("No major quality issues detected")
    return issues
```

## Visualization

```python
import matplotlib.pyplot as plt

def plot_sensorgrams(time, responses_matrix, concentrations,
                     fit_results=None, t_off=None):
    """Plot overlay of sensorgrams at multiple concentrations."""
    fig, ax = plt.subplots(figsize=(10, 6))

    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(concentrations)))

    for i, (conc, color) in enumerate(zip(concentrations, colors)):
        label = f'{conc*1e9:.1f} nM' if conc < 1e-6 else f'{conc*1e6:.1f} µM'
        ax.plot(time, responses_matrix[i], color=color, lw=1.5, label=label)

    if t_off:
        ax.axvline(x=t_off, color='grey', linestyle='--', alpha=0.5,
                   label='Dissociation start')

    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Response (RU)')
    ax.set_title('SPR Sensorgrams')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    return fig
```
