# Binding Kinetics Theory

## Core Definitions

```
Drug (D) + Target (T) ⇌ Drug-Target (DT)
              kon →
              ← koff

KD  = koff / kon          [M]       equilibrium dissociation constant
KA  = kon / koff          [M⁻¹]     association constant
RT  = 1 / koff            [s]        residence time
t½  = ln(2) / koff        [s]        complex half-life
```

```python
import numpy as np

def KD_from_rates(kon, koff):
    return koff / kon

def residence_time(koff):
    return 1.0 / koff

def t_half(koff):
    return np.log(2) / koff

def pKD(KD):
    return -np.log10(KD)

# Example: kon=1e6 M⁻¹s⁻¹, koff=1e-3 s⁻¹
# KD = 1e-9 M = 1 nM
# RT = 1000 s ≈ 17 min
# t½ = 693 s ≈ 12 min
```

## Diffusion Limit and kon

Diffusion-limited kon ≈ 10^8–10^9 M⁻¹s⁻¹.
Most drug-target pairs: kon = 10^4–10^7 M⁻¹s⁻¹ (slowed by conformational search).

**kon is usually less variable than koff** across a congeneric series.
→ Optimization primarily acts via koff modulation.

## Residence Time and In Vivo Efficacy

The **kinetic selectivity** concept (Copeland 2006):
- Drug plasma half-life (PK) may be << complex half-life (RT)
- If RT >> PK t½: drug stays on target after plasma clearance → prolonged effect
- If RT << PK t½: target driven by equilibrium (KD matters)

```python
def kinetic_selectivity_ratio(RT_target1, RT_target2):
    """
    Residence time selectivity: even if KD is equal for two targets,
    longer RT on target 1 gives functional selectivity.
    """
    return RT_target1 / RT_target2

# Example: RT_desired = 2h, RT_off-target = 10 min
# Kinetic selectivity = 12× despite KD selectivity = 1×
```

## Two-State Binding Models

### Simple 1:1 (Langmuir)
```
D + T ⇌ DT
```
Rate equations:
```
d[DT]/dt = kon[D][T] - koff[DT]
```

### Induced Fit / Two-State
```
D + T ⇌ DT* ⇌ DT**
       k1/k-1   k2/k-2
```
- DT* = initial encounter complex (fast)
- DT** = conformationally adapted complex (slow, tight)
- Observed koff,app = k-1 * k-2 / (k2 + k-2) — can be very slow

### Conformational Selection
```
T ⇌ T*    (rare open state)
D + T* → DT*
```
- Apparent kon,app depends on pre-equilibrium constant for T*
- Tight binders often exploit rare conformations (cryptic pockets)

## Thermodynamic Signature (ITC)

```
ΔG = ΔH - T·ΔS = RT_gas · ln(KD)    [kcal/mol]
```

| Signature | Interpretation | Risk |
|-----------|---------------|------|
| ΔH-driven (ΔH << 0, ΔS ≈ 0) | Specific H-bonds, electrostatics | Good — enthalpic binding robust |
| ΔS-driven (ΔH ≈ 0, ΔS > 0) | Hydrophobic burial | Promiscuous — ADMET risk |
| Mixed | Both contributions | Ideal |
| ΔH > 0, ΔS >> 0 | Entropy-driven only | Non-specific, aggregation suspect |

**Enthalpy-entropy compensation**: ΔH and -TΔS often correlate across a series.
Large ΔH gain is often cancelled by ΔS loss (rigidification).

```python
def thermodynamic_signature(dH, dS, T=298.15):
    """
    Decompose ΔG into enthalpic and entropic components.
    Returns fraction of ΔG from enthalpy.
    """
    dG = dH - T * dS
    frac_H = dH / dG if dG != 0 else None
    return {'dG': dG, 'dH': dH, 'TdS': T * dS,
            'frac_enthalpic': frac_H,
            'signature': 'enthalpic' if frac_H and frac_H > 0.7 else
                         'entropic' if frac_H and frac_H < 0.3 else 'mixed'}
```

## Kinetic Rate Map

Plot kon vs koff to visualize kinetic optimization space:

```python
import matplotlib.pyplot as plt
import numpy as np

def kinetic_rate_map(data_df):
    """
    Plot log(kon) vs log(koff) with KD iso-contours.
    data_df: columns kon [M⁻¹s⁻¹], koff [s⁻¹], label
    """
    fig, ax = plt.subplots(figsize=(8, 7))

    # KD iso-contours
    koff_range = np.logspace(-6, 0, 200)
    for kd_val, label in [(1e-9, '1 nM'), (1e-8, '10 nM'),
                           (1e-7, '100 nM'), (1e-6, '1 µM')]:
        kon_line = koff_range / kd_val
        ax.plot(np.log10(koff_range), np.log10(kon_line), '--',
                alpha=0.4, color='grey')
        ax.text(np.log10(koff_range[-10]), np.log10(kon_line[-10]),
                f'KD={label}', fontsize=8, color='grey')

    # RT iso-contours (horizontal lines — koff only)
    for rt_val, label in [(60, '1 min'), (3600, '1 h'), (86400, '1 day')]:
        koff_rt = 1.0 / rt_val
        ax.axvline(x=np.log10(koff_rt), linestyle=':', alpha=0.5, color='blue')
        ax.text(np.log10(koff_rt), 4.2, f'RT={label}',
                fontsize=8, color='blue', rotation=90)

    # Data points
    ax.scatter(np.log10(data_df['koff']), np.log10(data_df['kon']),
               s=80, zorder=3)
    if 'label' in data_df.columns:
        for _, row in data_df.iterrows():
            ax.annotate(row['label'],
                        (np.log10(row['koff']), np.log10(row['kon'])),
                        textcoords='offset points', xytext=(5, 5), fontsize=8)

    ax.set_xlabel('log₁₀(koff) [s⁻¹]', fontsize=12)
    ax.set_ylabel('log₁₀(kon) [M⁻¹s⁻¹]', fontsize=12)
    ax.set_title('Kinetic Rate Map', fontsize=13)
    ax.set_xlim(-6, 0)
    ax.set_ylim(3, 8)
    return fig
```

## Selectivity: Equilibrium vs Kinetic

```python
def compare_selectivity(compounds_df):
    """
    Compare equilibrium (KD) vs kinetic (RT) selectivity
    between two targets.
    Requires columns: KD_target1, KD_target2, RT_target1, RT_target2
    """
    df = compounds_df.copy()
    df['KD_selectivity'] = df['KD_target2'] / df['KD_target1']  # >1 = selective for T1
    df['RT_selectivity'] = df['RT_target1'] / df['RT_target2']   # >1 = longer on T1
    df['kinetic_gain'] = df['RT_selectivity'] / df['KD_selectivity']
    # kinetic_gain > 1: kinetics add selectivity beyond equilibrium
    return df[['smiles', 'KD_selectivity', 'RT_selectivity', 'kinetic_gain']]
```

## Kinetic Efficiency Metric (KEI)

Analogous to LE but for koff:
```
KEI = pKD / HAC      (standard LE using KD)
koff_LE = -log10(koff) / HAC   (how slow koff is per heavy atom)
```

```python
def koff_LE(koff, HAC):
    """Ligand efficiency based on koff (residence time efficiency)."""
    return -np.log10(koff) / HAC

# Interpretation: higher = longer residence time per atom added
# koff_LE > 0.25 → good kinetic efficiency
```
