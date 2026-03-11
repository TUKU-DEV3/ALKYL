# Ligand Efficiency Metrics

## Core Metrics

### Ligand Efficiency (LE)
```
LE = ΔG / HAC = 1.37 × pIC50 / HAC    [kcal/mol per heavy atom]
```
```python
def LE(pIC50, HAC):
    return 1.37 * pIC50 / HAC

# Thresholds:
# ≥ 0.40 → excellent fragment
# ≥ 0.30 → acceptable lead
# ≥ 0.25 → marginal
# < 0.25 → poor
```

### Lipophilic Ligand Efficiency (LLE / LipE)
```
LLE = pIC50 - cLogP
```
```python
def LLE(pIC50, cLogP):
    return pIC50 - cLogP

# Thresholds:
# ≥ 5 → excellent
# ≥ 3 → acceptable
# < 2 → lipophilicity-driven → ADMET risk
```

### Ligand-Lipophilicity Efficiency (LLEAT — Keseru & Makara)
Normalizes LLE by HAC to penalize large lipophilic compounds:
```
LLEAT = LLE / HAC^0.3    (or LLE / HAC for simpler version)
```
```python
def LLEAT(pIC50, cLogP, HAC):
    return (pIC50 - cLogP) / (HAC ** 0.3)

# LLEAT ≥ 6 → excellent
```

### Binding Efficiency Index (BEI)
```
BEI = pIC50 / (MW / 1000)    [per kDa]
```
```python
def BEI(pIC50, MW):
    return pIC50 / (MW / 1000)

# BEI ≥ 20 → good
```

### Surface Efficiency Index (SEI)
```
SEI = pIC50 / (PSA / 100)    [per 100 Å²]
```
```python
def SEI(pIC50, PSA):
    return pIC50 / (PSA / 100) if PSA > 0 else None

# SEI ≥ 5 → good
```

### Ligand-Efficiency-Dependent Lipophilicity (LELP)
```
LELP = cLogP / LE
```
Acceptable range: -10 ≤ LELP ≤ 10.
High LELP = potency driven by lipophilicity (bad).

```python
def LELP(cLogP, pIC50, HAC):
    le = LE(pIC50, HAC)
    return cLogP / le if le != 0 else None
```

### Group Efficiency (GE)
Efficiency gained by adding a specific substituent group:
```
GE = ΔpIC50 / Δ HAC
```
```python
def GE(pic50_new, pic50_parent, hac_new, hac_parent):
    delta_pic50 = pic50_new - pic50_parent
    delta_hac = hac_new - hac_parent
    return delta_pic50 / delta_hac if delta_hac > 0 else None

# GE ≥ 0.3 → growing is efficient
# GE < 0.1 → substituent adds atoms without proportional gain
```

## Computing All Metrics

```python
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Crippen
import pandas as pd

def compute_efficiency_metrics(smiles, pIC50):
    """Compute all LE metrics for a compound-activity pair."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    HAC = mol.GetNumHeavyAtoms()
    MW = Descriptors.ExactMolWt(mol)
    cLogP = Crippen.MolLogP(mol)
    PSA = rdMolDescriptors.CalcTPSA(mol)

    le = LE(pIC50, HAC)
    lle = LLE(pIC50, cLogP)

    return {
        'smiles': smiles,
        'pIC50': pIC50,
        'HAC': HAC,
        'MW': round(MW, 2),
        'cLogP': round(cLogP, 2),
        'PSA': round(PSA, 1),
        'LE': round(le, 3),
        'LLE': round(lle, 2),
        'LLEAT': round(LLEAT(pIC50, cLogP, HAC), 2),
        'BEI': round(BEI(pIC50, MW), 2),
        'SEI': round(SEI(pIC50, PSA), 2) if PSA > 0 else None,
        'LELP': round(LELP(cLogP, pIC50, HAC), 1),
        'LE_flag': 'good' if le >= 0.30 else ('marginal' if le >= 0.25 else 'poor'),
        'LLE_flag': 'good' if lle >= 3 else ('marginal' if lle >= 2 else 'poor'),
    }

def metrics_dataframe(smiles_list, pic50_list):
    records = [compute_efficiency_metrics(s, p)
               for s, p in zip(smiles_list, pic50_list)
               if compute_efficiency_metrics(s, p)]
    return pd.DataFrame(records)
```

## Efficiency Evolution Plots

Track how LE and LLE evolve across FBDD optimization rounds:

```python
import matplotlib.pyplot as plt
import numpy as np

def plot_efficiency_evolution(df, round_col='round'):
    """
    Plot LE and LLE across optimization rounds.
    df must have columns: round, LE, LLE, pIC50, HAC.
    """
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    for r in sorted(df[round_col].unique()):
        sub = df[df[round_col] == r]
        axes[0].scatter(sub['HAC'], sub['pIC50'], label=f'Round {r}',
                        alpha=0.7, s=60)
        axes[1].scatter(sub['HAC'], sub['LE'], alpha=0.7, s=60)
        axes[2].scatter(sub['cLogP'], sub['pIC50'], alpha=0.7, s=60)

    # LE iso-lines (constant LE contours on HAC vs pIC50 plot)
    hac_range = np.linspace(5, 35, 100)
    for le_val in [0.25, 0.30, 0.35, 0.40]:
        pic50_line = le_val * hac_range / 1.37
        axes[0].plot(hac_range, pic50_line, '--', alpha=0.4,
                     label=f'LE={le_val}')

    # LLE iso-lines on cLogP vs pIC50 plot
    logp_range = np.linspace(-2, 7, 100)
    for lle_val in [2, 3, 5]:
        axes[2].plot(logp_range, logp_range + lle_val, '--', alpha=0.4,
                     label=f'LLE={lle_val}')

    axes[0].set_xlabel('HAC'); axes[0].set_ylabel('pIC50')
    axes[0].set_title('LE plot (Abad-Zapatero)')
    axes[0].legend(fontsize=7)

    axes[1].set_xlabel('HAC'); axes[1].set_ylabel('LE (kcal/mol/HA)')
    axes[1].axhline(y=0.30, color='r', linestyle='--', label='LE=0.30')
    axes[1].axhline(y=0.25, color='orange', linestyle='--', label='LE=0.25')
    axes[1].set_title('LE evolution')
    axes[1].legend(fontsize=8)

    axes[2].set_xlabel('cLogP'); axes[2].set_ylabel('pIC50')
    axes[2].set_title('LLE plot')
    axes[2].legend(fontsize=7)

    plt.tight_layout()
    return fig
```

## Abad-Zapatero Plot

The canonical FBDD efficiency plot: pIC50 vs HAC with LE iso-contours.

```python
def abad_zapatero_plot(df, color_by='round', annotate=False):
    """
    Abad-Zapatero plot: pIC50 vs HAC with LE contours.
    Essential for tracking FBDD optimization.
    """
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import numpy as np

    fig, ax = plt.subplots(figsize=(9, 7))

    # LE contour lines
    hac_vals = np.linspace(5, 40, 200)
    for le_val, style, alpha in [(0.20, ':', 0.3), (0.25, '--', 0.4),
                                   (0.30, '--', 0.6), (0.35, '-', 0.5),
                                   (0.40, '-', 0.7)]:
        pic50_vals = le_val * hac_vals / 1.37
        ax.plot(hac_vals, pic50_vals, style, color='grey', alpha=alpha)
        ax.text(38, le_val * 38 / 1.37, f'LE={le_val}',
                fontsize=8, color='grey', va='center')

    # Data points
    if color_by and color_by in df.columns:
        groups = df[color_by].unique()
        colors = cm.viridis(np.linspace(0, 1, len(groups)))
        for g, c in zip(sorted(groups), colors):
            sub = df[df[color_by] == g]
            ax.scatter(sub['HAC'], sub['pIC50'], color=c, s=80,
                       alpha=0.8, label=str(g), zorder=3)
    else:
        ax.scatter(df['HAC'], df['pIC50'], s=80, alpha=0.8, zorder=3)

    if annotate and 'label' in df.columns:
        for _, row in df.iterrows():
            ax.annotate(row['label'], (row['HAC'], row['pIC50']),
                        textcoords='offset points', xytext=(4, 4),
                        fontsize=7)

    # Reference zones
    ax.axvspan(8, 20, alpha=0.05, color='green')
    ax.text(14, 0.5, 'Fragment zone', ha='center', color='green',
            fontsize=8, alpha=0.6)
    ax.axvspan(20, 30, alpha=0.05, color='blue')
    ax.text(25, 0.5, 'Lead-like', ha='center', color='blue',
            fontsize=8, alpha=0.6)

    ax.set_xlabel('Heavy Atom Count (HAC)', fontsize=12)
    ax.set_ylabel('pIC50', fontsize=12)
    ax.set_title('Abad-Zapatero Plot — LE Evolution', fontsize=13)
    ax.legend(title=color_by, fontsize=9)
    ax.set_xlim(5, 40)
    ax.set_ylim(0, 12)
    return fig
```

## Metric-Driven SAR Decision Table

| Observation | Interpretation | Action |
|-------------|----------------|--------|
| LE drops, pIC50 gains | Adding "bulk" lipophilic groups | Switch strategy: polar substituents, ring replace |
| LLE < 2 | cLogP growing faster than potency | Remove lipophilic groups, add polar |
| GE < 0.1 for last group | Substituent inefficient | Remove it; try different vector |
| LE maintained, pIC50 +1 | Efficient growing | Continue in same direction |
| LLE improves ≥ 1 unit | Lipophilicity-neutral potency gain | Excellent — prioritize this scaffold |
| LELP > 10 | Lipophilicity-dominated | Flag for ADMET risk |

## Quick Report Generator

```python
def fbdd_round_report(df, round_idx, parent_smiles=None, parent_pic50=None):
    """Print a FBDD round summary."""
    print(f"\n{'='*60}")
    print(f"FBDD Round {round_idx} — {len(df)} compounds")
    print(f"{'='*60}")

    # Best by LE
    best_le = df.loc[df['LE'].idxmax()]
    print(f"\nBest LE:   {best_le['LE']:.3f}  | {best_le['smiles']}")
    print(f"           pIC50={best_le['pIC50']:.2f}, HAC={best_le['HAC']}, "
          f"LLE={best_le['LLE']:.2f}")

    # Best by LLE
    best_lle = df.loc[df['LLE'].idxmax()]
    print(f"Best LLE:  {best_lle['LLE']:.2f}  | {best_lle['smiles']}")

    # Best by pIC50
    best_pot = df.loc[df['pIC50'].idxmax()]
    print(f"Best pot:  pIC50={best_pot['pIC50']:.2f}  | {best_pot['smiles']}")

    if parent_smiles and parent_pic50:
        parent_mol = Chem.MolFromSmiles(parent_smiles)
        parent_hac = parent_mol.GetNumHeavyAtoms()
        parent_le = LE(parent_pic50, parent_hac)
        print(f"\nParent:    pIC50={parent_pic50:.2f}, LE={parent_le:.3f}, HAC={parent_hac}")
        print(f"Best GE:   {df['pIC50'].max() - parent_pic50:.2f} log units "
              f"in {df['HAC'].iloc[df['pIC50'].idxmax()] - parent_hac} HA")

    print(f"\nLE distribution: ≥0.35: {(df['LE']>=0.35).sum()}, "
          f"0.30-0.35: {((df['LE']>=0.30)&(df['LE']<0.35)).sum()}, "
          f"<0.30: {(df['LE']<0.30).sum()}")
    print(f"LLE distribution: ≥5: {(df['LLE']>=5).sum()}, "
          f"3-5: {((df['LLE']>=3)&(df['LLE']<5)).sum()}, "
          f"<3: {(df['LLE']<3).sum()}")
```
