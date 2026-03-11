# Kinetic QSAR and Structure-Kinetics Relationships

## Structure-Kinetics Relationships (SKR)

Unlike equilibrium SAR, kinetic SAR studies how chemical modifications
affect kon and koff independently. Key principle:
- **koff** is sensitive to ligand-protein complementarity (shape/H-bond fit)
- **kon** is sensitive to charge, desolvation, conformational flexibility
- Modifying a H-bond donor often affects koff more than kon

## Feature Computation for Kinetic Models

```python
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Crippen, AllChem
import numpy as np
import pandas as pd

def kinetic_features(smiles):
    """
    Compute features relevant to koff/kon prediction.
    Returns dict of descriptors.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    # Standard descriptors
    features = {
        'MW': Descriptors.ExactMolWt(mol),
        'HAC': mol.GetNumHeavyAtoms(),
        'cLogP': Crippen.MolLogP(mol),
        'HBD': rdMolDescriptors.CalcNumHBD(mol),
        'HBA': rdMolDescriptors.CalcNumHBA(mol),
        'RotBonds': rdMolDescriptors.CalcNumRotatableBonds(mol),
        'PSA': rdMolDescriptors.CalcTPSA(mol),
        'Fsp3': rdMolDescriptors.CalcFractionCSP3(mol),
        'ArRings': rdMolDescriptors.CalcNumAromaticRings(mol),
        'Rings': rdMolDescriptors.CalcNumRings(mol),
        'Stereocenters': rdMolDescriptors.CalcNumAtomStereoCenters(mol),
    }

    # Kinetics-specific: rigidity (fewer conformations → slower association)
    n_rot = features['RotBonds']
    features['rigidity'] = 1.0 / (1.0 + n_rot)

    # Charge state descriptors
    features['formal_charge'] = Chem.GetFormalCharge(mol)
    pos_atoms = sum(1 for a in mol.GetAtoms() if a.GetFormalCharge() > 0)
    neg_atoms = sum(1 for a in mol.GetAtoms() if a.GetFormalCharge() < 0)
    features['n_pos_charge'] = pos_atoms
    features['n_neg_charge'] = neg_atoms

    # Molecular complexity (Bertz)
    from rdkit.Chem.GraphDescriptors import BertzCT
    features['bertz_complexity'] = BertzCT(mol)

    # Morgan fingerprint (for similarity-based models)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=1024)
    features['ecfp4'] = np.array(fp)

    return features

def compute_kinetic_dataset(smiles_list, kon_list=None, koff_list=None):
    """Build feature matrix for kinetic QSAR."""
    records = []
    for i, smi in enumerate(smiles_list):
        feat = kinetic_features(smi)
        if feat is None:
            continue
        row = {k: v for k, v in feat.items() if k != 'ecfp4'}
        row['smiles'] = smi
        if kon_list: row['log_kon'] = np.log10(kon_list[i])
        if koff_list: row['log_koff'] = np.log10(koff_list[i])
        records.append(row)
    return pd.DataFrame(records)
```

## QSAR Models for koff

```python
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.model_selection import cross_val_score, LeaveOneOut
from sklearn.preprocessing import StandardScaler
import numpy as np

def build_koff_model(df, target='log_koff', method='rf'):
    """
    Build koff QSAR model.
    LOOCV recommended for small kinetic datasets (n=20–100).
    """
    feature_cols = ['MW', 'HAC', 'cLogP', 'HBD', 'HBA', 'RotBonds',
                    'PSA', 'Fsp3', 'ArRings', 'rigidity', 'bertz_complexity',
                    'formal_charge', 'n_pos_charge', 'n_neg_charge']

    X = df[feature_cols].values
    y = df[target].values

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    if method == 'rf':
        model = RandomForestRegressor(n_estimators=200, min_samples_leaf=2,
                                       random_state=42)
    elif method == 'gbm':
        model = GradientBoostingRegressor(n_estimators=200, max_depth=3,
                                           random_state=42)

    # LOOCV (for small n)
    loo = LeaveOneOut()
    loo_scores = cross_val_score(model, X_scaled, y, cv=loo, scoring='r2')
    print(f"LOOCV R² = {loo_scores.mean():.3f} ± {loo_scores.std():.3f}")

    # Final model on all data
    model.fit(X_scaled, y)

    return model, scaler, feature_cols, {'loocv_r2': loo_scores.mean()}

def predict_koff(model, scaler, feature_cols, new_smiles_list):
    """Predict log_koff for new compounds."""
    records = []
    for smi in new_smiles_list:
        feat = kinetic_features(smi)
        if feat is None:
            records.append(None)
            continue
        X = np.array([[feat.get(c, 0) for c in feature_cols]])
        X_scaled = scaler.transform(X)
        log_koff_pred = model.predict(X_scaled)[0]
        records.append({
            'smiles': smi,
            'log_koff_pred': log_koff_pred,
            'koff_pred': 10**log_koff_pred,
            'RT_pred_s': 1.0 / 10**log_koff_pred
        })
    return pd.DataFrame(r for r in records if r)
```

## Fingerprint-Based Kinetic Model

For larger datasets (n > 100), use ECFP4 fingerprints:

```python
from rdkit.Chem import AllChem
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import PairwiseKernel

def tanimoto_kernel(X1, X2):
    dot = X1 @ X2.T
    n1 = X1.sum(axis=1, keepdims=True)
    n2 = X2.sum(axis=1, keepdims=True)
    return dot / (n1 + n2.T - dot + 1e-9)

def gp_koff_model(smiles_list, log_koff_list):
    """GP with Tanimoto kernel — uncertainty-aware koff prediction."""
    X = np.array([
        AllChem.GetMorganFingerprintAsBitVect(
            Chem.MolFromSmiles(s), 2, 2048
        ) for s in smiles_list
    ])
    y = np.array(log_koff_list)

    kernel = PairwiseKernel(metric=tanimoto_kernel)
    gp = GaussianProcessRegressor(kernel=kernel, alpha=0.2,
                                   normalize_y=True, random_state=42)
    gp.fit(X, y)
    return gp

# Predict with uncertainty for active learning (see active-learning skill)
def gp_koff_predict(gp, new_smiles):
    X_new = np.array([
        AllChem.GetMorganFingerprintAsBitVect(
            Chem.MolFromSmiles(s), 2, 2048
        ) for s in new_smiles
    ])
    mu, sigma = gp.predict(X_new, return_std=True)
    return pd.DataFrame({
        'smiles': new_smiles,
        'log_koff_pred': mu,
        'log_koff_std': sigma,
        'koff_pred': 10**mu,
        'RT_pred_s': 1.0 / 10**mu
    })
```

## Kinetic Cliffs Detection

Analogous to activity cliffs but for koff:

```python
from rdkit.Chem import DataStructs

def detect_kinetic_cliffs(df, similarity_threshold=0.7, cliff_threshold=2.0):
    """
    Find pairs with high structural similarity but large Δlog_koff.
    cliff_threshold: log units difference (100× koff difference)
    """
    fps = [AllChem.GetMorganFingerprintAsBitVect(
               Chem.MolFromSmiles(s), 2, 2048)
           for s in df.smiles]

    cliffs = []
    n = len(df)
    for i in range(n):
        for j in range(i+1, n):
            sim = DataStructs.TanimotoSimilarity(fps[i], fps[j])
            if sim >= similarity_threshold:
                delta_koff = abs(df.iloc[i]['log_koff'] - df.iloc[j]['log_koff'])
                if delta_koff >= cliff_threshold:
                    cliffs.append({
                        'compound_i': df.iloc[i]['smiles'],
                        'compound_j': df.iloc[j]['smiles'],
                        'tanimoto': sim,
                        'delta_log_koff': delta_koff,
                        'koff_ratio': 10**delta_koff,
                        'RT_i': 1.0 / df.iloc[i]['koff'],
                        'RT_j': 1.0 / df.iloc[j]['koff'],
                    })
    return pd.DataFrame(cliffs).sort_values('delta_log_koff', ascending=False)
```

## Kinetic Map: LE vs koff_LE

Track optimization efficiency in kinetic space:

```python
def kinetic_optimization_map(df):
    """
    2D map: LE (equilibrium efficiency) vs koff_LE (kinetic efficiency).
    Compounds in upper-right: both thermodynamically efficient and slow dissociating.
    """
    import matplotlib.pyplot as plt

    df = df.copy()
    df['LE'] = 1.37 * df['pKD'] / df['HAC']
    df['koff_LE'] = -np.log10(df['koff']) / df['HAC']

    fig, ax = plt.subplots(figsize=(8, 7))
    sc = ax.scatter(df['LE'], df['koff_LE'],
                    c=df['log_koff'], cmap='RdYlGn_r',
                    s=80, alpha=0.8, zorder=3)
    plt.colorbar(sc, ax=ax, label='log(koff)')

    ax.axhline(y=0.20, color='orange', linestyle='--', alpha=0.6,
               label='koff_LE=0.20 (min)')
    ax.axvline(x=0.30, color='blue', linestyle='--', alpha=0.6,
               label='LE=0.30 (min)')

    if 'label' in df.columns:
        for _, row in df.iterrows():
            ax.annotate(row['label'], (row['LE'], row['koff_LE']),
                        textcoords='offset points', xytext=(5, 3), fontsize=8)

    ax.set_xlabel('Ligand Efficiency (LE)', fontsize=12)
    ax.set_ylabel('Kinetic LE (-log koff / HAC)', fontsize=12)
    ax.set_title('Kinetic Optimization Map', fontsize=13)
    ax.legend(fontsize=9)
    return fig
```

## MMPA for Kinetic SAR

Identify transforms that improve koff without affecting kon (= RT improvement):

```python
def mmpa_kinetic_analysis(df, koff_col='log_koff', kon_col='log_kon'):
    """
    Use MMPA to find transforms that improve koff.
    Requires mmpdb-indexed dataset.
    Result: transforms ranked by Δlog_koff with context.
    """
    import subprocess

    # Write to mmpdb format
    df[['smiles', 'compound_id', koff_col]].to_csv(
        '/tmp/kinetic_props.csv', index=False
    )

    # Run mmpdb analysis (see mmpa skill)
    result = subprocess.run([
        'mmpdb', 'loadprops',
        '--properties', '/tmp/kinetic_props.csv',
        'fragments.mmpdb'
    ], capture_output=True, text=True)

    result2 = subprocess.run([
        'mmpdb', 'analyze',
        '--property', koff_col,
        '--min-count', '3',
        'fragments.mmpdb'
    ], capture_output=True, text=True)

    return result2.stdout

# Key insight: transforms that slow koff often add:
# - H-bond donors pointing into specific pockets
# - Rigid linkers reducing entropic cost of dissociation
# - Bulky substituents trapping the ligand (kinetic trap)
```

## Reporting Kinetic Data

```python
def kinetic_summary_table(df):
    """Standard table for kinetic SAR reporting."""
    cols = ['smiles', 'kon', 'koff', 'KD', 'KD_nM', 'RT_s', 't_half_min',
            'LE', 'koff_LE', 'dH', 'dG', 'TdS']
    available = [c for c in cols if c in df.columns]
    summary = df[available].copy()

    if 'kon' in summary:
        summary['kon'] = summary['kon'].apply(lambda x: f'{x:.2e}')
    if 'koff' in summary:
        summary['koff'] = summary['koff'].apply(lambda x: f'{x:.2e}')
    if 't_half_min' not in summary and 'koff' in df:
        summary['t_half_min'] = (np.log(2) / df['koff'] / 60).round(1)

    return summary.sort_values('koff')   # slowest koff first
```
