# SAR Delta Analysis — Property Statistics, Activity Cliffs, Bioisosteres, Visualization

## Aggregating SAR Rules from MMP Pairs

After finding all MMP pairs, aggregate ΔP by transform to extract statistically supported SAR rules:

```python
import pandas as pd
import numpy as np
from scipy import stats

def aggregate_sar_rules(mmp_df: pd.DataFrame, property_col: str,
                        min_pairs: int = 3) -> pd.DataFrame:
    """
    Aggregate ΔProperty per transform into SAR rules.

    Args:
        mmp_df: DataFrame from build_mmp_dataframe() or mmpdb CSV output
        property_col: name of delta column, e.g. "delta_pIC50"
        min_pairs: minimum pairs per rule (default: 3)
    """
    rules = mmp_df.groupby("transform")[property_col].agg(
        n_pairs="count",
        mean_delta="mean",
        std_delta="std",
        median_delta="median",
        min_delta="min",
        max_delta="max",
    ).reset_index()

    rules = rules[rules["n_pairs"] >= min_pairs]
    rules["sem"] = rules["std_delta"] / np.sqrt(rules["n_pairs"])
    rules["ci_95"] = 1.96 * rules["sem"]

    # Sort by |mean_delta| — most impactful transforms first
    rules = rules.reindex(rules["mean_delta"].abs().sort_values(ascending=False).index)

    return rules.reset_index(drop=True)


def top_rules(rules_df: pd.DataFrame, n: int = 20, direction: str = "both") -> pd.DataFrame:
    """Filter to top N most impactful rules."""
    if direction == "positive":
        return rules_df[rules_df["mean_delta"] > 0].head(n)
    elif direction == "negative":
        return rules_df[rules_df["mean_delta"] < 0].head(n)
    else:
        return rules_df.head(n)
```

## Visualizing ΔProperty Distributions

### Distribution of ΔLogP for all H→F transforms

```python
import matplotlib.pyplot as plt
import seaborn as sns

def plot_delta_distribution(mmp_df: pd.DataFrame, transform_query: str,
                             property_col: str = "delta_logp"):
    """Plot ΔP distribution for a specific transform."""
    subset = mmp_df[mmp_df["transform"].str.contains(transform_query, regex=False)]
    deltas = subset[property_col].dropna()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

    # Histogram + KDE
    sns.histplot(deltas, kde=True, ax=ax1, color="steelblue")
    ax1.axvline(deltas.mean(), color="red", ls="--", label=f"Mean={deltas.mean():.2f}")
    ax1.axvline(0, color="gray", ls=":", alpha=0.5)
    ax1.set_xlabel(f"Δ{property_col.replace('delta_', '')}")
    ax1.set_title(f"Transform: {transform_query}\nN={len(deltas)}, σ={deltas.std():.2f}")
    ax1.legend()

    # Box plot
    ax2.boxplot(deltas, vert=True)
    ax2.axhline(0, color="gray", ls=":", alpha=0.5)
    ax2.set_ylabel(f"Δ{property_col}")
    ax2.set_title(f"N={len(deltas)}")

    plt.tight_layout()
    return fig
```

### SAR heatmap: transform vs. property

```python
def plot_sar_heatmap(rules_df: pd.DataFrame, top_n: int = 25,
                     value_col: str = "mean_delta", figsize=(12, 8)):
    """
    Heatmap of top transforms × properties (when multiple properties available).
    """
    # Pivot: index=transform, columns=property, values=mean_delta
    # Assumes rules_df has multiple property_col rows per transform
    pivot = rules_df.pivot_table(
        index="transform", columns="property", values="mean_delta",
        aggfunc="mean"
    ).head(top_n)

    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(
        pivot, annot=True, fmt=".2f", center=0,
        cmap="RdBu_r", linewidths=0.5, ax=ax,
        cbar_kws={"label": "Mean ΔProperty"}
    )
    ax.set_title(f"SAR Rules Heatmap (Top {top_n} transforms)")
    plt.tight_layout()
    return fig
```

## Activity Cliff Detection

```python
from rdkit.Chem import DataStructs, rdMolDescriptors

def find_activity_cliffs(
    smiles_dict: dict,
    activity_dict: dict,      # {id: pIC50 value}
    similarity_threshold: float = 0.70,
    activity_threshold: float = 1.0,   # in log units
) -> pd.DataFrame:
    """
    Identify activity cliffs: high structural similarity but large activity difference.
    Uses Morgan2 Tanimoto + SALI score.
    """
    ids = list(smiles_dict.keys())
    mols = {id_: Chem.MolFromSmiles(smi) for id_, smi in smiles_dict.items()
            if Chem.MolFromSmiles(smi) is not None}
    fps = {id_: rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, 2048)
           for id_, mol in mols.items()}

    cliffs = []
    ids_with_activity = [id_ for id_ in ids if id_ in activity_dict and id_ in fps]

    for i in range(len(ids_with_activity)):
        for j in range(i + 1, len(ids_with_activity)):
            id_a = ids_with_activity[i]
            id_b = ids_with_activity[j]

            sim = DataStructs.TanimotoSimilarity(fps[id_a], fps[id_b])
            if sim < similarity_threshold:
                continue

            act_a = activity_dict[id_a]
            act_b = activity_dict[id_b]
            delta_act = abs(act_a - act_b)

            if delta_act >= activity_threshold:
                sali = delta_act / (1 - sim) if sim < 1.0 else float("inf")
                cliffs.append({
                    "id_a": id_a, "id_b": id_b,
                    "smiles_a": smiles_dict[id_a], "smiles_b": smiles_dict[id_b],
                    "activity_a": act_a, "activity_b": act_b,
                    "delta_activity": act_b - act_a,
                    "tanimoto": sim,
                    "sali": sali,
                })

    df = pd.DataFrame(cliffs)
    if not df.empty:
        df = df.sort_values("sali", ascending=False).reset_index(drop=True)
    return df
```

### Visualize activity landscape

```python
from rdkit.Chem import Draw
from mols2grid import MolsGrid  # pip install mols2grid

def visualize_cliff_pairs(cliffs_df: pd.DataFrame, top_n: int = 10):
    """Show top activity cliff pairs as side-by-side molecule images."""
    rows = []
    for _, row in cliffs_df.head(top_n).iterrows():
        rows.append({
            "SMILES": row["smiles_a"],
            "label": f"{row['id_a']}\npIC50={row['activity_a']:.2f}"
        })
        rows.append({
            "SMILES": row["smiles_b"],
            "label": f"{row['id_b']}\npIC50={row['activity_b']:.2f}\nΔ={row['delta_activity']:+.2f}"
        })
    return MolsGrid.from_list(rows, smiles_col="SMILES", subset=["img", "label"]).display()
```

## Bioisostere Table Construction

```python
def build_bioisostere_table(
    mmp_df: pd.DataFrame,
    delta_col: str = "delta_pIC50",
    min_pairs: int = 5,
    max_std: float = 0.8,  # exclude high-variance (context-dependent) transforms
) -> pd.DataFrame:
    """
    Build a bioisostere table: transforms with small ΔActivity and large ΔLogP.
    These are useful bioisosteres: same potency, improved properties.
    """
    rules = aggregate_sar_rules(mmp_df, delta_col, min_pairs=min_pairs)

    # Add LogP delta if available
    if "delta_logp" in mmp_df.columns:
        logp_rules = aggregate_sar_rules(mmp_df, "delta_logp", min_pairs=min_pairs)
        logp_rules = logp_rules.rename(columns={
            "mean_delta": "mean_delta_logp",
            "std_delta": "std_delta_logp",
            "n_pairs": "n_pairs_logp",
        })
        rules = rules.merge(
            logp_rules[["transform", "mean_delta_logp", "n_pairs_logp"]],
            on="transform", how="left"
        )

    # Filter: low ΔActivity + acceptable variance
    bioisosteres = rules[
        (rules["mean_delta"].abs() < 0.3) &   # < 2-fold activity change
        (rules["std_delta"] < max_std) &
        (rules["n_pairs"] >= min_pairs)
    ].copy()

    if "mean_delta_logp" in bioisosteres.columns:
        bioisosteres = bioisosteres.sort_values(
            "mean_delta_logp", ascending=False   # prefer LogP-improving bioisosteres
        )

    return bioisosteres.reset_index(drop=True)
```

## Statistical Tests for SAR Rules

```python
from scipy import stats

def test_transform_significance(delta_values: list, alpha: float = 0.05) -> dict:
    """
    Test if a transform has statistically significant effect.
    Uses one-sample t-test (H0: mean=0) and Wilcoxon signed-rank test.
    """
    arr = np.array(delta_values)
    n = len(arr)

    if n < 3:
        return {"significant": False, "reason": "N < 3"}

    # One-sample t-test vs 0
    t_stat, p_ttest = stats.ttest_1samp(arr, 0)

    # Wilcoxon signed-rank (non-parametric)
    if n >= 5:
        w_stat, p_wilcoxon = stats.wilcoxon(arr)
    else:
        p_wilcoxon = None

    return {
        "n": n,
        "mean": float(arr.mean()),
        "std": float(arr.std()),
        "t_statistic": float(t_stat),
        "p_ttest": float(p_ttest),
        "p_wilcoxon": p_wilcoxon,
        "significant_ttest": p_ttest < alpha,
        "significant_both": (p_wilcoxon is not None) and (p_ttest < alpha) and (p_wilcoxon < alpha),
        "effect_direction": "positive" if arr.mean() > 0 else "negative",
    }
```

## R-Group Contribution Heatmap

When analyzing a congeneric series (fixed scaffold, varying R-groups):

```python
from rdkit.Chem import RGroupDecomposition

def rgroup_activity_heatmap(smiles_list: list, activities: list, core_smarts: str):
    """
    Decompose molecules into core + R-groups, then plot R-group vs activity.
    """
    from rdkit.Chem.rdRGroupDecomposition import RGroupDecompose, RGroupDecompositionParameters

    core = Chem.MolFromSmarts(core_smarts)
    mols = [Chem.MolFromSmiles(s) for s in smiles_list]
    mols = [m for m in mols if m is not None]

    params = RGroupDecompositionParameters()
    params.removeHydrogensPostMatch = True

    groups, _ = RGroupDecompose([core], mols, asSmiles=True, asRows=True)

    df = pd.DataFrame(groups)
    df["pIC50"] = activities[:len(df)]

    # Pivot: R1 position vs pIC50 mean
    if "R1" in df.columns:
        r1_activity = df.groupby("R1")["pIC50"].agg(["mean", "count"]).reset_index()
        r1_activity.columns = ["R1_group", "mean_pIC50", "count"]

        fig, ax = plt.subplots(figsize=(10, 5))
        bars = ax.bar(range(len(r1_activity)),
                      r1_activity["mean_pIC50"],
                      color=["green" if v > df["pIC50"].mean() else "red"
                             for v in r1_activity["mean_pIC50"]])
        ax.set_xticks(range(len(r1_activity)))
        ax.set_xticklabels(r1_activity["R1_group"], rotation=45, ha="right")
        ax.set_ylabel("Mean pIC50")
        ax.set_title(f"R-group contribution (N={len(df)})")
        plt.tight_layout()
        return fig, r1_activity

    return None, df
```

## Summary Report

```python
def mmpa_report(mmp_df: pd.DataFrame, property_col: str = "delta_pIC50",
                output_path: str = "mmpa_report.html"):
    """Generate an HTML summary report of MMPA results."""
    rules = aggregate_sar_rules(mmp_df, property_col, min_pairs=3)
    cliffs = find_activity_cliffs(
        dict(zip(mmp_df["id_a"].tolist() + mmp_df["id_b"].tolist(),
                 mmp_df["smiles_a"].tolist() + mmp_df["smiles_b"].tolist())),
        activity_dict={},  # populate from your data
    )

    report = f"""
    <h2>MMPA Report</h2>
    <p>Total MMP pairs: {len(mmp_df)}</p>
    <p>Unique transforms: {mmp_df['transform'].nunique()}</p>
    <p>Rules (N≥3): {len(rules)}</p>
    <h3>Top 10 Positive Rules (improves {property_col})</h3>
    {rules[rules['mean_delta'] > 0].head(10).to_html(index=False)}
    <h3>Top 10 Negative Rules (worsens {property_col})</h3>
    {rules[rules['mean_delta'] < 0].head(10).to_html(index=False)}
    """

    with open(output_path, "w") as f:
        f.write(f"<html><body>{report}</body></html>")
    print(f"Report saved: {output_path}")
```

## Key Thresholds Summary

| Criterion | Threshold | Rationale |
|-----------|-----------|-----------|
| Min pairs per rule | N ≥ 3 | Statistical minimum; prefer N ≥ 10 |
| Activity cliff Δ | |ΔpIC50| ≥ 1.0 | = 10-fold potency difference |
| Structural similarity | Tanimoto ≥ 0.70 | Close analogues |
| SALI cliff threshold | SALI > 5 | Meaningful landscape discontinuity |
| Bioisostere ΔActivity | |ΔpIC50| < 0.3 | < 2-fold potency change |
| High context dependence | σ(ΔP) > 0.8 | Unreliable rule — context matters |
