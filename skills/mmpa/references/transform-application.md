# Transform Application — Analogue Generation, Library Design, Integration

## Applying MMPA Transforms to a Query Molecule

The primary application of MMPA: take a lead compound, apply all database-supported transforms, and generate predicted analogues with ΔProperty estimates.

### Via mmpdb CLI

```bash
# Apply all known transforms to a query molecule, predict ΔpIC50
mmpdb transform \
    --smiles "c1ccc(cc1)C(=O)NC2CCNCC2" \
    mmpdb.db \
    --property pIC50 \
    --min-pairs 3 \
    --score avg \
    --radius 0 \           # 0 = no context constraint (global rules)
    -o query_transforms.csv

# Output CSV columns:
# start,to_smiles,rule,property,count,avg,std,...
```

```bash
# Restrict to chemically similar contexts (radius=1 = nearest-neighbor context)
mmpdb transform \
    --smiles "c1ccc(cc1)C(=O)NC2CCNCC2" \
    mmpdb.db \
    --property pIC50 \
    --radius 1 \           # require context match at 1 hop from attachment point
    --min-pairs 3 \
    -o context_aware_transforms.csv
```

### Via Python

```python
import pandas as pd
from mmpdblib import cli

def predict_analogues(
    query_smiles: str,
    db_path: str,
    property_name: str = "pIC50",
    min_pairs: int = 3,
    radius: int = 0,
) -> pd.DataFrame:
    """
    Generate analogues from a query SMILES using MMPA transforms.
    Returns DataFrame sorted by predicted improvement.
    """
    import tempfile, os

    out_csv = tempfile.mktemp(suffix=".csv")
    cli.main([
        "transform",
        "--smiles", query_smiles,
        db_path,
        "--property", property_name,
        "--min-pairs", str(min_pairs),
        "--radius", str(radius),
        "--score", "avg",
        "-o", out_csv,
    ])

    if not os.path.exists(out_csv):
        return pd.DataFrame()

    df = pd.read_csv(out_csv)
    # Rename columns for clarity
    df = df.rename(columns={
        "to_smiles": "analogue_smiles",
        "avg": f"pred_delta_{property_name}",
        "count": "n_supporting_pairs",
        "std": f"std_delta_{property_name}",
    })
    # Sort: most positive predicted change first
    if f"pred_delta_{property_name}" in df.columns:
        df = df.sort_values(f"pred_delta_{property_name}", ascending=False)

    return df.reset_index(drop=True)
```

## Focused Library Generation Pipeline

Full pipeline: query compound → MMPA transforms → property prediction → filter → diversity → output.

```python
from rdkit import Chem
from rdkit.Chem import Descriptors, QED, DataStructs, rdMolDescriptors
import pandas as pd
import numpy as np

def generate_focused_library(
    lead_smiles: str,
    db_path: str,
    property_name: str = "pIC50",
    min_delta: float = 0.0,           # only keep improvements
    max_mw: float = 550.0,
    max_logp: float = 5.5,
    min_qed: float = 0.4,
    max_sa: float = 5.0,
    diversity_threshold: float = 0.70, # Tanimoto cutoff for MaxMin diversity
    n_diverse: int = 50,
) -> pd.DataFrame:
    """
    Generate a focused, filtered, diverse analogue library from MMPA transforms.
    """
    # Step 1: Get transforms
    df = predict_analogues(lead_smiles, db_path, property_name, min_pairs=3)
    if df.empty:
        return df

    delta_col = f"pred_delta_{property_name}"
    df = df[df[delta_col] >= min_delta]

    # Step 2: Compute drug-likeness filters
    def get_props(smi):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return None
        try:
            from sascorer import calculateScore
            sa = calculateScore(mol)
        except ImportError:
            sa = 3.0  # fallback
        return {
            "mw": Descriptors.ExactMolWt(mol),
            "logp": Descriptors.MolLogP(mol),
            "qed": QED.qed(mol),
            "sa": sa,
        }

    props_list = df["analogue_smiles"].apply(get_props)
    df = df[props_list.notna()].copy()
    props_df = pd.DataFrame(props_list.dropna().tolist())
    df = pd.concat([df.reset_index(drop=True), props_df.reset_index(drop=True)], axis=1)

    # Step 3: Drug-likeness filter
    df = df[
        (df["mw"] <= max_mw) &
        (df["logp"] <= max_logp) &
        (df["qed"] >= min_qed) &
        (df["sa"] <= max_sa)
    ]

    if df.empty:
        return df

    # Step 4: MaxMin diversity selection
    mols = [Chem.MolFromSmiles(s) for s in df["analogue_smiles"]]
    fps = [rdMolDescriptors.GetMorganFingerprintAsBitVect(m, 2, 2048)
           for m in mols if m is not None]

    if len(fps) <= n_diverse:
        return df.reset_index(drop=True)

    # MaxMin: seed with highest predicted delta
    selected_idx = [0]
    for _ in range(n_diverse - 1):
        max_min_dist = -1
        best_idx = -1
        for i in range(len(fps)):
            if i in selected_idx:
                continue
            min_sim = min(DataStructs.TanimotoSimilarity(fps[i], fps[j])
                         for j in selected_idx)
            if min_sim > max_min_dist:
                max_min_dist = min_sim
                best_idx = i
        if best_idx == -1:
            break
        selected_idx.append(best_idx)

    return df.iloc[selected_idx].reset_index(drop=True)
```

## Multi-Property Optimization

Apply transforms and rank by multiple property improvements simultaneously:

```python
def multi_property_transform(
    lead_smiles: str,
    db_path: str,
    properties: list[str],      # e.g. ["pIC50", "LogD", "clearance"]
    weights: list[float] = None,
    min_pairs: int = 3,
) -> pd.DataFrame:
    """
    Apply MMPA transforms with multi-property scoring.
    """
    if weights is None:
        weights = [1.0] * len(properties)

    # Get transforms per property
    dfs = {}
    for prop in properties:
        df = predict_analogues(lead_smiles, db_path, prop, min_pairs=min_pairs)
        if not df.empty and "analogue_smiles" in df.columns:
            dfs[prop] = df.set_index("analogue_smiles")[f"pred_delta_{prop}"]

    if not dfs:
        return pd.DataFrame()

    # Combine into one DataFrame
    combined = pd.DataFrame(dfs).reset_index()
    combined.columns = ["analogue_smiles"] + [f"delta_{p}" for p in dfs.keys()]
    combined = combined.dropna()

    # Normalize each delta column to [0, 1] range, then weighted sum
    for prop in dfs.keys():
        col = f"delta_{prop}"
        col_min, col_max = combined[col].min(), combined[col].max()
        if col_max > col_min:
            combined[f"{col}_norm"] = (combined[col] - col_min) / (col_max - col_min)
        else:
            combined[f"{col}_norm"] = 0.5

    combined["composite_score"] = sum(
        w * combined[f"delta_{p}_norm"]
        for w, p in zip(weights, dfs.keys())
    )

    return combined.sort_values("composite_score", ascending=False).reset_index(drop=True)
```

## Integration with REINVENT 4 (MMPA-Seeded RL)

Use MMPA-derived transforms to seed the REINVENT prior or define an analogue-focused scoring function:

```python
def mmpa_seeded_smiles_file(lead_smiles: str, db_path: str,
                             property_name: str = "pIC50",
                             output_smi: str = "seed_analogues.smi"):
    """
    Generate MMPA analogues as seed library for REINVENT transfer learning.
    """
    analogues = predict_analogues(lead_smiles, db_path, property_name, min_pairs=3)
    analogues = analogues[analogues[f"pred_delta_{property_name}"] > 0]

    valid = []
    for i, row in analogues.iterrows():
        smi = row["analogue_smiles"]
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            valid.append(smi)

    # Write as SMILES file for REINVENT transfer_learning
    with open(output_smi, "w") as f:
        for smi in valid:
            f.write(f"{smi}\n")

    print(f"Wrote {len(valid)} seed analogues to {output_smi}")
    return output_smi
```

## Integration with Docking (MMPA → Vina Scoring)

```python
import subprocess

def dock_mmpa_library(analogues_df: pd.DataFrame, receptor_pdbqt: str,
                      center: tuple, box_size: tuple = (20, 20, 20),
                      n_top: int = 20) -> pd.DataFrame:
    """
    Dock top MMPA analogues with AutoDock Vina.
    Requires: meeko (ligand prep), vina binary.
    """
    from vina import Vina
    import tempfile, os
    from rdkit.Chem import AllChem

    v = Vina(sf_name="vina")
    v.set_receptor(receptor_pdbqt)
    v.compute_vina_maps(center=list(center), box_size=list(box_size))

    results = []
    for i, row in analogues_df.head(n_top).iterrows():
        smi = row["analogue_smiles"]
        mol = Chem.MolFromSmiles(smi)
        if mol is None: continue

        # Generate 3D conformer
        mol_3d = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol_3d, AllChem.ETKDGv3())
        AllChem.MMFFOptimizeMolecule(mol_3d)

        # Write SDF
        sdf_path = f"/tmp/lig_{i}.sdf"
        with Chem.SDWriter(sdf_path) as w:
            w.write(mol_3d)

        # Convert to PDBQT (meeko)
        pdbqt_path = sdf_path.replace(".sdf", ".pdbqt")
        subprocess.run(["mk_prepare_ligand.py", "-i", sdf_path, "-o", pdbqt_path],
                       check=True, capture_output=True)

        # Dock
        v.set_ligand_from_file(pdbqt_path)
        v.dock(exhaustiveness=8, n_poses=1)
        score = v.score()[0]

        results.append({**row.to_dict(), "vina_score": score})

    df_docked = pd.DataFrame(results)
    if "vina_score" in df_docked.columns:
        df_docked = df_docked.sort_values("vina_score")  # lower = better
    return df_docked
```

## MMPA Transform Novelty Check

Verify generated analogues are not already in the training database:

```python
def novelty_check(analogues_df: pd.DataFrame, training_smiles: list,
                  similarity_threshold: float = 0.85) -> pd.DataFrame:
    """
    Flag analogues already known in training set (Tanimoto > threshold).
    """
    from rdkit.Chem import DataStructs, rdMolDescriptors

    train_fps = [
        rdMolDescriptors.GetMorganFingerprintAsBitVect(
            Chem.MolFromSmiles(s), 2, 2048
        )
        for s in training_smiles
        if Chem.MolFromSmiles(s) is not None
    ]

    def max_sim_to_train(smi):
        mol = Chem.MolFromSmiles(smi)
        if mol is None: return 1.0
        fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, 2048)
        sims = DataStructs.BulkTanimotoSimilarity(fp, train_fps)
        return max(sims) if sims else 0.0

    analogues_df = analogues_df.copy()
    analogues_df["max_train_similarity"] = analogues_df["analogue_smiles"].apply(max_sim_to_train)
    analogues_df["is_novel"] = analogues_df["max_train_similarity"] < similarity_threshold

    n_novel = analogues_df["is_novel"].sum()
    print(f"Novel analogues: {n_novel}/{len(analogues_df)} "
          f"({100*n_novel/len(analogues_df):.1f}%)")
    return analogues_df


def synthesis_feasibility_screen(analogues_df: pd.DataFrame,
                                  max_sa: float = 4.0) -> pd.DataFrame:
    """Filter analogues by synthetic accessibility (SA ≤ 4 = synthesizable)."""
    from sascorer import calculateScore

    def sa(smi):
        mol = Chem.MolFromSmiles(smi)
        return calculateScore(mol) if mol else 10.0

    analogues_df = analogues_df.copy()
    analogues_df["sa_score"] = analogues_df["analogue_smiles"].apply(sa)
    return analogues_df[analogues_df["sa_score"] <= max_sa].reset_index(drop=True)
```

## Typical MMPA Design Cycle

```
1. Build mmpdb from project dataset (compounds + measured properties)
        mmpdb fragment → index → loadprops

2. Analyze SAR rules
        mmpdb analyze → top transforms table

3. Apply to lead compound
        mmpdb transform → predicted analogues with ΔP

4. Filter and rank
        Drug-likeness (Lipinski) + SA score + novelty check

5. MaxMin diversity selection
        Pick 20-50 diverse representatives

6. Orthogonal validation
        a. Docking score (Vina/Gnina) for binding
        b. QikProp or SwissADME for ADMET

7. Wet-lab synthesis and testing
        Feed results back → new mmpdb entries → iterate
```

## Key Pitfalls

- **Transform direction matters**: A→B gives +ΔP; B→A gives −ΔP — always check sign when applying
- **`--radius 0` (global rules)**: ignores context — fast but may not apply to your scaffold; use `--radius 1` for context-aware predictions
- **Chained transforms**: mmpdb applies one transform per call; for multi-step modifications, call iteratively and re-filter after each step
- **SA score cap**: MMPA can suggest transforms that are known to improve potency but are hard to synthesize; always screen with SA ≤ 4
- **Property units consistency**: ensure all measured values use consistent units (pIC50, not IC50 nM) before loading into mmpdb
- **Circular transforms**: if context is too small, a transform might regenerate the starting material; add `smiles_a != smiles_b` check
