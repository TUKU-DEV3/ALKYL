# mmpdb 4 — Full CLI Workflow and Python API

## Installation

```bash
pip install mmpdb
# Verify
mmpdb --version    # 4.x
python -c "import mmpdblib; print('OK')"
```

## Input Format

SMILES file (`.smi`) — one molecule per line, `smiles id` tab-separated:

```
c1ccccc1C(=O)O	benzoic_acid
c1ccc(cc1)C(=O)O	4_aminobenzoic_acid
Cc1ccccc1C(=O)O	2_methylbenzoic_acid
CC(=O)Oc1ccccc1C(=O)O	aspirin
```

Properties CSV — must match molecule IDs:

```csv
ID,logD,pIC50,clearance_mL_min_kg
benzoic_acid,1.87,,
4_aminobenzoic_acid,0.42,6.23,12.4
aspirin,1.19,5.87,8.1
```

## Step 1: Fragment

```bash
# Basic fragmentation (default: max_heavies=10, max_cuts=3)
mmpdb fragment compounds.smi -o fragments.h5

# Options:
mmpdb fragment compounds.smi \
    --max-heavies 13 \          # allow larger variable parts (default: 10)
    --max-rotatable-bonds 10 \  # skip highly flexible compounds
    --cut-smarts "[#6;X4:1]-[!#1:2]" \  # custom cut SMARTS (default: C-heteroatom and ring exits)
    --num-jobs 4 \              # parallel fragmentation
    -o fragments.h5

# Check output
mmpdb fraginfo fragments.h5
# Reports: total molecules, fragments, unique variable parts
```

## Step 2: Index (Build MMP Database)

```bash
# Build MMP index from fragments
mmpdb index fragments.h5 -o mmpdb.db

# Options:
mmpdb index fragments.h5 \
    --symmetric \              # include both A→B and B→A transforms
    --max-variable-heavies 13 \
    --min-variable-heavies 0 \
    --title "My Dataset" \
    -o mmpdb.db

# Verify index
mmpdb info mmpdb.db
# Shows: #compounds, #pairs, #rules (transforms)
```

## Step 3: Load Properties

```bash
# Load experimental properties
mmpdb loadprops mmpdb.db properties.csv

# CSV format: first column = molecule ID (matches .smi IDs)
# All remaining columns = property names (numeric)

# Verify
mmpdb propinfo mmpdb.db
# Shows loaded properties and statistics
```

## Step 4: Analyze SAR Rules

```bash
# Extract all SAR rules for a property
mmpdb analyze mmpdb.db \
    --property pIC50 \
    --min-pairs 3 \            # minimum N for a rule to appear
    -o sar_rules_pIC50.csv

# Output columns:
# rule, smirks, from_smiles, to_smiles, count, avg, std, min, max, ...
```

Example output:

| from_smiles | to_smiles | count | avg_pIC50 | std_pIC50 |
|-------------|-----------|-------|-----------|-----------|
| `[*:1]c1ccccc1` | `[*:1]c1ccncc1` | 14 | +0.45 | 0.32 |
| `[*:1]Cl` | `[*:1]F` | 28 | -0.12 | 0.21 |
| `[*:1]OC` | `[*:1]NC` | 9 | +0.38 | 0.55 |

## Step 5: Transform a Query Molecule

```bash
# Apply all known transforms to a query SMILES
mmpdb transform \
    --smiles "c1ccc(cc1)C(=O)OC" \    # methyl benzoate
    mmpdb.db \
    --property pIC50 \
    --min-pairs 3 \
    --score avg \                       # or: max, min, std
    -o transformations.csv

# Output: proposed analogues + predicted ΔpIC50 from MMPA rules
```

## Step 6: Apply Transform Directly

```bash
# Apply a specific transform to a molecule
mmpdb transform \
    --smiles "c1ccccc1C(=O)O" \
    --transform "[*:1]c1ccccc1>>[*:1]c1ccncc1" \   # benzene → pyridine
    mmpdb.db \
    --property pIC50 \
    -o specific_transform.csv
```

## Python API (mmpdblib)

### Building the database programmatically

```python
from mmpdblib import do_fragment, do_index, do_loadprops

import subprocess
import pathlib

def build_mmpdb(smiles_file: str, props_csv: str, db_path: str):
    """Build mmpdb database from SMILES + properties."""
    frag_file = db_path.replace(".db", "_fragments.h5")

    # Step 1: fragment
    subprocess.run([
        "mmpdb", "fragment", smiles_file,
        "--max-heavies", "13",
        "-o", frag_file,
    ], check=True)

    # Step 2: index
    subprocess.run([
        "mmpdb", "index", frag_file,
        "--symmetric",
        "-o", db_path,
    ], check=True)

    # Step 3: loadprops
    subprocess.run([
        "mmpdb", "loadprops", db_path, props_csv,
    ], check=True)

    return db_path
```

### Querying MMPs programmatically

```python
import mmpdblib
from mmpdblib import cli

import sqlite3
import pandas as pd

def query_mmp_pairs(db_path: str, property_name: str, min_pairs: int = 3) -> pd.DataFrame:
    """Extract SAR rules as a DataFrame."""
    result_csv = "/tmp/sar_rules.csv"
    cli.main([
        "analyze", db_path,
        "--property", property_name,
        "--min-pairs", str(min_pairs),
        "-o", result_csv,
    ])
    return pd.read_csv(result_csv)


def transform_smiles(db_path: str, query_smiles: str, property_name: str,
                     min_pairs: int = 3) -> pd.DataFrame:
    """Get analogues with predicted property change for a query molecule."""
    result_csv = "/tmp/transforms.csv"
    cli.main([
        "transform",
        "--smiles", query_smiles,
        db_path,
        "--property", property_name,
        "--min-pairs", str(min_pairs),
        "--score", "avg",
        "-o", result_csv,
    ])
    return pd.read_csv(result_csv)
```

### Direct SQL queries on mmpdb SQLite database

mmpdb stores data in a SQLite file — query directly for custom analyses:

```python
import sqlite3
import pandas as pd

def get_pairs_for_smiles(db_path: str, smiles: str) -> pd.DataFrame:
    """Find all MMP partners of a given SMILES."""
    conn = sqlite3.connect(db_path)

    query = """
    SELECT
        c1.smiles AS smiles_A,
        c2.smiles AS smiles_B,
        rule.from_smiles AS variable_A,
        rule.to_smiles   AS variable_B,
        env.smiles       AS context
    FROM pair
    JOIN compound c1 ON pair.compound1_id = c1.id
    JOIN compound c2 ON pair.compound2_id = c2.id
    JOIN rule     ON pair.rule_id = rule.id
    JOIN environment env ON pair.environment_id = env.id
    WHERE c1.smiles = ?
    LIMIT 200
    """
    df = pd.read_sql_query(query, conn, params=(smiles,))
    conn.close()
    return df


def get_rules_with_properties(db_path: str, property_name: str, min_n: int = 5) -> pd.DataFrame:
    """Extract SAR rules with aggregated property statistics."""
    conn = sqlite3.connect(db_path)

    # Get property id
    prop_id = pd.read_sql_query(
        "SELECT id FROM property_name WHERE name = ?", conn, params=(property_name,)
    ).iloc[0]["id"]

    query = """
    SELECT
        r.from_smiles,
        r.to_smiles,
        COUNT(*) AS n_pairs,
        AVG(delta.value) AS mean_delta,
        group_concat(delta.value) AS all_deltas
    FROM rule r
    JOIN rule_environment re ON re.rule_id = r.id
    JOIN pair_property delta ON delta.pair_environment_id = re.id
        AND delta.property_name_id = ?
    GROUP BY r.id
    HAVING COUNT(*) >= ?
    ORDER BY ABS(AVG(delta.value)) DESC
    """
    df = pd.read_sql_query(query, conn, params=(int(prop_id), min_n))
    conn.close()

    # Parse all_deltas string into std
    import numpy as np
    df["std_delta"] = df["all_deltas"].apply(
        lambda s: float(np.std([float(x) for x in s.split(",")]))
    )
    return df[["from_smiles", "to_smiles", "n_pairs", "mean_delta", "std_delta"]]
```

## Generating a Large MMP Dataset from ChEMBL

```python
import pandas as pd
from mmpdblib import cli
import subprocess

def fetch_chembl_dataset(target_chembl_id: str, output_smi: str, output_props: str):
    """
    Download ChEMBL bioactivity data and prepare for mmpdb.
    Uses the ChEMBL MCP tool or rdkit-chembl script.
    """
    # Use ChEMBL MCP → get_bioactivity(target_chembl_id=..., standard_type="IC50")
    # Then save to .smi + .csv format

    # After fetching:
    # df columns: smiles, molecule_chembl_id, pIC50

    df = pd.read_csv("chembl_data.csv")

    # Save SMILES file
    with open(output_smi, "w") as f:
        for _, row in df.iterrows():
            if pd.notna(row["smiles"]):
                f.write(f"{row['smiles']}\t{row['molecule_chembl_id']}\n")

    # Save properties CSV
    props = df[["molecule_chembl_id", "pIC50"]].rename(columns={"molecule_chembl_id": "ID"})
    props.to_csv(output_props, index=False)
```

## Useful mmpdb CLI Reference

| Command | Purpose |
|---------|---------|
| `mmpdb fragment` | Fragment SMILES → fragments.h5 |
| `mmpdb index` | Build MMP pair index → .db |
| `mmpdb loadprops` | Load property CSV into .db |
| `mmpdb propinfo` | List loaded properties + statistics |
| `mmpdb info` | Database statistics |
| `mmpdb fraginfo` | Fragment file statistics |
| `mmpdb analyze` | Extract SAR rules for a property |
| `mmpdb transform` | Apply transforms to query molecule |
| `mmpdb predict` | Predict property for a SMILES |
| `mmpdb rgroup2smarts` | Convert R-group SMILES to SMARTS |
| `mmpdb help` | Full command listing |

## Key Pitfalls

- **mmpdb 2 → 4 migration**: `mmpdb make_fragments` became `mmpdb fragment`; `mmpdb make_index` became `mmpdb index`; Python API changed completely
- **Fragment file format**: `.h5` for mmpdb ≥ 3.x; older versions used `.fragdb`
- **Property CSV column name**: must be `ID` (exact case) for the identifier column
- **`--symmetric` flag**: without it, only A→B transforms are indexed; with it, both directions indexed (doubles pair count)
- **Large datasets**: >100k compounds may need >8 GB RAM for indexing; use `--num-jobs` for fragmentation
- **Custom cut SMARTS**: default cuts only single non-ring bonds; for ring opening add `--cut-smarts "[R:1]-[!R:2]"`
