# Docking — Virtual Screening Pipeline

High-throughput docking: from compound library to prioritized hit list.

---

## VS Pipeline Overview

```
Compound library (SMILES/SDF)
        │
        ▼
[1] Pre-filters (drug-likeness, PAINS, MW)      → chem_filter.py
        │
        ▼
[2] 3D conformer generation                     → chem_3d.py / obabel
        │
        ▼
[3] Ligand PDBQT preparation                    → meeko / obabel
        │
        ▼
[4] Docking (Vina batch)                         → vina (parallel)
        │
        ▼
[5] Score-based pre-selection (top 5–10%)
        │
        ▼
[6] CNN rescoring (Gnina)                        → gnina --score_only
        │
        ▼
[7] Interaction fingerprint filtering            → ProLIF
        │
        ▼
[8] Clustering + diversity selection             → chem_diversity.py
        │
        ▼
[9] Final ranked hit list + visualization        → py3Dmol
```

---

## Step 1 — Pre-filtering the Library

Always filter before docking — reduces compute cost dramatically.

```python
import subprocess
import json
from pathlib import Path

def prefilter_library(input_smi: str, output_smi: str) -> dict:
    """Apply drug-likeness + PAINS filters via chem_filter.py."""
    result = subprocess.run(
        ["python", "/path/to/scripts/chem_filter.py",
         "--input", input_smi,
         "--lipinski", "--pains", "--veber",
         "--out", output_smi],
        capture_output=True, text=True
    )
    stats = json.loads(result.stdout)
    return stats

stats = prefilter_library("zinc_subset.smi", "filtered.smi")
print(f"Kept {stats['passed']} / {stats['total']} compounds")
```

**Typical filter cascade for VS:**

| Filter | Threshold | Typical rejection rate |
|--------|-----------|----------------------|
| MW | 150–500 Da | ~30% of ZINC |
| cLogP | −2 to 5 | ~20% |
| PAINS | fail | ~5–15% |
| Structural alerts | fail | ~10% |
| Rotatable bonds | ≤ 10 | ~15% |

---

## Step 2 & 3 — 3D Generation + PDBQT Prep

```python
import subprocess
from pathlib import Path

def prepare_vs_library(filtered_smi: str, pdbqt_dir: str,
                       n_workers: int = 8) -> None:
    """
    Generate 3D conformers and convert to PDBQT.
    Uses obabel for speed (meeko for accuracy).
    """
    Path(pdbqt_dir).mkdir(exist_ok=True)

    # obabel: SMILES → individual PDBQT with embedded 3D
    subprocess.run([
        "obabel", filtered_smi,
        "-O", str(Path(pdbqt_dir) / "lig.pdbqt"),
        "--gen3d", "--fastest",   # 'best' for accuracy, 'fastest' for VS
        "-h",
        "-m",                     # split into separate files
        "--errorlevel", "1"
    ], check=True)

    n = len(list(Path(pdbqt_dir).glob("*.pdbqt")))
    print(f"Prepared {n} ligand PDBQT files")

prepare_vs_library("filtered.smi", "pdbqt_library/")
```

---

## Step 4 — Parallel Docking

```python
import subprocess, csv
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

def dock_ligand(args: tuple) -> dict:
    receptor, lig_pdbqt, box, out_dir, exhaustiveness = args
    name = Path(lig_pdbqt).stem
    out_pdbqt = str(Path(out_dir) / f"{name}.pdbqt")

    cmd = ["vina",
           "--receptor", receptor,
           "--ligand",   lig_pdbqt,
           "--out",      out_pdbqt,
           "--exhaustiveness", str(exhaustiveness),
           "--num_modes", "1"]
    for k, v in box.items():
        cmd += [f"--{k}", str(v)]

    proc = subprocess.run(cmd, capture_output=True, text=True)

    score = float("inf")
    for line in proc.stdout.splitlines():
        s = line.strip()
        if s.startswith("1 "):
            try:
                score = float(s.split()[1])
            except (ValueError, IndexError):
                pass
            break

    return {"name": name, "score": score,
            "pdbqt": out_pdbqt, "ok": proc.returncode == 0}


def run_vs(receptor: str, pdbqt_dir: str, box: dict,
           out_dir: str, exhaustiveness: int = 8,
           workers: int = 8) -> list:
    Path(out_dir).mkdir(exist_ok=True)
    ligands = list(Path(pdbqt_dir).glob("*.pdbqt"))
    print(f"Docking {len(ligands)} compounds ...")

    args = [(receptor, str(lig), box, out_dir, exhaustiveness)
            for lig in ligands]

    results = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for res in as_completed(ex.submit(dock_ligand, a) for a in args):
            results.append(res.result())

    results.sort(key=lambda x: x["score"])

    with open("vs_results.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["name", "score", "pdbqt", "ok"])
        w.writeheader()
        w.writerows(results)

    return results

box = {"center_x": 10.5, "center_y": -2.3, "center_z": 14.1,
       "size_x": 22, "size_y": 20, "size_z": 18}

results = run_vs("receptor.pdbqt", "pdbqt_library/", box,
                 "vs_output/", exhaustiveness=8, workers=8)

print(f"Top hit: {results[0]['name']} ({results[0]['score']:.2f} kcal/mol)")
```

---

## Step 5–6 — CNN Rescoring of Top Hits

```python
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd

def select_top_fraction(results: list, fraction: float = 0.05) -> list:
    """Select top N% by Vina score."""
    n = max(1, int(len(results) * fraction))
    return [r for r in results if r["ok"]][:n]

def convert_pdbqt_to_sdf(pdbqt_files: list, out_sdf: str) -> None:
    """Merge multiple PDBQT poses into one SDF for Gnina rescoring."""
    mols = []
    for pdbqt in pdbqt_files:
        result = subprocess.run(
            ["obabel", pdbqt, "-O", "/dev/stdout", "-osdf"],
            capture_output=True, text=True
        )
        # crude parse — obabel writes SDF to stdout
        mols.append(result.stdout)
    with open(out_sdf, "w") as f:
        f.write("\n".join(mols))

def gnina_rescore_batch(receptor: str, poses_sdf: str,
                        out_sdf: str) -> pd.DataFrame:
    subprocess.run([
        "gnina", "--receptor", receptor,
        "--ligand", poses_sdf,
        "--score_only", "--out", out_sdf
    ], check=True, capture_output=True)

    rows = []
    for mol in Chem.SDMolSupplier(out_sdf, removeHs=False):
        if mol is None:
            continue
        p = mol.GetPropsAsDict()
        rows.append({
            "name":         mol.GetProp("_Name") if mol.HasProp("_Name") else "?",
            "cnn_score":    float(p.get("CNNscore", 0)),
            "cnn_affinity": float(p.get("CNNaffinity", 0)),
            "vina_score":   float(p.get("minimizedAffinity", 0)),
        })

    df = pd.DataFrame(rows).sort_values("cnn_score", ascending=False)
    return df

top5pct = select_top_fraction(results, fraction=0.05)
convert_pdbqt_to_sdf([t["pdbqt"] for t in top5pct], "top5pct.sdf")
df = gnina_rescore_batch("receptor.pdbqt", "top5pct.sdf", "top5pct_gnina.sdf")
df.to_csv("rescored_hits.csv", index=False)
```

---

## VS Validation Metrics

Use these metrics when you have a set of known actives + decoys (e.g. DUD-E).

```python
import numpy as np
from sklearn.metrics import roc_auc_score

def enrichment_factor(scores: list, labels: list,
                      fraction: float = 0.01) -> float:
    """
    EF_x% = (actives in top x%) / (actives_total * x%)
    scores: docking scores (lower = better for Vina, reverse)
    labels: 1=active, 0=decoy
    """
    n = len(scores)
    n_actives = sum(labels)
    n_top = max(1, int(n * fraction))

    # Sort by score ascending (lower Vina score = better)
    sorted_labels = [l for _, l in sorted(zip(scores, labels))]
    actives_in_top = sum(sorted_labels[:n_top])

    random_expected = n_actives * fraction
    return actives_in_top / random_expected if random_expected > 0 else 0.0


def bedroc(scores: list, labels: list, alpha: float = 20.0) -> float:
    """
    BEDROC: Boltzmann-Enhanced Discrimination of ROC.
    alpha=20 weights early recognition heavily (default in literature).
    Range 0–1, random=0.5 approximately.
    """
    n = len(scores)
    n_actives = sum(labels)
    if n_actives == 0 or n_actives == n:
        return float("nan")

    ra = n_actives / n
    # Sort by score ascending (Vina: lower is better)
    sorted_labels = [l for _, l in sorted(zip(scores, labels))]

    risum = sum(
        np.exp(-alpha * i / n) * label
        for i, label in enumerate(sorted_labels, 1)
    )
    factor = (ra * np.sinh(alpha / 2)) / (np.cosh(alpha / 2) - np.cosh(alpha / 2 - alpha * ra))
    bedroc_val = risum * factor / n + 1 / (1 - np.exp(alpha * (1 - ra)))

    return float(bedroc_val)


# Example evaluation
vina_scores = [r["score"] for r in results]
labels = [1 if r["name"] in known_actives else 0 for r in results]

auc = roc_auc_score([-s for s in vina_scores], labels)  # negate: higher=better
ef1 = enrichment_factor(vina_scores, labels, fraction=0.01)
ef5 = enrichment_factor(vina_scores, labels, fraction=0.05)
bdr = bedroc(vina_scores, labels, alpha=20.0)

print(f"AUC-ROC: {auc:.3f}")
print(f"EF1%:    {ef1:.1f}x")
print(f"EF5%:    {ef5:.1f}x")
print(f"BEDROC:  {bdr:.3f}")
```

**Interpretation:**
- AUC > 0.70 : acceptable enrichment
- EF1% > 10 : strong early enrichment (drug discovery relevant)
- BEDROC > 0.70 : very good (alpha=20, weights top ~5%)

---

## Step 8 — Diversity Selection of Final Hits

Avoid reporting 10 near-identical compounds from the same scaffold.

```python
import subprocess, json

def diverse_hits(hits_smi: str, n_diverse: int = 50,
                 fp: str = "morgan") -> list:
    """Select diverse subset from top hits using MaxMin."""
    result = subprocess.run(
        ["python", "/path/to/scripts/chem_diversity.py",
         "--input", hits_smi,
         "--n", str(n_diverse),
         "--fingerprint", fp],
        capture_output=True, text=True
    )
    return json.loads(result.stdout)["selected"]

# Get SMILES for top 200 Vina hits
top200 = results[:200]
with open("top200.smi", "w") as f:
    for r in top200:
        f.write(f"{r['smiles']} {r['name']}\n")

diverse_50 = diverse_hits("top200.smi", n_diverse=50)
print(f"Selected {len(diverse_50)} diverse hits for synthesis")
```

---

## Common Issues in Virtual Screening

| Issue | Effect | Fix |
|-------|--------|-----|
| No 3D conformer | Vina error or flat molecule | Use `--gen3d best` (slower but correct) |
| Protonation not set | Wrong tautomer docked | Run at pH 7.4 via `obabel -p 7.4` |
| Box too small | Ligand extends outside → bad score | Add ≥ 8 Å padding around reference |
| Same scaffold floods results | Apparent diversity illusion | MaxMin selection (step 8) |
| Metal binder series | PDBQT atom types wrong | Use AutoDock Tools for metal-containing systems |
| Covalent warheads in library | Non-covalent docking invalid | Flag cysteine/serine reactive groups; use CovDock |
