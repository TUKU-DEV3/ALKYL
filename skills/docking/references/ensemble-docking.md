# Docking — Ensemble Docking & Protein Flexibility

Account for receptor conformational flexibility by docking against multiple snapshots.

---

## Why Ensemble Docking

Single-structure docking ignores protein dynamics. Reasons to use ensemble:
- **Induced fit**: ligand reshapes binding site (common in kinases, GPCRs)
- **Cryptic pockets**: only open in some conformations
- **Flexible loops**: gating residues (Asp-Phe-Gly loop in kinases)
- **Allosteric sites**: require conformer from apo/allosteric-bound state

**When to skip:** rigid lock-and-key binding sites (highly buried, rigid scaffolds like beta-sheets). Ensemble docking adds noise if protein is genuinely rigid.

---

## Step 1 — Generate Conformational Ensemble from MD

```python
# Using MDAnalysis — extract snapshots from MD trajectory
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
import numpy as np
from pathlib import Path

def extract_md_ensemble(topology: str, trajectory: str,
                        selection: str = "protein and backbone",
                        n_clusters: int = 10,
                        start: int = 0, stop: int = None,
                        step: int = 10,
                        out_dir: str = "ensemble") -> list:
    """
    Extract representative conformers from MD trajectory via clustering.
    Returns list of PDB file paths.
    """
    Path(out_dir).mkdir(exist_ok=True)

    u = mda.Universe(topology, trajectory)
    protein = u.select_atoms("protein")

    # Align to first frame
    aligner = align.AlignTraj(u, u, select=selection,
                               in_memory=True)
    aligner.run()

    # Collect Cα coordinates for clustering
    ca_sel = u.select_atoms("protein and name CA")
    coords = []
    frames = []

    for ts in u.trajectory[start:stop:step]:
        coords.append(ca_sel.positions.flatten())
        frames.append(ts.frame)

    coords = np.array(coords)

    # K-Means clustering
    from sklearn.cluster import KMeans
    km = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
    labels = km.fit_predict(coords)

    # Pick frame closest to each cluster center
    pdb_files = []
    for cluster_id in range(n_clusters):
        cluster_mask = labels == cluster_id
        cluster_coords = coords[cluster_mask]
        cluster_frames = [frames[i] for i, m in enumerate(cluster_mask) if m]

        center = km.cluster_centers_[cluster_id]
        dists = np.linalg.norm(cluster_coords - center, axis=1)
        best_idx = np.argmin(dists)
        best_frame = cluster_frames[best_idx]

        u.trajectory[best_frame]
        out_pdb = str(Path(out_dir) / f"conf_{cluster_id:02d}.pdb")
        protein.write(out_pdb)
        pdb_files.append(out_pdb)
        print(f"Cluster {cluster_id}: frame {best_frame} → {out_pdb}")

    return pdb_files

# Generate 10 representative conformers
conformers = extract_md_ensemble(
    "system.gro", "md_prod.xtc",
    n_clusters=10, step=10, out_dir="ensemble"
)
```

---

## Step 2 — Prepare Each Conformer for Docking

```python
import subprocess
from pathlib import Path
from pdbfixer import PDBFixer
from openmm.app import PDBFile

def prep_conformer(pdb_in: str, pdb_out: str, ph: float = 7.4) -> None:
    """Fix and protonate a single conformer."""
    fixer = PDBFixer(filename=pdb_in)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(ph)
    with open(pdb_out, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)


def prep_ensemble(conformers: list, out_dir: str) -> list:
    """Prep all conformers → PDBQT."""
    Path(out_dir).mkdir(exist_ok=True)
    pdbqt_files = []

    for pdb in conformers:
        stem = Path(pdb).stem
        fixed = str(Path(out_dir) / f"{stem}_fixed.pdb")
        pdbqt  = str(Path(out_dir) / f"{stem}.pdbqt")

        prep_conformer(pdb, fixed)
        subprocess.run([
            "obabel", fixed, "-O", pdbqt, "-xr"
        ], check=True, capture_output=True)
        pdbqt_files.append(pdbqt)

    return pdbqt_files

pdbqt_ensemble = prep_ensemble(conformers, "ensemble_pdbqt")
```

---

## Step 3 — Dock Against Each Conformer

```python
import subprocess, json
from concurrent.futures import ProcessPoolExecutor, as_completed

def dock_to_conformer(receptor_pdbqt: str, ligand_pdbqt: str,
                      box: dict, exhaustiveness: int = 16) -> dict:
    """Dock one ligand against one conformer."""
    import tempfile
    with tempfile.NamedTemporaryFile(suffix=".pdbqt",
                                     delete=False) as tmp:
        out = tmp.name

    cmd = ["vina",
           "--receptor", receptor_pdbqt,
           "--ligand",   ligand_pdbqt,
           "--out",      out,
           "--exhaustiveness", str(exhaustiveness),
           "--num_modes", "3"]
    for k, v in box.items():
        cmd += [f"--{k}", str(v)]

    proc = subprocess.run(cmd, capture_output=True, text=True)

    scores = []
    for line in proc.stdout.splitlines():
        s = line.strip()
        if s and s[0].isdigit():
            try:
                scores.append(float(s.split()[1]))
            except (IndexError, ValueError):
                pass

    return {
        "receptor":  receptor_pdbqt,
        "best_score": scores[0] if scores else float("inf"),
        "all_scores": scores,
        "out_pdbqt":  out,
    }


def ensemble_dock_ligand(ligand_pdbqt: str,
                          receptor_pdbqts: list,
                          box: dict,
                          exhaustiveness: int = 16,
                          workers: int = 4) -> dict:
    """
    Dock one ligand against all conformers.
    Aggregation: min (most favorable), mean, % below threshold.
    """
    args = [(r, ligand_pdbqt, box, exhaustiveness)
            for r in receptor_pdbqts]

    all_results = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futures = {
            ex.submit(dock_to_conformer, r, ligand_pdbqt, box, exhaustiveness): r
            for r in receptor_pdbqts
        }
        for fut in as_completed(futures):
            all_results.append(fut.result())

    best_scores = [r["best_score"] for r in all_results
                   if r["best_score"] < 0]

    if not best_scores:
        return {"min": float("inf"), "mean": float("inf"),
                "pct_below_minus7": 0.0, "results": all_results}

    return {
        "min":  min(best_scores),
        "mean": sum(best_scores) / len(best_scores),
        "pct_below_minus7": sum(1 for s in best_scores if s < -7) / len(best_scores),
        "results": all_results,
    }
```

---

## Step 4 — Aggregate Ensemble Scores

Multiple strategies for combining scores across conformers:

```python
import pandas as pd
import numpy as np

def aggregate_ensemble_results(ligand_results: list,
                                strategy: str = "min") -> pd.DataFrame:
    """
    ligand_results: list of {name, ensemble_scores: list[float]}
    strategy: 'min', 'mean', 'boltzmann', 'pct_active'
    """
    rows = []
    for lr in ligand_results:
        scores = [s for s in lr["ensemble_scores"] if s < 0]
        if not scores:
            rows.append({"name": lr["name"], "score": 0.0})
            continue

        if strategy == "min":
            agg = min(scores)
        elif strategy == "mean":
            agg = np.mean(scores)
        elif strategy == "boltzmann":
            # Boltzmann-weighted average (T=300K, kBT ~ 0.593 kcal/mol)
            kBT = 0.593
            weights = np.exp(-np.array(scores) / kBT)
            agg = float(np.dot(weights, scores) / weights.sum())
        elif strategy == "pct_active":
            threshold = -7.0
            agg = sum(1 for s in scores if s <= threshold) / len(lr["ensemble_scores"])

        rows.append({"name": lr["name"], "score": agg,
                     "n_conformers": len(lr["ensemble_scores"])})

    df = pd.DataFrame(rows).sort_values("score")
    return df

# Boltzmann weighting is theoretically most rigorous
# 'min' is most practical and widely used
```

**Strategy comparison:**

| Strategy | Pros | Cons |
|----------|------|------|
| `min` | Simple, reproducible, highlights strong binders | Prone to outliers / artefacts |
| `mean` | Smooth, less sensitive to outliers | Can dilute true binders |
| `boltzmann` | Thermodynamically grounded | Sensitive to score units calibration |
| `pct_active` | Robust, threshold-based | Arbitrary threshold |

---

## Alternative: Normal Mode Analysis Ensemble

For fast ensemble without MD:

```python
from ase import io
from ase.vibrations import Vibrations
from ase.calculators.emt import EMT  # replace with real calc
import numpy as np

def nma_ensemble(structure_file: str, n_modes: int = 6,
                 scale: float = 1.5, n_conformers: int = 10,
                 out_dir: str = "nma_ensemble") -> list:
    """
    Generate conformers by displacing along low-frequency normal modes.
    """
    from pathlib import Path
    Path(out_dir).mkdir(exist_ok=True)

    atoms = io.read(structure_file)
    # Attach real calculator here (e.g. xTB for proteins is impractical;
    # use this pattern for small systems or ligands only)
    atoms.calc = EMT()

    vib = Vibrations(atoms)
    vib.run()

    pdb_files = [structure_file]  # include original

    for mode_idx in range(6, 6 + n_modes):  # skip translations/rotations
        for direction in [+1, -1]:
            displaced = atoms.copy()
            mode_vec = np.real(vib.get_mode(mode_idx))
            mode_vec /= np.linalg.norm(mode_vec)
            displaced.positions += direction * scale * mode_vec

            out = str(Path(out_dir) / f"nma_mode{mode_idx}_{'+' if direction > 0 else '-'}.pdb")
            io.write(out, displaced)
            pdb_files.append(out)

    return pdb_files

# Note: NMA ensemble is better suited for ligands or small proteins
# For full receptor flexibility, MD-based ensemble is preferred
```

---

## Relaxed Complex Scheme (RCS)

Standard protocol for ensemble docking in drug discovery:

```
1. Run MD of apo receptor (or holo with reference ligand)
   → 10–50 ns, step 10 ps → 1000–5000 frames

2. Cluster trajectory → 10–20 representative structures
   (MDAnalysis K-Means on Cα, binding site residues only for clustering)

3. Prepare each conformer (pdbfixer → PDBQT)

4. For each conformer: dock full library (Vina, exhaustiveness=8)

5. Aggregate: take min score per ligand across conformers

6. Gnina CNN rescore top 5% of aggregated results
   (pick best CNNscore pose from all conformers)

7. ProLIF interaction filter:
   mandatory: retain key pharmacophoric interactions from known actives
```

```python
# Binding site-focused clustering (more relevant than full Cα)
def cluster_binding_site(topology: str, trajectory: str,
                          pocket_residues: list,  # e.g. [83, 145, 166]
                          n_clusters: int = 10) -> list:
    u = mda.Universe(topology, trajectory)
    resid_sel = " or ".join(f"resid {r}" for r in pocket_residues)
    pocket = u.select_atoms(f"protein and ({resid_sel}) and backbone")

    coords = []
    for ts in u.trajectory[::10]:
        coords.append(pocket.positions.flatten())

    from sklearn.cluster import KMeans
    km = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
    labels = km.fit_predict(np.array(coords))

    # Return frame indices of cluster centers
    centers = []
    for cid in range(n_clusters):
        mask = labels == cid
        cluster_coords = np.array(coords)[mask]
        dists = np.linalg.norm(cluster_coords - km.cluster_centers_[cid], axis=1)
        best = np.where(mask)[0][np.argmin(dists)]
        centers.append(best * 10)  # multiply by step

    return centers
```

---

## Induced Fit Docking (Conceptual)

Full induced fit (IFD) requires iterative docking + side chain refinement. Without commercial tools (Glide IFD), approximate with:

1. **Soft-core Vina**: increase `--weight_gauss1` to reduce steric penalty
2. **Gnina flexible residues**: `--flexres A:Leu83,Phe89` (up to 5 residues)
3. **ADFR Suite**: AutoDockFR supports flexible receptor side chains

```bash
# Gnina with flexible residues
gnina \
  --receptor receptor_rigid.pdbqt \
  --flex     receptor_flex.pdbqt \
  --ligand   ligand.sdf \
  --autobox_ligand ref_ligand.sdf \
  --exhaustiveness 32 \
  --out docked_flex.sdf

# Generate flex PDBQT with AutoDockTools:
# pythonsh prepare_flexreceptor4.py \
#   -r receptor.pdbqt \
#   -s receptor:A:LEU83_PHE89_ASP145
```
