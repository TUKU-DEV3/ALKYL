# Docking — AutoDock Vina & Gnina

Docking engines: AutoDock Vina 1.2 (empirical) and Gnina (CNN scoring).

---

## AutoDock Vina 1.2

### Python API (vina package)

```python
from vina import Vina

v = Vina(sf_name='vina')   # scoring function: 'vina' or 'vinardo'

v.set_receptor('receptor.pdbqt')
v.set_ligand_from_file('ligand.pdbqt')

v.compute_vina_maps(
    center=[10.5, -2.3, 14.1],
    box_size=[20, 20, 20]
)

# Dock
v.dock(exhaustiveness=16, n_poses=9)

# Results
energies = v.energies(n_poses=9)
print("Best score:", energies[0][0], "kcal/mol")

v.write_poses('docked.pdbqt', n_poses=9, overwrite=True)
```

### CLI (batch use)

```bash
# Single ligand
vina \
  --receptor receptor.pdbqt \
  --ligand    ligand.pdbqt \
  --config    box.conf \
  --exhaustiveness 16 \
  --num_modes 9 \
  --out       docked.pdbqt \
  --log       vina.log

# Config file (box.conf)
# center_x = 10.5
# center_y = -2.3
# center_z = 14.1
# size_x   = 20
# size_y   = 20
# size_z   = 20
```

### Batch docking — subprocess wrapper

```python
import subprocess
import json
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

def dock_one(receptor_pdbqt: str, ligand_pdbqt: str,
             box: dict, out_dir: str,
             exhaustiveness: int = 8) -> dict:
    """Dock a single ligand. Returns {name, score, out_file}."""
    name = Path(ligand_pdbqt).stem
    out_pdbqt = str(Path(out_dir) / f"{name}_docked.pdbqt")

    cmd = [
        "vina",
        "--receptor", receptor_pdbqt,
        "--ligand",   ligand_pdbqt,
        "--out",      out_pdbqt,
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", "1",    # only best pose for VS
    ]
    for k, v in box.items():
        cmd += [f"--{k}", str(v)]

    result = subprocess.run(cmd, capture_output=True, text=True)

    score = float("inf")
    for line in result.stdout.splitlines():
        stripped = line.strip()
        if stripped.startswith("1 "):
            try:
                score = float(stripped.split()[1])
            except (IndexError, ValueError):
                pass
            break

    return {"name": name, "score": score, "out": out_pdbqt,
            "success": result.returncode == 0}


def batch_dock(receptor_pdbqt: str, ligand_dir: str,
               box: dict, out_dir: str,
               exhaustiveness: int = 8, workers: int = 4) -> list:
    """Batch dock all PDBQT files in ligand_dir."""
    Path(out_dir).mkdir(exist_ok=True)
    ligands = list(Path(ligand_dir).glob("*.pdbqt"))

    results = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futures = {
            ex.submit(dock_one, receptor_pdbqt, str(lig),
                      box, out_dir, exhaustiveness): lig
            for lig in ligands
        }
        for fut in as_completed(futures):
            results.append(fut.result())

    results.sort(key=lambda x: x["score"])
    return results

# Usage
box = {"center_x": 10.5, "center_y": -2.3, "center_z": 14.1,
       "size_x": 20, "size_y": 20, "size_z": 20}

hits = batch_dock("receptor.pdbqt", "ligands/", box, "docked/",
                  exhaustiveness=8, workers=8)

# Save results
import csv
with open("screening_results.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["name", "score", "out", "success"])
    writer.writeheader()
    writer.writerows(hits)

print(f"Top 10 hits:")
for h in hits[:10]:
    print(f"  {h['name']}: {h['score']:.2f} kcal/mol")
```

### Key Parameters

| Parameter | Default | Notes |
|-----------|---------|-------|
| `exhaustiveness` | 8 | Higher = better sampling, slower. 16–32 for binding mode, 8 for VS |
| `num_modes` | 9 | Number of poses to output |
| `energy_range` | 3 | Max energy difference from best pose (kcal/mol) |
| `sf_name` | `vina` | `vina` (standard) or `vinardo` (alternative) |
| `min_rmsd` | 1.0 | Minimum RMSD between reported poses |

**Score thresholds (rough guide):**
- < −9 kcal/mol : strong predicted binder (for drug-like molecules)
- −7 to −9 : moderate, promising hits
- −5 to −7 : weak, borderline
- > −5 : typically noise / non-binder

> These are relative within a target/system. Never compare Vina scores across different receptors.

---

## Gnina — CNN Scoring

Gnina re-scores or docks using a convolutional neural network trained on PDBbind.

### Installation

```bash
# Prebuilt Linux binary (recommended)
wget https://github.com/gnina/gnina/releases/latest/download/gnina
chmod +x gnina
sudo mv gnina /usr/local/bin/

# Docker
docker pull gnina/gnina
```

### Basic Gnina docking

```bash
# Dock with autobox from reference ligand
gnina \
  --receptor receptor.pdbqt \
  --ligand    ligand.sdf \
  --autobox_ligand ref_ligand.sdf \
  --autobox_add 4 \
  --exhaustiveness 16 \
  --num_modes 9 \
  --out docked_gnina.sdf \
  --log gnina.log

# Output SDF contains:
#   minimizedAffinity  (Vina-like score, kcal/mol)
#   CNNscore           (0–1, higher is better — CNN predicted activity)
#   CNNaffinity        (CNN-predicted ΔG, kcal/mol)
```

### Gnina for rescoring Vina poses

```bash
# Score existing poses without redocking
gnina \
  --receptor receptor.pdbqt \
  --ligand    vina_poses.sdf \
  --score_only \
  --out gnina_rescored.sdf

# Then sort by CNNscore (better discriminator than minimizedAffinity)
```

```python
from rdkit import Chem
import subprocess

def gnina_rescore(receptor: str, poses_sdf: str,
                  out_sdf: str) -> list:
    """Rescore Vina poses with Gnina CNN."""
    subprocess.run([
        "gnina",
        "--receptor", receptor,
        "--ligand",   poses_sdf,
        "--score_only",
        "--out",      out_sdf
    ], check=True, capture_output=True)

    results = []
    supp = Chem.SDMolSupplier(out_sdf, removeHs=False)
    for mol in supp:
        if mol is None:
            continue
        props = mol.GetPropsAsDict()
        results.append({
            "name":        mol.GetProp("_Name") if mol.HasProp("_Name") else "?",
            "cnn_score":   float(props.get("CNNscore", 0.0)),
            "cnn_affinity": float(props.get("CNNaffinity", 0.0)),
            "vina_score":  float(props.get("minimizedAffinity", 0.0)),
        })

    results.sort(key=lambda x: -x["cnn_score"])
    return results

hits = gnina_rescore("receptor.pdbqt", "vina_poses.sdf", "gnina_out.sdf")
for h in hits[:5]:
    print(f"{h['name']}: CNN={h['cnn_score']:.3f}, "
          f"CNNaff={h['cnn_affinity']:.2f}, Vina={h['vina_score']:.2f}")
```

### Gnina vs Vina — When to Use Which

| Criterion | Vina | Gnina CNN |
|-----------|------|-----------|
| Speed | ★★★★★ | ★★★ |
| Pose prediction accuracy | ★★★ | ★★★★ |
| Virtual screening enrichment | ★★★ | ★★★★ |
| Out-of-domain molecules | ★★★ | ★★ (training set bias) |
| Fragment-sized ligands | ★★★ | ★★ (trained on drug-like) |
| No GPU required | ✓ | ✓ (CPU mode, slower) |

**Recommended workflow:**
1. Vina for large VS (speed, no bias)
2. Gnina CNN rescore top 1–5% from Vina
3. MM-GB/SA on final top hits (→ `free-energy` skill)

---

## Ligand Preparation for Docking

### meeko (recommended for Vina)

```python
# meeko: proper rotatable bond detection + PDBQT prep
import subprocess

def prep_ligand_meeko(sdf_file: str, out_pdbqt: str) -> None:
    """Prepare ligand PDBQT with meeko."""
    subprocess.run([
        "mk_prepare_ligand.py",
        "-i", sdf_file,
        "-o", out_pdbqt,
    ], check=True)

# Batch
from pathlib import Path

def prep_library_meeko(sdf_dir: str, pdbqt_dir: str) -> None:
    Path(pdbqt_dir).mkdir(exist_ok=True)
    for sdf in Path(sdf_dir).glob("*.sdf"):
        out = Path(pdbqt_dir) / (sdf.stem + ".pdbqt")
        try:
            prep_ligand_meeko(str(sdf), str(out))
        except subprocess.CalledProcessError:
            print(f"Failed: {sdf.name}")
```

### OpenBabel fallback

```bash
# Single SDF → PDBQT
obabel ligand.sdf -O ligand.pdbqt -h --gen3d

# Batch: SMILES file → individual PDBQT files
obabel library.smi -O ligands/lig.pdbqt --gen3d -h -m
# -m : split into separate files (lig0001.pdbqt, lig0002.pdbqt, ...)
```

---

## Scoring Function Theory (brief)

**Vina scoring function:**
```
ΔG = Σ w_i · f_i(d)
```
Terms: steric (Gauss1, Gauss2, repulsion), hydrophobic, hydrogen bond

**Gnina CNN scoring:**
- Input: 3D grid around binding site, atom type channels
- Architecture: ResNet-like, trained on PDBbind v2016
- `CNNscore`: binary active/inactive classifier (AUROC ~0.9 on DUD-E)
- `CNNaffinity`: regression on pKi/pKd

**Vinardo (alternative empirical):**
- More transferable across targets than Vina
- Use with `sf_name='vinardo'` in vina Python API
