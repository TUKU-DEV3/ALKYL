# Fragment Docking

## Special Challenges

Standard docking parameters are optimized for drug-like molecules (MW 350–500).
For fragments (MW < 300, HAC 8–17):
- **Underdocking**: too few poses generated for small ligand
- **Scoring noise**: small interaction surface → score variance dominates
- **False positives**: fragments fit anywhere (low selectivity of score)
- **Ligand desolation**: scoring functions underestimate desolvation cost for polar frags

## Vina Fragment Settings

```python
from vina import Vina

def dock_fragment(smiles, receptor_pdbqt, center, box_size,
                  exhaustiveness=32, n_poses=20):
    """
    Dock a fragment with higher exhaustiveness to compensate small size.
    exhaustiveness=32 (vs 8 for drug-like) for better sampling.
    """
    from meeko import MoleculePreparation
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    AllChem.MMFFOptimizeMolecule(mol)

    prep = MoleculePreparation()
    mol_setups = prep.prepare(mol)
    pdbqt_str = mol_setups[0].write_pdbqt_string()

    v = Vina(sf_name='vina', cpu=4, verbosity=0)
    v.set_receptor(receptor_pdbqt)
    v.compute_vina_maps(center=center, box_size=box_size)
    v.set_ligand_from_string(pdbqt_str)
    v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)

    energies = v.energies(n_poses=n_poses)
    return v, energies   # (Vina object, array of scores)
```

## Smina for Fragment Docking (rescoring)

Smina supports custom scoring functions and is useful for re-ranking fragment poses:

```bash
# Smina fragment docking with larger exhaustiveness
smina --receptor receptor.pdbqt \
      --ligand fragment.pdbqt \
      --center_x X --center_y Y --center_z Z \
      --size_x 20 --size_y 20 --size_z 20 \
      --exhaustiveness 32 \
      --num_modes 20 \
      --scoring vinardo \   # Vinardo often better for fragments
      --out fragment_poses.sdf
```

## Gnina CNN Rescoring for Fragments

```python
import subprocess
from pathlib import Path

def gnina_rescore_fragments(receptor_pdb, poses_sdf, n_poses=5):
    """CNN rescoring of fragment poses — often improves rank."""
    out_sdf = str(Path(poses_sdf).with_suffix('')) + '_gnina.sdf'
    cmd = [
        'gnina',
        '--receptor', receptor_pdb,
        '--ligand', poses_sdf,
        '--score_only',
        '--out', out_sdf,
        '--num_modes', str(n_poses),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)

    # Parse CNN scores from output
    scores = []
    for line in result.stdout.split('\n'):
        if 'CNNscore' in line:
            try:
                scores.append(float(line.split()[-1]))
            except ValueError:
                pass
    return scores, out_sdf
```

## Fragment Docking Pipeline (Batch)

```python
import pandas as pd
from pathlib import Path
import concurrent.futures

def batch_fragment_dock(fragment_df, receptor_pdbqt, center, box_size,
                         out_dir='fragment_docking', n_workers=4):
    """
    Dock all fragments in parallel.
    Returns DataFrame with best Vina score per fragment.
    """
    Path(out_dir).mkdir(exist_ok=True)
    results = []

    def dock_one(row):
        try:
            v, energies = dock_fragment(
                row.smiles, receptor_pdbqt, center, box_size,
                exhaustiveness=32, n_poses=20
            )
            best_score = energies[0][0]   # best pose Vina score
            return {
                'smiles': row.smiles,
                'vina_score': best_score,
                'status': 'ok'
            }
        except Exception as e:
            return {'smiles': row.smiles, 'vina_score': 0.0, 'status': str(e)}

    with concurrent.futures.ThreadPoolExecutor(max_workers=n_workers) as pool:
        futures = {pool.submit(dock_one, row): row
                   for _, row in fragment_df.iterrows()}
        for future in concurrent.futures.as_completed(futures):
            results.append(future.result())

    df = pd.DataFrame(results)
    df = df.sort_values('vina_score')   # most negative = best
    return df
```

## Pose Clustering (Remove Redundant Poses)

```python
from rdkit import Chem
from rdkit.Chem import rdMolAlign
import numpy as np

def cluster_poses(poses_sdf, rmsd_cutoff=1.5):
    """
    Cluster docking poses by RMSD.
    Returns representative pose (lowest energy) per cluster.
    """
    supplier = list(Chem.SDMolSupplier(poses_sdf, removeHs=True))
    n = len(supplier)
    if n == 0:
        return []

    # Compute pairwise RMSD
    rmsds = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            try:
                rmsd = rdMolAlign.CalcRMS(supplier[i], supplier[j])
            except Exception:
                rmsd = 99.0
            rmsds[i, j] = rmsds[j, i] = rmsd

    # Greedy clustering
    assigned = [False] * n
    clusters = []
    for i in range(n):   # assumes poses sorted by energy (best first)
        if assigned[i]:
            continue
        cluster = [i]
        assigned[i] = True
        for j in range(i+1, n):
            if not assigned[j] and rmsds[i, j] <= rmsd_cutoff:
                cluster.append(j)
                assigned[j] = True
        clusters.append(cluster)

    # Representative = first (lowest energy) in each cluster
    representatives = [cluster[0] for cluster in clusters]
    return [supplier[i] for i in representatives]
```

## Hotspot Validation

After docking, confirm fragment hits via predicted interaction with known hotspot:

```python
import prolif as plf
from rdkit import Chem
import MDAnalysis as mda

def validate_hotspot_interactions(receptor_pdb, poses_sdf,
                                   hotspot_residues=None):
    """
    Check if docked fragment contacts known hotspot residues.
    hotspot_residues: list of residue numbers (e.g. [145, 168, 202])
    """
    u = mda.Universe(receptor_pdb)
    ligands = Chem.SDMolSupplier(poses_sdf, removeHs=False)

    results = []
    for mol in ligands:
        if mol is None:
            continue
        fp = plf.Fingerprint()
        fp.run_from_mol(mol, u.select_atoms('protein'))
        ifp = fp.ifp[0]   # interaction fingerprint dict

        contacted = set()
        for (lig_atom, prot_res), interactions in ifp.items():
            resi = int(prot_res.split(':')[1].split('.')[0])
            contacted.add(resi)

        if hotspot_residues:
            overlap = contacted & set(hotspot_residues)
            results.append({
                'mol': mol,
                'contacted_residues': contacted,
                'hotspot_contacts': overlap,
                'hotspot_hit': len(overlap) > 0
            })
        else:
            results.append({'mol': mol, 'contacted_residues': contacted})

    return results
```

## Score Interpretation for Fragments

**Fragment Vina scores are systematically less negative than drug-like molecules.**
Do not compare across different MW classes.

Normalize by HAC:
```python
def vina_le(vina_score, n_heavy_atoms):
    """
    Vina Ligand Efficiency (kcal/mol per heavy atom).
    Vina score is already in kcal/mol.
    """
    return vina_score / n_heavy_atoms

# Fragment hit criteria:
# Vina_LE < -0.3 kcal/mol/HA → reasonable fragment hit
# Vina_LE < -0.4 kcal/mol/HA → good fragment hit
```

## Docking vs Experimental Fragment Screening

| Method | False positive rate | Throughput | Cost |
|--------|-------------------|------------|------|
| Vina docking | High (20–40%) | 10k/day | Low |
| SPR (Biacore) | Medium (5–15%) | 500/week | $$ |
| NMR WaterLOGSY | Low (5–10%) | 200/week | $$$ |
| X-ray (XChem) | Very low (<5%) | 1000/beamtime | $$$$ |

Docking best used for:
1. Pre-filtering large libraries before expensive biophysical assays
2. Prioritizing poses for elaboration after experimental hit confirmation
3. Guided growing: propose elaboration directions based on pose
