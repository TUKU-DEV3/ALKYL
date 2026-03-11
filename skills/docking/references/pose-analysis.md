# Docking — Pose Analysis & Interaction Fingerprints

After docking: selecting the best pose, extracting interactions, clustering, visualizing.

---

## ProLIF — Protein-Ligand Interaction Fingerprints

ProLIF encodes protein-ligand contacts as a binary fingerprint, enabling comparison and filtering of docking poses.

```python
import MDAnalysis as mda
import prolif as plf
import pandas as pd

def analyze_interactions(receptor_pdb: str, ligand_sdf: str,
                          protein_selection: str = "protein") -> pd.DataFrame:
    """
    Compute interaction fingerprint for one pose.
    Returns DataFrame with interaction types per residue.
    """
    u = mda.Universe(receptor_pdb)
    lig_u = mda.Universe(ligand_sdf)

    protein = u.select_atoms(protein_selection)
    ligand  = lig_u.select_atoms("all")

    fp = plf.Fingerprint([
        "HBDonor", "HBAcceptor",
        "Hydrophobic",
        "PiStacking", "PiCation",
        "Anionic", "Cationic",
        "XBAcceptor",           # halogen bond acceptor
        "MetalDonor",
    ])
    fp.run_from_iterable([ligand], protein)

    df = fp.to_dataframe()
    return df


# Usage
df = analyze_interactions("receptor.pdb", "best_pose.sdf")
print(df.T[df.any(axis=0)])  # show only active interactions
```

### Batch interaction fingerprints for VS hits

```python
from pathlib import Path
import prolif as plf
import MDAnalysis as mda
import numpy as np
import pandas as pd

def batch_prolif(receptor_pdb: str, poses_sdf_dir: str,
                 mandatory_interactions: list = None) -> pd.DataFrame:
    """
    Compute IFP for all poses. Optional: filter by mandatory interactions.
    mandatory_interactions: e.g. ['HBDonor_ASP123', 'Hydrophobic_LEU83']
    """
    u = mda.Universe(receptor_pdb)
    protein = u.select_atoms("protein")

    fp = plf.Fingerprint()
    results = []

    for sdf in Path(poses_sdf_dir).glob("*.sdf"):
        try:
            lig_u = mda.Universe(str(sdf))
            ligand = lig_u.select_atoms("all")
            fp.run_from_iterable([ligand], protein)
            row = fp.to_dataframe().iloc[0].to_dict()
            row["name"] = sdf.stem
            results.append(row)
        except Exception as e:
            print(f"Failed {sdf.name}: {e}")

    df = pd.DataFrame(results).fillna(False)

    if mandatory_interactions:
        mask = pd.Series([True] * len(df))
        for interaction in mandatory_interactions:
            if interaction in df.columns:
                mask &= df[interaction].astype(bool)
        df = df[mask]
        print(f"After mandatory filter: {len(df)} poses")

    return df

df = batch_prolif("receptor.pdb", "vs_output/",
                  mandatory_interactions=["HBAcceptor_ASP145A"])
```

---

## Interaction Analysis — Manual with RDKit + MDAnalysis

```python
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

def get_hbonds(receptor_pdb: str, ligand_sdf: str,
               cutoff_dist: float = 3.5,
               cutoff_angle: float = 120.0) -> list:
    """
    Simple geometric H-bond detection.
    Returns list of (donor_res, acceptor_res, distance, angle).
    """
    import MDAnalysis as mda
    from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis

    combined = mda.Universe(receptor_pdb)
    # Load ligand as separate universe, merge
    lig = mda.Universe(ligand_sdf)

    # MDAnalysis HBondAnalysis requires topology with H
    hbond = HydrogenBondAnalysis(
        universe=combined,
        donors_sel="protein and (name N* O*)",
        hydrogens_sel="protein and name H*",
        acceptors_sel="resname LIG and (name N* O* F*)",
        d_a_cutoff=cutoff_dist,
        d_h_a_angle_cutoff=cutoff_angle,
    )
    hbond.run()
    return hbond.results.hbonds


def get_hydrophobic_contacts(receptor_pdb: str, ligand_sdf: str,
                              cutoff: float = 4.5) -> list:
    """
    Identify hydrophobic contacts (C–C within cutoff).
    """
    import MDAnalysis as mda
    u = mda.Universe(receptor_pdb)

    lig = Chem.SDMolSupplier(ligand_sdf, removeHs=False)[0]
    if lig is None:
        return []

    conf = lig.GetConformer()
    lig_carbons = [
        np.array(conf.GetAtomPosition(i))
        for i in range(lig.GetNumAtoms())
        if lig.GetAtomWithIdx(i).GetSymbol() == "C"
    ]

    protein_carbons = u.select_atoms("protein and name C*")
    contacts = []
    for pc in protein_carbons:
        for lc in lig_carbons:
            dist = np.linalg.norm(pc.position - lc)
            if dist <= cutoff:
                contacts.append({
                    "protein_atom": f"{pc.resname}{pc.resid}_{pc.name}",
                    "dist": round(dist, 2)
                })

    # Deduplicate by residue
    seen = set()
    unique = []
    for c in contacts:
        res = c["protein_atom"].split("_")[0]
        if res not in seen:
            seen.add(res)
            unique.append(c)

    return sorted(unique, key=lambda x: x["dist"])
```

---

## Pose Clustering

Cluster multiple poses from one ligand (or similar ligands from VS) by RMSD.

```python
from rdkit import Chem
from rdkit.Chem import rdMolAlign
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster

def cluster_poses_rmsd(poses_sdf: str,
                       rmsd_cutoff: float = 2.0) -> dict:
    """
    Cluster docking poses by RMSD. Return cluster centers.
    """
    mols = [m for m in Chem.SDMolSupplier(poses_sdf, removeHs=True)
            if m is not None]
    n = len(mols)

    # RMSD matrix
    rmsd_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            try:
                r = rdMolAlign.CalcRMS(mols[i], mols[j])
            except Exception:
                r = 999.9
            rmsd_matrix[i, j] = rmsd_matrix[j, i] = r

    # Hierarchical clustering
    from scipy.spatial.distance import squareform
    condensed = squareform(rmsd_matrix)
    Z = linkage(condensed, method="complete")
    labels = fcluster(Z, t=rmsd_cutoff, criterion="distance")

    clusters = {}
    for i, label in enumerate(labels):
        clusters.setdefault(label, []).append(i)

    # Pick representative (lowest energy if scores available)
    centers = {label: idxs[0] for label, idxs in clusters.items()}
    print(f"{n} poses → {len(clusters)} clusters (cutoff {rmsd_cutoff} Å)")
    return {"clusters": clusters, "centers": centers, "mols": mols}


result = cluster_poses_rmsd("docked_poses.sdf", rmsd_cutoff=2.0)

# Write cluster centers
writer = Chem.SDWriter("clustered_poses.sdf")
for label, center_idx in sorted(result["centers"].items()):
    mol = result["mols"][center_idx]
    mol.SetProp("cluster_id", str(label))
    mol.SetProp("cluster_size", str(len(result["clusters"][label])))
    writer.write(mol)
writer.close()
```

---

## py3Dmol — Inline 3D Visualization

```python
import py3Dmol

def view_pose(receptor_pdb: str, ligand_sdf: str,
              width: int = 800, height: int = 500) -> None:
    """Display protein-ligand pose in Jupyter."""
    view = py3Dmol.view(width=width, height=height)

    # Protein
    with open(receptor_pdb) as f:
        pdb_str = f.read()
    view.addModel(pdb_str, "pdb")
    view.setStyle({"model": 0}, {
        "cartoon": {"color": "lightgrey", "opacity": 0.7}
    })

    # Binding site surface
    view.addSurface(py3Dmol.VDW, {
        "opacity": 0.3, "color": "white"
    }, {"model": 0})

    # Ligand
    with open(ligand_sdf) as f:
        sdf_str = f.read()
    view.addModel(sdf_str, "sdf")
    view.setStyle({"model": 1}, {
        "stick": {"colorscheme": "grayCarbon", "radius": 0.2},
        "sphere": {"colorscheme": "grayCarbon", "radius": 0.3}
    })

    view.zoomTo({"model": 1})
    view.show()


def view_multiple_poses(receptor_pdb: str, poses_sdf: str,
                        n_poses: int = 5) -> None:
    """Show N poses as sticks with different colors."""
    colors = ["cyan", "orange", "magenta", "lime", "yellow"]
    view = py3Dmol.view(width=800, height=500)

    with open(receptor_pdb) as f:
        view.addModel(f.read(), "pdb")
    view.setStyle({"model": 0}, {"cartoon": {"color": "lightgrey"}})

    from rdkit import Chem
    mols = [m for m in Chem.SDMolSupplier(poses_sdf) if m is not None]
    for i, mol in enumerate(mols[:n_poses]):
        sdf_block = Chem.MolToMolBlock(mol)
        view.addModel(sdf_block, "sdf")
        view.setStyle({"model": i + 1}, {
            "stick": {"color": colors[i % len(colors)], "radius": 0.2}
        })

    view.zoomTo()
    view.show()
```

---

## Interaction Summary Report

```python
def pose_report(receptor_pdb: str, ligand_sdf: str,
                vina_score: float = None) -> dict:
    """
    Generate a concise interaction report for one pose.
    """
    import MDAnalysis as mda
    import prolif as plf

    u  = mda.Universe(receptor_pdb)
    lu = mda.Universe(ligand_sdf)
    protein = u.select_atoms("protein")
    ligand  = lu.select_atoms("all")

    fp = plf.Fingerprint()
    fp.run_from_iterable([ligand], protein)
    df = fp.to_dataframe()

    active_cols = [col for col in df.columns if df[col].iloc[0]]

    interactions = {
        "hbond_donor":    [],
        "hbond_acceptor": [],
        "hydrophobic":    [],
        "pi_stacking":    [],
        "other":          [],
    }

    for col in active_cols:
        if isinstance(col, tuple):
            interaction_type = col[1] if len(col) > 1 else str(col)
            residue = str(col[0]) if len(col) > 0 else "?"
        else:
            parts = str(col).split("_")
            interaction_type = parts[-1] if parts else str(col)
            residue = "_".join(parts[:-1])

        t = interaction_type.lower()
        label = f"{residue} ({interaction_type})"
        if "donor" in t:
            interactions["hbond_donor"].append(label)
        elif "acceptor" in t and "hb" in t:
            interactions["hbond_acceptor"].append(label)
        elif "hydro" in t:
            interactions["hydrophobic"].append(label)
        elif "pi" in t:
            interactions["pi_stacking"].append(label)
        else:
            interactions["other"].append(label)

    report = {
        "vina_score": vina_score,
        "n_interactions": len(active_cols),
        **interactions,
    }

    print(f"Vina score: {vina_score} kcal/mol")
    print(f"Total interactions: {len(active_cols)}")
    for k, v in interactions.items():
        if v:
            print(f"  {k}: {', '.join(v)}")

    return report
```
