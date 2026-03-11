# Structure-Based Pharmacophore

Derive pharmacophore directly from a protein-ligand co-crystal structure. Interaction analysis reveals which features of the ligand contact the protein — these become the pharmacophore query.

---

## Workflow Overview

```
1. Load co-crystal (PDB) → extract ligand + protein
2. Detect protein-ligand interactions (H-bonds, hydrophobic, ionic, π-π)
3. Map interactions to pharmacophore features on ligand atoms
4. Set tolerance spheres + exclusion volumes from protein pocket
5. Validate: re-dock known actives and check retrieval
```

---

## Step 1 — Structure Preparation

```python
from rdkit import Chem
from rdkit.Chem import AllChem
import MDAnalysis as mda

# Load with MDAnalysis (handles PDB/mmCIF)
u = mda.Universe('complex.pdb')

# Extract ligand by residue name
lig_atoms = u.select_atoms('resname LIG')
prot_atoms = u.select_atoms('protein')

# Write separate files
lig_atoms.write('ligand.pdb')
prot_atoms.write('receptor.pdb')

# Convert ligand to RDKit
ligand = Chem.MolFromPDBFile('ligand.pdb', removeHs=False)
ligand = Chem.RemoveHs(ligand)
# Assign bond orders from SMILES if needed (PDB lacks bond orders)
from rdkit.Chem.AllChem import AssignBondOrdersFromTemplate
template = Chem.MolFromSmiles('known_smiles')
ligand_fixed = AllChem.AssignBondOrdersFromTemplate(template, ligand)
```

---

## Step 2 — Interaction Detection with ProLIF

```python
import prolif as plf
import MDAnalysis as mda

u = mda.Universe('complex.pdb')

# Select ligand and protein
lig = u.select_atoms('resname LIG')
prot = u.select_atoms('protein')

# Fingerprint (single frame)
fp = plf.Fingerprint()
fp.run_from_atomgroup(lig, prot)
df = fp.to_dataframe()

# Interactions present
interactions = df.columns.tolist()  # ('LIG1', 'SER100', 'HBDonor'), etc.
print(df.T[df.any()].index.tolist())  # active interactions only
```

### ProLIF Interaction Types → Pharmacophore Feature Mapping

```python
PROLIF_TO_PHARM = {
    'HBDonor':       'HBDonor',       # ligand donates H
    'HBAcceptor':    'HBAcceptor',    # ligand accepts H
    'Hydrophobic':   'Hydrophobe',
    'PiStacking':    'Aromatic',
    'PiCation':      'PosIonizable',  # ligand cation + π
    'CationPi':      'Aromatic',      # ligand π + protein cation
    'Anionic':       'NegIonizable',
    'Cationic':      'PosIonizable',
    'VdWContact':    None,            # too weak; skip
    'MetalAcceptor': 'NegIonizable',
    'FaceToFace':    'Aromatic',
    'EdgeToFace':    'Aromatic',
}
```

---

## Step 3 — Build Pharmacophore from Interactions

```python
import os
import numpy as np
from rdkit import RDConfig, Chem
from rdkit.Chem import AllChem, rdMolChemicalFeatures

fdefName = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
factory = rdMolChemicalFeatures.BuildFeatureFactory(fdefName)

def build_structure_based_pharmacophore(ligand_mol, prolif_df,
                                         tolerance=1.5):
    """
    Build pharmacophore query from ProLIF interaction dataframe.
    Returns list of (family, position, radius) tuples.
    """
    # Get 3D features on ligand
    AllChem.EmbedMolecule(ligand_mol, AllChem.ETKDGv3())
    feats = factory.GetFeaturesForMol(ligand_mol)

    # Map feature families to 3D positions
    feat_map = {}
    for feat in feats:
        fam = feat.GetFamily()
        if fam not in feat_map:
            feat_map[fam] = []
        pos = feat.GetPos()
        feat_map[fam].append({
            'pos': np.array([pos.x, pos.y, pos.z]),
            'atoms': list(feat.GetAtomIds()),
        })

    # Filter to features that have ProLIF-confirmed interactions
    pharmacophore_points = []
    for col in prolif_df.columns:
        _, residue, interaction = col
        pharm_family = PROLIF_TO_PHARM.get(interaction)
        if pharm_family is None or pharm_family not in feat_map:
            continue

        # Take the first matching feature (closest to interaction)
        for feat_entry in feat_map[pharm_family]:
            pharmacophore_points.append({
                'family': pharm_family,
                'position': feat_entry['pos'],
                'radius': tolerance,
                'residue': residue,
                'interaction': interaction,
            })
            break  # one point per interaction type (merge later)

    return pharmacophore_points

pharm_points = build_structure_based_pharmacophore(ligand_fixed, df)
for p in pharm_points:
    print(f"{p['family']:15s} via {p['residue']:8s} {p['interaction']}")
```

---

## Step 4 — PLIP Alternative (Python Library)

PLIP (Protein-Ligand Interaction Profiler) provides detailed interaction geometry.

```python
from plip.structure.preparation import PDBComplex

complex = PDBComplex()
complex.load_pdb('complex.pdb')
complex.analyze()

# Get interactions for first binding site
site_key = list(complex.interaction_sets.keys())[0]
interactions = complex.interaction_sets[site_key]

# H-bonds where ligand is donor
for hb in interactions.hbonds_ldon:
    print(f"HBD: ligand atom {hb.d_orig_idx} → {hb.restype}{hb.resnr}")
    print(f"  distance: {hb.distance_ad:.2f} Å, angle: {hb.angle:.1f}°")
    print(f"  donor pos: {hb.d.coords}")

# H-bonds where ligand is acceptor
for hb in interactions.hbonds_pdon:
    print(f"HBA: {hb.restype}{hb.resnr} → ligand atom {hb.a_orig_idx}")

# Hydrophobic contacts
for hc in interactions.hydrophobic_contacts:
    print(f"HYD: ligand atom {hc.ligatom_orig_idx}, dist {hc.distance:.2f} Å")

# π–π stacking
for pi in interactions.pistacking:
    print(f"π-π: {pi.type} (offset {pi.offset:.2f} Å, angle {pi.angle:.1f}°)")
```

---

## Step 5 — Add Exclusion Volumes

```python
import numpy as np
from rdkit import Chem

def add_exclusion_volumes(prot_pdb, ligand_pos, cutoff=5.0, bump_radius=1.5):
    """
    Find protein atoms within cutoff Å of any ligand atom.
    Each becomes an exclusion sphere.
    """
    protein = Chem.MolFromPDBFile(prot_pdb, removeHs=True)
    prot_conf = protein.GetConformer()
    prot_pos = prot_conf.GetPositions()

    ev_spheres = []
    for i, ppos in enumerate(prot_pos):
        # Minimum distance to any ligand atom
        dists = np.linalg.norm(ligand_pos - ppos, axis=1)
        min_dist = dists.min()

        if min_dist < cutoff:
            atom = protein.GetAtomWithIdx(i)
            vdW = {'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80}
            r = vdW.get(atom.GetSymbol(), 1.70) + bump_radius
            ev_spheres.append({
                'pos': ppos, 'radius': r,
                'atom_idx': i,
                'symbol': atom.GetSymbol()
            })

    return ev_spheres
```

---

## Pharmacophore Export (LigandScout LSFP / Pharmer JSON)

```python
import json

def export_pharmer_json(pharm_points, ev_spheres, outfile='pharma.json'):
    """Export pharmacophore in Pharmer/Pharmit JSON format."""
    features = []

    type_map = {
        'HBDonor': 'HydrogenDonor',
        'HBAcceptor': 'HydrogenAcceptor',
        'Aromatic': 'Aromatic',
        'Hydrophobe': 'Hydrophobic',
        'PosIonizable': 'PositiveIon',
        'NegIonizable': 'NegativeIon',
    }

    for p in pharm_points:
        ftype = type_map.get(p['family'], p['family'])
        features.append({
            'name': ftype,
            'x': float(p['position'][0]),
            'y': float(p['position'][1]),
            'z': float(p['position'][2]),
            'radius': p['radius'],
            'enabled': True,
            'vector_on': 0,
        })

    # Exclusion volumes
    for ev in ev_spheres:
        features.append({
            'name': 'ExclusionSphere',
            'x': float(ev['pos'][0]),
            'y': float(ev['pos'][1]),
            'z': float(ev['pos'][2]),
            'radius': ev['radius'],
            'enabled': True,
        })

    with open(outfile, 'w') as f:
        json.dump({'points': features}, f, indent=2)

    print(f"Exported {len(pharm_points)} features + {len(ev_spheres)} exclusion spheres")
    return outfile
```
