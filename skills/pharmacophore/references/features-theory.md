# Pharmacophore Features — Theory and FDEF Format

## Core Concepts

A pharmacophore model consists of:
1. **Feature points** — type + 3D centroid coordinates
2. **Tolerance spheres** — radius around each centroid (typically 1-2 Å)
3. **Distance/angle constraints** — distances between feature pairs (optional)
4. **Exclusion volumes** — spheres where no atoms are allowed (from protein pocket)

### Pharmacophore vs. Shape
- **Pharmacophore**: abstract chemical features, allows scaffold hopping
- **Shape**: 3D atom positions, encodes exact topology → complementary methods

---

## Feature Families

### H-Bond Donor (HBD)
- Atom with transferable H attached to electronegative atom (N-H, O-H)
- SMARTS: `[$([N;!H0;v3,v4&+1]),$([n;H1+0]),$([O,S;H1;+0]),$([nH1])]`
- Direction: along N-H or O-H bond → directional feature
- Radius: 1.5-2.5 Å (directional tolerance 30-60°)

### H-Bond Acceptor (HBA)
- Atom with lone pair available for H-bond
- SMARTS: `[$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N;v3;!$(N-*=!@[O,N,P,S])]),$([o,s;+0;!$([o,s]:n);!$([o,s]:c:n)])]`
- Common acceptors: carbonyl O, ether O, pyridine N, F
- Direction: along lone-pair axis

### Aromatic (AR)
- Aromatic ring centroid + normal vector
- SMARTS: `[a]` (any aromatic atom)
- Centroid = geometric center of ring atoms
- Angle: ring plane normal defines direction (edge-face vs. face-face stacking)

### Hydrophobic (HYD)
- Non-polar region; grease-grease interactions
- SMARTS: `[c,s,w,S&H0&v2,$([D3&!Ring1](~[CH2]~[CH2])~[CH2]),$([v3,v4;#6,Si]~[!#1;!#7;!#8;!#9;!#15;!#16;!#17;!#35;!#53])]`
- Simplified: aliphatic C, S, halogens (Cl, Br, I)
- Broad spheres (2.0-3.0 Å) due to non-directionality

### Positive Ionizable (POS)
- Groups protonated at physiological pH
- SMARTS: `[$([NH2+0]C(=NH)N),$([NH2+0]c1nccc1),$([n+0;H0]1cccc1),$([NH2+0]C(=[NH+0])N),$([nH+0]1ccc[nH+0]1),…]`
- Common: amines, guanidines, amidines, imidazoles (His)

### Negative Ionizable (NEG)
- Groups deprotonated at physiological pH
- SMARTS: `[$([C,S](=[O,S,P])-[OH]),$([P,S]([OH])(=[O]))]`
- Common: carboxylic acids, sulfonamides, phosphates, sulfonates

### Exclusion Volume (XV)
- Derived from protein pocket atoms within van der Waals contact
- Not a chemical feature of the ligand — blocks poses that clash with protein
- Sphere radius = vdW radius of protein atom (C: 1.7 Å, N/O: 1.5 Å)

---

## RDKit FDEF Format

FDEF files define feature families using AtomType + SMARTS patterns. RDKit ships `BaseFeatures.fdef`.

```
# Custom FDEF example (save as custom.fdef)
AtomType NDonor [N;!H0;v3,v4&+1]
AtomType NaroN [n;H1+0]
AtomType OHDonor [O,S;H1;+0]
AtomType NHaroN [nH1]

# Define feature family
DefineFeature SingleAtomDonor [$NDonor,$NaroN,$OHDonor,$NHaroN]
  Family HBDonor
  Weights 1.0
EndFeature

AtomType NAcceptor [n;+0;!$([n]:a:[!n]);!$([n]:1:[!n]:[!n]:[!n]1)]
AtomType OAcceptor [$([OH0;v2]),$([O;-]),$([o;+0])]

DefineFeature SingleAtomAcceptor [$NAcceptor,$OAcceptor]
  Family HBAcceptor
  Weights 1.0
EndFeature

DefineFeature AromRing a1aaaaa1,a1aaaa1
  Family Aromatic
  Weights 1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0
EndFeature

DefineFeature Hydrophobe c,s,[CH2v4;!R]
  Family Hydrophobe
  Weights 1.0
EndFeature
```

```python
# Load custom FDEF
from rdkit.Chem import rdMolChemicalFeatures
factory = rdMolChemicalFeatures.BuildFeatureFactory('custom.fdef')
```

---

## Feature Weights and Scoring

When scoring how well a molecule fits a pharmacophore:

```python
# Score = Σ_i w_i · fit_i
# fit_i = 1 if feature i within tolerance sphere
#       = Gaussian decay if outside
# fit_i(d) = exp(-d² / (2σ²))

# σ = tolerance sphere radius (Å)
# Typical: HBD/HBA σ = 1.5-2.0 Å (directional)
#          HYD σ = 2.5-3.0 Å (non-directional)
#          AR σ = 1.5 Å (centroid) + 30° plane angle

# Score threshold: typically 0.7-0.9 (out of 1.0)
```

---

## Pharmacophore vs. Fingerprint

| | Pharmacophore model | Pharmacophore fingerprint (Pharm2D) |
|-|---------------------|--------------------------------------|
| Input | 3D model + conformer | 2D SMILES (topological) |
| Search | 3D matching | Bit vector similarity (Tanimoto) |
| Speed | Slow (conformer gen needed) | Fast (no 3D) |
| Scaffold hopping | Excellent | Moderate |
| False positive rate | Low | Higher |
| Use case | VS with known 3D features | LBVS, clustering, diversity |

---

## Exclusion Volume Construction

```python
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

def get_exclusion_volumes(protein_pdb, ligand_mol, radius_cutoff=6.0, ev_radius=1.5):
    """
    Generate exclusion volume spheres from protein atoms
    within radius_cutoff Å of the ligand centroid.
    """
    protein = Chem.MolFromPDBFile(protein_pdb, removeHs=True)
    if protein is None:
        raise ValueError("Could not parse protein PDB")

    # Ligand centroid
    conf = ligand_mol.GetConformer()
    lig_pos = conf.GetPositions()
    centroid = lig_pos.mean(axis=0)

    # Protein atoms within cutoff
    prot_conf = protein.GetConformer()
    prot_pos = prot_conf.GetPositions()

    ev_spheres = []
    for i, pos in enumerate(prot_pos):
        dist = np.linalg.norm(pos - centroid)
        if dist < radius_cutoff:
            atom = protein.GetAtomWithIdx(i)
            # vdW radii (approximate)
            vdW = {'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80,
                   'H': 1.20, 'P': 1.80, 'F': 1.47}
            r = vdW.get(atom.GetSymbol(), 1.70) + ev_radius  # Å
            ev_spheres.append({'pos': pos, 'radius': r, 'atom': atom.GetSymbol()})

    return ev_spheres

ev = get_exclusion_volumes('protein.pdb', ligand_mol)
print(f"{len(ev)} exclusion volume spheres")
```
