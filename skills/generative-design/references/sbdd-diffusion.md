# Structure-Based 3D Generation — DiffSBDD, TargetDiff, Pocket2Mol, Linker Design

## Why Structure-Based Generation?

Condition molecule generation on a protein pocket → directly encode binding site geometry into generated structures. Higher relevance for tight binders vs. ligand-based methods.

**Challenge**: 3D generation is ~10-100× harder than 1D/2D — requires SE(3)-equivariant architectures or diffusion in Cartesian space.

## Protein Pocket Preparation

Quality of generated molecules depends critically on pocket representation. Use `homology-modeling` → `structure-prep` reference for full protocol.

```python
from pdbfixer import PDBFixer
from openmm.app import PDBFile

# 1. Fix receptor (remove water, add missing atoms, H at pH 7.4)
fixer = PDBFixer("receptor.pdb")
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(7.4)
PDBFile.writeFile(fixer.topology, fixer.positions, open("receptor_prepped.pdb", "w"))

# 2. Extract pocket residues within 8 Å of ligand (if co-crystal available)
import MDAnalysis as mda
import numpy as np

u = mda.Universe("complex.pdb")
ligand = u.select_atoms("resname LIG")
pocket = u.select_atoms("protein and around 8.0 resname LIG")
pocket.write("pocket.pdb")

# 3. Define pocket center and box (for docking-style tools)
center = ligand.center_of_geometry()
box_size = ligand.positions.max(0) - ligand.positions.min(0) + 10.0  # 5 Å padding each side
```

## DiffSBDD

Equivariant diffusion model conditioning molecule generation on protein pocket.

```bash
# GitHub: arneschneuing/DiffSBDD
git clone https://github.com/arneschneuing/DiffSBDD
cd DiffSBDD
pip install -e .
```

### Inference (generate molecules for a pocket)

```bash
# Using pre-trained CrossDocked2020 model
python generate_ligands.py \
    --checkpoint checkpoints/crossdocked_fullatom_cond.ckpt \
    --pdbfile receptor_prepped.pdb \
    --ref_ligand ligand.sdf \           # defines pocket center
    --n_samples 100 \
    --all_frags \                       # include incomplete fragments
    --sanitize \                        # filter invalid molecules
    --relax                             # MMFF relaxation of output
```

### Python API

```python
import torch
from src.lightning_modules import LigandPocketDDPM

# Load model (inference mode — use torch.no_grad() context below)
model = LigandPocketDDPM.load_from_checkpoint(
    "checkpoints/crossdocked_fullatom_cond.ckpt",
    map_location="cpu"
)
model.training = False  # disable dropout/batch norm training behavior

# Prepare pocket representation
from src.utils import prepare_context
pocket_context = prepare_context(
    pdb_path="receptor_prepped.pdb",
    ref_ligand_sdf="ligand.sdf",
    n_nodes=100,         # max pocket atoms
    include_radius=8.0,
)

# Generate
with torch.no_grad():
    molecules = model.generate(
        pocket_context,
        n_samples=50,
        ddim_steps=100,    # fewer = faster but lower quality
    )

# molecules: list of rdkit Mol objects
from rdkit import Chem
for mol in molecules:
    if mol is not None:
        print(Chem.MolToSmiles(mol))
```

## TargetDiff

SE(3)-equivariant diffusion model; improved sampling quality over DiffSBDD.

```bash
# GitHub: guanjq/targetdiff
git clone https://github.com/guanjq/targetdiff
pip install -r requirements.txt
```

```bash
# Generate
python scripts/sample_diffusion.py \
    configs/sampling.yml \
    --outdir results/ \
    --protein_path receptor.pdb \\
    --ligand_path ref_ligand.sdf \
    --num_samples 100 \
    --batch_size 8 \
    --steps 1000
```

```python
# Post-processing: extract valid SMILES with drug-likeness scoring
from rdkit import Chem
from rdkit.Chem import Descriptors, QED

def filter_generated_mols(mol_list, max_mw=600, max_logp=5, max_sa=5):
    """Filter DiffSBDD/TargetDiff output to drug-like subset."""
    from sascorer import calculateScore
    valid = []
    for mol in mol_list:
        if mol is None: continue
        try:
            smi = Chem.MolToSmiles(mol)
            mol_clean = Chem.MolFromSmiles(smi)
            if mol_clean is None: continue
            mw = Descriptors.ExactMolWt(mol_clean)
            lp = Descriptors.MolLogP(mol_clean)
            sa = calculateScore(mol_clean)
            qed = QED.qed(mol_clean)
            if mw <= max_mw and lp <= max_logp and sa <= max_sa:
                valid.append((smi, qed, sa))
        except Exception:
            continue
    return sorted(valid, key=lambda x: -x[1])  # sort by QED
```

## Pocket2Mol

Autoregressive 3D molecule generation conditioned on pocket (atom-by-atom placement).

```bash
# GitHub: pengxingang/Pocket2Mol
git clone https://github.com/pengxingang/Pocket2Mol
pip install -r requirements.txt
```

```bash
python sample.py \
    --config configs/sample.yml \
    --pdb_path receptor_prepped.pdb \
    --ligand_path ref_ligand.sdf \
    --result_path results/ \
    --num_samples 100
```

**When to prefer Pocket2Mol over DiffSBDD**:
- Better atom-level positioning (AR model places atoms one-by-one)
- Slower per molecule but higher interpretability
- Better for very constrained pockets

## Fragment Growing and Linker Design

### DiffLinker — 3D linker design

Connect two fragments in 3D via diffusion-based linker generation.

```bash
# GitHub: igashov/DiffLinker
git clone https://github.com/igashov/DiffLinker
pip install -r requirements.txt

# Generate linkers between two fragments in pocket
python generate.py \
    --fragments frag_a.sdf frag_b.sdf \
    --protein receptor.pdb \
    --n_samples 50 \
    --linker_size 5          # target linker heavy atom count (optional)
```

```python
# Fragment + linker assembly validation
from rdkit import Chem

def validate_linked_molecule(full_mol_smi, frag_a_smi, frag_b_smi):
    """Check that linked molecule contains both fragment substructures."""
    mol = Chem.MolFromSmiles(full_mol_smi)
    if mol is None: return False

    frag_a = Chem.MolFromSmiles(frag_a_smi)
    frag_b = Chem.MolFromSmiles(frag_b_smi)
    return mol.HasSubstructMatch(frag_a) and mol.HasSubstructMatch(frag_b)
```

### DeLinker — 2D graph-based linker design

When 3D structure not available:

```bash
# GitHub: oxpig/DeLinker
pip install dgl torch  # DeLinker uses DGL graph library
python generate.py --fragment1 "c1ccccc1" --fragment2 "C1CCNCC1" --n_samples 50
```

## Scoring Generated Molecules from SBDD

```python
from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS
import numpy as np

def score_sbdd_results(generated_mols, ref_ligand_mol):
    """
    Score SBDD output molecules for drug-likeness and novelty.
    Returns list of dicts with metrics per molecule.
    """
    from sascorer import calculateScore
    from rdkit.Chem import Descriptors, QED, DataStructs, rdMolDescriptors

    fp_ref = rdMolDescriptors.GetMorganFingerprintAsBitVect(ref_ligand_mol, 2, 2048)
    results = []

    for mol in generated_mols:
        if mol is None: continue
        smi = Chem.MolToSmiles(mol)
        mol2 = Chem.MolFromSmiles(smi)  # re-sanitize
        if mol2 is None: continue

        qed_val = QED.qed(mol2)
        sa_val = calculateScore(mol2)
        mw = Descriptors.ExactMolWt(mol2)
        logp = Descriptors.MolLogP(mol2)

        fp_gen = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol2, 2, 2048)
        tanimoto = DataStructs.TanimotoSimilarity(fp_gen, fp_ref)

        results.append({
            "smiles": smi, "qed": qed_val, "sa": sa_val,
            "mw": mw, "logp": logp, "tanimoto_to_ref": tanimoto
        })

    return sorted(results, key=lambda x: -x["qed"])

# RMSD from reference pose (if both are 3D)
def rmsd_3d_mcs(mol_gen, mol_ref):
    """RMSD between generated and reference ligand poses (after MCS alignment)."""
    mcs = rdFMCS.FindMCS([mol_gen, mol_ref], timeout=5)
    if mcs.numAtoms < 5: return None

    patt = Chem.MolFromSmarts(mcs.smartsString)
    idx_gen = mol_gen.GetSubstructMatch(patt)
    idx_ref = mol_ref.GetSubstructMatch(patt)

    pos_gen = np.array([mol_gen.GetConformer().GetAtomPosition(i) for i in idx_gen])
    pos_ref = np.array([mol_ref.GetConformer().GetAtomPosition(i) for i in idx_ref])

    return float(np.sqrt(np.mean(np.sum((pos_gen - pos_ref)**2, axis=1))))
```

## ProLIF Interaction Analysis for SBDD Output

```python
import prolif
import MDAnalysis as mda

def check_pocket_interactions(ligand_sdf, receptor_pdb, required_interactions=None):
    """
    Verify generated molecule recapitulates key binding interactions.
    required_interactions: list of (resid, interaction_type) tuples
      e.g. [("LYS203", "HBAcceptor"), ("PHE382", "Hydrophobic")]
    """
    u = mda.Universe(receptor_pdb, ligand_sdf)
    protein = u.select_atoms("protein")
    ligand = u.select_atoms("resname LIG")

    fp = prolif.Fingerprint()
    fp.run_from_atomgroup(ligand, protein)

    df = fp.to_dataframe()
    interactions_found = {}
    for col in df.columns:
        resid, int_type = col[0], col[1]
        if df[col].any():
            interactions_found[f"{resid}_{int_type}"] = True

    if required_interactions:
        satisfied = [f"{r}_{t}" in interactions_found
                     for r, t in required_interactions]
        return all(satisfied), interactions_found
    return interactions_found
```

## Clash Detection and 3D Quality

```python
from rdkit.Chem import AllChem
import numpy as np

def check_3d_clash(mol, clash_threshold=1.5):
    """Check 3D molecule for intra-molecular steric clashes."""
    if mol.GetNumConformers() == 0:
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        AllChem.MMFFOptimizeMolecule(mol)

    conf = mol.GetConformer()
    pos = conf.GetPositions()
    bonds_set = {(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()}

    for i in range(len(pos)):
        for j in range(i+2, len(pos)):
            if (i, j) in bonds_set or (j, i) in bonds_set:
                continue
            dist = float(np.linalg.norm(pos[i] - pos[j]))
            if dist < clash_threshold:
                return False, f"Clash atoms {i},{j} dist={dist:.2f}"
    return True, "OK"
```

## Model Comparison

| Model | Architecture | Speed | Validity | Pocket conditioning | Best for |
|-------|-------------|-------|----------|--------------------|----|
| **DiffSBDD** | E3-equivariant diffusion | Medium | ~70-80% | Strong | General SBDD |
| **TargetDiff** | SE(3) diffusion | Slow | ~75-85% | Strong | Tight binding sites |
| **Pocket2Mol** | Autoregressive | Slow | ~85-90% | Strong | Constrained pockets |
| **DiffLinker** | Diffusion (linker only) | Fast | ~80% | Optional | Fragment linking |
| **DeLinker** | GNN (2D graph) | Very fast | ~90% | None | No structure available |

## Key Pitfalls

- **Pocket quality is critical**: SBDD models trained on crystal structures; MD-derived pockets need careful selection
- **3D validity ≠ chemical validity**: re-sanitize via `Chem.MolFromSmiles(Chem.MolToSmiles(mol))`
- **SA score still needed**: 3D models can generate geometrically valid but synthetically intractable molecules
- **Generated 3D coordinates are approximate**: use as starting poses for docking refinement, not final structures
- **DiffSBDD needs CUDA**: inference is very slow on CPU; minimum GPU: RTX 3060 (~10 min for 100 mols)
- **Pocket radius matters**: 6 Å too tight; 10 Å includes non-interacting residues; 8 Å is standard
- **Clash check is not automatic**: always run `check_3d_clash` and optionally MMFF-relax generated poses
