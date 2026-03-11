# Fragment Growing, Linking, and Merging

## Growing via R-Group Enumeration

The most common computational strategy: enumerate substituted analogs
at identified exit vectors from the fragment pose.

```python
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
import pandas as pd

def find_exit_vectors(fragment_smiles):
    """
    Identify attachment points (Hs on non-NH positions)
    suitable for growing into adjacent pockets.
    Returns list of atom indices.
    """
    mol = Chem.MolFromSmiles(fragment_smiles)
    mol_h = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol_h, AllChem.ETKDGv3())

    vectors = []
    for atom in mol.GetAtoms():
        # Skip atoms already substituted or in rings
        if atom.GetDegree() >= 3:
            continue
        # Carbon, N, O with available H
        if atom.GetSymbol() in ('C', 'N', 'O') and atom.GetTotalNumHs() > 0:
            vectors.append(atom.GetIdx())
    return vectors

def enumerate_grown_library(fragment_smiles, r_groups, attachment_idx,
                             max_mw=400, max_hac=25):
    """
    Enumerate fragment + R-group at specified attachment point.
    Uses RDKit reaction: [*:1][H].[*:1]R >> [*:1]R
    """
    from rdkit.Chem import Descriptors, Crippen

    # Build attachment reaction
    rxn = AllChem.ReactionFromSmarts(
        '[cH1:1].[*:2]>>[c:1]-[*:2]'   # aryl C-H functionalization example
    )
    # More general: handle explicit attachment
    products = []
    frag_mol = Chem.MolFromSmiles(fragment_smiles)

    for rg_smiles in r_groups:
        rg_mol = Chem.MolFromSmiles(rg_smiles)
        if rg_mol is None:
            continue
        try:
            rxn_products = rxn.RunReactants((frag_mol, rg_mol))
            for prod_tuple in rxn_products:
                prod = prod_tuple[0]
                try:
                    Chem.SanitizeMol(prod)
                    smi = Chem.MolToSmiles(prod)
                    mw = Descriptors.ExactMolWt(prod)
                    hac = prod.GetNumHeavyAtoms()
                    if mw <= max_mw and hac <= max_hac:
                        products.append({'smiles': smi, 'MW': mw, 'HAC': hac,
                                         'parent': fragment_smiles, 'r_group': rg_smiles})
                except Exception:
                    continue
        except Exception:
            continue

    return pd.DataFrame(products).drop_duplicates('smiles')
```

## Growing with Explicit SMARTS Reaction

```python
from rdkit import Chem
from rdkit.Chem import AllChem

# Common growing reactions (use chem_react.py patterns)
GROWING_REACTIONS = {
    'amide_coupling': '[C:1](=O)[OH:2].[N:3][H]>>[C:1](=O)[N:3]',
    'reductive_amination': '[C:1]=O.[N:2][H]>>[C:1][N:2]',
    'buchwald': '[c:1][Br].[N:2][H]>>[c:1][N:2]',
    'suzuki': '[c:1][Br].[c:2][B]([OH])[OH]>>[c:1][c:2]',
    'click_triazole': '[C:1]#N=[N+]=[N-].[C:2]#[C:3]>>[c:1]1[n:2][n][n][c:3]1',
    'n_alkylation': '[N:1][H].[Br:2][C:3]>>[N:1][C:3]',
    'o_alkylation': '[O:1][H].[Br:2][C:3]>>[O:1][C:3]',
}

def apply_growing_reaction(fragment_smiles, reagent_smiles, reaction_name):
    """Apply a named growing reaction to a fragment + reagent."""
    if reaction_name not in GROWING_REACTIONS:
        raise ValueError(f"Unknown reaction: {reaction_name}")

    rxn = AllChem.ReactionFromSmarts(GROWING_REACTIONS[reaction_name])
    frag = Chem.MolFromSmiles(fragment_smiles)
    reagent = Chem.MolFromSmiles(reagent_smiles)

    products = []
    for prod_tuple in rxn.RunReactants((frag, reagent)):
        prod = prod_tuple[0]
        try:
            Chem.SanitizeMol(prod)
            products.append(Chem.MolToSmiles(prod))
        except Exception:
            continue
    return list(set(products))
```

## Fragment Merging via MCS

Merge two overlapping fragment hits into a single molecule:

```python
from rdkit import Chem
from rdkit.Chem import rdFMCS, AllChem

def merge_fragments(smiles_a, smiles_b, min_atoms=4):
    """
    Merge two fragment hits by finding their MCS and combining
    the unique portions of each fragment.

    Returns candidate merged SMILES or None if no suitable MCS found.
    """
    mol_a = Chem.MolFromSmiles(smiles_a)
    mol_b = Chem.MolFromSmiles(smiles_b)

    # Find Maximum Common Substructure
    mcs = rdFMCS.FindMCS(
        [mol_a, mol_b],
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=rdFMCS.BondCompare.CompareOrder,
        ringMatchesRingOnly=True,
        completeRingsOnly=True,
        timeout=10
    )

    if mcs.numAtoms < min_atoms:
        return None   # insufficient overlap

    mcs_mol = Chem.MolFromSmarts(mcs.smartsString)
    print(f"MCS: {mcs.smartsString} ({mcs.numAtoms} atoms, {mcs.numBonds} bonds)")

    # Match MCS in both fragments
    match_a = mol_a.GetSubstructMatch(mcs_mol)
    match_b = mol_b.GetSubstructMatch(mcs_mol)

    if not match_a or not match_b:
        return None

    # Find unique atoms in B (not in MCS)
    unique_b = [i for i in range(mol_b.GetNumAtoms()) if i not in match_b]
    if not unique_b:
        return smiles_a   # B is a substructure of A

    # Atoms in B connected to MCS (attachment points)
    attach_in_b = []
    attach_in_mcs = []
    for bond in mol_b.GetBonds():
        bi, bj = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if bi in match_b and bj in unique_b:
            attach_in_b.append(bj)
            attach_in_mcs.append(match_a[match_b.index(bi)])
        elif bj in match_b and bi in unique_b:
            attach_in_b.append(bi)
            attach_in_mcs.append(match_a[match_b.index(bj)])

    return {
        'mcs_smarts': mcs.smartsString,
        'mcs_size': mcs.numAtoms,
        'unique_b_atoms': len(unique_b),
        'attachment_points': list(zip(attach_in_mcs, attach_in_b))
    }
```

## Fragment Linking with Linker Design

```python
def design_linker(smiles_a, smiles_b, max_linker_atoms=4):
    """
    Design minimal linkers between two fragment hits.
    Simple approach: enumerate common linker motifs.
    """
    LINKER_SMARTS = {
        'direct_bond': '',
        'methylene': 'C',
        'ethylene': 'CC',
        'ethylamine': 'CN',
        'piperazine': 'C1CNCCN1',
        'propylene': 'CCC',
        'vinyl': 'C=C',
        'amide': 'C(=O)N',
        'ether': 'COC',
        'piperidine': 'C1CCNCC1',
    }

    # Fragment A exit: replace terminal H with attachment
    # Fragment B entry: add second attachment
    # → combine: A_smiles[*:1]-linker-[*:2]B_smiles

    results = []
    for linker_name, linker_smi in LINKER_SMARTS.items():
        if linker_smi == '':
            # Direct bond
            candidate = f'{smiles_a}.{smiles_b}'  # need actual SMILES manipulation
        else:
            candidate = f'[{smiles_a}]-{linker_smi}-[{smiles_b}]'
        results.append({'linker': linker_name, 'smiles_hint': candidate})
    return results
```

## Generative Growing with REINVENT

REINVENT with fragment constraint is the most powerful computational growing tool:

```toml
# REINVENT 4 config for scaffold growing
[parameters]
prior = "reinvent4_prior.prior"
agent = "agent.prior"

[[stage]]
type = "reinforcement-learning"
max_steps = 500

# Scaffold constraint: keep fragment core fixed
[stage.inception]
smiles = ["c1ccncc1C(=O)N"]   # fragment SMILES — used as seed

[stage.scoring]
type = "custom_product"

# Must contain fragment substructure
[[stage.scoring.component]]
[stage.scoring.component.matching_substructure]
smarts = ["c1ccncc1C(=O)N"]   # fragment SMARTS
weight = 1.0

# Drug-like expansion
[[stage.scoring.component]]
[stage.scoring.component.qed]
weight = 0.5

# SA feasibility
[[stage.scoring.component]]
[stage.scoring.component.sa_score]
weight = 0.3

# MW target: fragment + ~150 Da elaboration
[[stage.scoring.component]]
[stage.scoring.component.mol_weight]
weight = 0.2
[stage.scoring.component.mol_weight.params]
low = 200
high = 400

# Diversity filter
[stage.diversity_filter]
type = "IdenticalMurckoScaffold"
min_score = 0.4
bucket_size = 25
```

## Synthesis Feasibility Filter for Grown Compounds

```python
from rdkit.Chem import RDConfig
import os

def sa_score_filter(smiles_list, max_sa=4.0):
    """
    Filter grown compounds by Synthetic Accessibility score.
    SA ≤ 3: easily synthesizable
    SA ≤ 4: acceptable for fragment growing
    SA > 4: challenging — avoid unless strong justification
    """
    import sys
    sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
    import sascorer

    results = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        sa = sascorer.calculateScore(mol)
        results.append({'smiles': smi, 'SA_score': sa,
                        'synthesizable': sa <= max_sa})
    return pd.DataFrame(results)
```

## Elaboration Scoring: Combined Metric

After growing, score candidates by combined docking + efficiency:

```python
def score_elaborated_compounds(grown_df, parent_fragment,
                                parent_pIC50_pred, receptor_pdbqt,
                                center, box_size):
    """
    Score grown compounds by:
    1. Docking score improvement vs parent
    2. LE maintained
    3. LLE maintained or improved
    4. SA feasibility
    """
    from rdkit.Chem import Descriptors, Crippen
    import sascorer

    # Parent properties
    parent_mol = Chem.MolFromSmiles(parent_fragment)
    parent_hac = parent_mol.GetNumHeavyAtoms()
    parent_le = 1.37 * parent_pIC50_pred / parent_hac

    records = []
    for _, row in grown_df.iterrows():
        mol = Chem.MolFromSmiles(row.smiles)
        if mol is None:
            continue
        hac = mol.GetNumHeavyAtoms()
        mw = Descriptors.ExactMolWt(mol)
        logp = Crippen.MolLogP(mol)
        sa = sascorer.calculateScore(mol)

        # Dock
        _, energies = dock_fragment(row.smiles, receptor_pdbqt,
                                     center, box_size, exhaustiveness=16)
        vina_score = energies[0][0] if len(energies) > 0 else 0.0

        # Predicted pIC50 via scaling (rough estimate from Vina)
        pred_pic50 = -vina_score / 1.37  # inverse of LE formula (rough)
        le = 1.37 * pred_pic50 / hac
        lle = pred_pic50 - logp

        records.append({
            'smiles': row.smiles,
            'HAC': hac, 'MW': mw, 'cLogP': logp,
            'vina_score': vina_score,
            'pred_pIC50': pred_pic50,
            'LE': le, 'LLE': lle, 'SA': sa,
            'LE_vs_parent': le - parent_le,
            'HAC_added': hac - parent_hac,
            'keep': (le >= 0.28 and lle >= 2.0 and sa <= 4.0
                     and vina_score < -4.0)
        })

    return pd.DataFrame(records).sort_values('vina_score')
```

## Decision Tree: Which Strategy to Use

```
Fragment hit confirmed by biophysics?
├── Yes, have co-crystal structure
│   ├── Single fragment → GROWING (R-group enumeration at exit vectors)
│   └── Two fragments in adjacent pockets → LINKING or MERGING
└── Yes, only docking pose (no crystal)
    ├── High confidence pose → GROWING (verify with docking)
    └── Low confidence pose → MERGING (MCS-based, less pose-dependent)

Elaboration goal?
├── Improve potency (stay on target) → GROWING into hinge/solvent pocket
├── Improve selectivity → GROWING into selectivity pocket
├── Reduce cLogP → GROWING with polar substituents, avoid adding rings
└── Improve PK → REINVENT RL with ADMET scoring components
```
