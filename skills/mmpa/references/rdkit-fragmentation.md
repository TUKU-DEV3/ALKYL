# RDKit Programmatic MMP Generation

## Why RDKit-Based MMP?

When mmpdb is overkill (small dataset, custom fragmentation rules, or integration into an existing RDKit pipeline), MMPs can be generated directly using RDKit's bond manipulation and SMARTS matching.

## Core: BRICS Fragmentation

BRICS (Breaking Retrosynthetically Interesting Chemical Substructures) identifies medicinal-chemistry-relevant bond cuts:

```python
from rdkit import Chem
from rdkit.Chem import BRICS, AllChem

mol = Chem.MolFromSmiles("c1ccc(CC(=O)NC2CCNCC2)cc1")

# Get BRICS fragments (all possible cuts)
fragments = BRICS.BRICSDecompose(mol)
# Returns a set of SMILES with [1*], [2*], ... attachment points

# Get BRICS bonds (which bonds get cut)
bonds = list(BRICS.FindBRICSBonds(mol))
# Returns: [(bond_idx, (atom1_idx, atom2_idx)), ...]

# Enumerate BRICS recombinations (for analogue generation)
all_frags = list(BRICS.BRICSDecompose(mol))
new_mols = list(BRICS.BRICSBuild(list(BRICS.BRICSDecompose(mol))))
valid_mols = [m for m in new_mols if m is not None]
```

**BRICS attachment point labels**: `[1*]`=CH3, `[4*]`=ring, `[7*]`=aromatic ring, etc. — 16 types total based on bond environment.

## Single-Bond Cut MMP Generation

The core algorithm for MMPA: cut every single non-ring bond, enumerate all (context, variable) pairs across a dataset:

```python
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from collections import defaultdict
from itertools import combinations

def get_single_cut_fragments(mol, max_variable_heavies=10):
    """
    Fragment a molecule at each single non-ring bond.
    Returns list of (variable_smiles, context_smiles) tuples.
    Each uses [*:1] as attachment point.
    """
    results = []
    bonds = mol.GetBonds()

    for bond in bonds:
        if bond.IsInRing():
            continue
        if bond.GetBondTypeAsDouble() != 1.0:
            continue  # only single bonds

        begin = bond.GetBeginAtomIdx()
        end   = bond.GetEndAtomIdx()

        # Break this bond and get two fragments
        em = Chem.RWMol(mol)
        em.RemoveBond(begin, end)

        # Add attachment point dummy atoms ([*:1]) at cut sites
        idx_a = em.AddAtom(Chem.Atom(0))  # atom 0 = dummy *
        em.GetAtomWithIdx(idx_a).SetAtomMapNum(1)
        em.AddBond(begin, idx_a, Chem.BondType.SINGLE)

        idx_b = em.AddAtom(Chem.Atom(0))
        em.GetAtomWithIdx(idx_b).SetAtomMapNum(1)
        em.AddBond(end, idx_b, Chem.BondType.SINGLE)

        try:
            frags = Chem.GetMolFrags(em.GetMol(), asMols=True)
        except Exception:
            continue

        if len(frags) != 2:
            continue

        for i, frag in enumerate(frags):
            other = frags[1 - i]
            # frag = variable part; other = context
            n_heavy = sum(1 for a in frag.GetAtoms() if a.GetAtomicNum() > 0)
            if n_heavy > max_variable_heavies:
                continue
            var_smi = Chem.MolToSmiles(frag, rootedAtAtom=-1)
            ctx_smi = Chem.MolToSmiles(other, rootedAtAtom=-1)
            results.append((var_smi, ctx_smi))

    return results


def find_mmps_rdkit(smiles_dict: dict, max_variable_heavies=10) -> list:
    """
    Find all MMPs across a dictionary {id: smiles}.
    Returns list of (id_a, id_b, var_a, var_b, context) tuples.
    """
    # Step 1: fragment all molecules
    context_index = defaultdict(list)  # context_smiles → [(id, variable_smiles)]

    for mol_id, smiles in smiles_dict.items():
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        frags = get_single_cut_fragments(mol, max_variable_heavies=max_variable_heavies)
        for var_smi, ctx_smi in frags:
            context_index[ctx_smi].append((mol_id, var_smi))

    # Step 2: find pairs — same context, different variable
    pairs = []
    for ctx_smi, entries in context_index.items():
        if len(entries) < 2:
            continue
        seen_vars = {}
        for mol_id, var_smi in entries:
            if var_smi not in seen_vars:
                seen_vars[var_smi] = mol_id
        # All pairs with same context but different variable part
        var_list = list(seen_vars.items())  # [(var_smi, mol_id)]
        for (var_a, id_a), (var_b, id_b) in combinations(var_list, 2):
            if var_a != var_b:
                pairs.append((id_a, id_b, var_a, var_b, ctx_smi))

    return pairs
```

## Attachment Point Canonicalization

SMILES with `[*:1]` can have multiple canonical forms — normalize to avoid duplicates:

```python
def canonicalize_fragment(smiles_with_attachment: str) -> str:
    """Canonicalize a fragment SMILES containing [*:1] attachment point."""
    mol = Chem.MolFromSmiles(smiles_with_attachment)
    if mol is None:
        return smiles_with_attachment
    return Chem.MolToSmiles(mol)  # RDKit canonical form

# Example:
# "[*:1]CC" and "CC[*:1]" → both canonicalize to same string
```

## SMARTS-Based Custom Bond Cutting

For targeted cuts (e.g., cut only amide bonds, or only C-N bonds):

```python
from rdkit.Chem import AllChem

def cut_bonds_by_smarts(mol, cut_smarts: str, max_variable_heavies=10):
    """
    Fragment molecule only at bonds matching a SMARTS pattern.
    cut_smarts example: "[CX3](=O)-[NX3]" for amide bonds
    """
    pattern = Chem.MolFromSmarts(cut_smarts)
    if pattern is None:
        raise ValueError(f"Invalid SMARTS: {cut_smarts}")

    matches = mol.GetSubstructMatches(pattern)
    # For amide SMARTS "[CX3](=O)-[NX3]", match[0]=C, match[2]=N (bond between them)

    results = []
    for match in matches:
        # Find the bond between first and last matched atoms (the "cut bond")
        a1, a2 = match[0], match[-1]
        bond = mol.GetBondBetweenAtoms(a1, a2)
        if bond is None:
            # Try first two atoms
            bond = mol.GetBondBetweenAtoms(match[0], match[1])
        if bond is None:
            continue

        em = Chem.RWMol(mol)
        em.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

        idx_dummy_a = em.AddAtom(Chem.Atom(0))
        em.GetAtomWithIdx(idx_dummy_a).SetAtomMapNum(1)
        em.AddBond(bond.GetBeginAtomIdx(), idx_dummy_a, Chem.BondType.SINGLE)

        idx_dummy_b = em.AddAtom(Chem.Atom(0))
        em.GetAtomWithIdx(idx_dummy_b).SetAtomMapNum(1)
        em.AddBond(bond.GetEndAtomIdx(), idx_dummy_b, Chem.BondType.SINGLE)

        try:
            frags = Chem.GetMolFrags(em.GetMol(), asMols=True)
        except Exception:
            continue

        if len(frags) == 2:
            for i in range(2):
                var = frags[i]
                ctx = frags[1 - i]
                n_heavy = sum(1 for a in var.GetAtoms() if a.GetAtomicNum() > 0)
                if n_heavy <= max_variable_heavies:
                    results.append((
                        Chem.MolToSmiles(var),
                        Chem.MolToSmiles(ctx),
                    ))
    return results
```

## Reconstructing the Full Molecule from Fragments

Applying a transform = replacing variable part R₁ with R₂ in the context:

```python
from rdkit.Chem import AllChem

def apply_transform(context_smi: str, new_var_smi: str) -> str | None:
    """
    Replace [*:1] in context with [*:1] from new variable part.
    Returns SMILES of the new full molecule.
    """
    # Use RDKit reaction to combine context + new variable
    rxn_smarts = f"[*:1].[*:1]>>[*:1]"  # generic attachment point join

    context_mol = Chem.MolFromSmiles(context_smi)
    new_var_mol  = Chem.MolFromSmiles(new_var_smi)
    if context_mol is None or new_var_mol is None:
        return None

    rxn = AllChem.ReactionFromSmarts("[*:1][*:2].[*:1][*:3]>>[*:2][*:3]")
    products = rxn.RunReactants((context_mol, new_var_mol))

    valid_smiles = []
    for prod_tuple in products:
        for prod in prod_tuple:
            try:
                smi = Chem.MolToSmiles(Chem.RemoveHs(prod))
                if Chem.MolFromSmiles(smi) is not None:
                    valid_smiles.append(smi)
            except Exception:
                pass

    return valid_smiles[0] if valid_smiles else None


def reconstruct_molecule(context_smi: str, variable_smi: str) -> str | None:
    """
    Reconnect context [*:1] + variable [*:1] into full molecule.
    Attachment points are matched by atom map number 1.
    """
    ctx = Chem.MolFromSmiles(context_smi)
    var = Chem.MolFromSmiles(variable_smi)
    if ctx is None or var is None:
        return None

    # Combine and remove dummy atoms
    combo = Chem.CombineMols(ctx, var)
    combo_rw = Chem.RWMol(combo)

    # Find the two dummy atoms ([*:1])
    dummy_atoms = [
        a.GetIdx() for a in combo_rw.GetAtoms()
        if a.GetAtomicNum() == 0 and a.GetAtomMapNum() == 1
    ]
    if len(dummy_atoms) != 2:
        return None

    d1, d2 = dummy_atoms
    # Find their heavy atom neighbors
    nbr1 = [n.GetIdx() for n in combo_rw.GetAtomWithIdx(d1).GetNeighbors() if n.GetAtomicNum() > 0]
    nbr2 = [n.GetIdx() for n in combo_rw.GetAtomWithIdx(d2).GetNeighbors() if n.GetAtomicNum() > 0]
    if not nbr1 or not nbr2:
        return None

    # Add bond between neighbors, remove dummies
    combo_rw.AddBond(nbr1[0], nbr2[0], Chem.BondType.SINGLE)
    for d in sorted([d1, d2], reverse=True):
        combo_rw.RemoveAtom(d)

    try:
        Chem.SanitizeMol(combo_rw)
        return Chem.MolToSmiles(combo_rw.GetMol())
    except Exception:
        return None
```

## Building a MMP Property DataFrame

```python
import pandas as pd
from rdkit.Chem import Descriptors, QED

def build_mmp_dataframe(smiles_dict: dict, properties: dict) -> pd.DataFrame:
    """
    Find MMPs and compute ΔProperty for each pair.

    Args:
        smiles_dict: {id: smiles}
        properties: {id: {prop_name: value}}

    Returns:
        DataFrame with columns: id_a, id_b, transform, context, delta_logp, delta_qed, ...
    """
    pairs = find_mmps_rdkit(smiles_dict, max_variable_heavies=10)

    rows = []
    for id_a, id_b, var_a, var_b, ctx in pairs:
        if id_a not in properties or id_b not in properties:
            continue
        row = {
            "id_a": id_a, "id_b": id_b,
            "smiles_a": smiles_dict[id_a], "smiles_b": smiles_dict[id_b],
            "variable_a": var_a, "variable_b": var_b,
            "context": ctx,
            "transform": f"{var_a}>>{var_b}",
        }
        props_a = properties[id_a]
        props_b = properties[id_b]
        for prop in set(list(props_a.keys()) + list(props_b.keys())):
            if prop in props_a and prop in props_b:
                row[f"delta_{prop}"] = props_b[prop] - props_a[prop]
        rows.append(row)

    return pd.DataFrame(rows)
```

## Handling Stereochemistry

```python
from rdkit.Chem import EnumerateStereoisomers

def fragment_with_stereo(smiles: str) -> list:
    """Fragment and preserve stereochemistry annotation in fragments."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None: return []

    # Include Chem.SmilesWriteParams for stereo
    params = Chem.SmilesWriteParams()
    params.isomericSmiles = True

    frags = get_single_cut_fragments(mol)
    # When using mmpdb, use --stereo flag for chirality-aware MMPA
    return frags
```

## Key Pitfalls

- **Fragment canonicalization is critical**: same variable part from different molecules must map to identical SMILES string — always pass through `Chem.MolToSmiles(Chem.MolFromSmiles(smi))`
- **Ring bonds**: never cut ring bonds for standard MMPA; add `if bond.IsInRing(): continue`
- **Double bonds**: default is single bond only; for vinyl substitutions, extend to include double bonds
- **Implicit vs explicit H**: `[H]` attached to heteroatom can appear as variable part (H→R swap) — filter out `[H][*:1]` if unwanted
- **Context size explosion**: large molecules with many single bonds → O(n²) pair explosion; limit variable part to ≤13 heavy atoms and context must have ≥5 heavy atoms
- **mmpdb vs. custom**: for datasets >10k molecules, use mmpdb (optimized C++); RDKit approach is fine for <5k molecules
