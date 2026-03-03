#!/usr/bin/env python3
"""ALKYL — chem_pka.py: estimate protonation states at a given pH."""

import argparse
import json
import sys

# (type, SMARTS, pKa_min, pKa_max, charge_when_ionized, is_acid)
PKA_RULES = [
    ("carboxylic_acid",  "[CX3](=O)[OX2H1]",                    3.5,  5.0,  -1, True),
    ("phenol",           "[OX2H1]c",                              8.0, 11.0,  -1, True),
    ("thiol",            "[SX2H1]",                               8.0, 11.0,  -1, True),
    ("sulfonamide_NH",   "[#16X4](=[OX1])(=[OX1])N[H]",          9.0, 11.0,  -1, True),
    ("phosphate",        "[P](=O)([OH])([OH])",                   1.0,  3.0,  -1, True),
    ("aniline",          "[NX3;H2]c",                             3.5,  5.5,  +1, False),
    ("imidazole",        "c1cnc[nH]1",                            6.0,  7.5,  +1, False),
    ("pyridine",         "n1ccccc1",                              4.5,  6.5,  +1, False),
    ("guanidine",        "NC(=N)N",                              12.0, 14.0,  +1, False),
    ("amine_primary",    "[NX3;H2;!$(NC=O);!$(Nc)]",              9.0, 11.0,  +1, False),
    ("amine_secondary",  "[NX3;H1;!$(NC=O);!$(Nc)]",              9.0, 11.0,  +1, False),
    ("amine_tertiary",   "[NX3;H0;!$(NC=O);!$([n]);!$(Nc)]",     8.0, 10.5,  +1, False),
]


def load_mol(args):
    from rdkit import Chem
    if args.smiles:
        mol = Chem.MolFromSmiles(args.smiles)
        if mol is None:
            print(f"Invalid SMILES: {args.smiles}", file=sys.stderr)
            sys.exit(1)
        return mol
    suppl = Chem.SDMolSupplier(args.sdf)
    mol = next((m for m in suppl if m is not None), None)
    if mol is None:
        print(f"No valid molecule in: {args.sdf}", file=sys.stderr)
        sys.exit(1)
    return mol


def estimate_pka(mol, ph: float):
    """Apply SMARTS-based pKa rules and return groups + net charge."""
    from rdkit import Chem
    groups = []
    net_charge = 0

    for (gtype, smarts, pka_min, pka_max, charge_delta, is_acid) in PKA_RULES:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue
        matches = mol.GetSubstructMatches(pattern)
        if not matches:
            continue

        if is_acid:
            # acid loses proton (ionizes) when pH >> pKa
            state = ("ionized" if ph > pka_max
                     else "mixed" if ph > pka_min
                     else "neutral")
        else:
            # base: protonated (ionized) when pH << pKaH
            state = ("ionized" if ph < pka_min
                     else "mixed" if ph < pka_max
                     else "neutral")

        charge = charge_delta if state == "ionized" else 0
        net_charge += charge * len(matches)
        groups.append({
            "type": gtype,
            "pka_range": [pka_min, pka_max],
            "count": len(matches),
            "state_at_ph": state,
            "charge_contribution": charge * len(matches),
        })

    return groups, net_charge


def build_dominant_form(mol, groups: list) -> str:
    """Return SMILES of dominant form by editing formal charges based on ionization state.

    For each ionized group, applies the appropriate formal charge to the ionizable atom.
    Acid groups lose H and gain negative charge; base groups (when ionized at low pH)
    gain positive charge. Uses SetNoImplicit to avoid valence errors.
    """
    from rdkit import Chem
    from rdkit.Chem import RWMol
    from rdkit.Chem.MolStandardize import rdMolStandardize

    # Map group type → (SMARTS, index_in_match_of_ionizable_atom, charge_delta)
    # Index refers to position within GetSubstructMatches tuple.
    # For acids: ionizable atom is the H-bearing heteroatom.
    # For bases: ionizable atom is the N.
    TYPE_TO_SMARTS = {
        # Acids (ionized = deprotonated = negative)
        "carboxylic_acid": ("[CX3](=O)[OX2H1]",                    2, -1),
        "phenol":          ("[OX2H1]c",                             0, -1),
        "thiol":           ("[SX2H1]",                              0, -1),
        "sulfonamide_NH":  ("[#16X4](=[OX1])(=[OX1])[NH1]",        3, -1),
        "phosphate":       ("[PX4](=O)([OX2H1])([OX2H1])",         2, -1),
        # Bases (ionized = protonated = positive; only at very low pH)
        "aniline":         ("[NX3;H2]c",                            0, +1),
        "imidazole":       ("c1cnc[nH]1",                           4, +1),
        "pyridine":        ("n1ccccc1",                             0, +1),
        "guanidine":       ("NC(=N)N",                              0, +1),
        "amine_primary":   ("[NX3;H2;!$(NC=O);!$(Nc)]",            0, +1),
        "amine_secondary": ("[NX3;H1;!$(NC=O);!$(Nc)]",            0, +1),
        "amine_tertiary":  ("[NX3;H0;!$(NC=O);!$([n]);!$(Nc)]",    0, +1),
    }

    state_map = {g["type"]: g["state_at_ph"] for g in groups}
    rw = RWMol(mol)

    for gtype, state in state_map.items():
        if state != "ionized":
            continue
        if gtype not in TYPE_TO_SMARTS:
            continue
        smarts, match_pos, charge_delta = TYPE_TO_SMARTS[gtype]
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue
        matches = mol.GetSubstructMatches(pattern)
        for match in matches:
            atom_idx = match[match_pos]
            atom = rw.GetAtomWithIdx(atom_idx)
            atom.SetFormalCharge(charge_delta)
            # Acids: remove the implicit H by forcing no-implicit
            if charge_delta < 0:
                atom.SetNoImplicit(True)
                atom.SetNumExplicitHs(0)

    try:
        Chem.SanitizeMol(rw)
        return Chem.MolToSmiles(rw)
    except Exception:
        # Fallback: Uncharger for net-neutral, Reionizer otherwise
        net = sum(g["charge_contribution"] for g in groups)
        if net == 0:
            uncharger = rdMolStandardize.Uncharger()
            return Chem.MolToSmiles(uncharger.uncharge(mol))
        reionizer = rdMolStandardize.Reionizer()
        return Chem.MolToSmiles(reionizer.reionize(mol))


def main():
    parser = argparse.ArgumentParser(
        description="Estimate protonation states at a given pH.")
    inp = parser.add_mutually_exclusive_group(required=True)
    inp.add_argument("--smiles", help="Input SMILES")
    inp.add_argument("--sdf", help="Input SDF file")
    parser.add_argument("--ph", type=float, default=7.4,
                        help="Target pH (default: 7.4)")
    parser.add_argument("--out", help="Write JSON to file")
    args = parser.parse_args()

    from rdkit import Chem

    mol = load_mol(args)
    input_smi = Chem.MolToSmiles(mol)
    ph = args.ph

    groups, net_charge = estimate_pka(mol, ph)
    dominant_form = build_dominant_form(mol, groups)

    result = {
        "input_smiles": input_smi,
        "ph": ph,
        "dominant_form": dominant_form,
        "ionization_groups": groups,
        "net_charge_at_ph": net_charge,
    }

    output = json.dumps(result, indent=2)
    if args.out:
        with open(args.out, "w") as f:
            f.write(output)
    else:
        print(output)


if __name__ == "__main__":
    main()
