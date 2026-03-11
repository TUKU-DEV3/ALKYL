#!/usr/bin/env python3
"""ALKYL — chem_rgroup.py: R-group decomposition of a molecular library around a SMARTS core."""

import argparse
import json
import sys
from pathlib import Path


def read_smi(path: str):
    entries = []
    for i, line in enumerate(Path(path).read_text().splitlines()):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        smiles = parts[0]
        name = parts[1] if len(parts) > 1 else f"mol_{i}"
        entries.append((smiles, name))
    return entries


def main():
    parser = argparse.ArgumentParser(
        description="R-group decomposition of a .smi library around a SMARTS core.")
    parser.add_argument("--input", required=True, help="Input .smi file")
    parser.add_argument("--core", required=True,
                        help="SMARTS pattern for the core scaffold (use [*:1], [*:2], ... for R-group attachment points)")
    parser.add_argument("--out", help="Write JSON to file instead of stdout")
    args = parser.parse_args()

    from rdkit import Chem
    from rdkit.Chem.rdRGroupDecomposition import RGroupDecompose, RGroupDecompositionParameters

    core_mol = Chem.MolFromSmarts(args.core)
    if core_mol is None:
        print(f"Invalid core SMARTS: {args.core}", file=sys.stderr)
        sys.exit(1)

    raw = read_smi(args.input)
    if not raw:
        print("No molecules found in input.", file=sys.stderr)
        sys.exit(1)

    mols, names, original_smiles = [], [], []
    for smi, name in raw:
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            mols.append(mol)
            names.append(name)
            original_smiles.append(Chem.MolToSmiles(mol))

    if not mols:
        print("No valid molecules after parsing.", file=sys.stderr)
        sys.exit(1)

    params = RGroupDecompositionParameters()
    rgroup_rows, unmatched_indices = RGroupDecompose(
        [core_mol], mols, asSmiles=True, options=params
    )

    # Collect R-group labels (exclude "Core")
    rgroup_labels = []
    if rgroup_rows:
        rgroup_labels = sorted(
            k for k in rgroup_rows[0].keys() if k != "Core"
        )

    # Build matched results
    matched_set = set()
    decomposition = []
    for row_dict in rgroup_rows:
        # row_dict keys: "Core", "R1", "R2", ...
        # We need to figure out which molecule index this corresponds to
        # RGroupDecompose returns rows in the same order as input mols,
        # excluding unmatched ones — but unmatched_indices gives the indices of unmatched
        pass

    # Better: rebuild by iterating with unmatched knowledge
    unmatched_set = set(unmatched_indices)
    matched_indices = [i for i in range(len(mols)) if i not in unmatched_set]

    decomposition = []
    for row_idx, mol_idx in enumerate(matched_indices):
        row = rgroup_rows[row_idx]
        entry = {
            "name": names[mol_idx],
            "smiles": original_smiles[mol_idx],
            "core_match": row.get("Core", ""),
        }
        for label in rgroup_labels:
            entry[label] = row.get(label, "")
        decomposition.append(entry)

    unmatched_out = [
        {"name": names[i], "smiles": original_smiles[i]}
        for i in unmatched_indices
    ]

    result = {
        "core_smarts": args.core,
        "n_input": len(mols),
        "n_matched": len(decomposition),
        "n_unmatched": len(unmatched_out),
        "rgroup_labels": rgroup_labels,
        "decomposition": decomposition,
        "unmatched": unmatched_out,
    }

    output = json.dumps(result, indent=2)
    if args.out:
        Path(args.out).write_text(output)
    else:
        print(output)


if __name__ == "__main__":
    main()
