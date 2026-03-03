#!/usr/bin/env python3
"""ALKYL — chem_diversity.py: MaxMin diversity selection from a molecular library."""

import argparse
import json
import sys
from pathlib import Path


def read_smi(path: str):
    """Parse a .smi file: 'SMILES [name]' per line. Returns list of (smiles, name)."""
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


def compute_fingerprints(mols_smiles, fp_type: str):
    """Return (valid_entries, fingerprints) filtering out None mols."""
    from rdkit import Chem
    from rdkit.Chem import AllChem, MACCSkeys

    valid = []
    fps = []
    for idx, (smi, name) in enumerate(mols_smiles):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        if fp_type == "morgan":
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
        elif fp_type == "maccs":
            fp = MACCSkeys.GenMACCSKeys(mol)
        else:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
        valid.append((idx, smi, name, mol))
        fps.append(fp)
    return valid, fps


def maxmin_select(fps: list, n: int) -> list[int]:
    """Return list indices (into fps) of n maximally diverse molecules."""
    from rdkit import DataStructs

    if n >= len(fps):
        return list(range(len(fps)))

    selected = [0]
    remaining = list(range(1, len(fps)))

    while len(selected) < n and remaining:
        max_min_dist = -1.0
        best_pos = 0
        for pos, idx in enumerate(remaining):
            min_sim = min(
                DataStructs.TanimotoSimilarity(fps[idx], fps[s])
                for s in selected
            )
            dist = 1.0 - min_sim
            if dist > max_min_dist:
                max_min_dist = dist
                best_pos = pos
        selected.append(remaining[best_pos])
        remaining.pop(best_pos)

    return selected


def main():
    parser = argparse.ArgumentParser(
        description="MaxMin diversity selection from a .smi library.")
    parser.add_argument("--input", required=True, help="Input .smi file")
    parser.add_argument("--n", type=int, required=True,
                        help="Number of diverse molecules to select")
    parser.add_argument("--fingerprint", choices=["morgan", "maccs"],
                        default="morgan")
    parser.add_argument("--out")
    args = parser.parse_args()

    from rdkit import Chem

    raw = read_smi(args.input)
    if not raw:
        print("No molecules found in input.", file=sys.stderr)
        sys.exit(1)

    valid, fps = compute_fingerprints(raw, args.fingerprint)
    if not valid:
        print("No valid molecules after parsing.", file=sys.stderr)
        sys.exit(1)

    selected_positions = maxmin_select(fps, args.n)

    selected_out = []
    for pos in selected_positions:
        orig_idx, smi, name, mol = valid[pos]
        canonical = Chem.MolToSmiles(mol)
        selected_out.append({"idx": orig_idx, "name": name, "smiles": canonical})

    result = {
        "n_requested": args.n,
        "n_selected": len(selected_out),
        "n_library": len(valid),
        "fingerprint": args.fingerprint,
        "selected": selected_out,
    }

    output = json.dumps(result, indent=2)
    if args.out:
        with open(args.out, "w") as f:
            f.write(output)
    else:
        print(output)


if __name__ == "__main__":
    main()
