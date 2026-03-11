#!/usr/bin/env python3
"""ALKYL — chem_highlight.py: SVG/PNG of a molecule with substructure highlighted via SMARTS."""

import argparse
import json
import sys
from pathlib import Path


def highlight_mol(mol, smarts: str | None, width: int, height: int, fmt: str) -> bytes | str:
    from rdkit.Chem.Draw import rdMolDraw2D

    highlight_atoms = []
    highlight_bonds = []

    if smarts:
        from rdkit import Chem
        patt = Chem.MolFromSmarts(smarts)
        if patt is None:
            print(f"Invalid SMARTS: {smarts}", file=sys.stderr)
            sys.exit(1)
        matches = mol.GetSubstructMatches(patt)
        if not matches:
            print(f"SMARTS pattern found no matches.", file=sys.stderr)
            # Continue — draw without highlight rather than exit
        else:
            atom_set = set()
            for match in matches:
                atom_set.update(match)
            highlight_atoms = list(atom_set)

            # Collect bonds between highlighted atoms
            bond_set = set()
            for bond in mol.GetBonds():
                a1, a2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                if a1 in atom_set and a2 in atom_set:
                    bond_set.add(bond.GetIdx())
            highlight_bonds = list(bond_set)

    if fmt == "png":
        drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
    else:
        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)

    drawer.drawOptions().addAtomIndices = False
    drawer.DrawMolecule(
        mol,
        highlightAtoms=highlight_atoms if highlight_atoms else [],
        highlightBonds=highlight_bonds if highlight_bonds else [],
    )
    drawer.FinishDrawing()

    if fmt == "png":
        return drawer.GetDrawingText()  # bytes
    return drawer.GetDrawingText()      # str


def main():
    parser = argparse.ArgumentParser(
        description="Draw a molecule as SVG/PNG with an optional SMARTS substructure highlighted.")
    parser.add_argument("--smiles", required=True, help="Input SMILES")
    parser.add_argument("--smarts", help="SMARTS pattern to highlight (optional)")
    parser.add_argument("--out", help="Output file (.svg or .png). Prints SVG to stdout if omitted.")
    parser.add_argument("--width", type=int, default=400)
    parser.add_argument("--height", type=int, default=300)
    args = parser.parse_args()

    from rdkit import Chem

    mol = Chem.MolFromSmiles(args.smiles)
    if mol is None:
        print(f"Invalid SMILES: {args.smiles}", file=sys.stderr)
        sys.exit(1)

    fmt = "svg"
    if args.out and args.out.lower().endswith(".png"):
        fmt = "png"

    drawing = highlight_mol(mol, args.smarts, args.width, args.height, fmt)

    if args.out:
        out_path = Path(args.out)
        if fmt == "png":
            out_path.write_bytes(drawing)
        else:
            out_path.write_text(drawing)
        meta = {
            "smiles": Chem.MolToSmiles(mol),
            "smarts": args.smarts,
            "format": fmt,
            "output": args.out,
            "width": args.width,
            "height": args.height,
        }
        print(json.dumps(meta, indent=2))
    else:
        # SVG to stdout
        if fmt == "png":
            print("Use --out file.png to write PNG. Stdout only supports SVG.", file=sys.stderr)
            sys.exit(1)
        print(drawing)


if __name__ == "__main__":
    main()
