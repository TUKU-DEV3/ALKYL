#!/usr/bin/env python3
"""ALKYL — chem_convert.py: inter-format molecular I/O."""

import argparse
import json
import sys


def parse_input(args):
    from rdkit import Chem
    if args.smiles:
        mol = Chem.MolFromSmiles(args.smiles)
        if mol is None:
            print(f"Invalid SMILES: {args.smiles}", file=sys.stderr)
            sys.exit(1)
        return mol
    if args.sdf:
        suppl = Chem.SDMolSupplier(args.sdf)
        mol = next((m for m in suppl if m is not None), None)
        if mol is None:
            print(f"No valid molecule in: {args.sdf}", file=sys.stderr)
            sys.exit(1)
        return mol
    if args.inchi:
        from rdkit.Chem.inchi import MolFromInchi
        mol = MolFromInchi(args.inchi)
        if mol is None:
            print(f"Invalid InChI: {args.inchi}", file=sys.stderr)
            sys.exit(1)
        return mol
    print("No input provided. Use --smiles, --sdf, or --inchi.", file=sys.stderr)
    sys.exit(1)


def convert(mol, target: str) -> dict:
    from rdkit import Chem
    if target == "smiles":
        return {"smiles": Chem.MolToSmiles(mol)}
    if target == "inchi":
        from rdkit.Chem.inchi import MolToInchi, InchiToInchiKey
        inchi = MolToInchi(mol)
        return {"inchi": inchi, "inchikey": InchiToInchiKey(inchi)}
    if target == "sdf":
        return {"sdf": Chem.MolToMolBlock(mol)}
    if target == "svg":
        from rdkit.Chem.Draw import rdMolDraw2D
        drawer = rdMolDraw2D.MolDraw2DSVG(300, 200)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        return {"svg": drawer.GetDrawingText()}
    print(f"Unsupported target format: {target}", file=sys.stderr)
    sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description="Convert molecular formats.")
    inp = parser.add_mutually_exclusive_group(required=True)
    inp.add_argument("--smiles", help="Input SMILES string")
    inp.add_argument("--sdf", help="Input SDF file path")
    inp.add_argument("--inchi", help="Input InChI string")
    parser.add_argument("--to", required=True,
                        choices=["smiles", "inchi", "sdf", "svg"],
                        help="Target format")
    parser.add_argument("--out", help="Write output to file instead of stdout")
    args = parser.parse_args()

    mol = parse_input(args)
    result = convert(mol, args.to)

    output = json.dumps(result, indent=2)
    if args.out:
        with open(args.out, "w") as f:
            f.write(output)
    else:
        print(output)


if __name__ == "__main__":
    main()
