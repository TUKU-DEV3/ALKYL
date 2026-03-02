#!/usr/bin/env python3
"""ALKYL — chem_props.py: physico-chemical descriptors, drug-likeness, fingerprints."""

import argparse
import json
import sys

DESCRIPTOR_MAP = {
    "mw":       ("rdkit.Chem.Descriptors", "MolWt"),      # average MW (Lipinski convention)
    "exact_mw": ("rdkit.Chem.Descriptors", "ExactMolWt"), # monoisotopic MW
    "logp":     ("rdkit.Chem.Descriptors", "MolLogP"),
    "hbd":      ("rdkit.Chem.Descriptors", "NumHDonors"),
    "hba":      ("rdkit.Chem.Descriptors", "NumHAcceptors"),
    "tpsa":     ("rdkit.Chem.Descriptors", "TPSA"),
    "rotbonds": ("rdkit.Chem.Descriptors", "NumRotatableBonds"),
    "rings":    ("rdkit.Chem.Descriptors", "RingCount"),
    "fsp3":     ("rdkit.Chem.Descriptors", "FractionCSP3"),
}
ALL_DESCRIPTORS = list(DESCRIPTOR_MAP.keys())


def load_mol(args):
    from rdkit import Chem
    if args.smiles:
        mol = Chem.MolFromSmiles(args.smiles)
    elif args.sdf:
        suppl = Chem.SDMolSupplier(args.sdf)
        mol = next((m for m in suppl if m is not None), None)
    else:
        print("Provide --smiles or --sdf.", file=sys.stderr)
        sys.exit(1)
    if mol is None:
        print("Invalid molecule input.", file=sys.stderr)
        sys.exit(1)
    return mol


def calc_descriptors(mol, names: list[str]) -> dict:
    import importlib
    result = {}
    for name in names:
        if name not in DESCRIPTOR_MAP:
            continue
        module_path, func_name = DESCRIPTOR_MAP[name]
        module = importlib.import_module(module_path)
        result[name] = round(getattr(module, func_name)(mol), 4)
    return result


def check_lipinski(mol) -> dict:
    from rdkit.Chem import Descriptors
    violations = []
    if Descriptors.MolWt(mol) > 500:
        violations.append("MW > 500")
    if Descriptors.MolLogP(mol) > 5:
        violations.append("LogP > 5")
    if Descriptors.NumHDonors(mol) > 5:
        violations.append("HBD > 5")
    if Descriptors.NumHAcceptors(mol) > 10:
        violations.append("HBA > 10")
    return {"pass": len(violations) == 0, "violations": violations}


def check_pains(mol) -> dict:
    from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog(params)
    matches = catalog.GetMatches(mol)
    alerts = [m.GetDescription() for m in matches]
    return {"alerts": alerts, "clean": len(alerts) == 0}


def calc_fingerprint(mol, fp_type: str, radius: int, nbits: int) -> dict:
    from rdkit.Chem import MACCSkeys
    if fp_type == "morgan":
        from rdkit.Chem import rdFingerprintGenerator
        gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=nbits)
        fp = gen.GetFingerprint(mol)
    elif fp_type == "maccs":
        fp = MACCSkeys.GenMACCSKeys(mol)
    else:
        from rdkit.Chem.rdmolops import RDKFingerprint
        fp = RDKFingerprint(mol, fpSize=nbits)
    on_bits = list(fp.GetOnBits())
    return {"type": fp_type, "nbits": fp.GetNumBits(), "bits": on_bits,
            "density": round(len(on_bits) / fp.GetNumBits(), 4)}


def main():
    parser = argparse.ArgumentParser(description="Compute molecular properties.")
    inp = parser.add_mutually_exclusive_group(required=True)
    inp.add_argument("--smiles", help="Input SMILES")
    inp.add_argument("--sdf", help="Input SDF file")
    parser.add_argument("--descriptors",
                        help="Comma-separated descriptors or 'all'")
    parser.add_argument("--lipinski", action="store_true")
    parser.add_argument("--pains", action="store_true")
    parser.add_argument("--fingerprint", choices=["morgan", "maccs", "rdkit"])
    parser.add_argument("--radius", type=int, default=2)
    parser.add_argument("--nbits", type=int, default=2048)
    args = parser.parse_args()

    mol = load_mol(args)
    result = {}

    if args.descriptors:
        names = ALL_DESCRIPTORS if args.descriptors == "all" \
            else [d.strip() for d in args.descriptors.split(",")]
        result.update(calc_descriptors(mol, names))

    if args.lipinski:
        result["lipinski"] = check_lipinski(mol)

    if args.pains:
        result["pains"] = check_pains(mol)

    if args.fingerprint:
        result["fingerprint"] = calc_fingerprint(
            mol, args.fingerprint, args.radius, args.nbits)

    if not result:
        parser.error("Specify at least one of: --descriptors, --lipinski, "
                     "--pains, --fingerprint")

    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
