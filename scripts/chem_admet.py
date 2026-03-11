#!/usr/bin/env python3
"""ALKYL — chem_admet.py: ADMET property prediction (ESOL solubility, BBB, hERG, P-gp, PPB)."""

import argparse
import json
import sys
from pathlib import Path


# ── ESOL (Delaney 2004) ──────────────────────────────────────────────────────

def _esol(mol) -> dict:
    """Estimate aqueous solubility via the Delaney (2004) equation."""
    from rdkit.Chem import Descriptors, rdMolDescriptors

    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    rotbonds = Descriptors.NumRotatableBonds(mol)
    ar_rings = rdMolDescriptors.CalcNumAromaticRings(mol)

    log_s = 0.16 - 0.63 * logp - 0.0062 * mw + 0.066 * rotbonds - 0.74 * ar_rings

    if log_s > -1:
        sol_class = "Highly soluble"
    elif log_s > -2:
        sol_class = "Soluble"
    elif log_s > -4:
        sol_class = "Moderately soluble"
    elif log_s > -6:
        sol_class = "Poorly soluble"
    else:
        sol_class = "Insoluble"

    return {
        "logS_esol": round(log_s, 3),
        "class": sol_class,
        "method": "Delaney 2004",
        "inputs": {
            "MW": round(mw, 2),
            "LogP": round(logp, 4),
            "RotBonds": rotbonds,
            "ArRings": ar_rings,
        },
    }


# ── BBB penetration (heuristic) ───────────────────────────────────────────────

def _bbb(mol) -> dict:
    """Estimate blood-brain barrier permeability using a multi-parameter heuristic."""
    from rdkit.Chem import Descriptors

    logp = Descriptors.MolLogP(mol)
    mw = Descriptors.MolWt(mol)
    tpsa = Descriptors.TPSA(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)

    # Each criterion contributes +1 to score (out of 5)
    criteria = {
        "LogP_1_to_3": 1.0 <= logp <= 3.0,
        "MW_lt_400": mw < 400,
        "TPSA_lt_90": tpsa < 90,
        "HBD_le_3": hbd <= 3,
        "HBA_le_7": hba <= 7,
    }
    score = sum(criteria.values())

    if score >= 4:
        prediction = "likely+"
    elif score == 3:
        prediction = "uncertain"
    else:
        prediction = "likely-"

    return {
        "prediction": prediction,
        "score": f"{score}/5",
        "criteria": {k: bool(v) for k, v in criteria.items()},
        "note": "Heuristic — not a validated QSAR model",
    }


# ── hERG risk ─────────────────────────────────────────────────────────────────

_HERG_ALERTS = [
    ("basic_amine_aliphatic",
     "[NH2,NH1,NH0+0;!$(NC=O);!$(NS=O);!$(N[c,n]);!$(N[CX3]=[CX3])]"),
    ("basic_amine_aromatic_adj",
     "[NH2,NH1,NH0;!$(NC=O)]cc"),
    ("piperidine_or_piperazine",
     "[NX3;R1;!$(NC=O)]1CCNCC1"),
    ("tertiary_amine_near_aryl",
     "[NX3;!$(NC=O);!$(N[c,n])]([CX4])[CX4]c"),
    ("quaternary_N",
     "[N+;!$(NC=O)]"),
    ("long_aliphatic_chain",
     "[CH2][CH2][CH2][CH2][CH2]"),
]


def _herg(mol) -> dict:
    from rdkit import Chem
    from rdkit.Chem import Descriptors

    logp = Descriptors.MolLogP(mol)
    alerts_found = []
    for name, smarts in _HERG_ALERTS:
        patt = Chem.MolFromSmarts(smarts)
        if patt and mol.HasSubstructMatch(patt):
            alerts_found.append(name)

    has_basic_n = any("amine" in a or "N" in a for a in alerts_found)
    if len(alerts_found) >= 2 or (has_basic_n and logp > 3.0):
        risk = "high"
    elif len(alerts_found) == 1 or (has_basic_n and logp > 1.5):
        risk = "medium"
    else:
        risk = "low"

    return {
        "risk": risk,
        "alerts": alerts_found,
        "n_alerts": len(alerts_found),
        "logP": round(logp, 4),
        "note": "SMARTS-based heuristic — not a validated hERG QSAR model",
    }


# ── P-glycoprotein substrate ──────────────────────────────────────────────────

def _pgp(mol) -> dict:
    from rdkit.Chem import Descriptors

    mw = Descriptors.MolWt(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)
    logp = Descriptors.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)

    # Criteria from Ecker & Chiba-Falek (simplified)
    flags = {
        "MW_gt_400": mw > 400,
        "HBD_plus_HBA_gt_7": (hbd + hba) > 7,
        "LogP_gt_4": logp > 4.0,
        "TPSA_30_to_120": 30 <= tpsa <= 120,
    }
    n_flags = sum(flags.values())
    substrate_likely = n_flags >= 2

    return {
        "substrate_likely": substrate_likely,
        "n_criteria_met": n_flags,
        "criteria": {k: bool(v) for k, v in flags.items()},
        "note": "Heuristic — not a validated P-gp QSAR model",
    }


# ── Plasma protein binding ────────────────────────────────────────────────────

def _ppb(mol) -> dict:
    from rdkit.Chem import Descriptors

    logp = Descriptors.MolLogP(mol)
    mw = Descriptors.MolWt(mol)
    tpsa = Descriptors.TPSA(mol)

    flags = {
        "LogP_gt_2": logp > 2.0,
        "MW_gt_500": mw > 500,
        "TPSA_lt_75": tpsa < 75,
    }
    n_flags = sum(flags.values())
    high_ppb_likely = n_flags >= 2

    return {
        "high_ppb_likely": high_ppb_likely,
        "predicted_ppb": ">90%" if high_ppb_likely else "<90% (tentative)",
        "n_criteria_met": n_flags,
        "criteria": {k: bool(v) for k, v in flags.items()},
        "note": "Crude heuristic — not a validated PPB model",
    }


# ── Main ──────────────────────────────────────────────────────────────────────

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


def main():
    parser = argparse.ArgumentParser(
        description="ADMET property prediction: ESOL solubility, BBB, hERG, P-gp, PPB.")
    inp = parser.add_mutually_exclusive_group(required=True)
    inp.add_argument("--smiles", help="Input SMILES")
    inp.add_argument("--sdf", help="Input SDF file (first molecule used)")
    parser.add_argument("--out", help="Write JSON to file instead of stdout")
    args = parser.parse_args()

    from rdkit import Chem

    mol = load_mol(args)

    result = {
        "smiles": Chem.MolToSmiles(mol),
        "solubility": _esol(mol),
        "bbb_penetration": _bbb(mol),
        "herg_risk": _herg(mol),
        "pgp_substrate": _pgp(mol),
        "ppb": _ppb(mol),
        "disclaimer": (
            "All predictions are heuristic/SMARTS-based. "
            "Use validated QSAR models (SwissADME, pkCSM, ADMETlab) for serious decisions."
        ),
    }

    output = json.dumps(result, indent=2)
    if args.out:
        Path(args.out).write_text(output)
    else:
        print(output)


if __name__ == "__main__":
    main()
