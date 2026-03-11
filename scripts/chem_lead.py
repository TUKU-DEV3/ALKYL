#!/usr/bin/env python3
"""ALKYL — chem_lead.py: Ligand efficiency metrics (LE, LLE, BEI, LELP) for lead optimization."""

import argparse
import csv
import json
import math
import sys
from pathlib import Path

RT_LN10 = 1.364  # kcal/mol at 298 K (RT × ln10)


def _pic50(activity: float, unit: str) -> float:
    """Convert activity value to pIC50."""
    if unit == "pic50":
        return activity
    if unit == "nm":
        return 9.0 - math.log10(activity)
    if unit == "um":
        return 6.0 - math.log10(activity)
    if unit == "m":
        return -math.log10(activity)
    raise ValueError(f"Unknown unit: {unit}")


def _compute_metrics(mol, pic50: float) -> dict:
    from rdkit.Chem import Descriptors

    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)
    n_heavy = mol.GetNumHeavyAtoms()

    le = round(RT_LN10 * pic50 / n_heavy, 4) if n_heavy > 0 else None
    lle = round(pic50 - logp, 4)
    bei = round(pic50 / mw * 1000, 4) if mw > 0 else None
    lelp = round(logp / le, 4) if (le and le != 0) else None

    return {
        "pIC50": round(pic50, 4),
        "MW": round(mw, 2),
        "LogP": round(logp, 4),
        "TPSA": round(tpsa, 2),
        "n_heavy": n_heavy,
        "LE": le,
        "LLE": lle,
        "BEI": bei,
        "LELP": lelp,
    }


def _thresholds(m: dict) -> dict:
    """Flag metrics against standard thresholds."""
    flags = {}
    if m["LE"] is not None:
        flags["LE_acceptable"] = m["LE"] >= 0.30
    if m["LLE"] is not None:
        flags["LLE_acceptable"] = m["LLE"] >= 3.0
    if m["BEI"] is not None:
        flags["BEI_acceptable"] = m["BEI"] >= 12.0
    if m["LELP"] is not None:
        flags["LELP_acceptable"] = m["LELP"] <= 10.0
    return flags


def _round_summary(compounds: list) -> dict:
    """Compute per-round mean metrics if 'round' column is present."""
    from collections import defaultdict
    by_round = defaultdict(list)
    for c in compounds:
        r = c.get("round")
        if r is not None:
            by_round[r].append(c["metrics"])

    if not by_round:
        return {}

    summary = {}
    for r, metrics_list in sorted(by_round.items()):
        keys = ["pIC50", "LE", "LLE", "BEI", "LELP", "MW", "LogP"]
        avgs = {}
        for k in keys:
            vals = [m[k] for m in metrics_list if m.get(k) is not None]
            avgs[f"mean_{k}"] = round(sum(vals) / len(vals), 4) if vals else None
        avgs["n"] = len(metrics_list)
        summary[str(r)] = avgs

    # Compute deltas between consecutive rounds
    round_keys = sorted(summary.keys(), key=lambda x: (isinstance(x, str), x))
    for i in range(1, len(round_keys)):
        prev, curr = round_keys[i - 1], round_keys[i]
        for k in ["mean_pIC50", "mean_LE", "mean_LLE"]:
            pv = summary[prev].get(k)
            cv = summary[curr].get(k)
            if pv is not None and cv is not None:
                summary[curr][f"delta_{k[5:]}"] = round(cv - pv, 4)

    return summary


def main():
    parser = argparse.ArgumentParser(
        description="Compute ligand efficiency metrics (LE, LLE, BEI, LELP) from a CSV of structures and activities.")
    parser.add_argument("--csv", required=True, help="Input CSV with SMILES and activity columns")
    parser.add_argument("--smiles-col", default="smiles", help="Column name for SMILES (default: smiles)")
    parser.add_argument("--activity-col", default=None,
                        help="Column name for activity (auto-detected if omitted)")
    parser.add_argument("--unit", choices=["nm", "um", "m", "pic50"], default="nm",
                        help="Activity unit (default: nm). Use pic50 if already log-transformed.")
    parser.add_argument("--name-col", default=None, help="Column name for molecule name (optional)")
    parser.add_argument("--round-col", default=None, help="Column name for round/cycle (optional)")
    parser.add_argument("--out", help="Write JSON to file instead of stdout")
    args = parser.parse_args()

    from rdkit import Chem

    rows = list(csv.DictReader(Path(args.csv).open()))
    if not rows:
        print("Empty CSV.", file=sys.stderr)
        sys.exit(1)

    # Auto-detect activity column
    act_col = args.activity_col
    if act_col is None:
        candidates = ["pIC50", "pic50", "IC50", "ic50", "pKd", "pKi", "activity", "Ki", "Kd"]
        for c in candidates:
            if c in rows[0]:
                act_col = c
                break
        if act_col is None:
            print(f"Cannot auto-detect activity column. Available: {list(rows[0].keys())}",
                  file=sys.stderr)
            sys.exit(1)

    compounds = []
    errors = []
    for i, row in enumerate(rows):
        smi = row.get(args.smiles_col, "").strip()
        act_str = row.get(act_col, "").strip()
        name = row.get(args.name_col, f"mol_{i}") if args.name_col else f"mol_{i}"
        round_val = row.get(args.round_col) if args.round_col else None

        if not smi or not act_str:
            errors.append({"row": i, "reason": "missing SMILES or activity"})
            continue

        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            errors.append({"row": i, "name": name, "reason": "invalid SMILES"})
            continue

        try:
            activity = float(act_str)
        except ValueError:
            errors.append({"row": i, "name": name, "reason": f"non-numeric activity: {act_str}"})
            continue

        if activity <= 0:
            errors.append({"row": i, "name": name, "reason": "activity ≤ 0 (cannot log-transform)"})
            continue

        pic50 = _pic50(activity, args.unit)
        metrics = _compute_metrics(mol, pic50)
        thresholds = _thresholds(metrics)

        compounds.append({
            "name": name,
            "smiles": Chem.MolToSmiles(mol),
            "activity_raw": activity,
            "activity_unit": args.unit,
            "round": round_val,
            "metrics": metrics,
            "thresholds": thresholds,
        })

    round_summary = _round_summary(compounds)

    result = {
        "activity_column": act_col,
        "unit": args.unit,
        "n_input": len(rows),
        "n_computed": len(compounds),
        "n_errors": len(errors),
        "thresholds_reference": {
            "LE": "≥ 0.30 kcal/mol per heavy atom",
            "LLE": "≥ 3.0",
            "BEI": "≥ 12.0 (pIC50 per kDa)",
            "LELP": "≤ 10.0",
        },
        "compounds": compounds,
        "round_summary": round_summary,
        "errors": errors,
    }

    output = json.dumps(result, indent=2)
    if args.out:
        Path(args.out).write_text(output)
    else:
        print(output)


if __name__ == "__main__":
    main()
