#!/usr/bin/env python3
"""ALKYL — chem_cluster.py: Butina clustering of a molecular library by Tanimoto similarity."""

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


def compute_fingerprints(entries, fp_type: str):
    from rdkit import Chem
    from rdkit.Chem import AllChem, MACCSkeys

    valid, fps = [], []
    for i, (smi, name) in enumerate(entries):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        if fp_type == "maccs":
            fp = MACCSkeys.GenMACCSKeys(mol)
        else:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
        valid.append((i, smi, name))
        fps.append(fp)
    return valid, fps


def butina_cluster(fps, cutoff: float):
    """Return list of tuples (centroid_idx, member_idx, ...) — RDKit Butina format."""
    from rdkit import DataStructs
    from rdkit.ML.Cluster import Butina

    dists = []
    for i in range(1, len(fps)):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        dists.extend(1.0 - s for s in sims)

    return Butina.ClusterData(dists, len(fps), cutoff, isDistData=True)


def main():
    parser = argparse.ArgumentParser(
        description="Butina clustering of a .smi library by Tanimoto similarity.")
    parser.add_argument("--input", required=True, help="Input .smi file")
    parser.add_argument("--cutoff", type=float, default=0.4,
                        help="Tanimoto distance cutoff (default 0.4 → similarity ≥ 0.6)")
    parser.add_argument("--fingerprint", choices=["morgan", "maccs"], default="morgan")
    parser.add_argument("--out", help="Write JSON to file instead of stdout")
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

    clusters_raw = butina_cluster(fps, args.cutoff)

    clusters_out = []
    n_singletons = 0
    for cluster_id, cluster in enumerate(clusters_raw):
        centroid_pos = cluster[0]
        orig_idx, smi, name = valid[centroid_pos]
        mol = Chem.MolFromSmiles(smi)
        centroid = {
            "pos": centroid_pos,
            "idx": orig_idx,
            "name": name,
            "smiles": Chem.MolToSmiles(mol),
        }
        members = []
        for pos in cluster[1:]:
            mi, ms, mn = valid[pos]
            members.append({"pos": pos, "idx": mi, "name": mn,
                            "smiles": Chem.MolToSmiles(Chem.MolFromSmiles(ms))})
        if len(cluster) == 1:
            n_singletons += 1
        clusters_out.append({
            "cluster_id": cluster_id,
            "size": len(cluster),
            "centroid": centroid,
            "members": members,
        })

    result = {
        "n_molecules": len(valid),
        "n_clusters": len(clusters_out),
        "n_singletons": n_singletons,
        "cutoff_distance": args.cutoff,
        "similarity_threshold": round(1.0 - args.cutoff, 4),
        "fingerprint": args.fingerprint,
        "clusters": clusters_out,
    }

    output = json.dumps(result, indent=2)
    if args.out:
        Path(args.out).write_text(output)
    else:
        print(output)


if __name__ == "__main__":
    main()
