# Homology Modeling — Model Quality Assessment

Validate your model before using it for docking or MD. A bad model gives misleading results.

---

## Overview: Validation Hierarchy

```
1. DOPE score (MODELLER)   → overall energy, compare models of SAME target
2. pLDDT (AlphaFold2)      → per-residue AI confidence
3. Ramachandran plot        → backbone dihedral geometry
4. Clash score (MolProbity) → steric clashes (all-atom)
5. Rotamer outliers         → side-chain geometry
6. RMSD to template/exp.   → when experimental structure available
```

---

## DOPE Score (MODELLER)

DOPE (Discrete Optimized Protein Energy) is a statistical potential. **Lower = better.**

```python
from modeller import Environ
from modeller.scripts import complete_pdb

env = Environ()
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')

def score_model_dope(pdb_file: str) -> float:
    """Compute DOPE score for a model PDB."""
    mdl = complete_pdb(env, pdb_file)
    from modeller import Selection
    sel = Selection(mdl)
    score = sel.assess('DOPE')
    return score

# Compare multiple models of the same target
model_files = [
    'target.B99990001.pdb',
    'target.B99990002.pdb',
    'target.B99990003.pdb',
]
scores = [(score_model_dope(f), f) for f in model_files]
scores.sort()
best_dope, best_file = scores[0]
print(f"Best model: {best_file}")
print(f"DOPE: {best_dope:.2f}  (more negative = better)")
```

**Interpretation:**
- DOPE scores are only meaningful when **comparing models of the same protein**
- A DOPE score > 0 almost always indicates very poor geometry
- Typical good models: −40 000 to −90 000 (depends on sequence length)
- Per-residue DOPE > −0.1 → loop needing refinement

---

## pLDDT from B-Factor Column (AlphaFold2)

```python
from Bio.PDB import PDBParser
import numpy as np

def assess_plddt(pdb_file: str, pocket_residues: list = None) -> dict:
    """
    Compute pLDDT statistics from AF2 PDB.
    pocket_residues: list of residue numbers to check specifically.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("model", pdb_file)

    plddt_all    = []
    plddt_pocket = []

    for model in structure:
        for chain in model:
            for res in chain:
                if res.id[0] != ' ':
                    continue
                if 'CA' not in res:
                    continue
                val = res['CA'].get_bfactor()
                plddt_all.append(val)
                if pocket_residues and res.id[1] in pocket_residues:
                    plddt_pocket.append(val)
        break

    plddt_all = np.array(plddt_all)
    report = {
        'mean_plddt'      : float(plddt_all.mean()),
        'fraction_high'   : float((plddt_all > 90).mean()),
        'fraction_medium' : float(((plddt_all >= 70) & (plddt_all <= 90)).mean()),
        'fraction_low'    : float((plddt_all < 70).mean()),
        'n_disordered'    : int((plddt_all < 50).sum()),
        'docking_ready'   : bool(plddt_all.mean() > 70),
    }
    if plddt_pocket:
        arr = np.array(plddt_pocket)
        report['pocket_mean_plddt'] = float(arr.mean())
        report['pocket_docking_ok'] = bool(arr.mean() > 80)

    return report

report = assess_plddt("best_model.pdb", pocket_residues=[45, 67, 89, 120, 155])
for k, v in report.items():
    print(f"  {k}: {v}")
```

---

## Ramachandran Plot (Biopython + Matplotlib)

```python
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import PPBuilder, three_to_one
import matplotlib.pyplot as plt
import numpy as np

def ramachandran_plot(pdb_file: str, out_png: str = "ramachandran.png") -> dict:
    """
    Compute phi/psi dihedral angles and plot Ramachandran diagram.
    Returns statistics: fraction in favored, allowed, outlier regions.
    """
    parser  = PDBParser(QUIET=True)
    builder = PPBuilder()
    structure = parser.get_structure("model", pdb_file)

    phi_list, psi_list = [], []

    for model in structure:
        for chain in model:
            for pp in builder.build_peptides(chain):
                angles = pp.get_phi_psi_list()
                for phi, psi in angles:
                    if phi is not None and psi is not None:
                        phi_list.append(np.degrees(phi))
                        psi_list.append(np.degrees(psi))
        break

    phi = np.array(phi_list)
    psi = np.array(psi_list)

    # Favored regions (approximate Ramachandran cores)
    def in_helix(p, s):
        return (-100 < p < -30) and (-70 < s < -10)

    def in_sheet(p, s):
        return ((-160 < p < -50) and (90 < s < 180)) or \
               ((-160 < p < -50) and (-180 < s < -150))

    n_total    = len(phi)
    n_favored  = sum(1 for p, s in zip(phi, psi) if in_helix(p, s) or in_sheet(p, s))
    frac_fav   = n_favored / n_total * 100

    # Plot
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.scatter(phi, psi, s=4, alpha=0.6, c='steelblue')

    # Draw approximate favored regions
    helix_box = plt.Rectangle((-100, -70), 70, 60,
                               fill=True, facecolor='green', alpha=0.1,
                               label='Alpha-helix region')
    ax.add_patch(helix_box)
    beta_box  = plt.Rectangle((-160, 90),  110, 90,
                               fill=True, facecolor='blue', alpha=0.1,
                               label='Beta-sheet region')
    ax.add_patch(beta_box)

    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)
    ax.axhline(0, color='k', linewidth=0.5)
    ax.axvline(0, color='k', linewidth=0.5)
    ax.set_xlabel('Phi (°)')
    ax.set_ylabel('Psi (°)')
    ax.set_title(f'Ramachandran Plot — {n_total} residues, {frac_fav:.1f}% favored')
    ax.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    print(f"Saved: {out_png}")

    return {
        'n_residues'     : n_total,
        'pct_favored'    : round(frac_fav, 1),
        'pct_outlier'    : round((n_total - n_favored) / n_total * 100, 1),
    }

stats = ramachandran_plot("best_model.pdb", "rama.png")
print(stats)
# Good model: > 90% favored  |  Excellent: > 95%
```

**Thresholds (MolProbity standards):**
- Favored > 98%: excellent (crystal structure quality)
- Favored > 90%: acceptable for MD starting point
- Favored < 85%: poor — refine loops before use

---

## MolProbity / PROCHECK (CLI)

```bash
# MolProbity command-line (if installed via phenix or standalone)
molprobity.clashscore model.pdb
molprobity.ramalyze  model.pdb
molprobity.rotalyze  model.pdb

# Full MolProbity report
phenix.molprobity model.pdb

# Alternatively, upload to web server: https://molprobity.biochem.duke.edu/
```

```python
import subprocess, json, re

def run_molprobity(pdb_file: str) -> dict:
    """Run MolProbity clashscore (requires phenix or standalone install)."""
    try:
        result = subprocess.run(
            ['molprobity.clashscore', pdb_file],
            capture_output=True, text=True, timeout=120
        )
        output = result.stdout

        # Parse clashscore
        clash_match = re.search(r'clashscore\s*=\s*([\d.]+)', output)
        clash = float(clash_match.group(1)) if clash_match else None

        return {'clashscore': clash,
                'raw': output}
    except FileNotFoundError:
        print("MolProbity not installed. Use web server: molprobity.biochem.duke.edu")
        return {}

# Clashscore interpretation:
#   0–10  : excellent
#   10–20 : acceptable
#   > 30  : problematic (serious steric clashes)
```

---

## Sequence Identity and Template Coverage Rules

```python
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

def assess_template_reliability(target_seq: str, template_seq: str) -> dict:
    """
    Compute sequence identity and coverage. Rough reliability estimate.
    """
    alignments = pairwise2.align.globalms(
        target_seq, template_seq,
        2, -1, -3, -0.5,   # match, mismatch, gap open, gap extend
        penalize_end_gaps=False
    )
    if not alignments:
        return {'identity': 0.0, 'coverage': 0.0, 'zone': 'safe bet'}

    aln = alignments[0]
    seq1, seq2 = aln.seqA, aln.seqB

    matches   = sum(a == b and a != '-' for a, b in zip(seq1, seq2))
    aligned   = sum(1 for a, b in zip(seq1, seq2) if a != '-' and b != '-')
    identity  = matches / aligned * 100 if aligned > 0 else 0
    coverage  = sum(1 for c in seq1 if c != '-') / len(target_seq) * 100

    # Model reliability zones
    if identity > 50:
        zone = 'safe bet — accurate backbone expected'
    elif identity > 30:
        zone = 'reliable — backbone good, some side-chain uncertainty'
    elif identity > 25:
        zone = 'twilight zone — use multi-template + careful validation'
    else:
        zone = 'dark zone — AI prediction (AlphaFold2) preferred'

    return {
        'identity' : round(identity, 1),
        'coverage' : round(coverage, 1),
        'zone'     : zone,
    }

result = assess_template_reliability("MKTAYIAKQR...", "MKTAYIAKQR...")
print(result)
```

---

## RMSD to Experimental Structure (Biopython)

```python
from Bio.PDB import PDBParser, Superimposer
import numpy as np

def compute_backbone_rmsd(model_pdb: str, ref_pdb: str,
                           chain_id: str = 'A') -> float:
    """
    Superimpose model onto reference, return CA RMSD.
    Both structures must have same number of residues in chain_id.
    """
    parser = PDBParser(QUIET=True)
    model_struct = parser.get_structure("model", model_pdb)
    ref_struct   = parser.get_structure("ref",   ref_pdb)

    model_chain = model_struct[0][chain_id]
    ref_chain   = ref_struct[0][chain_id]

    # Collect CA atoms for common residues
    model_cas, ref_cas = [], []
    ref_residues = {r.id[1]: r for r in ref_chain if r.id[0] == ' '}

    for res in model_chain:
        if res.id[0] != ' ':
            continue
        resnum = res.id[1]
        if resnum in ref_residues and 'CA' in res and 'CA' in ref_residues[resnum]:
            model_cas.append(res['CA'])
            ref_cas.append(ref_residues[resnum]['CA'])

    if len(model_cas) < 10:
        raise ValueError(f"Too few common CA atoms: {len(model_cas)}")

    sup = Superimposer()
    sup.set_atoms(ref_cas, model_cas)
    sup.apply(model_struct.get_atoms())

    rmsd = sup.rms
    print(f"Backbone CA RMSD: {rmsd:.2f} Å  ({len(model_cas)} residues aligned)")
    return rmsd

rmsd = compute_backbone_rmsd("homology_model.pdb", "crystal_structure.pdb", "A")
# < 1.5 Å: excellent  |  1.5–3.0 Å: acceptable  |  > 3.0 Å: poor
```

---

## ProDy Structural Analysis

```python
import prody

# Parse and analyze
model = prody.parsePDB("best_model.pdb")
crystal = prody.parsePDB("reference.pdb")

# Match chains and compute RMSD
matches = prody.matchChains(model, crystal, subset='calpha', overlap=50)
if matches:
    mdl_matched, ref_matched = matches[0][0], matches[0][1]
    rmsd = prody.calcRMSD(mdl_matched, ref_matched)
    print(f"RMSD: {rmsd:.2f} Å")

# Superimpose
prody.superpose(mdl_matched, ref_matched)

# Write superimposed model
prody.writePDB("model_superimposed.pdb", model)

# B-factor / pLDDT profile
ca = model.select('calpha')
bfactors = ca.getBetas()
print(f"Mean pLDDT: {bfactors.mean():.1f}")
print(f"Min  pLDDT: {bfactors.min():.1f}")
```

---

## Validation Checklist

```
[ ] DOPE score < 0 (MODELLER models only)
[ ] Mean pLDDT > 70 (AF2/ESMFold models)
[ ] Pocket residues pLDDT > 80 (before docking)
[ ] Ramachandran favored > 90%
[ ] Clashscore < 20 (MolProbity)
[ ] Rotamer outliers < 5%
[ ] No chain breaks in pocket region
[ ] Backbone RMSD < 3.0 Å vs template (if available)
```
