# Protocol — Molecule Evaluation (Rigid)

Use when performing a **complete characterization** of one or more molecules.
Follow steps in order; document outputs at each stage.

---

## Stage 1 — Structure Validation ⚡
```bash
python chem_standardize.py --smiles "SMILES"
# → canonical SMILES, desalted, neutralized, changes list
```
If changes reported → use the standardized SMILES for all subsequent steps.

## Stage 2 — 2D Properties ⚡
```bash
python chem_props.py --smiles "SMILES" --lipinski --pains
python chem_analyze.py --smiles "SMILES"
```
Key outputs to record:
- MW, LogP, TPSA, HBD, HBA (Lipinski)
- QED (0–1, higher = more drug-like)
- SA score (1–10, > 6 = synthesis concern)
- Functional groups, ring count, stereocenters

## Stage 3 — Drug-Likeness Filters ⚡
```bash
python chem_filter.py --smiles "SMILES"
# Lipinski / Veber / Egan / Ghose / PAINS
```
Interpret: any FAIL is a flag, not a disqualifier. Record which rules fail and why.

## Stage 4 — Ionization at Physiological pH ⚡
```bash
python chem_pka.py --smiles "SMILES" --ph 7.4
```
Record: dominant_form, net_charge_at_ph, ionization groups.
Charged molecules: note impact on membrane permeability.

## Stage 5 — Metabolic Liabilities ⚡
```bash
python chem_metabolism.py --smiles "SMILES"
```
Flag: any "high" risk sites → structural modification candidates.

## Stage 6 — External Bioactivity Lookup ⚡ (requires ChEMBL MCP)
```
ChEMBL.compound_search(query="canonical SMILES or name")
→ get ChEMBL ID

ChEMBL.get_bioactivity(molecule_chembl_id="CHEMBLXXX")
→ known IC50/EC50/Ki against targets

ChEMBL.get_mechanism(molecule_chembl_id="CHEMBLXXX")
→ mechanism of action if approved drug
```

## Stage 7 — 3D Conformation ⚡⚡ (optional, if needed for docking/QM)
```bash
python chem_3d.py --smiles "SMILES" --conformers 10 --minimize
```

## Stage 8 — Summary Table

| Property | Value | Flag? |
|----------|-------|-------|
| MW | | |
| LogP | | |
| TPSA | | |
| QED | | |
| SA score | | |
| Lipinski | PASS/FAIL | |
| PAINS | PASS/FAIL | |
| net_charge pH 7.4 | | |
| CYP risk | none/low/medium/high | |
| Known activity | (ChEMBL) | |
