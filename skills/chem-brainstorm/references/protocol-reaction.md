# Protocol — Reaction Evaluation (Rigid)

Use when assessing a chemical reaction: feasibility, conditions, mechanism, or retrosynthesis.

---

## Stage 1 — Define the Transformation

Encode as reaction SMILES if possible:
```
reactants >> products
[reactant1].[reactant2]>>[product]
```

If conditions known: note solvent, base, temperature, catalyst.

## Stage 2 — Apply Reaction SMARTS ⚡

For a known transformation:
```bash
python chem_react.py --smiles "reactant_SMILES" --reaction "SMARTS"
# → products (deduplicated, sanitized)
# load daylight-theory references/smirks.md to write correct SMIRKS
```

## Stage 3 — Product Analysis ⚡

For each product:
```bash
python chem_analyze.py --smiles "product_SMILES"
# → SA score: if > 6, product itself is hard to isolate or purify
# → stereocenter count: uncontrolled stereocenters = problem

python chem_enum.py --smiles "product_SMILES"
# → if multiple stereocenters: enumerate all possible diastereomers
```

## Stage 4 — Tautomer / Ionization Check ⚡

```bash
python chem_tautomers.py --smiles "product_SMILES"
# → is the drawn tautomeric form the stable one?

python chem_pka.py --smiles "product_SMILES" --ph 7.4
# → will the product be charged under workup conditions?
```

## Stage 5 — Literature Precedent ⚡ (ChEMBL + bioRxiv)

```
ChEMBL.compound_search(query="product name or SMILES")
→ is this compound known? Any reported synthesis?

bioRxiv.search_preprints(query="[reaction type] conditions selectivity")
→ recent methodological developments
```

## Stage 6 — Retrosynthesis (if forward route unclear)

Use TorchDrug skill (`torchdrug`) for ML-based retrosynthesis:
```python
# load torchdrug skill → references/retrosynthesis.md
```

Or SynKit for rule-based approach:
```python
# load synkit skill → references/rule-dpo.md
# DPO rule extraction + SynReactor
```

## Stage 7 — QM Feasibility (for novel or strained systems) ⚡⚡⚡

```bash
python chem_qm.py --smiles "reactant" --engine orca --task opt --method B3LYP --basis 6-31G*
# → geometry optimization → confirm ground state structure

python chem_qm.py --smiles "product" --engine orca --task freq
python chem_qm.py --parse output.log --parse-ir
# → n_imaginary=0 confirms minimum; n_imaginary=1 = transition state
```

## Stage 8 — Reaction Summary

| Aspect | Assessment |
|--------|-----------|
| Feasibility | likely / uncertain / unlikely |
| Key concern | (selectivity / side reactions / stability) |
| Precedent | (known / analogous / novel) |
| SA score (product) | |
| Stereocenters controlled? | yes / no / N/A |
| Next step | (run / optimize / redesign) |
