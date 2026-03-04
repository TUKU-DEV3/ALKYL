# Protocol — Hypothesis Construction (Rigid)

Use when building mechanistic or SAR hypotheses from scratch or from data.

---

## Stage 1 — Define the Observation

State clearly:
- What is observed? (activity, selectivity, toxicity, property change)
- What is the chemical space? (scaffold, series, targets)
- What data exists? (IC50 table, selectivity profile, in vivo result)

## Stage 2 — Target & Pathway Context ⚡ (MCPs)

```
OpenTargets.search_entities(query="[gene/protein]")
→ genetic evidence score, disease associations, tractability

ChEMBL.target_search(query="[gene symbol]")
→ target class, known ligands, reference compounds

ChEMBL.get_bioactivity(target_chembl_id="CHEMBLXXX", limit=50)
→ activity landscape: what binds, what doesn't
```

## Stage 3 — Literature Scan ⚡ (bioRxiv MCP)

```
bioRxiv.search_preprints(query="[target] [mechanism] inhibitor", limit=10)
→ recent methods, negative results, competing hypotheses

bioRxiv.search_preprints(query="[scaffold] selectivity mechanism")
→ SAR trends in the field
```

## Stage 4 — Structural Analysis of Active Compounds ⚡

For each key active compound:
```bash
python chem_analyze.py --smiles "SMILES"
# → FG, rings, stereocenters, complexity

python chem_scaffold.py --smiles "SMILES"
# → Murcko scaffold (common pharmacophore?)

python chem_compare.py --smiles-a "active" --smiles-b "inactive"
# → MCS, Δ properties → what structural feature drives activity?
```

## Stage 5 — SMARTS-Based SAR Hypothesis

Formulate hypothesis as a SMARTS pattern (load `daylight-theory` skill):
```
"All actives contain [pattern X] → test with chem_search.py --mode substructure"
```

```bash
python chem_search.py --query "[SMARTS]" --library actives.sdf --mode substructure
# confirm pattern present in actives, absent in inactives
```

## Stage 6 — Generate Hypotheses

Write 2–3 hypotheses in the form:
```
H1: [Structural feature] is responsible for [observed effect] because [mechanism].
   Evidence: [data points]
   Test: [specific experiment or calculation]
   Falsifiable by: [counter-experiment]

H2: ...
H3: ...
```

## Stage 7 — Prioritize

Rank by:
- **Evidence strength** (literature > SAR data > structural analogy)
- **Testability** (can be tested with available ALKYL tools vs. requires wet lab)
- **Impact** (explains the most observations)

## Stage 8 — Clinical Relevance (if applicable)

```
ClinicalTrials.search_trials(condition="[disease]", intervention="[drug class]")
→ are there trials targeting this mechanism?

ClinicalTrials.analyze_endpoints(condition="[disease]")
→ what endpoints matter clinically?
```
