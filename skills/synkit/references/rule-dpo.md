# SynKit — Rule: DPO Graph Transformation & Reaction Application

## Overview

SynKit represents reactions as **DPO (Double Pushout) graph transformation rules** in GML format. A DPO rule encodes exactly what changes in a reaction: atoms/bonds in the LHS but not RHS are destroyed; atoms/bonds in the RHS but not LHS are created; atoms/bonds in both (the interface/context) are preserved.

---

## DPO Rule Structure

```
Rule = LHS ← Interface → RHS
         (left)  (context) (right)
```

In GML format:
```gml
rule [
  ruleID "ester_hydrolysis"
  left [
    edge [ source 1 target 2 label "=" ]   ← bond broken
  ]
  context [
    node [ id 1 label "C" ]
    node [ id 2 label "O" ]
    node [ id 3 label "O" ]
  ]
  right [
    edge [ source 1 target 2 label "-" ]   ← bond formed (single)
    edge [ source 3 target 4 label "-" ]
  ]
]
```

---

## Extracting Rules from Reactions

### From a single reaction SMILES

```python
from synkit.Rule import MoleculeRule, ReactorRule

rxn = "[CH3:1][C:2](=[O:3])[O:4][CH3:5].[OH2:6]>>[CH3:1][C:2](=[O:3])[OH:4].[CH3:5][OH:6]"

# Extract DPO rule as GML string
mol_rule = MoleculeRule(rxn)
gml = mol_rule.extract()
print(gml)
```

### From an ITS graph

```python
from synkit.Graph import ITSConstruction
from synkit.Rule import ReactorRule

its = ITSConstruction.from_reaction_smiles(rxn)
reactor_rule = ReactorRule.from_its(its)
gml = reactor_rule.to_gml()
```

---

## Applying Rules (Forward Prediction)

### SynReactor — NetworkX-based reactor

```python
from synkit.Synthesis import SynReactor
from rdkit import Chem

# Load a rule
reactor = SynReactor(gml_string)

# Apply to reactant molecule(s)
reactant_smiles = "CC(=O)OC"   # methyl acetate
products = reactor.apply(reactant_smiles)
# returns list of product SMILES strings

for prod in products:
    print(prod)
```

SynReactor features:
- Comprehensive atom-mapping output
- Native handling of implicit hydrogens
- Minimal dependencies (NetworkX only, no MØD required)
- DPO-like semantics (faithful graph rewriting)

### MODReactor — MØD-based reactor (optional)

```python
from synkit.Synthesis import MODReactor

# Requires MØD installation
reactor = MODReactor(gml_string)
products = reactor.apply(reactant_smiles)
```

---

## Rule Composition

Compose two reactions rules into a combined rule representing a two-step transformation.

```python
from synkit.Rule import RuleCompose

gml_step1 = "..."   # protection reaction
gml_step2 = "..."   # functionalization reaction

# Compose: step1 followed by step2
composed = RuleCompose(gml_step1, gml_step2).compose()
composed_gml = composed.to_gml()

# Apply composed rule
from synkit.Synthesis import SynReactor
reactor = SynReactor(composed_gml)
final_products = reactor.apply(starting_material)
```

---

## MODAAM — Atom Mapping via MØD

```python
from synkit.Synthesis import MODAAM

# Map atoms in a reaction using MØD's canonical mapper
mapped_rxn = MODAAM(unmapped_rxn).map()
print(mapped_rxn)   # atom-mapped SMILES
```

---

## Rule Management

### Load / save rules

```python
from synkit.Rule import MoleculeRule
import json

# Extract and serialize multiple rules
rules = {}
for rxn_id, rxn in reactions.items():
    try:
        gml = MoleculeRule(rxn).extract()
        rules[rxn_id] = gml
    except Exception:
        continue

with open("rules.json", "w") as f:
    json.dump(rules, f)

# Reload
with open("rules.json") as f:
    rules = json.load(f)
```

### Convert GML → NetworkX and back

```python
from synkit.IO import GMLToNX, NXToGML

lhs, ctx, rhs = GMLToNX(gml).convert()   # three NetworkX DiGraphs
gml_out = NXToGML(lhs, ctx, rhs).convert()
```

---

## Retrosynthesis via Rule Inversion

A DPO rule is invertible: swap LHS and RHS to get the retro-direction.

```python
from synkit.Rule import ReactorRule

rule = ReactorRule.from_gml(gml)
retro_rule = rule.invert()                  # swap LHS ↔ RHS
retro_gml  = retro_rule.to_gml()

# Apply in retro direction
from synkit.Synthesis import SynReactor
reactor = SynReactor(retro_gml)
precursors = reactor.apply(product_smiles)
```

---

## Rule Visualization

```python
from synkit.Vis import RuleVis

# Visualize LHS / context / RHS side by side
RuleVis(gml_string).draw(
    show_context=True,
    highlight_changes=True
)
```

---

## Key Concepts

### ITS subgraph = DPO rule

The reaction center subgraph of an ITS (edges with `change="formed"` or `change="broken"`, plus their incident nodes) is exactly the DPO rule's interface + LHS/RHS edges. SynKit exploits this bijection to move fluidly between the two representations.

### Why DPO (not SMARTS)?

| | DPO / GML | SMARTS |
|-|-----------|--------|
| Direction | Bidirectional (invert LHS↔RHS for retro) | Forward only |
| Composition | Algebraic rule composition | Manual |
| Graph theory foundation | Formal (Ehrig et al.) | Heuristic |
| Expressiveness | Full graph rewriting | Pattern matching |
| Human-readable | Less | More |

---

## Workflow: Build a Rule Library from USPTO-50k

```python
from synkit.Chem import CanonRSMI, AAMValidator
from synkit.Rule import MoleculeRule
from synkit.Graph import ITSConstruction, GraphCanonicaliser, BatchCluster
from collections import defaultdict
import pandas as pd

df = pd.read_csv("uspto_50k.csv")
rule_clusters = defaultdict(list)

for rxn in df["rxn_smiles"]:
    try:
        # Validate and canonicalize
        if not AAMValidator(rxn).is_valid():
            continue
        canon = CanonRSMI(rxn).standardize()

        # Build canonical ITS and extract rule
        its = ITSConstruction.from_reaction_smiles(canon)
        canon_its = GraphCanonicaliser(its).canonicalize(method="exact")
        gml = MoleculeRule(canon).extract()

        # Use WL hash as cluster key
        from synkit.Graph import CanonicalGraph
        key = CanonicalGraph(canon_its).wl_hash(n_iter=3)
        rule_clusters[key].append({"rxn": canon, "gml": gml})

    except Exception:
        continue

print(f"Extracted {len(rule_clusters)} unique reaction templates")
```

---

## Best Practices

- Extract rules only from **validated, balanced, atom-mapped** reactions
- Use `SynReactor` (no MØD) for portability; use `MODReactor` for performance in large-scale applications
- **Rule composition** is powerful but combinatorially explosive — limit to 2–3 step compositions
- Always **visualize** composed rules with `RuleVis` before deploying in a pipeline
- Store rules in GML format (human-readable, version-control friendly)
