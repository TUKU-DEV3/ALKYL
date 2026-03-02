# SynKit — Synthesis: Planning & Route Exploration

## Overview

The Synthesis module provides the reactive engine for applying DPO rules to molecules (forward prediction) and exploring synthetic routes. SynKit is **not a turnkey retrosynthesis planner** — it provides modular primitives for building one. Full planning requires external rule databases, search strategies, and feasibility scoring.

---

## SynReactor — Apply a Rule to Reactants

```python
from synkit.Synthesis import SynReactor

# Load a rule (GML string)
reactor = SynReactor(gml_string)

# Apply to a reactant SMILES
products = reactor.apply("CC(=O)OC")          # one SMILES string
# or with explicit atom mapping preserved:
products = reactor.apply("CC(=O)OC", keep_aam=True)

for prod in products:
    print(prod)   # SMILES of each product set
```

SynReactor uses NetworkX DPO-like graph rewriting:
- Handles implicit hydrogens natively
- Returns all valid application sites (not just the first)
- No MØD dependency

---

## MODReactor — MØD-Based Reactor (optional)

```python
from synkit.Synthesis import MODReactor

# Requires: pip install synkit[all] + MØD installation
reactor = MODReactor(gml_string)
products = reactor.apply("CC(=O)OC")
```

Differences from SynReactor:

| Feature | SynReactor | MODReactor |
|---------|-----------|------------|
| Dependency | NetworkX only | MØD required |
| Speed | Moderate | Faster at scale |
| Atom mapping output | Full | Full |
| Installation | Easy | Complex |

---

## MODAAM — Atom Mapping

```python
from synkit.Synthesis import MODAAM

# Map atoms in an unmapped reaction using MØD
unmapped = "CC(=O)OC.O>>CC(=O)O.CO"
mapped = MODAAM(unmapped).map()
print(mapped)   # atom-mapped SMILES
```

---

## Forward Prediction Pipeline

```python
from synkit.Synthesis import SynReactor
from synkit.Rule import MoleculeRule
from synkit.Chem import CanonRSMI

# 1. Build rule library from training reactions
rules = {}
for rxn in training_reactions:
    try:
        gml = MoleculeRule(rxn).extract()
        key = CanonRSMI(rxn).canonicalize()
        rules[key] = gml
    except Exception:
        continue

# 2. Apply all compatible rules to a query molecule
query = "CCOc1ccc(N)cc1"    # aromatic amine

all_products = []
for rule_key, gml in rules.items():
    reactor = SynReactor(gml)
    try:
        products = reactor.apply(query)
        all_products.extend(products)
    except Exception:
        continue

# 3. Deduplicate and filter
from rdkit import Chem
unique = {Chem.MolToSmiles(Chem.MolFromSmiles(p)) for p in all_products if Chem.MolFromSmiles(p)}
print(f"Generated {len(unique)} unique products")
```

---

## Single-Step Retrosynthesis

```python
from synkit.Rule import ReactorRule
from synkit.Synthesis import SynReactor

# Invert rule (LHS ↔ RHS swap)
rule = ReactorRule.from_gml(gml_string)
retro_rule = rule.invert()
retro_reactor = SynReactor(retro_rule.to_gml())

# Apply in retro direction
target_smiles = "CC(=O)Nc1ccc(O)cc1"   # paracetamol
precursors = retro_reactor.apply(target_smiles)
for prec in precursors:
    print("Precursor set:", prec)
```

---

## Multi-Step Synthesis Planning (BFS / DFS / MCTS)

SynKit does not include a built-in tree search, but provides the primitives. Example BFS:

```python
from synkit.Synthesis import SynReactor
from synkit.Rule import ReactorRule
from collections import deque

def is_commercial(smiles: str, building_blocks: set) -> bool:
    from rdkit import Chem
    canon = Chem.MolToSmiles(Chem.MolFromSmiles(smiles)) if smiles else None
    return canon in building_blocks

def bfs_retrosynthesis(target: str, retro_rules: list, building_blocks: set, max_depth: int = 5):
    """BFS synthesis tree search using SynKit retro reactors."""
    queue = deque([(target, [], 0)])
    visited = set()

    while queue:
        molecule, path, depth = queue.popleft()

        if is_commercial(molecule, building_blocks):
            return path + [molecule]   # route found!

        if depth >= max_depth or molecule in visited:
            continue
        visited.add(molecule)

        for reactor in retro_rules:
            try:
                precursor_sets = reactor.apply(molecule)
                for precursor_set in precursor_sets:
                    for prec in precursor_set.split("."):
                        queue.append((prec, path + [molecule], depth + 1))
            except Exception:
                continue

    return None   # no route found within max_depth

# Set up retro reactors from rule library
retro_reactors = [
    SynReactor(ReactorRule.from_gml(gml).invert().to_gml())
    for gml in rule_library.values()
]

route = bfs_retrosynthesis(
    "CC(=O)Nc1ccc(O)cc1",    # paracetamol
    retro_reactors,
    building_blocks,
    max_depth=4
)
```

---

## Rule Composition for Multi-Step Sequences

```python
from synkit.Rule import RuleCompose

# Combine two single-step rules into a two-step rule
rule_a = "..."   # step 1 GML
rule_b = "..."   # step 2 GML

composed = RuleCompose(rule_a, rule_b).compose()
reactor  = SynReactor(composed.to_gml())

# Apply the two-step transformation in a single pass
final_products = reactor.apply(starting_material)
```

---

## ReactionNetwork (Synthesis Module)

Build a local reaction network by exhaustive rule application:

```python
from synkit.Synthesis import ReactionNetwork

# Start from a set of molecules and apply rules iteratively
network = ReactionNetwork(
    initial_species=["CCO", "Cl", "CC(=O)Cl"],
    rules=list(rule_library.values()),
    max_depth=3
)

network.build()    # applies all compatible rules at each depth

# Inspect the generated network
print(f"Generated {network.num_species} species, {network.num_reactions} reactions")

# Export to SynKit CRN for analysis
from synkit.CRN import CRNNetwork
crn = network.to_crn()
```

---

## Integration with Other Tools

### With TorchDrug (ML retrosynthesis)

```python
# SynKit: rule-based retrosynthesis primitives
# TorchDrug: ML-based (GNN) retrosynthesis

# Use SynKit to validate TorchDrug predictions:
from synkit.Chem import AAMValidator, BalanceReactionCheck

for predicted_rxn in torchdrug_predictions:
    if (BalanceReactionCheck(predicted_rxn).is_balanced()
            and AAMValidator(predicted_rxn).is_valid()):
        # Confirmed prediction — extract rule for rule library
        from synkit.Rule import MoleculeRule
        gml = MoleculeRule(predicted_rxn).extract()
```

### With RDKit (SMARTS-based filtering)

```python
from rdkit import Chem
from rdkit.Chem import AllChem
from synkit.Synthesis import SynReactor

# Pre-filter molecules with SMARTS before applying DPO rules
smarts = Chem.MolFromSmarts("[NH2]c1ccccc1")   # aromatic amine

candidates = [mol for mol in molecule_library
              if mol.HasSubstructMatch(smarts)]

for mol in candidates:
    products = reactor.apply(Chem.MolToSmiles(mol))
```

---

## Workflow: Build a Mini Synthesis Tree

```python
from synkit.Rule import MoleculeRule, ReactorRule
from synkit.Synthesis import SynReactor
from synkit.Chem import AAMValidator
from rdkit import Chem

# 1. Build validated rule library
rules = {}
for rxn in curated_reactions:
    if AAMValidator(rxn).is_valid():
        try:
            gml = MoleculeRule(rxn).extract()
            rules[rxn] = gml
        except Exception:
            continue

# 2. Build retro reactors
retro_reactors = [
    SynReactor(ReactorRule.from_gml(gml).invert().to_gml())
    for gml in rules.values()
]

# 3. Single-step retrosynthesis of target
target = "O=C(Nc1ccc(F)cc1)c1ccc(Cl)cc1"
precursor_sets = []
for reactor in retro_reactors:
    try:
        sets = reactor.apply(target)
        precursor_sets.extend(sets)
    except Exception:
        continue

print(f"Found {len(precursor_sets)} 1-step disconnections")
for ps in precursor_sets[:5]:
    print(" +".join(ps.split(".")))
```

---

## Best Practices

- **Rule quality matters most**: derive rules only from validated, balanced, atom-mapped reactions
- **SynReactor** for portability; **MODReactor** for large-scale production (requires MØD)
- **Filter before applying**: use SMARTS or fingerprint similarity to pre-select compatible rules (avoids O(n×m) rule applications)
- **Rule composition**: limit to 2–3 step compositions — exponential blowup otherwise
- SynKit synthesis is **rule-based** (mechanistically grounded); TorchDrug is **ML-based** (data-driven) — use both in combination for best coverage
- Always **validate generated products** with RDKit before presenting to chemists
