# SynKit — CRN: Chemical Reaction Network Analysis

## Overview

The CRN module builds and analyzes **Chemical Reaction Networks (CRNs)** using graph theory and Feinberg's mathematical theory of reaction kinetics. It covers:
- CRN construction from reaction tables or SBML
- Structural properties (complexes, linkage classes, deficiency)
- Feinberg deficiency zero / one theorems
- Petri-net analysis (semiflows, siphons, traps)
- Autocatalysis detection
- Species-reaction graph properties

---

## Core Data Structures

```python
from synkit.CRN import CRNSpecies, CRNReaction, CRNNetwork

# Define species
glucose = CRNSpecies("glucose")
pyruvate = CRNSpecies("pyruvate")
atp = CRNSpecies("ATP")

# Define a reaction
glycolysis = CRNReaction(
    reactants={"glucose": 1, "ATP": 2},
    products={"pyruvate": 2, "ATP": 4},
    rate_constant=0.1,           # optional
    reversible=False
)

# Build network
crn = CRNNetwork()
crn.add_reaction(glycolysis)
```

---

## Build CRN from Reaction Table

```python
from synkit.CRN import crn_from_rxn_table
import pandas as pd

df = pd.DataFrame({
    "reactants": ["A + B", "C", "A + C"],
    "products":  ["C",     "D", "D + B"],
})

crn = crn_from_rxn_table(df, reactant_col="reactants", product_col="products")
print(crn.num_species)
print(crn.num_reactions)
```

## Build CRN from SBML

```python
from synkit.CRN import crn_from_sbml

crn = crn_from_sbml("pathway_model.xml")   # standard SBML format
```

---

## Structural Properties

```python
from synkit.CRN import CRNAnalyzer, CRNStructuralProperties

analyzer = CRNAnalyzer(crn)
props    = analyzer.compute_structural_properties()

# Core properties
print(props.num_species)          # |S|
print(props.num_complexes)        # |C|
print(props.num_reactions)        # |R|
print(props.num_linkage_classes)  # l
print(props.deficiency)           # δ = |C| - l - rank(stoichiometric matrix)
print(props.rank)                 # rank of stoichiometric matrix
print(props.reversible)           # all reactions reversible?
```

### Complexes and Complex Graph

```python
from synkit.CRN import compute_complexes, build_complex_graph

complexes = compute_complexes(crn)       # list of (species, stoichiometry) pairs
cx_graph  = build_complex_graph(crn)     # NetworkX DiGraph: complexes → reactions
```

---

## Feinberg Deficiency Theory

The deficiency δ = |C| - l - s (complexes - linkage classes - stoichiometric rank).

### Deficiency Zero Theorem

If δ = 0 and the network is weakly reversible → admits a unique positive steady state per stoichiometric compatibility class.

```python
from synkit.CRN import is_deficiency_zero_applicable

result = is_deficiency_zero_applicable(crn)
if result.applicable:
    print("Deficiency zero theorem applies — unique positive steady state exists")
    print(f"δ = {result.deficiency}, weakly reversible = {result.weakly_reversible}")
```

### Deficiency One Theorem

```python
from synkit.CRN import (
    is_deficiency_one_theorem_applicable,
    compute_linkage_class_deficiencies,
    run_deficiency_one_algorithm,
    DeficiencyOneAlgorithmResult
)

# Check applicability conditions
check = is_deficiency_one_theorem_applicable(crn)
print(check.applicable)
print(f"δ = {check.deficiency}")
print(f"Per-linkage-class deficiencies: {check.linkage_deficiencies}")
print(f"Regular network: {check.is_regular}")

# Full deficiency one analysis
result: DeficiencyOneAlgorithmResult = run_deficiency_one_algorithm(crn)
print(result.has_multiple_steady_states)   # can multiple positive SS exist?
print(result.bistability_possible)
```

### Linkage class deficiencies

```python
lc_deficiencies = compute_linkage_class_deficiencies(crn)
# List of per-linkage-class δ values; must all be ≤ 1 for deficiency one theorem
```

---

## Network Regularity

```python
from synkit.CRN import is_regular_network

print(is_regular_network(crn))   # required for deficiency one theorem
```

---

## Network Properties (Concordance, Endotacticity)

```python
from synkit.CRN import (
    is_concordant, is_accordant,
    is_endotactic, is_strongly_endotactic
)

print(is_concordant(crn))          # for injectivity / multistationarity analysis
print(is_accordant(crn))
print(is_endotactic(crn))          # for persistence analysis
print(is_strongly_endotactic(crn)) # stronger persistence guarantee
```

---

## Petri-Net Analysis

CRNs are bipartite Petri nets: species = places, reactions = transitions.

### Semiflows

```python
from synkit.CRN import compute_P_semiflows, compute_T_semiflows

# P-semiflows (place invariants): conservation laws
p_semiflows = compute_P_semiflows(crn)
for sf in p_semiflows:
    print(f"Conservation law: {sf}")   # e.g., [A] + [B] = const

# T-semiflows (transition invariants): cyclic pathways
t_semiflows = compute_T_semiflows(crn)
for sf in t_semiflows:
    print(f"Cyclic pathway: {sf}")
```

### Siphons and Traps

```python
from synkit.CRN import find_siphons, find_traps, check_persistence_sufficient

siphons = find_siphons(crn)    # sets of species that cannot be produced if depleted
traps   = find_traps(crn)      # sets of species that cannot be consumed if produced

# Sufficient condition for persistence (no species goes to zero)
persistent = check_persistence_sufficient(crn)
print(f"Sufficient condition for persistence: {persistent}")
```

---

## Autocatalysis Detection

```python
from synkit.CRN import (
    build_species_reaction_graph,
    is_autocatalytic,
    find_sr_graph_cycles,
    check_species_reaction_graph_conditions,
    is_SSD
)

sr_graph = build_species_reaction_graph(crn)
cycles   = find_sr_graph_cycles(sr_graph)

# Check if any species catalyzes its own production
autocatalytic = is_autocatalytic(crn)
print(f"Autocatalytic: {autocatalytic}")

# Check SSD (Single Species Domination)
ssd = is_SSD(crn)

# Full SR-graph condition check (for multistationarity analysis)
conditions = check_species_reaction_graph_conditions(sr_graph)
```

---

## Full Analysis Pipeline

```python
from synkit.CRN import CRNNetwork, CRNAnalyzer, crn_from_rxn_table

crn = crn_from_rxn_table(reaction_df)

analyzer = CRNAnalyzer(crn)
result   = analyzer.run_full_analysis()

# Summary report
print(f"Species:           {result.num_species}")
print(f"Complexes:         {result.num_complexes}")
print(f"Linkage classes:   {result.num_linkage_classes}")
print(f"Deficiency δ:      {result.deficiency}")
print(f"Weakly reversible: {result.weakly_reversible}")
print(f"Def-zero applies:  {result.deficiency_zero_applicable}")
print(f"Def-one applies:   {result.deficiency_one_applicable}")
print(f"Autocatalytic:     {result.autocatalytic}")
print(f"Persistent:        {result.persistence_sufficient}")
print(f"P-semiflows:       {len(result.p_semiflows)}")
```

---

## Key Concepts

### Deficiency δ

```
δ = |C| - l - rank(N)
```
- |C| = number of complexes (distinct LHS / RHS of reactions)
- l   = number of linkage classes (connected components of the complex graph)
- rank(N) = rank of the stoichiometric matrix

| δ | Implications |
|---|--------------|
| 0 | Simplest case — strong theorems apply (unique SS, no oscillation) |
| 1 | Deficiency one theorem may apply — conditional multistationarity possible |
| ≥ 2 | Less constrained — general analysis needed |

### Feinberg Theorems in a nutshell

| Theorem | Condition | Conclusion |
|---------|-----------|------------|
| Deficiency zero | δ = 0 + weakly reversible | Unique positive SS, no multistationarity, no oscillation |
| Deficiency one | δ = 1 + regular + δₗ ≤ 1 | Can determine if multistationarity possible |

---

## Best Practices

- Always call `compute_structural_properties()` first — it's fast and guides all subsequent analysis
- **Deficiency zero** is the most powerful theorem — check it before investing in deeper analysis
- Use **Petri-net semiflows** to find conservation laws automatically (useful for biochemical networks)
- **Siphons** are the primary persistence concern in biochemical networks — always check for minimal siphons
- For metabolic networks from SBML, use `crn_from_sbml` directly — it preserves stoichiometry and reversibility flags
