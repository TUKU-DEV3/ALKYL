# SynKit — Graph: ITS Construction, Canonicalization & Matching

## Overview

The Graph module provides graph-based representations of chemical reactions. The two core abstractions are:
- **ITS (Imaginary Transition State)**: merges reactants + products into one multigraph
- **MTG (Mechanistic Transition Graph)**: extends ITS to model bond-making/breaking sequences and transient intermediates

---

## ITS — Imaginary Transition State Graph

```
Reactant graph  +  Product graph  =  ITS graph
                                      (edges carry "formed" / "broken" / None)
```

The ITS is equivalent to the CGR (Condensed Graph of Reaction) used in other frameworks. It is the universal substrate for all downstream SynKit operations.

### Construct from reaction SMILES

```python
from synkit.Graph import ITSConstruction

rxn = "[CH3:1][OH:2].[Br:3][CH2:4][CH3:5]>>[CH3:1][O:2][CH2:4][CH3:5].[H:6][Br:3]"

# Standard construction
its = ITSConstruction.from_reaction_smiles(rxn)

# With explicit hydrogens at reaction centers (required for GML → SMILES reversion)
its_h = ITSConstruction.from_reaction_smiles(rxn, explicit_h=True)
```

### Inspect the ITS graph

```python
import networkx as nx

# Nodes (atoms)
for node, data in its.nodes(data=True):
    print(node, data["symbol"], data["charge"], data.get("in_rc", False))

# Edges (bonds)
for u, v, data in its.edges(data=True):
    bt   = data["bond_type"]                  # "SINGLE", "DOUBLE", etc.
    chng = data.get("change")                 # "formed", "broken", or None
    print(f"{u}-{v}: {bt} ({chng})")

# Reaction center: edges with change != None
rc_edges = [(u, v) for u, v, d in its.edges(data=True) if d.get("change")]
```

### Construct from separate reactant / product SMILES

```python
reactants = "CCO.BrCCl"
products  = "CCOCCl.Br"
its = ITSConstruction.from_smiles(reactants, products, atom_map=True)
```

---

## MTG — Mechanistic Transition Graph

Extends ITS by explicitly modeling the **sequence** of bond events and transient intermediates (novel in SynKit).

```python
from synkit.Graph import ITSExpand

# Expand ITS → MTG (sequence of elementary steps)
its = ITSConstruction.from_reaction_smiles(rxn)
mtg = ITSExpand(its).expand()

# MTG nodes include transient intermediates
# MTG edges carry step ordering (bond formation before/after cleavage)

for node, data in mtg.nodes(data=True):
    if data.get("intermediate"):
        print("Transient:", data["symbol"])
```

MTG enables:
- Full reaction pathway analysis
- Distinguishing concerted vs. stepwise mechanisms
- Modeling formal charge and radical state transitions

---

## GraphCanonicaliser — Canonical ITS

Produces a canonical form of the ITS graph independent of atom indexing/ordering.

### Exact canonicalization (recommended)

Uses Nauty/Bliss individualization-refinement algorithms.

```python
from synkit.Graph import GraphCanonicaliser

canon = GraphCanonicaliser(its)
canonical_its = canon.canonicalize(method="exact")

# Validate: two equivalent reactions should give identical canonical ITS
its_a = ITSConstruction.from_reaction_smiles(rxn_a)
its_b = ITSConstruction.from_reaction_smiles(rxn_b)
same = GraphCanonicaliser.are_isomorphic(its_a, its_b)
print(same)   # True if same reaction type
```

| Method | Accuracy | Speed | Notes |
|--------|----------|-------|-------|
| `"exact"` | 100% | Median 0.77 ms (worst case: seconds) | Uses Nauty/Bliss |
| `"wl"` / `"wlgh3"` | ~95% | Fast | 5% non-canonical edge cases |

### Approximate canonicalization (WL hash)

```python
from synkit.Graph import CanonicalGraph

# WL hash fingerprint for rapid clustering (bucket prefilter)
wl_hash = CanonicalGraph(its).wl_hash(n_iter=3)   # WLGH3
```

---

## Graph Isomorphism Testing

```python
from synkit.Graph import GraphMatcherEngine

matcher = GraphMatcherEngine()

# Full graph isomorphism
is_same = matcher.is_isomorphic(its_a, its_b)

# Node/edge attribute matching (atom type + bond type)
is_same_typed = matcher.is_isomorphic(
    its_a, its_b,
    node_match=lambda a, b: a["symbol"] == b["symbol"],
    edge_match=lambda a, b: a["bond_type"] == b["bond_type"]
)
```

---

## Subgraph Matching

Find all occurrences of a query subgraph inside target graphs.

```python
from synkit.Graph import SubgraphMatch, SubgraphSearchEngine

# Single match
query  = ITSConstruction.from_reaction_smiles(query_rxn)
target = ITSConstruction.from_reaction_smiles(target_rxn)

match = SubgraphMatch(query, target).find()
# Returns list of node mapping dicts: {query_node: target_node, ...}

# Bulk search across a database of ITS graphs
engine = SubgraphSearchEngine(query)
hits = engine.search(its_database)          # list of matching ITS graphs
```

### SING / TurboISO backends

SynKit exposes two subgraph isomorphism backends:

```python
from synkit.Graph import SING, TurboISO

# SING: general-purpose, exact
results_sing    = SING(query, target).run()

# TurboISO: optimized for large graphs
results_turbo   = TurboISO(query, target).run()
```

---

## Reaction Clustering

Cluster a reaction database by structural similarity using WL graph hashing.

```python
from synkit.Graph import GraphCluster, BatchCluster

reactions = [its_1, its_2, its_3, ...]   # list of ITS graphs

# Single-step clustering via WL hash buckets
clusters = GraphCluster(reactions).cluster(method="wl", n_iter=3)
# Returns dict: {hash_key: [ITS graphs in bucket]}

# Batch clustering with refinement (WL prefilter → exact isomorphism)
batch = BatchCluster(reactions)
fine_clusters = batch.run(prefilter="wl", refine="exact")
```

### Workflow for reaction database deduplication

```python
from synkit.Graph import BatchCluster
from synkit.Graph import ITSConstruction

# Build ITS graphs from database
its_list = [ITSConstruction.from_reaction_smiles(r) for r in rxn_smiles_list]

# Cluster
clusters = BatchCluster(its_list).run(prefilter="wl", refine="exact")

# Deduplicate: keep one representative per cluster
unique_its = [members[0] for members in clusters.values()]
print(f"Reduced {len(its_list)} reactions to {len(unique_its)} unique templates")
```

---

## Graph Utilities

```python
from synkit.Graph import (
    remove_wildcard_nodes,
    add_wildcard_subgraph_for_unmapped,
    clean_graph_keep_largest_component,
    has_wildcard_node,
    print_graph_attributes
)

# Remove wildcard (*) atoms from graph
clean_its = remove_wildcard_nodes(its)

# Add wildcard subgraph for unmapped atoms (for partial templates)
partial = add_wildcard_subgraph_for_unmapped(its)

# Check if wildcards are present
if has_wildcard_node(its):
    its = remove_wildcard_nodes(its)

# Keep only the largest connected component
largest = clean_graph_keep_largest_component(its)

# Debug: print all attributes
print_graph_attributes(its)
```

---

## Key Concepts: ITS vs. MTG

| Feature | ITS | MTG |
|---------|-----|-----|
| Bond events | Unordered set | Ordered sequence |
| Intermediates | Not represented | Explicit nodes |
| Concerted vs. stepwise | Indistinguishable | Distinguishable |
| Use for templates | Yes (standard) | Research / deep mechanistic analysis |
| Computational cost | Low | Higher |

---

## Best Practices

- Use **exact canonicalization** for deduplication in production databases
- Use **WL hash as prefilter** then exact isomorphism for large-scale clustering (two-stage)
- **Wildcards**: add them for unmapped atoms when building partial reaction templates; remove them before isomorphism tests
- When comparing ITS graphs, always match on both **node attributes** (symbol, charge) and **edge attributes** (bond type, change flag)
- MTG is computationally heavier — use only when mechanistic sequence analysis is needed
