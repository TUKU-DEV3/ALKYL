# SynKit — IO: Format Conversion & Data Loading

## Overview

The IO module is SynKit's adapter layer: it translates between reaction SMILES, ITS graphs (NetworkX), and DPO rules (GML). All conversions preserve atom-mapping and provenance.

---

## Conversion Pipeline

```
Reaction SMILES  ──────────► ITS Graph (NetworkX)  ──────────► DPO Rule (GML)
                 IO.convert        Graph.canonicalize       Rule.extract
                 ◄──────────────────────────────────────────────────────────
                          reverse conversion (lossless for balanced reactions)
```

**Performance on USPTO-50k:**
- SMILES → GML: 2.19 ± 2.15 ms / reaction
- GML → SMILES: 3.63 ± 2.69 ms / reaction

---

## SMILES → ITS Graph

```python
from synkit.IO import load_reaction_smiles  # high-level adapter
from synkit.Graph import ITSConstruction    # low-level control

# High-level: string → ITS in one call
rxn = "[CH3:1][OH:2].[Br:3][CH2:4][CH3:5]>>[CH3:1][O:2][CH2:4][CH3:5].[H:6][Br:3]"
its_graph = load_reaction_smiles(rxn)   # returns NetworkX DiGraph

# Low-level: fine-grained control over ITS construction
its_graph = ITSConstruction.from_reaction_smiles(
    rxn,
    sanitize=True,        # apply RDKit sanitization
    keep_mapping=True,    # preserve atom-map numbers as node attributes
    explicit_h=False      # include explicit H only at reaction centers
)

# Inspect result
print(its_graph.number_of_nodes())
print(its_graph.number_of_edges())
for u, v, data in its_graph.edges(data=True):
    print(u, v, data)  # data["bond_type"], data["change"] = "formed"|"broken"|None
```

---

## ITS Graph → SMILES

```python
from synkit.IO import its_to_smiles

rxn_smiles = its_to_smiles(its_graph)   # reconstructs atom-mapped reaction SMILES

# Or via the Graph module's reverse interface
from synkit.Graph import ITSConstruction
rxn_smiles = ITSConstruction.to_reaction_smiles(its_graph)
```

**Caveats for GML → SMILES reversion:**
- Requires explicit H at reaction centers in the GML rule
- Stereochemistry is not preserved (omitted in ITS representation)
- Only valid for balanced, fully-mapped reactions

---

## SMILES → GML (DPO Rule)

```python
from synkit.IO import smiles_to_gml

rxn = "[CH3:1][OH:2].[Br:3][CH2:4]>>[CH3:1][O:2][CH2:4].[H:6][Br:3]"
gml_string = smiles_to_gml(rxn)
print(gml_string)
# rule [
#   ruleID "..."
#   left  [ ... ]       ← LHS: reactant bonds/atoms
#   context [ ... ]     ← Interface: unchanged atoms
#   right [ ... ]       ← RHS: product bonds/atoms
# ]
```

---

## GML → SMILES

```python
from synkit.IO import gml_to_smiles

rxn_smiles = gml_to_smiles(gml_string)  # reconstructs reaction SMILES from DPO rule
```

---

## GML ↔ NetworkX

```python
from synkit.IO import GMLToNX, NXToGML

# GML string → NetworkX graph (three-part DPO rule)
lhs, context, rhs = GMLToNX(gml_string).convert()

# NetworkX graphs → GML string
gml_out = NXToGML(lhs=lhs, context=context, rhs=rhs).convert()
```

---

## Molecule ↔ Graph

```python
from synkit.IO import MolToGraph, GraphToMol
from rdkit import Chem

# RDKit Mol → NetworkX graph (single molecule)
mol = Chem.MolFromSmiles("CCO")
g = MolToGraph(mol).convert()       # nodes=atoms, edges=bonds

# NetworkX graph → RDKit Mol
mol2 = GraphToMol(g).convert()
print(Chem.MolToSmiles(mol2))       # "CCO"
```

---

## Batch / File I/O

```python
from synkit.IO import load_reactions_from_file, save_reactions_to_file

# Load from CSV / JSON with reaction SMILES column
reactions = load_reactions_from_file(
    "reactions.csv",
    smiles_col="rxn_smiles",
    sep=","
)

# Save ITS graphs or processed reactions
save_reactions_to_file(processed, "output.json")
```

---

## Visualization

```python
from synkit.Vis import RXNVis, RuleVis

# Visualize reaction from SMILES
RXNVis(rxn_smiles).draw()

# Visualize DPO rule (LHS / context / RHS side by side)
RuleVis(gml_string).draw()
```

---

## Debug Helper

```python
from synkit.IO import print_graph_attributes

# Pretty-print all node and edge attributes of an ITS graph
print_graph_attributes(its_graph)
```

---

## ITS Graph Node & Edge Attributes

### Node attributes (atoms)

| Attribute | Type | Description |
|-----------|------|-------------|
| `symbol` | str | Atomic symbol (`"C"`, `"O"`, `"N"`, …) |
| `charge` | int | Formal charge |
| `radical` | int | Number of radical electrons |
| `hybridization` | str | `"SP3"`, `"SP2"`, `"SP"`, `"AROMATIC"` |
| `hcount` | int | Explicit H count at reaction center |
| `atom_map` | int | Atom mapping number from SMILES |
| `in_rc` | bool | True if this atom is in the reaction center |

### Edge attributes (bonds)

| Attribute | Type | Values |
|-----------|------|--------|
| `bond_type` | str | `"SINGLE"`, `"DOUBLE"`, `"TRIPLE"`, `"AROMATIC"` |
| `change` | str / None | `"formed"` / `"broken"` / `None` (unchanged) |
| `order` | float | 1.0, 1.5, 2.0, 3.0 |

---

## Common Workflow: Reaction Database Processing

```python
from synkit.IO import load_reactions_from_file, smiles_to_gml
from synkit.Graph import ITSConstruction, GraphCanonicaliser
from synkit.Chem import AAMValidator
import pandas as pd

df = pd.read_csv("uspto_50k.csv")  # column "rxn_smiles"

results = []
for rxn in df["rxn_smiles"]:
    # Validate AAM
    if not AAMValidator(rxn).is_valid():
        continue

    # Build canonical ITS
    its = ITSConstruction.from_reaction_smiles(rxn)
    canon = GraphCanonicaliser(its).canonicalize()

    # Extract DPO rule
    gml = smiles_to_gml(rxn)

    results.append({"rxn": rxn, "its": its, "gml": gml})
```

---

## Best Practices

- Always **validate AAM first** (`AAMValidator`) before constructing ITS
- Use **high-level `load_reaction_smiles`** for pipelines; use `ITSConstruction` directly for research/debugging
- Keep `explicit_h=False` unless you need GML → SMILES reversion (which requires H at reaction centers)
- For large batches, process reactions in parallel — each conversion is stateless
- **Stereochemistry** is lost in ITS representation; keep original SMILES if stereochemistry matters downstream
