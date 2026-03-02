# SynKit — Chem: Standardization, Canonicalization & AAM Validation

## Overview

The Chem module extends RDKit's molecule protocols to full **transformational settings** (reactions). It handles reaction SMILES canonicalization, tautomer normalization, atom-mapping validation, and balance checking.

---

## CanonRSMI — Canonical Reaction SMILES

Produces a deterministic, canonical form of a reaction SMILES independent of atom ordering or RDKit's non-deterministic output.

```python
from synkit.Chem import CanonRSMI

rxn = "[CH3:1][OH:2].[Na:3][H:4]>>[CH3:1][O:2][Na:3].[H:4][H:5]"
canon = CanonRSMI(rxn)

# Canonicalize
canonical_smiles = canon.canonicalize()
print(canonical_smiles)   # deterministic regardless of input atom order

# Normalize tautomers (reduces false negatives in duplicate detection)
normalized = canon.normalize()

# Full pipeline: normalize then canonicalize
result = canon.standardize()
```

### When to use

| Use case | Method |
|----------|--------|
| Deduplication of reaction databases | `standardize()` |
| Consistent reaction SMILES for hashing | `canonicalize()` |
| Pre-processing before AAM comparison | `normalize()` |

---

## AAMValidator — Atom-to-Atom Mapping Validation

Tests whether a predicted atom-to-atom mapping is chemically correct by comparing the resulting ITS graph against the curated ground-truth ITS via graph isomorphism.

```python
from synkit.Chem import AAMValidator

rxn = "[CH3:1][OH:2].[Br:3][CH2:4]>>[CH3:1][O:2][CH2:4].[H:6][Br:3]"
validator = AAMValidator(rxn)

# Check validity (graph isomorphism test)
is_valid = validator.is_valid()          # True / False
print(is_valid)

# Detailed comparison with reference mapping
score = validator.compare(reference_rxn)  # returns similarity score
```

### How it works

1. Constructs ITS graph from the predicted mapping
2. Constructs ITS graph from the reference (curated) mapping
3. Tests graph isomorphism between the two ITS graphs
4. Returns True iff the predicted mapping reproduces the reference ITS

### Common use: evaluate RXNMapper / LocalMapper output

```python
from synkit.Chem import AAMValidator
from rxnmapper import RXNMapper

unmapped_rxn = "CCO.BrCCl>>CCOCCl.Br"
mapper = RXNMapper()
mapped = mapper.get_attention_guided_atom_maps([unmapped_rxn])[0]["mapped_rxn"]

reference = "[CH3:1][CH2:2][OH:3].[Br:4][CH2:5][Cl:6]>>[CH3:1][CH2:2][O:3][CH2:5][Cl:6].[H:7][Br:4]"
is_correct = AAMValidator(mapped).compare(reference)
print(f"Mapping correct: {is_correct}")
```

---

## Reaction — Reaction Representation & Analysis

```python
from synkit.Chem import Reaction

rxn = Reaction("[CH3:1][OH:2].[Cl:3][CH2:4]>>[CH3:1][O:2][CH2:4].[H:5][Cl:3]")

# Basic properties
print(rxn.reactants)        # list of SMILES
print(rxn.products)         # list of SMILES
print(rxn.atom_map)         # atom mapping dict

# Standardize the reaction
std = rxn.standardize()

# Check atom balance
balanced = rxn.is_balanced()

# Get reaction center atoms
rc_atoms = rxn.get_reaction_center()

# Fingerprint for clustering
fp = rxn.fingerprint(method="ITS")    # ITS-based fingerprint
```

---

## BalanceReactionCheck

Verify that a reaction conserves atoms and charges:

```python
from synkit.Chem import BalanceReactionCheck

check = BalanceReactionCheck("[CH3:1][OH:2]>>[CH3:1][H:2]")
print(check.is_balanced())           # False — oxygen atom lost
print(check.get_balance_report())    # detailed atom count diff
```

---

## Standardization Pipeline

Recommended pipeline before building a reaction database:

```python
from synkit.Chem import CanonRSMI, AAMValidator, BalanceReactionCheck
from synkit.Graph import ITSConstruction

def standardize_reaction(rxn_smiles: str) -> dict | None:
    """Validate, balance-check, and canonicalize a reaction SMILES."""

    # 1. Balance check
    if not BalanceReactionCheck(rxn_smiles).is_balanced():
        return None

    # 2. Canonicalize SMILES
    canon_rxn = CanonRSMI(rxn_smiles).standardize()

    # 3. Validate AAM (requires reference or self-consistency check)
    validator = AAMValidator(canon_rxn)
    if not validator.is_valid():
        return None

    # 4. Build ITS
    its = ITSConstruction.from_reaction_smiles(canon_rxn)

    return {"canon_smiles": canon_rxn, "its": its}
```

---

## FPCalculator — Reaction Fingerprints

```python
from synkit.Chem import FPCalculator

calc = FPCalculator(rxn_smiles)

# Difference fingerprint (products FP minus reactants FP)
diff_fp = calc.difference_fingerprint(radius=2, n_bits=2048)

# ITS-based fingerprint (from reaction center graph)
its_fp = calc.its_fingerprint()

# Use for similarity search
from rdkit import DataStructs
similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
```

---

## Standardize — Reaction Standardization

```python
from synkit.Chem import Standardize

std = Standardize(rxn_smiles)

# Full standardization pipeline
result = std.run()
# 1. Fragment removal (salts, spectator ions)
# 2. Charge normalization
# 3. Tautomer normalization
# 4. Canonical SMILES
# 5. Agent/reagent separation (>> syntax)
```

---

## Key Concepts

### Why tautomer normalization matters

Two atom-mapped reactions can have identical atom mappings but different tautomeric forms. Without normalization, graph isomorphism comparisons return false negatives. `CanonRSMI.normalize()` applies RDKit's tautomer enumerator to pick a canonical tautomeric form before comparison.

### AAM validation vs. canonicalization

| Operation | Class | Purpose |
|-----------|-------|---------|
| Make SMILES deterministic | `CanonRSMI` | Deduplication, hashing |
| Check if mapping is correct | `AAMValidator` | Quality control of AAM tools |
| Verify atom conservation | `BalanceReactionCheck` | Data cleaning |
| Full standardization | `Standardize` | Database building |

---

## Best Practices

- Run `BalanceReactionCheck` **before** `ITSConstruction` — unbalanced reactions produce malformed ITS graphs
- Use `CanonRSMI.standardize()` (not just `.canonicalize()`) for database deduplication — tautomers are frequent in real data
- `AAMValidator` is only as good as the reference mapping — use curated USPTO/Reaxys references
- For large datasets, vectorize: apply `CanonRSMI` and `BalanceReactionCheck` as row functions before building ITS graphs
