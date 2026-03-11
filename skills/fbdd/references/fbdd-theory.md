# FBDD Theory

## Why Fragments Work: The Hann Complexity Model

Hann et al. (2001) showed that as molecular complexity increases, the probability
of a complementary binding event falls exponentially. Fragments are "simpler"
than drug-like molecules → higher hit rate in screens (5–30% vs 0.01–0.1% for HTS).

```
P(binding) ∝ exp(-complexity)
complexity ∝ MW, n_rotbonds, n_stereocenters, n_pharmacophoric_features
```

Fragments bind weakly (Kd 10 µM – 10 mM) but:
1. Occupy hot spots efficiently (high LE)
2. Leave room for elaboration
3. Require biophysical assays (SPR, NMR, X-ray) not cell-based HTS

## Rule of 3 (Ro3) — Astex Pharmaceuticals (2003)

```
MW ≤ 300 Da
cLogP ≤ 3
HBD ≤ 3
HBA ≤ 3
```
Extensions often used:
```
HAC ≤ 20 (heavy atom count)
Rotatable bonds ≤ 3
PSA ≤ 60 Å²
```

## Ligand Efficiency (LE)

```
LE = ΔG_binding / HAC = (RT ln Kd) / HAC
   ≈ 1.37 × pIC50 / HAC    [kcal/mol per heavy atom, at 298K]
```

| LE value | Interpretation |
|----------|---------------|
| ≥ 0.4 | Excellent fragment |
| 0.3–0.4 | Good (acceptable for elaboration) |
| 0.25–0.3 | Marginal — check binding mode |
| < 0.25 | Poor — likely non-specific |

**The LE Paradox**: LE tends to decrease as MW increases.
You cannot maintain LE = 0.4 above MW ~400 (Snyder et al., 2009).
Target LE ≥ 0.3 for leads; ≥ 0.3 for drug candidates.

## Lipophilic Ligand Efficiency (LLE / LipE)

```
LLE = pIC50 - cLogP
```

Corrects for potency gained "for free" via lipophilicity.
- LLE ≥ 3: acceptable
- LLE ≥ 5: excellent
- LLE < 2: lipophilicity-driven binding → ADMET liabilities

## Fragment-to-Lead Strategies

### 1. Fragment Growing
Extend fragment into adjacent pockets by adding substituents.
Requires: co-crystal structure or reliable docking pose.

```
Fragment → Fragment + R-group → Lead
HAC: 10 → 18–25
LE: 0.40 → 0.32–0.36 (acceptable decline)
```

Key: grow into unexplored pockets, not just add lipophilicity.

### 2. Fragment Linking
Merge two fragments that bind to adjacent sub-pockets using a linker.
Theoretical enthalpy additivity: ΔG_combined ≈ ΔG_A + ΔG_B + ΔG_linker
In practice: rarely fully additive (entropy cost of linker).

```
Fragment A + Fragment B → Fragment A–linker–Fragment B
Linker design: short (2–4 atoms), rigid, no new HBA/HBD unless needed
```

### 3. Fragment Merging (most common)
Combine overlapping fragments sharing a common pharmacophoric core.
Find MCS, keep pharmacophoric features of both.

```
Fragment A: [core]-[feature_A]
Fragment B: [core]-[feature_B]
Merged:     [core]-[feature_A]-[feature_B]
```

### 4. Fragment Hopping
Use fragment as query for 3D similarity search → find scaffold with same pharmacophore
but better ADMET.

## FBDD vs HTS: Key Differences

| | FBDD | HTS |
|-|------|-----|
| Library size | 500–2000 cpds | 500k–5M cpds |
| Hit rate | 1–30% | 0.01–0.1% |
| Hit potency | mM–µM | µM–nM |
| LE of hits | 0.3–0.5 | 0.15–0.3 |
| Detection | NMR, SPR, X-ray, DSF | Biochemical/cell assay |
| Starting point quality | High LE, low MW | Variable, often lipophilic |

## Key Success Cases

| Drug | Target | Fragment LE | Final IC50 |
|------|--------|-------------|-----------|
| Venetoclax (ABT-199) | BCL-2 | 0.37 | 0.01 nM |
| Vemurafenib | BRAF V600E | ~0.35 | 31 nM |
| AT13387 | HSP90 | 0.4 | 0.7 nM |
| Navitoclax (ABT-263) | BCL-2/XL | 0.37 | 1 nM |

## Binding Hotspots and FTMap

Hotspot residues contribute >70% of binding energy (DeLano, 2002).
Computational mapping identifies them:

```bash
# FTMap server (online): submit receptor PDB → hotspot residues
# Local: fpocket for cavity detection
fpocket -f receptor.pdb
# → pockets ranked by druggability score

# Python: use MDAnalysis + volmap for pocket volume
```

Fragment X-ray screening detects hotspots empirically:
- XChem (Diamond Light Source) — 1000 fragments, semi-automated
- Pan-Dataset Density Analysis (PanDDA) — weak density event detection

## Solubility and Aggregation

Fragments must be soluble at screening concentration (1–10 mM in DMSO assays,
0.1–1 mM aqueous for NMR/SPR).

Common aggregation filter:
```python
# Colloidal aggregators: aromatic compounds, cLogP > 3, MW > 200
# Rishton PAINS (broader) already in chem_filter.py
# Additional check: colloidal aggregator fingerprint
```

Rule of thumb: if fragment shows non-linear dose-response → suspect aggregation.
Confirm with detergent (0.01% Triton X-100) counter-screen.
