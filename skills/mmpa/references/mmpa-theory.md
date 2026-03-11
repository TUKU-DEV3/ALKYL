# MMPA Theory — MMP Definition, Fragmentation, SMIRKS Transforms, SAR Framework

## The Matched Molecular Pair (MMP)

A **matched molecular pair** (A, B) is defined as two molecules that differ by a single structural transformation at exactly one site. Formally:

```
A = Context + Fragment₁
B = Context + Fragment₂

Where:
- Context = identical molecular environment (the "scaffold")
- Fragment₁, Fragment₂ = the "variable parts" at the site of change
- Transform = Fragment₁ → Fragment₂ (expressed as SMIRKS)
```

**Property change (transform effect)**:
```
ΔProperty = Property(B) - Property(A)
```

**Example** — H→F transformation:
```
A: c1ccccc1CH₂OH   (benzyl alcohol,  LogP = 1.10)
B: c1ccccc1CF      (benzyl fluoride, LogP = 2.14)

Transform: [CH2:1]O >> [CH2:1]F    ΔLogP = +1.04
Context:   c1ccccc1[*:1]
```

## Single-Cut vs. Double-Cut Fragmentation

### Single-cut (standard MMPA)

One bond is broken between context and variable part. The variable part is typically ≤10 heavy atoms.

```
       [context]─[variable]
                 ↑
           cut bond here
```

This captures: R-group changes, functional group swaps, ring opening/closing at a chain end.

### Double-cut (linker analysis)

Two bonds are broken; the central piece (linker) is the variable part and two contexts flank it.

```
[context₁]─[variable/linker]─[context₂]
            ↑              ↑
       cut bond 1      cut bond 2
```

Captures: linker changes between two pharmacophores; scaffold hops. Combinatorially larger — use with restraint.

## SMIRKS Representation of Transforms

A SMIRKS transform describes the variable parts:

```
[Fragment₁:1] >> [Fragment₂:1]

Examples:
  [CH3:1] >> [CF3:1]             # methyl → trifluoromethyl
  [OH:1] >> [NH2:1]              # hydroxyl → amine (bioisostere)
  [c:1][H:1] >> [c:1]F          # aromatic H → F (ring substitution)
  [CX4:1][Cl:1] >> [CX4:1]F     # Cl → F at sp3 carbon
  [N:1]([H:1]) >> [N:1]         # N-H → N-alkyl (if double cut captures alkyl)
  [*:1] >> [*:1]                 # generic placeholder (mmpdb uses [*:1])
```

**Canonical form**: mmpdb canonicalizes transforms so `A→B` and `B→A` are distinct (direction matters for ΔP sign).

## Property Cliffs

An **activity cliff** is an MMP where the structural change is small (Tanimoto > 0.7) but the property change is large (|ΔpIC50| > 1.0 log unit).

```
Activity cliff criteria (common thresholds):
- Structural similarity: Tanimoto(Morgan2, 2048) > 0.70
- Property difference:   |ΔpIC50| ≥ 1.0  (or |ΔLogD| ≥ 1.0)

SALI (Structure-Activity Landscape Index) per pair:
  SALI(A,B) = |Activity(A) - Activity(B)| / (1 - Similarity(A, B))
  SALI > 5 = cliff; > 10 = strong cliff
```

## Statistical SAR Framework

A **SAR rule** aggregates ΔP across all pairs sharing the same transform:

```
Transform R₁→R₂ paired with N instances:
  Mean(ΔP) = average effect of this transform
  Std(ΔP)  = variability (context dependence)
  N        = number of supporting pairs

Reliability thresholds:
  N ≥ 3   → minimum for reporting
  N ≥ 10  → reasonable confidence
  N ≥ 30  → high confidence, context-independent rule
```

**Context specificity**: some transforms are context-dependent.
- F→Cl on aromatic ring: ΔLogP = +0.70 ± 0.15 across 200 pairs (context-independent)
- OH→NHMe on heterocycle: ΔpIC50 = +0.3 ± 0.8 (N=12, context-dependent — large variance)

## Variable Part Size and Complexity

The variable part can range from a single atom (H→F) to a fragment (phenyl→pyridyl):

| Variable part size | Examples | Pairs found |
|--------------------|----------|-------------|
| 1 heavy atom | H→F, H→Cl, OH→NH₂ | Most abundant |
| 2-3 heavy atoms | CH₃→CF₃, NH₂→NHMe, COOH→CN | Common |
| 4-6 heavy atoms | phenyl→pyridyl, cyclohexyl→pip | Moderate |
| 7-10 heavy atoms | scaffold hops | Fewer, more specific |
| > 10 heavy atoms | mmpdb default cut-off | Not used by default |

## Fragmentation Algorithm (Hussain-Rea)

mmpdb uses the **Hussain-Rea** algorithm (JCIM 2010):

1. Identify all "heavy atom bonds" between non-ring atoms
2. For each single bond: break it to give variable and context parts
3. Represent variable part as SMILES with `[*:1]` attachment point
4. Represent context as SMILES with `[*:1]` where the bond was broken
5. Index all (context, variable_part) pairs across the database
6. Find MMPs: same context, different variable part

**Attachment point notation**:
```python
# Variable part for ethyl (at attachment point):
"CC[*:1]"     # context attaches at [*:1]

# Context for para-substituted benzene:
"c1ccc([*:1])cc1"  # substituent attaches at [*:1]
```

## R-Group Decomposition vs. MMPA

| Aspect | R-group decomposition | MMPA |
|--------|----------------------|------|
| Requires shared scaffold | Yes (explicit core) | No |
| Handles scaffold hops | No | Yes |
| Statistics | Positional | Transform-based |
| Scale | Tens of R-groups | Millions of pairs |
| Tool | RDKit RGroupDecomposition | mmpdb |

**Use MMPA when**: no obvious shared scaffold; diverse chemotypes; large dataset.
**Use R-group decomp when**: congeneric series; known core; small dataset.

## Common Bioisostere Transforms (literature-validated)

| Transform | ΔLogP | ΔpIC50 | Application |
|-----------|-------|--------|-------------|
| H → F | +0.14 ± 0.10 | variable | Block metabolism (CYP aromatic) |
| COOH → tetrazole | −1.0 ± 0.5 | ≈0 | Metabolic stability, oral BA |
| OH → NH₂ | −0.7 ± 0.4 | variable | H-bond donor change |
| Ph → pyridine | −0.7 ± 0.3 | variable | LogP reduction |
| CH₃ → CF₃ | +1.5 ± 0.3 | variable | LogP increase, metabolic block |
| CH₂CH₂ → CH=CH | −0.3 ± 0.2 | variable | Rigidification |
| NH → O (amide) | −0.5 ± 0.3 | variable | H-bond acceptor change |
| Cl → F | −0.2 ± 0.2 | variable | Metabolic block, reduced MW |

*Note: ΔpIC50 is highly context-dependent; use only as directional guidance.*

## MMPA vs. Free Energy Perturbation (FEP)

| Aspect | MMPA | FEP |
|--------|------|-----|
| Data needed | Experimental ΔP from database | Crystal structure + force field |
| Speed | Seconds | Hours-days |
| Accuracy | Statistical (Std ≈ 0.5-1.0 log) | ~1 kcal/mol |
| Novelty | Interpolation from known | Extrapolation possible |
| Scope | Any property with experimental data | Only binding ΔG |

**Complementarity**: use MMPA for rapid prioritization → confirm top analogues with FEP.
