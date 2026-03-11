# Electron Flow, Arrow Pushing & HSAB

---

## Arrow Conventions

| Arrow | Meaning |
|---|---|
| Curved full arrow (→) | Electron pair movement |
| Fish-hook (⇀) | Single electron (radical only) |
| Double-headed (↔) | Resonance (not a reaction step) |
| Straight → | Irreversible reaction |
| ⇌ | Equilibrium |
| ⟹ | Retrosynthesis disconnection |

**Direction rule**: arrowhead points to where electrons *go* (from Nu− to E+, from lone pair to electrophilic center, from π bond to adjacent E+).

---

## Valence Check After Each Step

After drawing each arrow, verify:
- C: 4 bonds (no lone pairs in neutral state)
- O: 2 bonds neutral; 3 bonds = oxonium (+1); 1 bond + 3 LP = alkoxide (−1)
- N: 3 bonds neutral; 4 bonds = ammonium (+1)
- Halogens: 1 bond neutral; 0 bonds = halide (−1)
- Total charge and electron count must be conserved across each step

---

## HSAB — 1,2 vs 1,4 Attack

Hard/Soft Acid-Base theory (Pearson) determines regioselectivity on conjugated systems:

### Hard Nucleophiles → 1,2 addition (attack at carbonyl carbon)
Examples: RO⁻, HO⁻, RMgX (Grignard), RLi, hydride (LiAlH₄, NaBH₄)

### Soft Nucleophiles → 1,4 addition (Michael addition, attack at β-carbon)
Examples: RS⁻, I⁻, R₂CuLi (always 1,4), enolates (carbanions), CN⁻ (borderline)

```
Hard Nu (RO⁻, OH⁻, RMgX)     →  attacks C=O (1,2)
Soft Nu (RS⁻, R₂CuLi, enolates) → attacks β-carbon (1,4)
```

**Organocuprates (R₂CuLi)**: ALWAYS 1,4. This is their defining feature vs Grignard.

**Under thermodynamic control (reversible conditions)**: 1,4 product often preferred (more stable).
**Under kinetic control (−78 °C, irreversible)**: hard Nu at 1,2 prevails.

---

## Intramolecular vs Intermolecular

Check whether Nu and E+ are on the **same molecule**:
- If yes → intramolecular reaction is possible and often preferred if a 5- or 6-membered ring forms
- Dilute conditions favor intramolecular (unimolecular rate independent of concentration)
- Ring sizes 5 and 6 are kinetically and thermodynamically preferred; 3 forms under forcing conditions (SN2 within molecule); 7+ rings require high dilution

**Intramolecular examples**: lactone formation, epoxide formation (KOH + halohydrin), Dieckmann condensation.

---

## Carbonyl Reactivity Hierarchy

More electrophilic at C=O → more reactive toward Nu:

```
Acid chloride > Anhydride > Aldehyde > Ketone > Ester > Amide > Carboxylate
```

Reasoning: electron-withdrawing attached groups (Cl, O-acyl) withdraw density from C=O; electron-donating groups (NR₂, OR) donate via resonance → stabilize C=O → less electrophilic.

---

## Common Mechanisms Step by Step

### Nucleophilic Addition to Carbonyl (Nu + C=O)
```
1. Nu− attacks carbonyl C → tetrahedral alkoxide intermediate
2. Protonate O (workup, H₃O⁺) → alcohol product
If LG on carbonyl (ester, acyl Cl): LG leaves after tetrahedral intermediate → substitution
```

### Enolate Formation + Alkylation
```
1. Base deprotonates alpha-H → enolate (Nu−)
2. Enolate attacks electrophile (alkyl halide, carbonyl) → C–C bond
3. Protonation of O during workup
```

### EAS (Electrophilic Aromatic Substitution)
```
1. Generate electrophile (E+) from reagent (HNO₃/H₂SO₄ → NO₂⁺, etc.)
2. Pi electrons of ring attack E+ → arenium (Wheland) intermediate
3. Deprotonation restores aromaticity
```

### SNAr (Aromatic Nucleophilic Substitution)
```
Requires: strong EWG at ortho/para to LG, good Nu−
1. Nu− adds to ipso carbon → Meisenheimer complex (anionic sigma complex)
2. LG⁻ departs → substituted aromatic
```
