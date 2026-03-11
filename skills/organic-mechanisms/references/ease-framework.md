# EASE Framework — Full Decision Logic

Source: AceOrganicChem.com E.A.S.E. method (2013), adapted.

---

## Step 1 — E: Identify Electrophile and Nucleophile

**Electrophile clues** (wants electrons, E+):
- Atom with full or partial positive charge
- Atom that gains positive charge via resonance (e.g., carbonyl carbon)
- Atom bearing a leaving group (creates electrophilic site on departure)
- Polarizable bond where one end is δ+ (e.g., C–Cl, C=O, CO₂, epoxide)

**Nucleophile clues** (has electrons, Nu−):
- Atom with full or partial negative charge
- Atom with lone pairs (O, N, S, halide)
- π system (alkene, enolate, aromatic ring with EDG)
- Negative charge stabilized by resonance = weaker Nu but still available

**If Nu is unclear**: draw all resonance structures. The most nucleophilic site will become obvious (highest electron density).

**Key signal**: label both E+ and Nu− on your drawing before proceeding.

---

## Step 2 — A: Acid/Base Check

**Rule**: If a strong acid or strong base is present, move the proton *before* doing anything else.

**Strong base present** → deprotonate the most acidic H available:
- Alpha-H of carbonyl (pKa ~20) → enolate
- Terminal alkyne H (pKa ~25)
- OH / NH (pKa 10–16)
- Depends on base strength (see `electrophiles-nucleophiles.md` pKa table)

**Strong acid present** → protonate the most basic site:
- Carbonyl oxygen (forms oxocarbenium, activates C=O as better E+)
- Alkene (Markovnikov addition of H+)
- Leaving group departure facilitated by protonation (e.g., OH → OH₂+ → H₂O)

**Lewis acid present** (AlCl₃, BF₃, ZnCl₂):
- Coordinates to most Lewis basic site (usually carbonyl O > ether O)
- Withdraws electron density → makes adjacent carbon more electrophilic
- "Cascading bully": activates 1,4 position in Michael acceptors

**Bulky base** (KOtBu, LDA, 2,6-lutidine, Et₃N):
- Cannot act as nucleophile → elimination only (E2), never SN2
- Deprotonates alpha-H; which alpha-H depends on approach geometry

**If no strong acid/base**: skip to Step 3.

---

## Step 3 — S: Sterics

**Where steric impedance occurs**:
1. On the electrophile (most common): blocked backside → SN1 or E
2. On the nucleophile: bulky base/Nu → shifts to base-only behavior

**Three outcomes when steric bulk is present**:
1. **Blocks reaction entirely** — e.g., norbornyl halide + NaCN (bridged system, backside inaccessible)
2. **Stabilizes carbocation intermediate** — tertiary C+ forms, SN1/E1 follows
3. **Overcome** — small Nu, high dilution, forced intramolecular geometry

**Bulky groups to flag immediately**:
- t-butyl, isopropyl (on electrophile → secondary/tertiary substrate)
- neo-pentyl (β-quaternary C, even though 1° substrate)
- Adamantane / bridged bicyclics (approach blocked from all sides)
- Bulky bases: KOtBu, LDA, 2,6-dimethylpyridine, Et₃N → elimination only

**FAQ**: t-Bu group = immediate flag for SN2 impossibility or E2 with bulky base.

---

## Step 4 — E: Electron Flow

**Golden rule**: electrons move from Nu− to E+. Arrow tip points toward destination.

**Arrow types**:
- Full curved arrow = electron pair movement
- Half-headed (fish-hook) = single electron (radical only)

**Bond counting check after each step**:
- C normally forms 4 bonds
- N: 3 bonds (neutral) or 4 bonds (positive charge)
- O: 2 bonds (neutral) or 3 bonds (positive charge, oxonium)
- Halogens: 1 bond (neutral)
- If E+ already has 4 bonds → must break one bond simultaneously (e.g., C=O → C–O−)

**Iterate**: after Step 4, look at the new intermediate. Does it still have reactive sites? Repeat E–A–S–E from Step 1 on the new intermediate.

**Termination**: reaction is complete when all charged sites are neutralized (or you've reached the target product).

---

## Worked Example: Aldol (from book)

```
Reagents: 2 × CH₃CHO + NaOH

E-Step 1: C=O is E+; nucleophile unclear
A-Step 2: NaOH is base; alpha-H is acid → deprotonate alpha carbon → enolate (Nu−)
S-Step 3: No steric issue
E-Step 4: Enolate (Nu−) attacks C=O (E+) intermolecularly → beta-hydroxyaldehyde
```

**Repeat cycle**: after the C–C bond forms, the alkoxide is protonated by solvent → product.
