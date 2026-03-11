# Electrophiles, Nucleophiles, pKa & Nucleophilicity

---

## pKa Table (aqueous, 25 °C)

| Acid | pKa | Conjugate Base |
|------|-----|----------------|
| HI | −10 | I⁻ |
| HBr | −9 | Br⁻ |
| HCl | −7 | Cl⁻ |
| H₂SO₄ | −3 | HSO₄⁻ |
| H₃O⁺ | −1.7 | H₂O |
| HF | 3.2 | F⁻ |
| RCOOH (acetic) | 4.75 | RCOO⁻ |
| PhOH | 10 | PhO⁻ |
| RSH (thiol) | 10–11 | RS⁻ |
| HCO₃⁻ | 10.3 | CO₃²⁻ |
| RNH₂ (ammonium) | ~10 | RNH₂ |
| H₂O | 15.7 | HO⁻ |
| ROH (MeOH) | 15.5 | MeO⁻ |
| ROH (EtOH) | 16 | EtO⁻ |
| ROH (t-BuOH) | 19 | t-BuO⁻ |
| Alpha-H of ketone | ~20 | enolate |
| Terminal alkyne | 25 | RC≡C⁻ |
| Alpha-H of ester | ~25 | ester enolate |
| Alpha-H of nitrile | ~25 | α-cyano carbanion |
| NH₃ | 38 | NH₂⁻ (NaNH₂) |
| RH (sp³ C–H) | ~50 | R⁻ (alkyl) |

**Rule**: A base deprotonates only if its conjugate acid pKa > substrate pKa (by ≥ 2 units for reliable reaction). Stronger base = higher pKa conjugate acid.

---

## Common Bases and Their Effective pKa Range

| Base | Effective deprotonation | Notes |
|------|------------------------|-------|
| NaOH / KOH | pKa < 15.7 | carboxylic acids, phenols |
| NaOMe / NaOEt | pKa < 16 | esters, carboxylic acids |
| K₂CO₃ | pKa < 10 | thiols, phenols, soft alkylation |
| KOtBu | pKa < 19 | alcohols, alpha-H of ketones; **bulky → E2** |
| NaNH₂ | pKa < 38 | terminal alkynes |
| LDA | pKa < 36 | ketone alpha-H kinetically (−78 °C); **bulky → E2 only** |
| n-BuLi | pKa < 50 | almost anything |

---

## Nucleophilicity

### Polar Protic Solvents (H₂O, ROH, RCOOH)
Solvation stabilizes small/hard anions more → nucleophilicity follows **polarizability**:

```
I⁻ > Br⁻ > Cl⁻ > F⁻
RS⁻ >> RO⁻
```

### Polar Aprotic Solvents (DMSO, DMF, acetone, MeCN, THF)
No H-bonding solvation → nucleophilicity follows **basicity**:

```
F⁻ > Cl⁻ > Br⁻ > I⁻
RO⁻ > RS⁻
```

**Key consequence**: F⁻, CN⁻, N₃⁻, RO⁻ are far better nucleophiles in DMSO/DMF than in water. SN2 reactions are dramatically faster in polar aprotic solvents.

### Nucleophilicity Ranking (general, polar aprotic)
```
RS⁻ ≈ I⁻ > CN⁻ ≈ N₃⁻ > Br⁻ > RO⁻ ≈ HO⁻ > Cl⁻ > F⁻ > H₂O ≈ ROH
R₂NH > RNH₂ > NH₃ > ArNH₂
```

### Nucleophile vs. Base
A nucleophile attacks carbon (C–X bond). A base attacks hydrogen (C–H or X–H bond).
Same species, different role — determined by substrate and steric context.

---

## Common Electrophile Classes

| Electrophile | E+ site | Reactivity |
|---|---|---|
| Alkyl halide (1°) | C bearing X | SN2 > E2 |
| Alkyl halide (2°) | C bearing X | SN2 / E2 / SN1 depends on conditions |
| Alkyl halide (3°) | C bearing X | SN1 / E1 / E2 |
| Carbonyl C=O | carbonyl C | 1,2-addition (hard Nu) |
| α,β-unsaturated carbonyl | β-C (soft), carbonyl C (hard) | 1,4 Michael (soft Nu) |
| Epoxide | C–O (less hindered end) | SN2 under basic; SN1-like under acid at more substituted C |
| Acyl chloride (RCOCl) | carbonyl C | very reactive, Nu must not be bulky |
| Aromatic ring + EWG | ipso (SNAr) or electrophilic site (EAS) | requires activation |
| Proton (H⁺) | H+ | protonates most basic lone pair / π bond |

---

## Leaving Group Ability

Good LG = stable after departure (weak base, pKa conjugate acid < 0):

```
Best  →  TfO⁻ > I⁻ > TsO⁻ > Br⁻ > Cl⁻ >> F⁻ > RO⁻ > HO⁻ > NH₂⁻
```

**Trick**: Convert bad LG (OH, NH₂) to good LG before SN2 — activate with SOCl₂, TsCl, or H⁺ (→ OH₂⁺).
