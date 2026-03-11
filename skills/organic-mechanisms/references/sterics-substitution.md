# Sterics, Substitution & Elimination — Decision Tree

---

## SN1 / SN2 / E1 / E2 Master Decision Tree

```
Substrate?
│
├─ Primary (1°)
│   ├─ Strong Nu, polar aprotic → SN2
│   ├─ Strong bulky base → E2 (only if alpha-H available)
│   └─ Weak base/Nu → no reaction (1° doesn't form carbocation)
│
├─ Secondary (2°)
│   ├─ Strong Nu, polar aprotic, small Nu → SN2
│   ├─ Strong bulky base → E2 (Hofmann product if very bulky)
│   ├─ Polar protic, weak Nu/base → SN1 / E1 mixture
│   └─ Strong non-bulky base, protic → E2 (Zaitsev)
│
└─ Tertiary (3°)
    ├─ Any Nu in polar protic → SN1 + E1 (SN1 dominates if Nu is good, T room)
    ├─ Heat + polar protic → E1 (Zaitsev)
    ├─ Strong base → E2 (100% elimination)
    └─ SN2 is essentially impossible (steric block)
```

---

## Conditions Summary Table

| Conditions | Prediction |
|---|---|
| 1° substrate + strong Nu + DMSO/DMF | SN2 |
| 1° substrate + KOtBu / LDA | E2 |
| 2° substrate + NaI, DMSO | SN2 |
| 2° substrate + KOtBu, EtOH | E2 (Zaitsev or Hofmann depending on base) |
| 2° substrate + H₂O, weak acid | SN1 / E1 |
| 3° substrate + EtOH, heat | E1 (Zaitsev) |
| 3° substrate + KOtBu | E2 |
| 3° substrate + H₂O/ROH, room T | SN1 |

---

## Regiochemistry of Elimination

### Zaitsev Rule (more substituted alkene)
- Forms when base is **not bulky** (NaOEt, NaOH, NaOMe)
- Thermodynamic product (more stable)

### Hofmann Rule (less substituted alkene)
- Forms when base is **very bulky** (KOtBu, LDA, 2,6-lutidine)
- Kinetic product — base approaches least hindered alpha-H

---

## Stereochemistry

| Mechanism | Stereochemical Outcome |
|---|---|
| SN2 | **Inversion** (Walden) — backside attack only |
| SN1 | **Racemization** (flat carbocation, attack from both faces) |
| E2 | **Anti-periplanar** required — trans H and LG must be coplanar |
| E1 | Mixed cis/trans (follows Zaitsev, no stereo constraint) |

**E2 consequence**: cyclohexane substrates need diaxial conformation of H and LG. If equatorial LG → must ring-flip first.

---

## Steric Groups to Flag Immediately

| Group | Effect |
|---|---|
| t-Bu (C(CH₃)₃) | Severe — blocks SN2 entirely; E2/E1 only |
| Isopropyl | Moderate — disfavors SN2 on 2° substrate |
| Neo-pentyl (CH₂C(CH₃)₃) | 1° substrate but β-quaternary C → SN2 very slow |
| Bridged bicyclics (norbornane, adamantane) | Backside attack geometrically impossible → no SN2 |

**Bulky Nu/base** → can only act as base:
- KOtBu, LDA, 2,6-lutidine, Et₃N → E2 exclusively
- Exception: very strained systems where even E2 is slow

---

## Solvent Effect on SN1 vs SN2

| Solvent type | Examples | Favors |
|---|---|---|
| Polar protic | H₂O, EtOH, AcOH | SN1 / E1 (stabilizes carbocation) |
| Polar aprotic | DMSO, DMF, acetone, MeCN | SN2 (desolvates Nu−) |
| Nonpolar | hexane, Et₂O | Rare — usually organolithium chemistry |
