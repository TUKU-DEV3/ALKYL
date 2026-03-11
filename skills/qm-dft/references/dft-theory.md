# DFT Theory — Functionals, Basis Sets, Dispersion

## Jacob's Ladder of DFT Functionals

```
Heaven (exact) ─────────────────────────────────────────────
Rung 5: Double hybrids   — B2-PLYP, ωB97M(2), DSD-PBEP86
         (DFT + MP2 correlation, O(N⁵), very accurate)
Rung 4: Hybrids          — B3LYP, PBE0, ωB97X-D, M06-2X
         (mix HF exchange, O(N⁴))
Rung 3: Meta-GGA         — TPSS, M06-L, r²SCAN, SCAN
         (uses kinetic energy density τ, O(N³))
Rung 2: GGA              — PBE, BLYP, BP86
         (gradient of density, O(N³))
Rung 1: LDA              — SVWN, VWN
         (density only, O(N³), overbinds)
Earth (Hartree-Fock)  ──────────────────────────────────────
```

---

## Functional Recommendations

### General Recommendations (drug-like organic molecules)

| Task | Functional | Basis set | Notes |
|------|-----------|-----------|-------|
| Geometry optimization | B3LYP-D3BJ | def2-SVP | Fast, reliable |
| Geometry (recommended) | r²SCAN-3c | def2-mTZVPP | Composite, includes D4 + gCP |
| Single-point energy | B3LYP-D3BJ | def2-TZVP | On opt geometry |
| Single-point (better) | ωB97X-D | def2-TZVP | Better for CT, barriers |
| Reaction barriers | ωB97X-D or M06-2X | def2-TZVP | Hybrids needed |
| Non-covalent interactions | B3LYP-D3BJ or BLYP-D3BJ | def2-TZVP | D3BJ critical |
| NMR shifts | B3LYP or PBE0 | def2-TZVP | Gauge-including AOs (GIAO) |
| UV-Vis (TD-DFT) | ωB97X-D or CAM-B3LYP | def2-TZVP | Range-separated for CT states |
| High-accuracy benchmark | DLPNO-CCSD(T) | def2-TZVPP or CBS | Single-point on DFT opt |
| Fast pre-screening | GFN2-xTB | — | Semi-empirical (see xtb ref) |

### Composite Methods (one keyword, implicit corrections)
```
r²SCAN-3c   = r²SCAN + D4 dispersion + gCP + modified basis (def2-mTZVPP)
B3LYP-3c    = B3LYP + D3BJ + gCP + def2-mSVP
PBEh-3c     = PBEh + D3BJ + gCP + def2-mSVP
HF-3c       = HF + three-body correction + STO-nG basis
```
→ `-3c` methods designed for geometry optimization; good accuracy/cost ratio.

---

## Basis Sets

### Pople (older, common in older literature)
```
6-31G         — minimal DZ
6-31G*        — + polarization on heavy atoms (= 6-31G(d))
6-31G**       — + polarization on H (= 6-31G(d,p))
6-31+G*       — + diffuse on heavy atoms (for anions)
6-311+G(2d,p) — TZ quality
```

### Ahlrichs def2 (recommended for ORCA/Gaussian)
```
def2-SV(P)    — DZ quality, fast
def2-SVP      — DZ + polarization (standard optimization)
def2-TZVP     — TZ + polarization (single-point energy)
def2-TZVPP    — TZ + double polarization (more accurate)
def2-QZVP     — QZ (benchmark)
```

### Dunning (correlated methods, CBS extrapolation)
```
cc-pVDZ       — DZ (Dunning)
cc-pVTZ       — TZ
cc-pVQZ       — QZ
aug-cc-pVDZ   — + diffuse functions (anions, weak interactions)
aug-cc-pVTZ   — standard for CCSD(T) single points
```

### Auxiliary Basis Sets (RI/DF approximation — speeds up DFT ~10×)
```
def2/J        — Coulomb fitting (RI-J): use with all DFT
def2-TZVP/C   — correlation fitting: MP2, CCSD
def2/JK       — exchange fitting: hybrid functionals
```
ORCA keyword: `! RI-J def2/J` or auto-enabled with `! RIJCOSX`

### Heavy Elements (relativistic)
```
def2-TZVP    — built-in ECP for elements Z > 36 (Kr)
SARC-DKH     — for very heavy elements (ORCA)
```

---

## Dispersion Corrections (mandatory for non-covalent)

Without dispersion, DFT severely underestimates π–π stacking, van der Waals, and hydrophobic interactions.

| Correction | Notes |
|-----------|-------|
| D3(0) | Becke-Johnson damping predecessor |
| **D3BJ** | Standard, recommended for most functionals |
| **D4** | More accurate, accounts for charge dependence |
| XDM | Exchange-dispersion model (Psi4) |
| VV10 | Non-local correlation (SCAN+rVV10) |

```
ORCA keyword:   ! B3LYP D3BJ def2-SVP
Gaussian:       # B3LYP/def2SVP EmpiricalDispersion=GD3BJ
PySCF:          mol.DFT().with_dftd3(xc='b3lyp', version='d3bj')
```

**D3BJ parameters are functional-dependent** — the program looks them up automatically when you specify the functional name.

---

## Solvation Models

| Model | Type | Cost | Notes |
|-------|------|------|-------|
| CPCM | Implicit | Negligible | Cavity + dielectric; standard |
| COSMO | Implicit | Negligible | Similar to CPCM |
| SMD | Implicit + non-elec | Negligible | Recommended for ΔG_solv |
| ALPB (xTB) | Analytical | Negligible | Fast, xTB-specific |
| GBSA (xTB) | Generalized Born | Negligible | xTB solvation |
| Explicit water | QM water | High | For specific H-bonding |

```
ORCA:    ! CPCM(water)
         or %cpcm smd true solvent "water" end
Gaussian: # B3LYP/def2SVP SCRF=(SMD,Solvent=Water)
PySCF:   mol.with_solvent.DDCOSMO('water')
```

---

## Understanding Key Numbers

### Energy units
```
1 Hartree (Eh) = 627.51 kcal/mol = 2625.5 kJ/mol = 27.211 eV
1 kcal/mol = 0.001594 Eh = 4.184 kJ/mol
kBT at 298K = 0.593 kcal/mol ≈ 1 kcal/mol (chemical accuracy threshold)
```

### Convergence thresholds (ORCA defaults)
```
TightSCF:   energy 1e-8 Eh; gradient 1e-5 Eh/Bohr
VeryTightSCF: 1e-9 Eh
TightOpt:   gradient RMS 3e-5; max force 1e-4; RMSD 2e-3
```

### When DFT fails
- **Charge transfer (CT) excitations** → range-separated functionals (ωB97X-D, CAM-B3LYP)
- **Multi-reference character** (diradicals, bond breaking) → CASSCF/NEVPT2 or DLPNO-CCSD(T)
- **Strong correlation** (transition metal d-electrons) → DFT+U or CASSCF
- **Van der Waals complexes without D3** → add dispersion correction always
- **Anions / diffuse density** → add diffuse basis functions (aug-/+)
