# Force Field Fundamentals

## The Molecular Mechanics Energy Function

Total potential energy is a sum of bonded and non-bonded contributions:

```
E_total = E_bonds + E_angles + E_torsions + E_impropers + E_vdW + E_electrostatics
```

All terms are classical approximations — no electrons, no bond breaking/forming.

---

## Bonded Terms

### Bond Stretching (harmonic)
```
E_bond = k_b (r - r₀)²
```
- `k_b` — force constant (kcal/mol/Å²)
- `r₀` — equilibrium bond length (Å)
- Morse potential (exponential) more accurate but rarely used in MD

### Angle Bending (harmonic)
```
E_angle = k_θ (θ - θ₀)²
```
- `θ₀` — equilibrium angle (degrees)
- CHARMM uses Urey-Bradley cross term: adds k_UB(r_1,3 - S₀)² between atoms 1 and 3

### Torsion (dihedral)
```
E_torsion = Σₙ Vₙ/2 · [1 + cos(nφ - δ)]
```
- `Vₙ` — barrier height (kcal/mol), `n` — periodicity, `δ` — phase angle
- Multiple n terms summed for accuracy
- Drives rotameric states, conformational equilibria

### Improper Dihedral (out-of-plane)
```
E_improper = k_ξ (ξ - ξ₀)²   [CHARMM-style]
or  E_improper = Vₙ/2 · [1 + cos(nφ - δ)]  [AMBER-style]
```
- Maintains planarity of sp² centers (aromatic rings, peptide bonds, carbonyls)

---

## Non-Bonded Terms

### van der Waals — Lennard-Jones 12-6
```
E_LJ = ε [(r_min/r)¹² - 2(r_min/r)⁶]
     = 4ε [(σ/r)¹² - (σ/r)⁶]
```
- `ε` — well depth (energy minimum), `r_min` = 2^(1/6) σ
- AMBER uses r_min; CHARMM uses r_min/2; OpenFF uses σ
- 12-6 is a compromise — repulsive (12) + attractive dispersion (6)
- Combining rules: AMBER uses Lorentz-Berthelot; OPLS uses geometric mean

### Electrostatics — Coulomb
```
E_elec = q_i · q_j / (4πε₀ · r_ij)
```
- Partial charges `q` assigned per atom (AM1-BCC, RESP, DDEC)
- Long-range (1/r) — requires PME for periodic systems
- 1-4 interactions scaled: AMBER 1/1.2 (charge) and 1/2.0 (LJ); CHARMM 1.0 with explicit pairs

---

## Force Field Families

### AMBER
- **ff14SB** — standard protein FF (2015); backbone φ/ψ recalibrated vs ff99SB
- **ff19SB** — improved backbone/sidechain torsions (2020); requires OPC water
- **GAFF2** — General Amber Force Field 2 for organic small molecules
- **lipid21** — phospholipids; **GLYCAM06** — carbohydrates
- Files: `.prmtop` / `.inpcrd` or `.parm7` / `.rst7`
- Parameter units: energy in kcal/mol, distance in Å

### CHARMM
- **CHARMM36m** — protein FF optimized for disordered regions (2017)
- **CGenFF** — General FF for drug-like molecules (CHARMM parametrization)
- Uses Urey-Bradley + CMAP torsion correction maps
- Files: `.psf` (topology) + `.dcd` (trajectory) + `.prm` (parameters)

### OPLS
- **OPLS-AA/M** — optimized for condensed phase (solvation free energies)
- **OPLS3e** — Schrödinger commercial, good for drug-like molecules
- All-atom, geometric combining rules

### SMIRNOFF / Open Force Field
- **Sage (openff-2.2.0)** — current recommended for drug-like molecules (2023)
- **Parsley (1.3.1)** — predecessor; Sage preferred
- **Espaloma** — ML-based typing (message-passing NN, less validated)
- SMARTS-based typing: no atom types, parameters attached to chemical patterns
- Files: `.offxml` format

### GROMOS
- United-atom model (aliphatic CH treated as single particle → faster)
- **54A7** / **53A6** — protein + lipid
- Mainly GROMACS ecosystem; less common in Python workflows

---

## Choosing a Force Field

| System | Recommended |
|--------|-------------|
| Standard protein in water | ff14SB + TIP3P |
| Protein with disordered regions | CHARMM36m + TIP3P-CHARMM |
| Drug-like small molecule | OpenFF Sage 2.2 |
| Drug-like (GROMACS workflow) | GAFF2 via acpype |
| Protein-ligand complex | ff14SB + OpenFF Sage (openmmforcefields) |
| Protein-ligand (CHARMM users) | CHARMM36m + CGenFF |
| Lipid bilayer | CHARMM36 + lipid21 |
| Ionic liquid / unusual molecule | CGenFF or custom GAFF2 |
| Materials / periodic solid | ReaxFF, UFF, DREIDING (via ASE/LAMMPS) |

---

## Key Concepts

### 1-2, 1-3, 1-4 Exclusions
- 1-2 (bonded) and 1-3 (angle) non-bonded interactions: **excluded** (= 0)
- 1-4 (across dihedral): **scaled** — reduces double-counting with torsion parameters

### Combining Rules
| Rule | σ_ij | ε_ij |
|------|------|------|
| Lorentz-Berthelot (AMBER) | (σ_i + σ_j)/2 | √(ε_i · ε_j) |
| Geometric (OPLS, GROMOS) | √(σ_i · σ_j) | √(ε_i · ε_j) |
| Waldman-Hagler | (σ_i³ + σ_j³)^(1/3)/2 | rare |

### Covalent Bond Model Limitation
Force fields assume fixed connectivity — no bond breaking or formation. For reactive MD, use ReaxFF or QM/MM.
