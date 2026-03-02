# RDKit — Properties & Descriptors

## Quick Reference: Key Descriptors

```python
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Crippen, QED

mol = Chem.MolFromSmiles('CC(=O)Oc1ccccc1C(=O)O')  # aspirin
```

| Property | Function | Value (aspirin) |
|----------|----------|-----------------|
| Molecular Weight | `Descriptors.MolWt(mol)` | 180.16 |
| Exact MW | `Descriptors.ExactMolWt(mol)` | 180.042 |
| Heavy Atom MW | `Descriptors.HeavyAtomMolWt(mol)` | 172.14 |
| LogP (Crippen) | `Descriptors.MolLogP(mol)` | 1.31 |
| TPSA | `Descriptors.TPSA(mol)` | 63.6 |
| HBD (Lipinski) | `rdMolDescriptors.CalcNumLipinskiHBD(mol)` | 1 |
| HBA (Lipinski) | `rdMolDescriptors.CalcNumLipinskiHBA(mol)` | 4 |
| HBD (Ertl) | `rdMolDescriptors.CalcNumHBD(mol)` | 1 |
| HBA (Ertl) | `rdMolDescriptors.CalcNumHBA(mol)` | 3 |
| Rotatable bonds | `rdMolDescriptors.CalcNumRotatableBonds(mol)` | 3 |
| Rings | `rdMolDescriptors.CalcNumRings(mol)` | 1 |
| Aromatic rings | `rdMolDescriptors.CalcNumAromaticRings(mol)` | 1 |
| Molecular formula | `rdMolDescriptors.CalcMolFormula(mol)` | C9H8O4 |
| Fsp3 | `rdMolDescriptors.CalcFractionCSP3(mol)` | 0.0 |
| QED | `QED.qed(mol)` | 0.55 |
| Amide bonds | `rdMolDescriptors.CalcNumAmideBonds(mol)` | 0 |
| Stereocenters | `rdMolDescriptors.CalcNumAtomStereoCenters(mol)` | 0 |

---

## Lipinski Rule of Five (Ro5)

```python
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

def lipinski_ro5(mol) -> dict:
    """Returns Ro5 pass/fail for each criterion + overall."""
    mw   = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd  = rdMolDescriptors.CalcNumLipinskiHBD(mol)
    hba  = rdMolDescriptors.CalcNumLipinskiHBA(mol)

    return {
        'MW':   (mw,   mw   <= 500),
        'LogP': (logp, logp <= 5.0),
        'HBD':  (hbd,  hbd  <= 5),
        'HBA':  (hba,  hba  <= 10),
        'pass': all([mw <= 500, logp <= 5.0, hbd <= 5, hba <= 10])
    }

mol = Chem.MolFromSmiles('CC(=O)Oc1ccccc1C(=O)O')
print(lipinski_ro5(mol))
```

---

## Drug-Likeness Filters

### Lipinski Ro5
| Rule | Threshold | Descriptor |
|------|-----------|------------|
| MW ≤ 500 Da | `MolWt(mol) <= 500` | `Descriptors.MolWt` |
| LogP ≤ 5 | `MolLogP(mol) <= 5` | `Descriptors.MolLogP` |
| HBD ≤ 5 | `CalcNumLipinskiHBD(mol) <= 5` | `rdMolDescriptors` |
| HBA ≤ 10 | `CalcNumLipinskiHBA(mol) <= 10` | `rdMolDescriptors` |

### Veber Rules (oral bioavailability)
```python
def veber(mol) -> bool:
    rotb = rdMolDescriptors.CalcNumRotatableBonds(mol)
    tpsa = Descriptors.TPSA(mol)
    return rotb <= 10 and tpsa <= 140
```

### Ghose Filter
```python
def ghose(mol) -> bool:
    mw   = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    mr   = Crippen.MolMR(mol)
    na   = mol.GetNumHeavyAtoms()
    return (160 <= mw <= 480 and
            -0.4 <= logp <= 5.6 and
            40 <= mr <= 130 and
            20 <= na <= 70)
```

### Lead-like (Oprea)
```python
def lead_like(mol) -> bool:
    mw   = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd  = rdMolDescriptors.CalcNumHBD(mol)
    hba  = rdMolDescriptors.CalcNumHBA(mol)
    rotb = rdMolDescriptors.CalcNumRotatableBonds(mol)
    return (mw <= 350 and logp <= 3.5 and
            hbd <= 3 and hba <= 7 and rotb <= 7)
```

### PAINS (Pan Assay Interference) via SMARTS
```python
from rdkit.Chem import FilterCatalog

params = FilterCatalog.FilterCatalogParams()
params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_A)
params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_B)
params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_C)
catalog = FilterCatalog.FilterCatalog(params)

def is_pains(mol) -> bool:
    return catalog.HasMatch(mol)

def get_pains_alerts(mol) -> list[str]:
    matches = catalog.GetMatches(mol)
    return [m.GetDescription() for m in matches]
```

---

## QED (Quantitative Estimate of Drug-likeness)

```python
from rdkit.Chem import QED

mol = Chem.MolFromSmiles('CC(=O)Oc1ccccc1C(=O)O')
qed_score = QED.qed(mol)          # float 0-1, higher = more drug-like

# Full QED properties breakdown
props = QED.properties(mol)
# QEDproperties(MW=180.157, ALOGP=1.312, HBA=4, HBD=1, PSA=63.6,
#               ROTB=3, AROM=1, ALERTS=0)
```

---

## Full Descriptor Set

```python
from rdkit.Chem import Descriptors

# Calculate all ~200 descriptors at once
mol = Chem.MolFromSmiles('c1ccccc1C(=O)O')
all_desc = Descriptors.CalcMolDescriptors(mol)
# → dict {'MolWt': 122.12, 'MolLogP': 1.36, 'TPSA': 37.3, ...}

# List all available descriptor names
names = [name for name, _ in Descriptors.descList]
print(f"{len(names)} descriptors available")

# Calculate as DataFrame (for a set of molecules)
import pandas as pd

def descriptors_df(mols: list) -> pd.DataFrame:
    rows = []
    for mol in mols:
        if mol is None:
            continue
        d = Descriptors.CalcMolDescriptors(mol)
        d['smiles'] = Chem.MolToSmiles(mol)
        rows.append(d)
    return pd.DataFrame(rows)
```

### Selected important descriptors from `Descriptors` module

| Category | Descriptor | Description |
|----------|-----------|-------------|
| **Weight** | `MolWt` | Average MW |
| | `ExactMolWt` | Exact monoisotopic MW |
| **Lipophilicity** | `MolLogP` | Crippen logP |
| | `MolMR` | Molar refractivity |
| **Polarity** | `TPSA` | Topological PSA |
| **H-bonding** | `NumHDonors` | H-bond donors (Ertl) |
| | `NumHAcceptors` | H-bond acceptors (Ertl) |
| **Flexibility** | `NumRotatableBonds` | Rotatable bonds |
| | `FractionCSP3` | sp3 C fraction |
| **Size** | `HeavyAtomCount` | Non-H atom count |
| **Rings** | `RingCount` | Total ring count |
| | `NumAromaticRings` | Aromatic ring count |
| | `NumAliphaticRings` | Aliphatic ring count |
| **Electronic** | `MaxPartialCharge` | Max Gasteiger charge |
| | `MinPartialCharge` | Min Gasteiger charge |
| | `MaxAbsPartialCharge` | Max absolute charge |
| **Connectivity** | `BertzCT` | Bertz complexity |
| | `Chi0`, `Chi1`, `Chi2v`… | Chi path indices |
| | `Kappa1`, `Kappa2`, `Kappa3` | Kappa shape indices |
| **Fingerprint** | `FpDensityMorgan1/2/3` | Morgan FP density |

---

## rdMolDescriptors — Fast C++ Descriptors

```python
from rdkit.Chem import rdMolDescriptors

mol = Chem.MolFromSmiles('CC(N)C(=O)O')   # alanine

# Core Lipinski / ADME
rdMolDescriptors.CalcNumHBD(mol)                     # H-bond donors (Ertl definition)
rdMolDescriptors.CalcNumHBA(mol)                     # H-bond acceptors (Ertl)
rdMolDescriptors.CalcNumLipinskiHBD(mol)             # Lipinski HBD
rdMolDescriptors.CalcNumLipinskiHBA(mol)             # Lipinski HBA
rdMolDescriptors.CalcTPSA(mol)                       # TPSA (polar surface area)
rdMolDescriptors.CalcTPSA(mol, includeSandP=True)    # include S and P in PSA

# Structural
rdMolDescriptors.CalcNumRotatableBonds(mol)
rdMolDescriptors.CalcNumRotatableBonds(mol, rdMolDescriptors.NumRotatableBondsOptions.Strict)
rdMolDescriptors.CalcNumRings(mol)
rdMolDescriptors.CalcNumAromaticRings(mol)
rdMolDescriptors.CalcNumHeterocycles(mol)
rdMolDescriptors.CalcNumAromaticHeterocycles(mol)
rdMolDescriptors.CalcNumAliphaticRings(mol)
rdMolDescriptors.CalcNumSaturatedRings(mol)
rdMolDescriptors.CalcNumBridgeheadAtoms(mol)
rdMolDescriptors.CalcNumSpiroAtoms(mol)
rdMolDescriptors.CalcNumAmideBonds(mol)

# Composition
rdMolDescriptors.CalcMolFormula(mol)                 # 'C3H7NO2'
rdMolDescriptors.CalcFractionCSP3(mol)               # Fsp3
rdMolDescriptors.CalcNumAtomStereoCenters(mol)
rdMolDescriptors.CalcNumUnspecifiedAtomStereoCenters(mol)

# Crippen
logp, mr = rdMolDescriptors.CalcCrippenDescriptors(mol)

# AUTOCORR descriptors (2D)
autocorr = rdMolDescriptors.CalcAUTOCORR2D(mol)      # 192-element vector
```

---

## Gasteiger Charges

```python
from rdkit.Chem import AllChem

mol = Chem.MolFromSmiles('CC(=O)O')
AllChem.ComputeGasteigerCharges(mol)

for atom in mol.GetAtoms():
    charge = float(atom.GetProp('_GasteigerCharge'))
    print(f"Atom {atom.GetIdx()} ({atom.GetSymbol()}): {charge:.4f}")
```

---

## 3D Descriptors (require conformer)

```python
from rdkit.Chem import rdMolDescriptors, AllChem

mol = Chem.MolFromSmiles('CC(N)C(=O)O')
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
AllChem.MMFFOptimizeMolecule(mol)

# Shape descriptors
rdMolDescriptors.CalcAsphericity(mol)          # 0 = sphere, 1 = rod
rdMolDescriptors.CalcEccentricity(mol)
rdMolDescriptors.CalcInertialShapeFactor(mol)  # spherical vs linear
rdMolDescriptors.CalcSpherocityIndex(mol)      # 0-1
rdMolDescriptors.CalcRadiusOfGyration(mol)     # Angstroms

# Principal moments of inertia
rdMolDescriptors.CalcPMI1(mol)
rdMolDescriptors.CalcPMI2(mol)
rdMolDescriptors.CalcPMI3(mol)

# Descriptor vectors (for QSAR models)
rdf   = rdMolDescriptors.CalcRDF(mol)           # 210-element vector
morse = rdMolDescriptors.CalcMORSE(mol)         # 224-element vector
whim  = rdMolDescriptors.CalcWHIM(mol)          # 114-element vector
getaway = rdMolDescriptors.CalcGETAWAY(mol)     # 273-element vector
```
