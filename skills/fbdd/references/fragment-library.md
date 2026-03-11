# Fragment Library Design

## Property Filters (Rule of 3 + extras)

```python
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Crippen
import numpy as np

def filter_fragment_library(smiles_list, strict=True):
    """
    Apply Rule of 3 + quality filters to a SMILES list.
    Returns DataFrame with pass/fail and property values.
    """
    import pandas as pd
    from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

    # PAINS catalog
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    pains_catalog = FilterCatalog(params)

    records = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue

        props = {
            'smiles': smi,
            'MW': Descriptors.ExactMolWt(mol),
            'HAC': mol.GetNumHeavyAtoms(),
            'cLogP': Crippen.MolLogP(mol),
            'HBD': rdMolDescriptors.CalcNumHBD(mol),
            'HBA': rdMolDescriptors.CalcNumHBA(mol),
            'RotBonds': rdMolDescriptors.CalcNumRotatableBonds(mol),
            'PSA': rdMolDescriptors.CalcTPSA(mol),
            'Rings': rdMolDescriptors.CalcNumRings(mol),
            'ArRings': rdMolDescriptors.CalcNumAromaticRings(mol),
        }

        # Sp3 character (fraction of sp3 carbons)
        props['Fsp3'] = rdMolDescriptors.CalcFractionCSP3(mol)

        # PAINS check
        props['PAINS'] = bool(pains_catalog.HasMatch(mol))

        # Rule of 3
        ro3 = (props['MW'] <= 300 and props['cLogP'] <= 3 and
               props['HBD'] <= 3 and props['HBA'] <= 3)
        props['Ro3'] = ro3

        # Extended fragment quality
        quality_ok = (
            ro3 and
            props['HAC'] >= 8 and props['HAC'] <= 20 and
            props['RotBonds'] <= 3 and
            props['PSA'] <= 60 and
            not props['PAINS']
        )

        if strict:
            # Additional quality filters
            quality_ok = quality_ok and _quality_checks(mol, props)

        props['pass'] = quality_ok
        records.append(props)

    return pd.DataFrame(records)

def _quality_checks(mol, props):
    """Additional quality checks beyond Ro3."""
    from rdkit.Chem import rdMolDescriptors

    # No reactive groups (check common problematic patterns)
    REACTIVE_SMARTS = [
        '[F,Cl,Br,I][C;!$(C=O)]',   # alkyl halides (electrophilic)
        'C(=O)Cl',                    # acid chlorides
        '[N+]=[N-]',                  # diazo
        'N=[N+]=[N-]',               # azide
        'C(=O)OC(=O)',               # anhydride
        '[SH]',                       # thiol (promiscuous)
        'C=CC=O',                     # Michael acceptor
        'c1cc[nH]cc1',               # pyrrole (too reactive)
    ]
    for smarts in REACTIVE_SMARTS:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(smarts)):
            return False

    # Not purely flat (some 3D character preferred)
    # Fsp3 == 0 means entirely aromatic (flat — often promiscuous)
    if props['Fsp3'] == 0 and props['ArRings'] >= 2:
        return False   # fused polycyclic aromatic → promiscuous

    # Molecular complexity: not too many stereocenters for a fragment
    n_stereo = rdMolDescriptors.CalcNumAtomStereoCenters(mol)
    if n_stereo > 1:
        return False  # fragments should be ≤ 1 stereocenter

    return True
```

## sp3 Character and 3D Fragments

Fragments with Fsp3 > 0.2 are preferred — 3D fragments access deeper pockets
and show better selectivity (Lovering 2009 "Escape from Flatland").

```python
# Flag flat fragments
def classify_dimensionality(mol):
    fsp3 = rdMolDescriptors.CalcFractionCSP3(mol)
    ar_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    if fsp3 == 0:
        return '2D (fully flat)'
    elif fsp3 < 0.2:
        return '2.5D (mostly flat)'
    else:
        return '3D (sp3-enriched)'
```

## Fragment Diversity Selection

```python
from rdkit.Chem import AllChem
import numpy as np

def select_diverse_fragments(df, n=500, fingerprint='morgan'):
    """MaxMin diversity selection from filtered fragment library."""
    import subprocess, sys

    # Write to SMI file
    df[df['pass']]['smiles'].to_csv('/tmp/frags_filtered.smi',
                                     index=False, header=False)

    # Use ALKYL chem_diversity.py
    result = subprocess.run(
        [sys.executable,
         '/path/to/ALKYL/scripts/chem_diversity.py',
         '/tmp/frags_filtered.smi',
         '--n', str(n),
         '--fingerprint', fingerprint,
         '--out', '/tmp/frags_diverse.smi'],
        capture_output=True, text=True
    )
    print(result.stdout)
    return pd.read_csv('/tmp/frags_diverse.smi', header=None, names=['smiles'])
```

## Scaffold Analysis

Ensure the library covers diverse scaffolds (Murcko):

```python
from rdkit.Chem.Scaffolds import MurckoScaffold
from collections import Counter

def analyze_scaffolds(smiles_list):
    scaffolds = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            scaf = MurckoScaffold.GetScaffoldForMol(mol)
            scaffolds.append(Chem.MolToSmiles(scaf))
    counter = Counter(scaffolds)
    n_unique = len(counter)
    n_singletons = sum(1 for v in counter.values() if v == 1)
    print(f"Unique scaffolds: {n_unique}/{len(smiles_list)}")
    print(f"Scaffold singletons: {n_singletons} "
          f"({100*n_singletons/n_unique:.0f}%)")
    return counter
```

## Commercial Fragment Sources

| Source | Library size | Price range |
|--------|-------------|-------------|
| Enamine REAL Fragments | ~200k | €0.50–2/cpd |
| Maybridge HitFinder | 14,400 | $$$ |
| Sigma-Aldrich Discovery Series | 10,240 | $$ |
| Life Chemicals | ~30k | $ |
| Asinex | ~20k | $ |
| Prestwick Chemical | ~1280 (known drugs, fragments) | $$$ |

**Recommended starting point**: Enamine REAL Fragments with Ro3 pre-filter
→ download SMILES → run filter_fragment_library → MaxMin select 500–1000.

## Fragment Library Quality Metrics

```python
def library_qc_report(df):
    """Summary statistics for fragment library."""
    passed = df[df['pass']]
    print(f"Total: {len(df)}, Passed: {len(passed)} "
          f"({100*len(passed)/len(df):.1f}%)")
    print(f"\nProperty distributions (passed only):")
    for col in ['MW', 'HAC', 'cLogP', 'HBD', 'HBA', 'RotBonds', 'PSA', 'Fsp3']:
        print(f"  {col:12s}: mean={passed[col].mean():.2f}, "
              f"median={passed[col].median():.2f}, "
              f"range=[{passed[col].min():.1f}, {passed[col].max():.1f}]")
    print(f"\n3D classification:")
    for cat, count in passed.apply(
        lambda r: classify_dimensionality(Chem.MolFromSmiles(r.smiles)),
        axis=1
    ).value_counts().items():
        print(f"  {cat}: {count}")
```

## Fragment Library for NMR Screening

NMR (ligand-observed: STD, WaterLOGSY) requires:
- Solubility: ≥ 1 mM in PBS + 5% DMSO
- No proton-less fragments (need NMR signal)
- No paramagnetic groups
- Deuterated solvent compatibility

Additional filters for NMR:
```python
# Exclude quaternary ammonium, metals, excessive fluorination
NMR_EXCLUSIONS = [
    '[F][F]',        # gem-difluoro (suppresses 1H signal nearby)
    '[#7+]',         # quaternary N (charged, often aggregates)
    '[B,Si,Ge,As,Se,Te,Sn,Pb]',  # unusual atoms
]
```

## Solubility Prediction Filter

```python
# Quick logS estimate via ESOL (Delaney equation)
def esol_logs(mol):
    """ESOL: estimated aqueous solubility (log mol/L)."""
    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    rb = rdMolDescriptors.CalcNumRotatableBonds(mol)
    ap = sum(1 for a in mol.GetAromaticAtoms()) / mol.GetNumHeavyAtoms()
    return 0.16 - 0.63*logp - 0.0062*mw + 0.066*rb - 0.74*ap

# Fragments with logS < -4 (< 0.1 mM) unlikely soluble at screening conc
df['logS_ESOL'] = df['smiles'].apply(
    lambda s: esol_logs(Chem.MolFromSmiles(s)) if Chem.MolFromSmiles(s) else None
)
df['soluble'] = df['logS_ESOL'] > -4  # > 0.1 mM
```
