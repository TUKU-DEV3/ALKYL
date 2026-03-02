# RDKit — Substructure Search, MCS, Rings & Stereochemistry

## Substructure Search (SMARTS)

```python
from rdkit import Chem

mol = Chem.MolFromSmiles('CC(=O)Oc1ccccc1C(=O)O')   # aspirin

# Query from SMARTS
pattern = Chem.MolFromSmarts('[C:1](=O)[OH]')    # carboxylic acid
query   = Chem.MolFromSmarts('c1ccccc1')          # benzene ring

# Check presence
mol.HasSubstructMatch(query)             # True

# First match (returns atom indices)
match = mol.GetSubstructMatch(query)     # (3, 4, 5, 6, 7, 8)

# All matches
all_matches = mol.GetSubstructMatches(query)
# → ((3,4,5,6,7,8), (3,8,7,6,5,4))  [both ring traversal directions]

# With chirality / stereo checking
mol.HasSubstructMatch(pattern, useChirality=True)
mol.GetSubstructMatches(pattern, useChirality=True, useQueryQueryMatches=False)
```

### Common SMARTS Patterns for Drug Discovery

```python
SMARTS_PATTERNS = {
    # Functional groups
    'carboxylic_acid':   '[C;!$(C=O)](=O)[OH]',
    'ester':             '[#6][C](=O)O[#6]',
    'amide':             'C(=O)[NX3;H2,H1;!$(NC=O)]',
    'primary_amine':     '[NX3;H2;!$(NC=O)]',
    'secondary_amine':   '[NX3;H1;!$(NC=O)]',
    'tertiary_amine':    '[NX3;H0;!$(NC=O);!$([n])]',
    'hydroxyl':          '[OX2H]',
    'thiol':             '[SX2H]',
    'sulfonamide':       'S(=O)(=O)[NX3;H2,H1]',
    'phosphate':         '[PX4](=O)([OX2H])[OX2H]',

    # Rings
    'aromatic_ring':     'c1ccccc1',
    'pyridine':          'n1ccccc1',
    'piperidine':        '[NH]1CCCCC1',
    'piperazine':        'N1CCNCC1',

    # Reactive groups / alerts
    'michael_acceptor':  'C=CC(=O)',
    'aldehyde':          '[CX3H1](=O)[#6]',
    'acyl_halide':       'C(=O)[F,Cl,Br,I]',
    'nitro':             '[$([NX3](=O)=O),$([NX3+](=O)[O-])]',

    # HBD/HBA
    'hbd':               '[#7,#8;!H0]',
    'hba':               '[#7,#8]',
}

def find_groups(mol, patterns: dict) -> dict:
    return {name: bool(mol.HasSubstructMatch(Chem.MolFromSmarts(sma)))
            for name, sma in patterns.items()}
```

### Highlight Substructure (for visualization)

```python
from rdkit.Chem.Draw import rdMolDraw2D

def highlight_smarts(mol, smarts: str, size=(400, 300)) -> bytes:
    pattern = Chem.MolFromSmarts(smarts)
    match   = mol.GetSubstructMatch(pattern)
    d = rdMolDraw2D.MolDraw2DCairo(*size)
    d.drawOptions().addAtomIndices = False
    d.DrawMolecule(mol, highlightAtoms=list(match))
    d.FinishDrawing()
    return d.GetDrawingText()
```

---

## Maximum Common Substructure (MCS)

```python
from rdkit.Chem import rdFMCS

mols = [
    Chem.MolFromSmiles('O=C(NCc1cc(OC)c(O)cc1)CCCC/C=C/C(C)C'),
    Chem.MolFromSmiles('CC(C)CCCCCC(=O)NCC1=CC(=C(C=C1)O)OC'),
    Chem.MolFromSmiles('c1(C=O)cc(OC)c(O)cc1'),
]

# Basic MCS
result = rdFMCS.FindMCS(mols)
print(f"MCS atoms: {result.numAtoms}, bonds: {result.numBonds}")
print(f"MCS SMARTS: {result.smartsString}")
print(f"Timed out: {result.canceled}")

mcs_mol = Chem.MolFromSmarts(result.smartsString)
```

### MCS Options

```python
# Compare by element type (default)
result = rdFMCS.FindMCS(mols, atomCompare=rdFMCS.AtomCompare.CompareElements)

# Any atom matches any atom (topology only)
result = rdFMCS.FindMCS(mols, atomCompare=rdFMCS.AtomCompare.CompareAny)

# Compare by isotope labels
result = rdFMCS.FindMCS(mols, atomCompare=rdFMCS.AtomCompare.CompareIsotopes)

# Bond comparison
result = rdFMCS.FindMCS(mols, bondCompare=rdFMCS.BondCompare.CompareAny)
result = rdFMCS.FindMCS(mols, bondCompare=rdFMCS.BondCompare.CompareOrder)      # default
result = rdFMCS.FindMCS(mols, bondCompare=rdFMCS.BondCompare.CompareOrderExact)

# Ring constraints
result = rdFMCS.FindMCS(mols, ringMatchesRingOnly=True)
result = rdFMCS.FindMCS(mols, completeRingsOnly=True)

# Valence matching
result = rdFMCS.FindMCS(mols, matchValences=True)

# Maximize bonds instead of atoms
result = rdFMCS.FindMCS(mols, maximizeBonds=True)

# Timeout
result = rdFMCS.FindMCS(mols, timeout=5)   # stop after 5 seconds
```

### MCS Visualization

```python
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import AllChem

result = rdFMCS.FindMCS(mols[:2])
mcs_smarts = result.smartsString

# Highlight MCS in each molecule
mcs_query = Chem.MolFromSmarts(mcs_smarts)
for mol in mols[:2]:
    match = mol.GetSubstructMatch(mcs_query)
    # ... use match indices to highlight
```

---

## Ring Analysis

```python
mol = Chem.MolFromSmiles('OC1C2C1CC2')

# Check if atom is in ring
mol.GetAtomWithIdx(0).IsInRing()           # False (O)
mol.GetAtomWithIdx(1).IsInRing()           # True

# Ring size membership
mol.GetAtomWithIdx(2).IsInRingSize(3)      # True
mol.GetAtomWithIdx(2).IsInRingSize(5)      # False

# SSSR (Smallest Set of Smallest Rings)
ssr = Chem.GetSSSR(mol)               # number of SSSR rings
sssr_list = list(Chem.GetSymmSSSR(mol))   # symmetrized rings
for ring in sssr_list:
    print(list(ring))  # atom indices

# RingInfo — efficient ring queries
ri = mol.GetRingInfo()
ri.NumAtomRings(2)              # how many rings atom 2 belongs to
ri.NumBondRings(1)              # how many rings bond 1 belongs to
ri.IsAtomInRingOfSize(2, 3)    # True: atom 2 is in a 3-membered ring
ri.IsBondInRingOfSize(1, 3)    # True: bond 1 is in a 3-membered ring
ri.AtomRingSizes(2)            # sizes of all rings containing atom 2
ri.BondRingSize(1)             # size of ring containing bond 1
```

### Ring Systems

```python
def get_ring_systems(mol) -> list[frozenset]:
    """Return each fused ring system as a frozenset of atom indices."""
    ri = mol.GetRingInfo()
    ring_atom_sets = [frozenset(r) for r in ri.AtomRings()]
    systems = []
    for ring in ring_atom_sets:
        merged = False
        for i, sys in enumerate(systems):
            if sys & ring:
                systems[i] = sys | ring
                merged = True
                break
        if not merged:
            systems.append(ring)
    return systems
```

---

## Stereochemistry

### Chiral Centers (R/S)

```python
from rdkit import Chem

mol = Chem.MolFromSmiles('CC[C@H](F)Cl')

# Assign CIP stereochemistry (R/S labels)
Chem.AssignStereochemistry(mol, cleanIt=True, force=True)

for atom in mol.GetAtoms():
    if atom.HasProp('_CIPCode'):
        print(f"Atom {atom.GetIdx()} ({atom.GetSymbol()}): {atom.GetProp('_CIPCode')}")
    tag = atom.GetChiralTag()
    # rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW
    # rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW
    # rdkit.Chem.rdchem.ChiralType.CHI_UNSPECIFIED

# Count stereocenters
from rdkit.Chem import rdMolDescriptors
n_centers = rdMolDescriptors.CalcNumAtomStereoCenters(mol)
n_unspec   = rdMolDescriptors.CalcNumUnspecifiedAtomStereoCenters(mol)
```

### E/Z Double Bond Stereo

```python
mol = Chem.MolFromSmiles(r'Cl/C=C/F')   # E-isomer
Chem.AssignStereochemistry(mol, cleanIt=True, force=True)

for bond in mol.GetBonds():
    stereo = bond.GetStereo()
    # rdkit.Chem.rdchem.BondStereo.STEREOANY
    # rdkit.Chem.rdchem.BondStereo.STEREOE
    # rdkit.Chem.rdchem.BondStereo.STEREOZ
    # rdkit.Chem.rdchem.BondStereo.STEREONONE
    if stereo not in (Chem.rdchem.BondStereo.STEREONONE,
                      Chem.rdchem.BondStereo.STEREOANY):
        print(f"Bond {bond.GetIdx()}: {stereo}")

# Substructure matching with chirality
m = Chem.MolFromSmiles('CC[C@H](F)Cl')
m.HasSubstructMatch(Chem.MolFromSmiles('C[C@H](F)Cl'), useChirality=True)   # True
m.HasSubstructMatch(Chem.MolFromSmiles('C[C@@H](F)Cl'), useChirality=True)  # False
```

### Enumerate Stereocenters

```python
from rdkit.Chem import EnumerateStereoisomers

mol = Chem.MolFromSmiles('CC(F)C(Cl)CC')   # 2 unspecified stereocenters

opts = EnumerateStereoisomers.StereoEnumerationOptions(
    unique=True,        # deduplicate
    onlyUnassigned=True # only vary unspecified centers
)
isomers = list(EnumerateStereoisomers.EnumerateStereoisomers(mol, options=opts))
print(f"{len(isomers)} stereoisomers")
for iso in isomers:
    print(Chem.MolToSmiles(iso))
```

---

## Atom and Bond Properties

```python
mol = Chem.MolFromSmiles('c1cccnc1')   # pyridine

for atom in mol.GetAtoms():
    print(
        atom.GetIdx(),
        atom.GetSymbol(),          # 'C' or 'N'
        atom.GetAtomicNum(),       # 6 or 7
        atom.GetIsAromatic(),      # True
        atom.GetDegree(),          # number of bonds (explicit)
        atom.GetTotalValence(),    # explicit + implicit
        atom.GetFormalCharge(),    # ionic charge
        atom.GetNumImplicitHs(),   # implicit H count
        atom.GetHybridization(),   # SP2, SP3, etc.
        atom.IsInRing(),           # ring membership
    )

for bond in mol.GetBonds():
    print(
        bond.GetIdx(),
        bond.GetBeginAtomIdx(),
        bond.GetEndAtomIdx(),
        bond.GetBondTypeAsDouble(),  # 1.0, 2.0, 1.5 (aromatic)
        bond.GetIsAromatic(),
        bond.IsInRing(),
    )

# Get bond between specific atoms
bond = mol.GetBondBetweenAtoms(0, 1)
```

---

## Chemical Features & Pharmacophores

```python
from rdkit.Chem import ChemicalFeatures
from rdkit import RDConfig
import os

# Load default feature definitions
fdef_path = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
factory   = ChemicalFeatures.BuildFeatureFactory(fdef_path)

mol   = Chem.MolFromSmiles('OCc1ccccc1CN')
feats = factory.GetFeaturesForMol(mol)

for feat in feats:
    print(feat.GetFamily(),      # 'Donor', 'Acceptor', 'Aromatic', 'Hydrophobe'
          feat.GetType(),        # e.g., 'SingleAtomDonor'
          feat.GetAtomIds())     # tuple of atom indices

# Feature families available in BaseFeatures.fdef
# Donor, Acceptor, NegIonizable, PosIonizable, ZnBinder, Aromatic, Hydrophobe, LumpedHydrophobe
```

### 2D Pharmacophore Fingerprints

```python
from rdkit.Chem.Pharm2D.SigFactory import SigFactory
from rdkit.Chem.Pharm2D import Generate

feat_factory = ChemicalFeatures.BuildFeatureFactory(fdef_path)
sig_factory  = SigFactory(feat_factory, minPointCount=2, maxPointCount=3)
sig_factory.SetBins([(0, 2), (2, 5), (5, 8)])
sig_factory.Init()

mol = Chem.MolFromSmiles('OCC(=O)CCCN')
fp  = Generate.Gen2DFingerprint(mol, sig_factory)

print(fp.GetNumOnBits())
print(sig_factory.GetBitDescription(list(fp.GetOnBits())[0]))
```
