# RDKit — Transformations, Reactions & Fragmentation

## Chemical Reactions via SMARTS

```python
from rdkit.Chem import AllChem, rdChemReactions
from rdkit import Chem

# Define reaction from SMARTS
rxn = AllChem.ReactionFromSmarts(
    '[C:1](=[O:2])-[OD1].[N!H0:3]>>[C:1](=[O:2])[N:3]'   # amide coupling
)

# Apply reaction to reactants
acid = Chem.MolFromSmiles('CC(=O)O')
amine = Chem.MolFromSmiles('NC')
products = rxn.RunReactants((acid, amine))

# products is a tuple of tuples: ((product1,), (product2,), ...)
unique_products = {}
for prods in products:
    for p in prods:
        smi = Chem.MolToSmiles(p)
        unique_products[smi] = p

print(list(unique_products.keys()))   # ['CNC(C)=O']
```

### Reaction from MDL RXN file

```python
rxn = AllChem.ReactionFromRxnFile('reaction.rxn')
```

### Validate and inspect reaction

```python
# Validate
nwarn, nerr, nagan, ndup = rxn.Validate()
if nerr > 0:
    raise ValueError(f"Reaction has {nerr} errors")

rxn.GetNumReactantTemplates()   # number of reactant patterns
rxn.GetNumProductTemplates()    # number of product patterns
rxn.GetNumAgentTemplates()      # catalysts / solvents

# Draw
from rdkit.Chem import Draw
d2d = Draw.MolDraw2DCairo(800, 300)
d2d.DrawReaction(rxn, highlightByReactant=True)
d2d.FinishDrawing()
png_data = d2d.GetDrawingText()
```

### Protect atoms from reacting

```python
# Protect amide nitrogens
amide_n = Chem.MolFromSmarts('[N;$(NC=[O,S])]')
mol_rw = Chem.RWMol(base_mol)
for match in base_mol.GetSubstructMatches(amide_n):
    base_mol.GetAtomWithIdx(match[0]).SetProp('_protected', '1')

# Run reaction — protected atoms are skipped
products = rxn.RunReactants((acid, base_mol))
```

### Reaction Fingerprints (compare reactions)

```python
from rdkit.Chem import rdChemReactions
from rdkit import DataStructs

rxn1 = rdChemReactions.ReactionFromSmarts('CCCO>>CCC=O')
rxn2 = rdChemReactions.ReactionFromSmarts('NCCO>>NCC=O')

fp1 = rdChemReactions.CreateDifferenceFingerprintForReaction(rxn1)
fp2 = rdChemReactions.CreateDifferenceFingerprintForReaction(rxn2)

sim = DataStructs.TanimotoSimilarity(fp1, fp2)
```

---

## Fragmentation

### Murcko Scaffolds

```python
from rdkit.Chem.Scaffolds import MurckoScaffold

# Extract Murcko scaffold (ring system + linkers, no side chains)
mol = Chem.MolFromSmiles('c1ccccc1C(=O)NC2CC2')
scaffold = MurckoScaffold.GetScaffoldForMol(mol)
print(Chem.MolToSmiles(scaffold))   # 'C1CC1Nc1ccccc1'

# Generic scaffold (all atoms become C, bonds become single)
generic = MurckoScaffold.MakeScaffoldGeneric(scaffold)

# Hash-based scaffold ID
from rdkit.Chem import rdMolHash
murcko_hash = rdMolHash.MolHash(mol, rdMolHash.HashFunction.MurckoScaffold)
```

### BRICS Decomposition

```python
from rdkit.Chem import BRICS

mol = Chem.MolFromSmiles('COc1ccc(CC(NC(=O)C2CCCO2)C(F)(F)F)cc1')

# Decompose into fragments
fragments = BRICS.BRICSDecompose(mol)
# Returns set of SMILES with attachment points: {'[3*]CC(...)...', ...}

# Process a whole library
all_frags = set()
with Chem.SDMolSupplier('library.sdf') as suppl:
    for mol in suppl:
        if mol is None:
            continue
        all_frags.update(BRICS.BRICSDecompose(mol))

# Recombine fragments to generate new molecules
fragmols = [Chem.MolFromSmiles(s) for s in sorted(all_frags)]
fragmols = [m for m in fragmols if m is not None]

builder = BRICS.BRICSBuild(fragmols)
new_mols = []
for i, m in enumerate(builder):
    if i >= 100:
        break
    m.UpdatePropertyCache(strict=False)
    Chem.SanitizeMol(m, catchErrors=True)
    new_mols.append(m)
```

### RECAP Decomposition (retrosynthetic)

```python
from rdkit.Chem import Recap

mol = Chem.MolFromSmiles('c1ccccc1OCCOC(=O)CC')
hierarchy = Recap.RecapDecompose(mol)

# Root molecule SMILES
print(hierarchy.smiles)

# Direct children (one bond cleavage)
for smi, node in hierarchy.children.items():
    print(smi)

# All leaf fragments (terminal decomposition)
leaves = hierarchy.GetLeaves()
for smi, node in leaves.items():
    print(smi, Chem.MolToSmiles(node.mol))
```

### Generic Bond Fragmentation

```python
# Break all ring-chain bonds
mol = Chem.MolFromSmiles('CC1CC(O)C1CCC1CC1')
ring_chain_bonds = mol.GetSubstructMatches(Chem.MolFromSmarts('[!R][R]'))
bond_ids = [mol.GetBondBetweenAtoms(a, b).GetIdx() for a, b in ring_chain_bonds]

fragmented = Chem.FragmentOnBonds(mol, bond_ids)
print(Chem.MolToSmiles(fragmented, True))
```

---

## Molecule Standardization

### MolStandardize (RDKit built-in)

```python
from rdkit.Chem.MolStandardize import rdMolStandardize

# Standard cleanup pipeline
def standardize(mol):
    """Standardize molecule: normalize, neutralize, remove fragments."""
    # Remove salts / largest fragment
    largest = rdMolStandardize.LargestFragmentChooser()
    mol = largest.choose(mol)

    # Normalize functional groups
    normalizer = rdMolStandardize.Normalizer()
    mol = normalizer.normalize(mol)

    # Neutralize charges
    uncharger = rdMolStandardize.Uncharger()
    mol = uncharger.uncharge(mol)

    # Reionize
    reionizer = rdMolStandardize.Reionizer()
    mol = reionizer.reionize(mol)

    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    return mol

mol = Chem.MolFromSmiles('CC(CNC[O-])[N+]([O-])=O')
mol_std = standardize(mol)
```

### Largest Fragment (salt stripping)

```python
from rdkit.Chem import rdmolops

# Get all disconnected fragments
mol = Chem.MolFromSmiles('c1ccccc1.CC(=O)O.[Na+].[Cl-]')
frags  = rdmolops.GetMolFrags(mol, asMols=True)
parent = max(frags, key=lambda m: m.GetNumAtoms())
print(Chem.MolToSmiles(parent))   # 'c1ccccc1'
```

### Charge Neutralization

```python
def neutralize_atoms(mol):
    """Neutralize charged atoms by adjusting H count."""
    pattern = Chem.MolFromSmarts(
        "[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]"
    )
    matches = [x[0] for x in mol.GetSubstructMatches(pattern)]
    rwmol = Chem.RWMol(mol)
    for idx in matches:
        atom = rwmol.GetAtomWithIdx(idx)
        chg  = atom.GetFormalCharge()
        nhyd = atom.GetTotalNumHs()
        atom.SetFormalCharge(0)
        atom.SetNumExplicitHs(nhyd - chg)
        atom.UpdatePropertyCache()
    return rwmol.GetMol()
```

---

## Tautomer Handling

```python
from rdkit.Chem.MolStandardize import rdMolStandardize

enumerator = rdMolStandardize.TautomerEnumerator()

mol = Chem.MolFromSmiles('Oc1ncc(OC(CC)C)cc1')

# Enumerate all tautomers
tautomers    = enumerator.Enumerate(mol)
tautomer_smi = [Chem.MolToSmiles(t) for t in tautomers]

# Get canonical (preferred) tautomer
canonical_taut = enumerator.Canonicalize(mol)
print(Chem.MolToSmiles(canonical_taut))

# Score tautomers (higher = more stable)
scored = [(Chem.MolToSmiles(t), enumerator.ScoreTautomer(t)) for t in tautomers]
scored.sort(key=lambda x: x[1], reverse=True)
```

---

## Substructure Replacement

```python
# Replace substructure with a new fragment
from rdkit.Chem import AllChem

mol = Chem.MolFromSmiles('c1ccccc1NC(=O)C')

# Replace amide with sulfonamide
repl = AllChem.ReplaceSubstructs(
    mol,
    Chem.MolFromSmarts('NC(=O)'),
    Chem.MolFromSmiles('NS(=O)(=O)'),
)
print([Chem.MolToSmiles(m) for m in repl])

# Delete a substructure
mol2 = AllChem.DeleteSubstructs(mol, Chem.MolFromSmarts('[CH3]'))

# Replace sidechains / core
core     = Chem.MolFromSmarts('c1ccccc1')
no_side  = AllChem.ReplaceSidechains(mol, core)   # keep only core
core_out = AllChem.ReplaceCore(mol, core)          # keep only sidechains
```

---

## MolHash — Canonical Identifiers

```python
from rdkit.Chem import rdMolHash

mol = Chem.MolFromSmiles('CC(C(C1=CC(=C(C=C1)O)O)O)N(C)C(=O)OCC2=CC=CC=C2')

# Available hash functions
hash_fns = rdMolHash.HashFunction.names   # dict of all available

# Common hashes
inchikey = rdMolHash.MolHash(mol, rdMolHash.HashFunction.InchiKey)
murcko   = rdMolHash.MolHash(mol, rdMolHash.HashFunction.MurckoScaffold)
net_charge = rdMolHash.MolHash(mol, rdMolHash.HashFunction.NetCharge)
regioisomer = rdMolHash.MolHash(mol, rdMolHash.HashFunction.Regioisomer)
anon_graph  = rdMolHash.MolHash(mol, rdMolHash.HashFunction.AnonymousGraph)

# Group molecules by scaffold
from collections import defaultdict
scaffold_groups = defaultdict(list)
for mol in molecules:
    scaffold_id = rdMolHash.MolHash(mol, rdMolHash.HashFunction.MurckoScaffold)
    scaffold_groups[scaffold_id].append(mol)
```

---

## SMILES Enumeration (Data Augmentation)

```python
from rdkit import Chem

mol = Chem.MolFromSmiles('CC(N)C(=O)O')

# Generate N random SMILES (for ML augmentation)
rand_smiles = Chem.MolToRandomSmilesVect(mol, 100, randomSeed=42)
unique_smiles = list(set(rand_smiles))
print(f"{len(unique_smiles)} unique SMILES representations")

# Custom random SMILES
smiles_set = set()
for i in range(200):
    s = Chem.MolToSmiles(mol, doRandom=True, canonical=False)
    smiles_set.add(s)
```
