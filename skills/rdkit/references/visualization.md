# RDKit — Visualization & Drawing

## 2D Coordinates

```python
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import rdMolDraw2D

mol = Chem.MolFromSmiles('CC(=O)Oc1ccccc1C(=O)O')   # aspirin

# Generate 2D layout (required before drawing/saving MOL blocks)
AllChem.Compute2DCoords(mol)
```

---

## Simple Drawing

```python
# PIL Image (Jupyter / file output)
img = Draw.MolToImage(mol, size=(400, 300))
img.save('mol.png')
img   # displays in Jupyter

# To SVG string
svg = Draw.MolToSVGString(mol, size=(400, 300))
with open('mol.svg', 'w') as f:
    f.write(svg)

# Multiple molecules in a grid
mols = [Chem.MolFromSmiles(s) for s in ['CC(=O)O', 'c1ccccc1', 'NC(=O)c1ccccc1']]
img  = Draw.MolsToGridImage(
    mols,
    molsPerRow=3,
    subImgSize=(300, 250),
    legends=['Acetic acid', 'Benzene', 'Benzamide']
)
img.save('grid.png')
```

---

## High-Quality Rendering (rdMolDraw2D)

### Cairo (PNG — production quality)

```python
d = rdMolDraw2D.MolDraw2DCairo(400, 300)
d.drawOptions().addStereoAnnotation = True
d.drawOptions().addAtomIndices       = False
d.DrawMolecule(mol)
d.FinishDrawing()

png_bytes = d.GetDrawingText()
with open('mol.png', 'wb') as f:
    f.write(png_bytes)
```

### SVG (scalable, web)

```python
d = rdMolDraw2D.MolDraw2DSVG(400, 300)
d.DrawMolecule(mol)
d.FinishDrawing()

svg_text = d.GetDrawingText()
```

### Drawing Options

```python
d = rdMolDraw2D.MolDraw2DCairo(400, 300)
opts = d.drawOptions()

opts.addAtomIndices        = True     # show atom indices
opts.addStereoAnnotation   = True     # show R/S, E/Z
opts.atomLabelFontSize     = 0.6      # relative font size
opts.bondLineWidth         = 1.5      # line width
opts.padding               = 0.1      # whitespace padding (fraction)
opts.useBWAtomPalette()               # grayscale
opts.useDefaultAtomPalette()          # color by element (default)
```

---

## Substructure Highlighting

### Highlight atoms and bonds

```python
mol = Chem.MolFromSmiles('CC(=O)Oc1ccccc1C(=O)O')
pattern = Chem.MolFromSmarts('C(=O)[OH]')
match   = mol.GetSubstructMatch(pattern)     # (7, 8, 9) or similar

d = rdMolDraw2D.MolDraw2DCairo(400, 300)
d.DrawMolecule(
    mol,
    highlightAtoms=list(match),
    highlightAtomColors={i: (0.9, 0.2, 0.2) for i in match},   # red
)
d.FinishDrawing()
```

### Multi-color highlighting

```python
from collections import defaultdict

mol    = Chem.MolFromSmiles('CC1=CC(=CC=C1)NC(=O)CCC2=CC=CC=C2')
colors = {
    'aromatic': (0.0, 0.0, 1.0, 0.3),   # blue, translucent
    'other':    (1.0, 0.0, 0.0, 0.3),   # red, translucent
}

atom_highlights = defaultdict(list)
atom_radii      = {}
for atom in mol.GetAtoms():
    aid = atom.GetIdx()
    color = colors['aromatic'] if atom.GetIsAromatic() else colors['other']
    atom_highlights[aid].append(color)
    atom_radii[aid] = 0.3

bond_highlights = defaultdict(list)
for bond in mol.GetBonds():
    bid = bond.GetIdx()
    a1, a2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
    color = colors['aromatic'] if bond.GetIsAromatic() else colors['other']
    bond_highlights[bid].append(color)

d = rdMolDraw2D.MolDraw2DCairo(400, 400)
d.DrawMoleculeWithHighlights(
    mol, "",
    dict(atom_highlights),
    dict(bond_highlights),
    atom_radii, {}
)
d.FinishDrawing()
```

---

## Atom Index Labels

```python
def mol_with_indices(mol):
    """Return copy of mol with atom map numbers set to atom index."""
    mol2 = Chem.RWMol(mol)
    for atom in mol2.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx())
    return mol2.GetMol()

mol_indexed = mol_with_indices(mol)
AllChem.Compute2DCoords(mol_indexed)
img = Draw.MolToImage(mol_indexed, size=(400, 300))
```

---

## Reaction Drawing

```python
rxn = AllChem.ReactionFromSmarts('[C:1](=[O:2])-[OD1].[N!H0:3]>>[C:1](=[O:2])[N:3]')

d = rdMolDraw2D.MolDraw2DCairo(800, 300)
d.DrawReaction(rxn)
d.FinishDrawing()
png = d.GetDrawingText()

# With reactant color coding
d2 = rdMolDraw2D.MolDraw2DCairo(800, 300)
d2.DrawReaction(rxn, highlightByReactant=True)
d2.FinishDrawing()
```

---

## Similarity Maps (atom contribution to similarity)

```python
from rdkit.Chem.Draw import SimilarityMaps

ref = Chem.MolFromSmiles('CCCN(CCCCN1CCN(c2ccccc2OC)CC1)Cc1ccc2ccccc2c1')
mol = Chem.MolFromSmiles('COc1cccc2cc(C(=O)NCCCCN3CCN(c4cccc5nccnc54)CC3)oc21')

d2d = rdMolDraw2D.MolDraw2DCairo(400, 400)
_, maxweight = SimilarityMaps.GetSimilarityMapForFingerprint(
    ref, mol,
    SimilarityMaps.GetMorganFingerprint,
    d2d
)
d2d.FinishDrawing()
png = d2d.GetDrawingText()
```

---

## Align Molecules to Template (2D)

```python
from rdkit.Chem import AllChem

template = Chem.MolFromSmiles('c1ccccc1C(=O)O')
AllChem.Compute2DCoords(template)

mols = [Chem.MolFromSmiles(s) for s in [
    'c1ccccc1C(=O)NC',
    'c1ccccc1C(=O)NCC(=O)O',
    'c1ccccc1C(=O)N',
]]

aligned = []
for mol in mols:
    AllChem.Compute2DCoords(mol)
    AllChem.GenerateDepictionMatching2DStructure(mol, template)
    aligned.append(mol)

img = Draw.MolsToGridImage(
    aligned,
    molsPerRow=3,
    subImgSize=(300, 250),
    legends=[Chem.MolToSmiles(m) for m in aligned]
)
```

---

## Jupyter / IPython Inline Display

```python
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw

IPythonConsole.ipython_useSVG = True   # SVG in notebooks

# Direct display
mol   # just output mol in Jupyter cell → image

# Grid of molecules
Draw.MolsToGridImage(mols, molsPerRow=4)
```

---

## CoordGen (Better 2D for Macrocycles)

```python
from rdkit.Chem import rdCoordGen

mol = Chem.MolFromSmiles('C1CCCCCCCCCCCCCCC1')   # large ring
rdCoordGen.AddCoords(mol)   # better ring layout than Compute2DCoords

img = Draw.MolToImage(mol, size=(400, 300))
```

---

## Drawing to File (SVG / PNG / PDF)

```python
from rdkit.Chem import Draw

# Single molecule
Draw.MolToFile(mol, 'mol.png', size=(400, 300))
Draw.MolToFile(mol, 'mol.svg', size=(400, 300))

# Grid
Draw.MolsToGridImage(mols, molsPerRow=4).save('grid.png')
```
