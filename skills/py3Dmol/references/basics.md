# py3Dmol Basics

## Installation
```bash
pip install py3Dmol
# For Jupyter Lab (optional widget extension):
pip install ipywidgets
jupyter labextension install jupyterlab_3dmol  # older JupyterLab only
```

## View Creation

```python
import py3Dmol

# Basic view
view = py3Dmol.view(width=800, height=500)

# With linked views (synchronize rotation)
view1, view2 = py3Dmol.view(width=400, height=400), py3Dmol.view(width=400, height=400)

# Inline Jupyter display
view.show()   # or just `view` on the last line of a cell
```

## Loading Structures

### From PDB file
```python
view = py3Dmol.view(width=800, height=500)
with open('protein.pdb') as f:
    pdb_str = f.read()
view.addModel(pdb_str, 'pdb')
view.setStyle({'cartoon': {'color': 'spectrum'}})
view.zoomTo()
view.show()
```

### From SDF file (ligand)
```python
with open('ligand.sdf') as f:
    sdf_str = f.read()
view.addModel(sdf_str, 'sdf')
view.setStyle({}, {'stick': {'colorscheme': 'cyanCarbon'}})
view.zoomTo()
```

### From SMILES via RDKit (generate 3D first)
```python
from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol

def smiles_to_3d_view(smiles, width=600, height=400):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    AllChem.MMFFOptimizeMolecule(mol)
    mol = Chem.RemoveHs(mol)
    sdf = Chem.MolToMolBlock(mol)

    view = py3Dmol.view(width=width, height=height)
    view.addModel(sdf, 'mol')
    view.setStyle({}, {'stick': {'colorscheme': 'cyanCarbon', 'radius': 0.15},
                       'sphere': {'colorscheme': 'cyanCarbon', 'radius': 0.3}})
    view.zoomTo()
    return view
```

### From PDB ID (fetch online)
```python
view = py3Dmol.view(query='pdb:1HSG', width=800, height=500)
view.setStyle({'cartoon': {'color': 'spectrum'}})
view.show()
```

### From URL
```python
view = py3Dmol.view(width=800, height=500)
view.setBackgroundColor('white')
# Load from URL directly
view.addModelFromUrl('https://files.rcsb.org/download/1HSG.pdb', 'pdb')
```

## Basic Style Types

| Style | Code | Use for |
|-------|------|---------|
| Cartoon | `{'cartoon': {'color': 'spectrum'}}` | Protein secondary structure |
| Stick | `{'stick': {'radius': 0.15}}` | Small molecules, ligands |
| Sphere | `{'sphere': {'scale': 0.3}}` | VDW representation |
| Line | `{'line': {}}` | Solvent, large systems |
| Cross | `{'cross': {'radius': 0.1}}` | Ligand minimalistic |

```python
# Cartoon: color options
view.setStyle({'cartoon': {'color': 'spectrum'}})        # rainbow N→C
view.setStyle({'cartoon': {'color': 'chain'}})           # by chain
view.setStyle({'cartoon': {'color': 'residue'}})         # by residue type
view.setStyle({'cartoon': {'color': '#4488cc'}})         # custom hex

# Stick + sphere combined (ball-and-stick)
view.setStyle({}, {
    'stick': {'radius': 0.15, 'colorscheme': 'cyanCarbon'},
    'sphere': {'radius': 0.3, 'colorscheme': 'cyanCarbon'}
})

# Lines for solvent
view.setStyle({'resn': 'HOH'}, {'line': {}})
```

## Background & Camera

```python
view.setBackgroundColor('white')           # white, black, or hex '#1a1a2e'
view.setBackgroundColor('0xffffff')        # hex with 0x prefix also works

view.zoomTo()                              # fit all models
view.zoomTo({'resn': 'LIG'})              # fit to ligand
view.rotate(90, 'y')                       # rotate 90° around y-axis
view.zoom(1.5)                             # zoom in
view.center({'chain': 'A'})               # center on chain A
```

## Multiple Models

```python
view = py3Dmol.view(width=800, height=500)

# Add protein
view.addModel(protein_pdb, 'pdb')
view.setStyle({'model': 0}, {'cartoon': {'color': 'grey', 'opacity': 0.8}})

# Add ligand
view.addModel(ligand_sdf, 'sdf')
view.setStyle({'model': 1}, {'stick': {'colorscheme': 'greenCarbon'}})

view.zoomTo({'model': 1})  # zoom to ligand
view.show()
```

## Removing Hydrogens from Display

```python
# Hide H atoms (keep in model, just don't style them)
view.setStyle({'elem': 'H'}, {})   # empty style = invisible
```

## Quick Protein View Function

```python
def show_protein(pdb_path, style='cartoon', color='spectrum', width=800, height=500):
    """One-liner protein viewer."""
    view = py3Dmol.view(width=width, height=height)
    view.addModel(open(pdb_path).read(), 'pdb')
    view.setStyle({style: {'color': color}})
    view.setBackgroundColor('white')
    view.zoomTo()
    return view

# Usage
show_protein('protein.pdb').show()
```
