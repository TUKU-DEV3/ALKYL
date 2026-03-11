# Selections, Styles, and Surfaces

## Selection Language

Selections are Python dicts passed to `setStyle`, `addSurface`, `zoomTo`, etc.

### Atom/Residue Properties

```python
# By chain
{'chain': 'A'}
{'chain': ['A', 'B']}   # multiple chains

# By residue number
{'resi': 42}
{'resi': [42, 43, 44]}          # list of residues
{'resi': '10-50'}               # range string

# By residue name
{'resn': 'LIG'}
{'resn': ['HIS', 'ASP', 'GLU']}  # multiple names

# By atom name
{'atom': 'CA'}
{'atom': ['N', 'CA', 'C', 'O']}  # backbone atoms

# By element
{'elem': 'N'}
{'elem': ['N', 'O']}

# By model index (when multiple models loaded)
{'model': 0}   # first model (0-indexed)
{'model': -1}  # last model
```

### Logical Operators

```python
# AND (dict merging = implicit AND)
{'chain': 'A', 'resi': 42}   # chain A AND residue 42

# OR (use list)
[{'chain': 'A'}, {'chain': 'B'}]   # chain A OR chain B
# Note: list syntax for OR is NOT standard in all py3Dmol versions
# More reliable: make separate setStyle calls
```

### Proximity Selection
```python
# Residues within N Å of selection
{'within': {'distance': 5.0, 'sel': {'resn': 'LIG'}}}

# NOT within (i.e., far from ligand)
{'byres': True, 'within': {'distance': 10, 'sel': {'resn': 'LIG'}}}
```

### Selecting Everything / Clearing
```python
# Select everything
view.setStyle({}, {'cartoon': {'color': 'grey'}})   # {} = all atoms

# Clear style for specific selection (make invisible)
view.setStyle({'resn': 'HOH'}, {})    # empty dict = remove style = invisible
```

## Color Schemes

### Built-in Carbon Color Schemes (for organic molecules)
```python
# Format: '{color}Carbon' — only carbons are colored, heteroatoms use standard CPK
'greenCarbon'    # green C, standard N/O/S/P
'cyanCarbon'     # cyan C
'magentaCarbon'  # magenta C
'yellowCarbon'   # yellow C
'whiteCarbon'    # white C (good on dark background)
'orangeCarbon'   # orange C
'purpleCarbon'   # purple C
'blueCarbon'     # blue C
```

### Standard Color Names
```python
'red', 'green', 'blue', 'cyan', 'magenta', 'yellow', 'orange',
'white', 'grey', 'lightgrey', 'black', 'purple', 'pink', 'brown'
```

### Hex Colors
```python
{'stick': {'color': '#FF6B6B', 'radius': 0.15}}   # custom red
```

### Property-Based Color Schemes
```python
# Spectrum (N→C direction)
{'cartoon': {'color': 'spectrum'}}

# By secondary structure
{'cartoon': {'color': 'ssPyMol'}}    # helix=red, sheet=yellow, coil=green
{'cartoon': {'color': 'ssJmol'}}     # helix=pink, sheet=yellow, coil=white

# By chain
{'cartoon': {'color': 'chain'}}

# By residue type (Jmol scheme: polar=green, nonpolar=white, etc.)
{'cartoon': {'color': 'residue'}}

# By B-factor / custom property
{
    'cartoon': {
        'colorscheme': {
            'prop': 'b',              # use B-factor column
            'gradient': 'rwb',        # red-white-blue
            'min': 0, 'max': 100
        }
    }
}
# gradient options: 'rwb', 'roygb', 'sinebow', custom min/max colors
```

## Style Options

### Cartoon
```python
{'cartoon': {
    'color': 'spectrum',
    'opacity': 0.8,        # 0.0–1.0
    'style': 'rectangle',  # 'oval', 'parabola', 'rectangle', 'edgedRectangle'
    'arrows': True,        # show strand arrows
    'tubes': True,         # show loop as tube
    'thickness': 0.4,      # ribbon thickness
    'width': 1.5,          # sheet width
}}
```

### Stick
```python
{'stick': {
    'radius': 0.15,                     # default 0.25
    'colorscheme': 'greenCarbon',
    'color': '#4488cc',                 # override all atoms
    'singleBonds': False,               # show bond orders
    'opacity': 0.9,
}}
```

### Sphere
```python
{'sphere': {
    'scale': 0.3,                       # fraction of VDW radius (default 1.0)
    'radius': 0.4,                      # absolute radius (overrides scale)
    'colorscheme': 'cyanCarbon',
    'opacity': 0.8,
}}
```

### Line
```python
{'line': {
    'linewidth': 1.0,
    'colorscheme': 'whiteCarbon',
    'opacity': 0.5,
}}
```

## Molecular Surfaces

```python
# Surface types (use py3Dmol constants)
py3Dmol.SES   # Solvent Excluded Surface (default, slow)
py3Dmol.SAS   # Solvent Accessible Surface
py3Dmol.VDW   # Van der Waals surface (fast)
py3Dmol.MS    # Molecular surface

# Basic surface
view.addSurface(py3Dmol.SES, {'opacity': 0.7, 'color': 'white'})

# Colored by electrostatics (B-factor proxy)
view.addSurface(py3Dmol.SES, {
    'opacity': 0.8,
    'colorscheme': {
        'prop': 'b',
        'gradient': 'rwb',
        'min': -20, 'max': 20
    }
})

# Surface for subset only
view.addSurface(
    py3Dmol.VDW,
    {'opacity': 0.5, 'color': 'lightblue'},
    {'within': {'distance': 8, 'sel': {'resn': 'LIG'}}}  # selection
)

# Remove surfaces
view.removeAllSurfaces()
view.removeSurface(surface_id)   # surface_id returned by addSurface()
```

## Labels

```python
# Atom label
view.addLabel('LIG', {
    'position': {'x': 0, 'y': 0, 'z': 0},  # or use selection
    'fontColor': 'white',
    'fontSize': 14,
    'backgroundColor': 'black',
    'backgroundOpacity': 0.5,
    'borderColor': 'grey',
    'borderThickness': 1,
})

# Residue labels (label all residues in selection)
view.addResLabels(
    {'chain': 'A', 'resi': [42, 45, 67]},
    {'fontColor': 'yellow', 'fontSize': 12, 'showBackground': True}
)

# Atom labels for ligand
view.addPropertyLabels(
    'atom',                          # property to label
    {'resn': 'LIG'},                 # selection
    {'fontColor': 'white', 'fontSize': 10}
)

# Remove labels
view.removeAllLabels()
```

## Bonds / Interaction Lines

```python
# Draw a dashed line (e.g., for H-bond visualization)
view.addCylinder({
    'start': {'x': x1, 'y': y1, 'z': z1},
    'end': {'x': x2, 'y': y2, 'z': z2},
    'color': 'yellow',
    'radius': 0.05,
    'dashed': True,
    'dashLength': 0.2,
    'gapLength': 0.15,
})

# Sphere at a point (e.g., pharmacophore feature)
view.addSphere({
    'center': {'x': cx, 'y': cy, 'z': cz},
    'radius': 1.5,
    'color': 'blue',
    'opacity': 0.4,
})
```

## Pharmacophore Feature Overlay

```python
FEATURE_COLORS = {
    'HBD': 'blue',
    'HBA': 'red',
    'AR': 'orange',
    'HYD': 'yellow',
    'POS': 'cyan',
    'NEG': 'magenta',
}

def add_pharmacophore_features(view, features):
    """
    features: list of dicts with keys: type, x, y, z, [radius]
    """
    for feat in features:
        color = FEATURE_COLORS.get(feat['type'], 'grey')
        r = feat.get('radius', 1.2)
        view.addSphere({
            'center': {'x': feat['x'], 'y': feat['y'], 'z': feat['z']},
            'radius': r, 'color': color, 'opacity': 0.4
        })
        view.addLabel(feat['type'], {
            'position': {'x': feat['x'], 'y': feat['y'], 'z': feat['z']},
            'fontColor': color, 'fontSize': 11, 'showBackground': False
        })
```
