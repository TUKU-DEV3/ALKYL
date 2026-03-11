# Jupyter Patterns, Saving, and NGLview

## Jupyter Notebook Integration

```python
import py3Dmol

# In a Jupyter cell — just return the view object or call .show()
view = py3Dmol.view(width=800, height=500)
view.addModel(pdb_str, 'pdb')
view.setStyle({'cartoon': {}})
view.zoomTo()
view   # ← display inline (or view.show())
```

## Saving as PNG

```python
# Get base64-encoded PNG
png_data = view.png()  # returns 'data:image/png;base64,...'

# Save to file
import base64
header, data = png_data.split(',', 1)
with open('structure.png', 'wb') as f:
    f.write(base64.b64decode(data))

# Or use PIL
from PIL import Image
import io
img_bytes = base64.b64decode(png_data.split(',')[1])
img = Image.open(io.BytesIO(img_bytes))
img.save('structure.png', dpi=(300, 300))
```

## Saving as HTML (interactive, shareable)

```python
# Write standalone interactive HTML
html_str = view._make_html()
with open('structure.html', 'w') as f:
    f.write(html_str)
```

## ipywidgets: Pose Slider

```python
import ipywidgets as widgets
from IPython.display import display
import py3Dmol
from rdkit import Chem

def interactive_pose_viewer(receptor_pdb, poses_sdf):
    """Slider to flip through docking poses."""
    supplier = list(Chem.SDMolSupplier(poses_sdf, removeHs=True))
    n_poses = len(supplier)

    view = py3Dmol.view(width=850, height=550)

    def update(pose_idx):
        view.removeAllModels()
        view.removeAllSurfaces()

        # Receptor
        view.addModel(open(receptor_pdb).read(), 'pdb')
        view.setStyle({'model': 0}, {'cartoon': {'color': 'lightgrey', 'opacity': 0.7}})
        view.setStyle({'model': 0, 'resn': 'HOH'}, {})

        # Selected pose
        mol = supplier[pose_idx]
        if mol:
            sdf = Chem.MolToMolBlock(mol)
            view.addModel(sdf, 'sdf')
            view.setStyle({'model': 1},
                          {'stick': {'colorscheme': 'greenCarbon', 'radius': 0.2},
                           'sphere': {'colorscheme': 'greenCarbon', 'radius': 0.35}})
            # Get score from mol props
            score = mol.GetProp('minimizedAffinity') if mol.HasProp('minimizedAffinity') else '?'
            view.addLabel(f'Pose {pose_idx+1} | Score: {score}',
                          {'position': {'x': 0, 'y': 0, 'z': 0},
                           'fontColor': 'white', 'fontSize': 12,
                           'backgroundColor': 'black', 'backgroundOpacity': 0.5})
        view.zoomTo({'model': 1})
        view.render()

    slider = widgets.IntSlider(min=0, max=n_poses-1, step=1, value=0,
                                description='Pose:')
    widgets.interact(update, pose_idx=slider)
    display(view.show())

interactive_pose_viewer('receptor.pdb', 'poses.sdf')
```

## Animation (Conformer Loop)

```python
import py3Dmol
from rdkit import Chem

def animate_conformers(mol_or_sdf, interval=200):
    """
    Animate conformer ensemble.
    mol_or_sdf: RDKit Mol with multiple conformers, or SDF path.
    interval: ms between frames.
    """
    if isinstance(mol_or_sdf, str):
        mol = Chem.SDMolSupplier(mol_or_sdf, removeHs=False)[0]
    else:
        mol = mol_or_sdf

    view = py3Dmol.view(width=600, height=450)
    view.setBackgroundColor('white')

    # Add all conformers as models
    for conf_id in range(mol.GetNumConformers()):
        sdf = Chem.MolToMolBlock(mol, confId=conf_id)
        view.addModel(sdf, 'sdf')

    view.setStyle({}, {
        'stick': {'colorscheme': 'cyanCarbon', 'radius': 0.15},
        'sphere': {'colorscheme': 'cyanCarbon', 'radius': 0.28}
    })
    view.animate({'loop': 'forward', 'interval': interval})
    view.zoomTo()
    return view
```

## Batch Grid View (multiple molecules)

```python
import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem
import math

def grid_view(smiles_list, labels=None, cols=3, width=300, height=250):
    """
    Display molecules in a grid.
    Each molecule gets its own view embedded in HTML.
    """
    from IPython.display import HTML
    import base64

    rows = math.ceil(len(smiles_list) / cols)
    html_parts = ['<div style="display:flex;flex-wrap:wrap">']

    for i, smi in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        AllChem.MMFFOptimizeMolecule(mol)
        mol = Chem.RemoveHs(mol)
        sdf = Chem.MolToMolBlock(mol)

        view = py3Dmol.view(width=width, height=height)
        view.addModel(sdf, 'mol')
        view.setStyle({}, {'stick': {'colorscheme': 'cyanCarbon', 'radius': 0.12},
                           'sphere': {'colorscheme': 'cyanCarbon', 'radius': 0.25}})
        view.setBackgroundColor('white')
        view.zoomTo()

        label = labels[i] if labels and i < len(labels) else f'#{i+1}'
        html_parts.append(
            f'<div style="text-align:center;margin:5px">'
            f'<div>{label}</div>'
            f'{view._make_html()}'
            f'</div>'
        )

    html_parts.append('</div>')
    return HTML('\n'.join(html_parts))

# Usage
grid_view(['CCO', 'c1ccccc1', 'CC(=O)O'],
          labels=['Ethanol', 'Benzene', 'Acetic acid'])
```

## NGLview (Alternative — Better for MD Trajectories)

```python
# Install: pip install nglview; jupyter nbextension enable nglview --py --sys-prefix
import nglview as nv
import MDAnalysis as mda

# Show PDB
view = nv.show_file('protein.pdb')
view

# Show MDAnalysis Universe (trajectory-aware)
u = mda.Universe('topology.prmtop', 'traj.dcd')
view = nv.show_mdanalysis(u)
view.add_cartoon(selection='protein', color_scheme='residueindex')
view.add_licorice(selection='resname LIG', color_scheme='element')
view

# NGLview vs py3Dmol:
# py3Dmol: simpler API, pure HTML/JS, easier saving, good for static poses
# nglview: better trajectory support, smoother animation, more selection ops
```

## RDKit → py3Dmol Integration

```python
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import py3Dmol

def rdkit_mol_to_view(mol, remove_hs=True, width=600, height=400,
                       colorscheme='cyanCarbon'):
    """Convert RDKit Mol (with 3D conformer) to py3Dmol view."""
    if remove_hs:
        mol = Chem.RemoveHs(mol)
    sdf = Chem.MolToMolBlock(mol)
    view = py3Dmol.view(width=width, height=height)
    view.addModel(sdf, 'mol')
    view.setStyle({}, {
        'stick': {'colorscheme': colorscheme, 'radius': 0.15},
        'sphere': {'colorscheme': colorscheme, 'radius': 0.28}
    })
    view.setBackgroundColor('white')
    view.zoomTo()
    return view

# Load SDF and display all molecules
def view_sdf(sdf_path, n_max=10):
    supplier = list(Chem.SDMolSupplier(sdf_path, removeHs=False))[:n_max]
    views = []
    for mol in supplier:
        if mol and mol.GetNumConformers() > 0:
            views.append(rdkit_mol_to_view(mol))
    return views
```

## py3Dmol in Non-Notebook Scripts (offline PNG)

```python
import py3Dmol
import base64
from pathlib import Path

def save_png(pdb_path, out_png, width=800, height=600, style='cartoon'):
    """Render structure to PNG (requires display/Xvfb in headless env)."""
    view = py3Dmol.view(width=width, height=height)
    view.addModel(open(pdb_path).read(), 'pdb')
    if style == 'cartoon':
        view.setStyle({'cartoon': {'color': 'spectrum'}})
    elif style == 'stick':
        view.setStyle({}, {'stick': {'colorscheme': 'cyanCarbon'}})
    view.setBackgroundColor('white')
    view.zoomTo()

    # Note: .png() only works in Jupyter context
    # For headless: use write_html → capture with playwright/selenium
    html = view._make_html()
    html_path = str(out_png).replace('.png', '.html')
    Path(html_path).write_text(html)
    print(f"Saved HTML: {html_path}")
    print("Convert to PNG with: python -m playwright screenshot "
          f"--viewport 800,600 {html_path} {out_png}")

# Headless PNG via playwright (install: pip install playwright && playwright install)
# from playwright.sync_api import sync_playwright
# with sync_playwright() as p:
#     browser = p.chromium.launch()
#     page = browser.new_page(viewport={'width': 800, 'height': 600})
#     page.goto(f'file://{html_path}')
#     page.wait_for_timeout(2000)  # wait for WebGL render
#     page.screenshot(path=out_png)
#     browser.close()
```

## Common Pitfalls

| Problem | Cause | Fix |
|---------|-------|-----|
| View renders blank | Model not added before setStyle | Always addModel before setStyle |
| Wrong ligand highlighted | Wrong resn (LIG vs UNL vs MOL) | Check: `grep '^HETATM' file.pdb \| awk '{print $4}' \| sort -u` |
| `view.show()` does nothing | Not in Jupyter | Use `view._make_html()` or save HTML |
| Surface very slow | SES on full protein | Use VDW or restrict to pocket selection |
| setStyle clears previous style | `setStyle` on `{}` resets all | Use specific selection + `addStyle` or ordered calls |
| Labels overlap | Too many addLabel calls | Use `addResLabels` for selective residues |
