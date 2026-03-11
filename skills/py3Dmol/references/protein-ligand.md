# Protein-Ligand Visualization

## Standard Protein-Ligand Complex

```python
import py3Dmol

def view_complex(pdb_path, ligand_resn='LIG', width=900, height=550):
    """
    Standard protein-ligand visualization.
    ligand_resn: residue name of ligand in PDB (LIG, UNL, MOL, etc.)
    """
    view = py3Dmol.view(width=width, height=height)
    view.addModel(open(pdb_path).read(), 'pdb')
    view.setBackgroundColor('white')

    # Protein: cartoon, grey
    view.setStyle({'cartoon': {'color': 'lightgrey', 'opacity': 0.8}})

    # Binding site residues: stick (within 5 Å of ligand)
    view.setStyle(
        {'within': {'distance': 5, 'sel': {'resn': ligand_resn}}},
        {'stick': {'colorscheme': 'whiteCarbon', 'radius': 0.15}}
    )

    # Ligand: colored stick
    view.setStyle(
        {'resn': ligand_resn},
        {'stick': {'colorscheme': 'greenCarbon', 'radius': 0.2},
         'sphere': {'colorscheme': 'greenCarbon', 'radius': 0.35}}
    )

    # Hide water
    view.setStyle({'resn': 'HOH'}, {})

    view.zoomTo({'resn': ligand_resn})
    return view

view_complex('docked_complex.pdb', ligand_resn='LIG').show()
```

## Binding Pocket Surface

```python
def view_binding_pocket(pdb_path, ligand_resn='LIG', pocket_radius=8.0):
    """Show molecular surface of binding pocket around ligand."""
    view = py3Dmol.view(width=900, height=550)
    view.addModel(open(pdb_path).read(), 'pdb')
    view.setBackgroundColor('#1a1a2e')   # dark background for surface

    # Transparent cartoon
    view.setStyle({'cartoon': {'color': 'grey', 'opacity': 0.3}})

    # Pocket surface (within N Å of ligand)
    view.addSurface(
        py3Dmol.SES,
        {'opacity': 0.6, 'color': 'white'},
        {'within': {'distance': pocket_radius, 'sel': {'resn': ligand_resn}}}
    )

    # Ligand on top
    view.setStyle(
        {'resn': ligand_resn},
        {'stick': {'colorscheme': 'greenCarbon', 'radius': 0.25},
         'sphere': {'colorscheme': 'greenCarbon', 'radius': 0.4}}
    )

    view.zoomTo({'resn': ligand_resn})
    return view
```

## Docking Pose Batch Viewer

Inspect multiple poses from Vina/Gnina output:

```python
from rdkit import Chem
import py3Dmol

def view_docking_poses(receptor_pdb, poses_sdf, n_poses=5,
                       ligand_resn=None, width=900, height=550):
    """
    Show top N docking poses overlaid on receptor.
    poses_sdf: SDF file with multiple conformers (Vina/Gnina output)
    """
    # Color palette for poses
    colors = ['green', 'cyan', 'magenta', 'orange', 'yellow']

    view = py3Dmol.view(width=width, height=height)
    view.setBackgroundColor('white')

    # Receptor
    view.addModel(open(receptor_pdb).read(), 'pdb')
    view.setStyle({'cartoon': {'color': 'lightgrey', 'opacity': 0.7}})

    # Hide water/ions
    view.setStyle({'resn': 'HOH'}, {})

    # Binding site (use first pose as reference for zoom)
    supplier = Chem.SDMolSupplier(poses_sdf, removeHs=False)
    first_pose_sdf = Chem.MolToMolBlock(next(iter(supplier)))

    # Add each pose
    supplier = Chem.SDMolSupplier(poses_sdf, removeHs=True)
    for i, mol in enumerate(supplier):
        if i >= n_poses or mol is None:
            break
        sdf_str = Chem.MolToMolBlock(mol)
        model_idx = i + 1   # model 0 = receptor
        view.addModel(sdf_str, 'sdf')
        color = colors[i % len(colors)]
        view.setStyle(
            {'model': model_idx},
            {'stick': {'color': color, 'radius': 0.15},
             'sphere': {'color': color, 'radius': 0.28}}
        )

    # Zoom to first pose
    view.addModel(first_pose_sdf, 'sdf')
    view.setStyle({'model': n_poses + 1}, {})   # invisible, just for zoom
    view.zoomTo({'model': 1})
    return view

view_docking_poses('receptor.pdb', 'vina_poses.sdf', n_poses=5).show()
```

## Two-Structure Overlay (e.g., before/after minimization)

```python
def view_overlay(pdb1, pdb2, label1='Structure 1', label2='Structure 2'):
    """Overlay two structures, colored differently."""
    view = py3Dmol.view(width=900, height=500)
    view.setBackgroundColor('white')

    view.addModel(open(pdb1).read(), 'pdb')
    view.setStyle({'model': 0}, {'cartoon': {'color': '#4488cc', 'opacity': 0.8}})

    view.addModel(open(pdb2).read(), 'pdb')
    view.setStyle({'model': 1}, {'cartoon': {'color': '#cc4444', 'opacity': 0.8}})

    # Add text labels
    view.addLabel(label1, {'position': {'x': 0, 'y': 0, 'z': 0},
                            'fontColor': '#4488cc', 'fontSize': 14})
    view.addLabel(label2, {'position': {'x': 5, 'y': 0, 'z': 0},
                            'fontColor': '#cc4444', 'fontSize': 14})
    view.zoomTo()
    return view
```

## Highlighting Residues from ProLIF Analysis

```python
def view_prolif_interactions(pdb_path, hbond_residues, hydrophobic_residues,
                              ligand_resn='LIG'):
    """
    Color binding site residues by interaction type from ProLIF analysis.
    hbond_residues: list of residue numbers (int) forming H-bonds
    hydrophobic_residues: list of residue numbers forming hydrophobic contacts
    """
    view = py3Dmol.view(width=900, height=600)
    view.addModel(open(pdb_path).read(), 'pdb')
    view.setBackgroundColor('white')

    # Protein backbone cartoon (grey)
    view.setStyle({'cartoon': {'color': 'lightgrey', 'opacity': 0.6}})

    # H-bond residues: blue sticks
    for resi in hbond_residues:
        view.setStyle(
            {'resi': resi},
            {'stick': {'color': '#4488ff', 'radius': 0.2},
             'cartoon': {'color': '#4488ff'}}
        )

    # Hydrophobic residues: orange sticks
    for resi in hydrophobic_residues:
        view.setStyle(
            {'resi': resi},
            {'stick': {'color': '#ff8844', 'radius': 0.2},
             'cartoon': {'color': '#ff8844'}}
        )

    # Ligand: green
    view.setStyle(
        {'resn': ligand_resn},
        {'stick': {'colorscheme': 'greenCarbon', 'radius': 0.22},
         'sphere': {'colorscheme': 'greenCarbon', 'radius': 0.38}}
    )

    view.zoomTo({'resn': ligand_resn})
    return view
```

## MD Snapshot Viewer

```python
import MDAnalysis as mda

def view_md_snapshot(topology, trajectory, frame_idx=0,
                     protein_sel='protein', ligand_sel='resname LIG'):
    """Visualize a single MD frame."""
    u = mda.Universe(topology, trajectory)
    u.trajectory[frame_idx]

    # Write frame to PDB string in memory
    import io
    protein = u.select_atoms(f'{protein_sel} or {ligand_sel}')
    with mda.Writer('/tmp/_snapshot.pdb', n_atoms=len(protein)) as w:
        w.write(protein)

    return view_complex('/tmp/_snapshot.pdb',
                        ligand_resn=ligand_sel.split()[-1])
```

## Coloring by B-factor (e.g., pLDDT from AlphaFold)

```python
def view_by_bfactor(pdb_path, label='pLDDT', vmin=0, vmax=100):
    """
    Color structure by B-factor column.
    For AlphaFold: B-factor = pLDDT (0–100).
    """
    view = py3Dmol.view(width=800, height=500)
    view.addModel(open(pdb_path).read(), 'pdb')
    view.setBackgroundColor('white')

    # Use built-in b-factor coloring
    view.setStyle({
        'cartoon': {
            'colorscheme': {
                'prop': 'b',
                'gradient': 'roygb',   # red=low, blue=high
                'min': vmin,
                'max': vmax
            }
        }
    })
    view.zoomTo()
    return view

# AlphaFold pLDDT: >90 dark blue (very high), 70-90 cyan (confident),
#                  50-70 yellow (low), <50 orange (very low)
```
