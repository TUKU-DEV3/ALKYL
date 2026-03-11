# MARTINI 3 Membrane Simulations

## MARTINI 3 Lipid Library

MARTINI 3 includes pre-parameterized lipids. Key ones:

| Lipid | MARTINI name | Notes |
|-------|-------------|-------|
| DPPC | DPPC | Saturated PC (dipalmitoyl) |
| POPC | POPC | Unsaturated PC (most common in models) |
| POPE | POPE | PE headgroup, unsaturated |
| POPS | POPS | PS headgroup (anionic) |
| CHOL | CHOL | Cholesterol |
| POPG | POPG | PG headgroup (bacterial) |
| DPPE | DPPE | PE saturated |
| GM3 | GM3 | Ganglioside |
| PIP2 | POP2 | PIP2 (anionic signaling) |
| Cer | CER | Ceramide (sphingolipid) |

Download: https://cgmartini.nl → Downloads → Lipids

## insane.py — Membrane Builder

The standard MARTINI membrane builder. Creates bilayer + water + ions.

```bash
# Install
pip install insane

# Simple POPC bilayer
insane \
  -l POPC:1 \              # inner leaflet composition
  -u POPC:1 \              # outer leaflet composition
  -d 0 \                   # distance between leaflets (0 = auto)
  -sol W \                 # solvent (W = MARTINI water)
  -salt 0.15 \             # salt concentration (M)
  -x 15 -y 15 -z 12 \      # box dimensions (nm)
  -o membrane.gro \
  -p membrane.top

# Mixed bilayer (e.g. raft-like composition)
insane \
  -l DPPC:0.4,DIPC:0.4,CHOL:0.2 \
  -u DPPC:0.4,DIPC:0.4,CHOL:0.2 \
  -sol W \
  -salt 0.15 \
  -x 20 -y 20 -z 15 \
  -o raft_membrane.gro \
  -p raft_membrane.top

# Asymmetric bilayer (plasma membrane model)
insane \
  -l POPE:0.4,POPS:0.15,POPC:0.4,CHOL:0.05 \   # inner leaflet
  -u POPC:0.7,CHOL:0.3 \                          # outer leaflet
  -sol W \
  -salt 0.15 \
  -x 20 -y 20 -z 15 \
  -o asymm_bilayer.gro \
  -p asymm_bilayer.top
```

### insane Python API

```python
import subprocess

def build_membrane(lipids_inner, lipids_outer, box=(15, 15, 12),
                    salt=0.15, out_prefix='membrane'):
    """
    Build MARTINI bilayer with insane.
    lipids_inner/outer: dict {lipid_name: fraction}
    box: (x, y, z) in nm
    """
    def format_composition(comp_dict):
        return ','.join(f'{k}:{v}' for k, v in comp_dict.items())

    cmd = [
        'insane',
        '-l', format_composition(lipids_inner),
        '-u', format_composition(lipids_outer),
        '-sol', 'W',
        '-salt', str(salt),
        '-x', str(box[0]), '-y', str(box[1]), '-z', str(box[2]),
        '-o', f'{out_prefix}.gro',
        '-p', f'{out_prefix}.top'
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"insane failed:\n{result.stderr}")
    print(result.stdout)
    return f'{out_prefix}.gro', f'{out_prefix}.top'
```

## CHARMM-GUI MARTINI Maker (Recommended for Complex Systems)

For production membranes, CHARMM-GUI Membrane Builder → MARTINI Bilayer Builder
gives better-validated starting structures with proper water padding.

Outputs: GROMACS .gro + .top + .mdp files ready to run.

## Protein-Membrane Embedding

```bash
# Method 1: insane with protein
insane \
  -f protein_cg.pdb \      # CG protein
  -l POPC:1 \
  -u POPC:1 \
  -sol W \
  -dm 0 \                  # protein insertion depth (nm, 0=auto)
  -o protein_membrane.gro \
  -p protein_membrane.top

# Method 2: embed into existing bilayer (MemProtMD approach)
# 1. Build bilayer without protein (insane)
# 2. Center protein at bilayer COM
# 3. Remove overlapping lipids (gmx editconf + gmx select)
```

```python
import MDAnalysis as mda
import numpy as np

def embed_protein_in_membrane(protein_cg_pdb, membrane_gro,
                               membrane_top, output_prefix='complex',
                               overlap_cutoff=0.4):
    """
    Embed CG protein into pre-built MARTINI membrane.
    Removes lipids within overlap_cutoff nm of protein beads.
    """
    # Load systems
    u_prot = mda.Universe(protein_cg_pdb)
    u_mem = mda.Universe(membrane_gro)

    # Center protein at bilayer center (z = box_z/2)
    box_z = u_mem.dimensions[2]
    prot_com_z = u_prot.atoms.center_of_geometry()[2]
    u_prot.atoms.translate([0, 0, box_z/2 - prot_com_z])

    # Find overlapping lipid atoms
    prot_pos = u_prot.atoms.positions
    lip_atoms = u_mem.select_atoms('not resname W NA CL')
    lip_pos = lip_atoms.positions

    # Distance matrix
    from scipy.spatial.distance import cdist
    dists = cdist(prot_pos, lip_pos)
    min_dists = dists.min(axis=0)

    # Remove overlapping lipid residues
    overlap_mask = min_dists < overlap_cutoff * 10  # Å
    overlap_resids = set(lip_atoms[overlap_mask].resids)
    keep_atoms = u_mem.select_atoms(
        f'resname W NA CL or (not resid {" ".join(map(str, overlap_resids))})'
    )

    # Merge
    from MDAnalysis.coordinates.memory import MemoryReader
    merged = mda.Merge(u_prot.atoms, keep_atoms)
    merged.atoms.write(f'{output_prefix}.pdb')
    print(f"Removed {len(overlap_resids)} lipid residues")
    return merged
```

## Membrane Topology File

```python
def update_topology(top_template, n_lipids_inner, n_lipids_outer,
                     lipid_composition, n_water, n_na, n_cl,
                     protein=False, output='system.top'):
    """Generate GROMACS topology for CG membrane system."""
    lines = [
        '#include "martini_v3.0.0.itp"',
        '#include "martini_v3.0.0_solvents_v1.itp"',
        '#include "martini_v3.0.0_ions_v1.itp"',
    ]
    # Lipid includes
    for lipid in set(lipid_composition.keys()):
        lines.append(f'#include "martini_v3.0.0_{lipid.lower()}.itp"')

    if protein:
        lines.append('#include "protein.itp"')

    lines += [
        '',
        '[ system ]',
        'MARTINI Membrane System',
        '',
        '[ molecules ]',
    ]
    if protein:
        lines.append('Protein    1')

    for lipid, count in lipid_composition.items():
        total = int((n_lipids_inner + n_lipids_outer) * count)
        lines.append(f'{lipid:<10} {total}')

    lines += [f'W          {n_water}',
              f'NA         {n_na}',
              f'CL         {n_cl}']

    with open(output, 'w') as f:
        f.write('\n'.join(lines))
    return output
```

## Membrane Equilibration Protocol

```
Step 1: Energy minimization (10,000 steps, steep)
Step 2: NVT, dt=0.005 ps, 10,000 steps (with position restraints on lipid headgroups)
Step 3: NPT, dt=0.010 ps, 10,000 steps (release restraints gradually)
Step 4: NPT production, dt=0.020 ps (20 fs)
```

```bash
# Self-assembly test: start from random lipid positions
insane -l POPC:1 -u POPC:1 -sol W -random -o random_start.gro -p random.top
# Run for 100 ns CG (~400 ns effective) → bilayer forms spontaneously
```
