# MDAnalysis — Universe, Selections & I/O

## Universe Creation

```python
import MDAnalysis as mda

# Topology only (single frame)
u = mda.Universe('protein.pdb')
u = mda.Universe('system.gro')

# Topology + trajectory
u = mda.Universe('system.prmtop', 'traj.nc')       # AMBER
u = mda.Universe('system.psf', 'traj.dcd')         # CHARMM/NAMD
u = mda.Universe('system.tpr', 'traj.xtc')         # GROMACS
u = mda.Universe('system.tpr', ['run1.xtc', 'run2.xtc'])  # multiple trajs

# Universe info
print(u)                     # <Universe with N atoms>
print(u.atoms.n_atoms)       # total atoms
print(u.residues.n_residues) # total residues
print(len(u.trajectory))     # frames
print(u.trajectory.dt)       # time step, ps
print(u.trajectory.totaltime) # total time, ps
```

## Atom Attributes

```python
atoms = u.atoms

atoms.positions      # (N, 3) array, Å
atoms.names          # atom names: ['N', 'CA', 'C', ...]
atoms.resnames       # residue names: ['ALA', 'GLY', ...]
atoms.resids         # residue numbers: [1, 1, 1, 2, ...]
atoms.masses         # atomic masses, amu
atoms.charges        # partial charges (if in topology)
atoms.types          # atom types
atoms.elements       # element symbols (if available)
atoms.ids            # 1-based serial indices
atoms.indices        # 0-based indices

# Via residues
u.residues.resnames
u.residues.resids
u.residues.segids
```

## Selection Language

```python
ag = u.select_atoms('selection string')
```

### Standard keywords

```python
# Structural
u.select_atoms('protein')                   # all protein atoms
u.select_atoms('backbone')                  # N, CA, C, O
u.select_atoms('nucleic')                   # nucleic acid atoms
u.select_atoms('water')                     # HOH, TIP3, etc.

# By name/type
u.select_atoms('name CA')                   # alpha-carbons
u.select_atoms('name CA CB')               # multiple names
u.select_atoms('name CA*')                  # wildcard: CA, CAY, ...
u.select_atoms('type C')                    # by force field atom type

# By residue
u.select_atoms('resname ALA')
u.select_atoms('resname ALA GLY')          # multiple resnames
u.select_atoms('resid 10')
u.select_atoms('resid 1:50')               # range (inclusive)
u.select_atoms('resid 10 20 30')           # specific list
u.select_atoms('resnum 42')               # canonical PDB numbering

# By chain / segment
u.select_atoms('segid A')
u.select_atoms('chainID B')

# By element / charge
u.select_atoms('element N O')
u.select_atoms('formalcharge -1')
```

### Boolean operators

```python
u.select_atoms('protein and name CA')
u.select_atoms('resname ASP or resname GLU')
u.select_atoms('not water')
u.select_atoms('backbone and not name O')
u.select_atoms('(resname GLU or resname HIS) and name CA and resid 1:100')
```

### Geometric selections

```python
# Sphere around a selection
u.select_atoms('around 5.0 resname LIG')       # within 5 Å of ligand

# Point in space
u.select_atoms('point 10.5 20.0 5.0 8.0')     # within 8 Å of point

# Property (coordinate-based)
u.select_atoms('prop z > 10.0')                # atoms with z > 10 Å
u.select_atoms('prop abs z < 5.0')

# Spherical layer
u.select_atoms('sphlayer 5.0 10.0 resname LIG')   # 5–10 Å from ligand
```

### Expanding selections

```python
# Same residue as matched atoms
u.select_atoms('same residue as around 5.0 resname LIG')

# Whole residue
u.select_atoms('byres name CA')    # all atoms of residues containing CA

# Bonded neighbors
u.select_atoms('bonded resname LIG')
```

### Updating selections (dynamic — recalculate each frame)

```python
# Static: evaluated ONCE
static  = u.select_atoms('around 5.0 resname LIG')

# Updating: re-evaluated EVERY frame — essential for binding analysis
dynamic = u.select_atoms('around 5.0 resname LIG', updating=True)

for ts in u.trajectory:
    print(dynamic.n_atoms)   # changes as atoms move in/out
```

### SMARTS selection (requires RDKit)

```python
u.select_atoms('smarts [NH2]')    # primary amine
```

---

## Trajectory Iteration

```python
# Full iteration
for ts in u.trajectory:
    t     = ts.time          # ps
    frame = ts.frame         # 0-indexed
    pos   = u.atoms.positions  # current (N,3) array

# Sliced iteration
for ts in u.trajectory[10:50:2]:   # frames 10–50 every 2
    pass

# Jump to specific frame
u.trajectory[100]            # jump to frame 100
u.trajectory[-1]             # last frame

# Current frame info
ts = u.trajectory.ts
print(ts.frame, ts.time, ts.dt)
```

---

## Writing Trajectories and Structures

```python
# Single frame
selection = u.select_atoms('protein')
selection.write('protein.pdb')
selection.write('protein.gro')
selection.write('protein.mol2')

# Full trajectory (subset of atoms)
ca = u.select_atoms('name CA')
with mda.Writer('ca_traj.xtc', ca.n_atoms) as w:
    for ts in u.trajectory:
        w.write(ca)

# Slice trajectory
with mda.Writer('last100.xtc', u.atoms.n_atoms) as w:
    for ts in u.trajectory[-100:]:
        w.write(u.atoms)
```

---

## Adding Topology Attributes

```python
# Add B-factor column (for visualization)
u.add_TopologyAttr('tempfactors')
protein = u.select_atoms('protein')
protein.atoms.tempfactors = 0.0   # set to 0 initially

# Assign per-residue values
for residue, value in zip(protein.residues, rmsf_values):
    residue.atoms.tempfactors = value

protein.write('colored.pdb')   # B-factor column = RMSF
```

---

## Useful AtomGroup Properties

```python
ag = u.select_atoms('protein')

ag.center_of_mass()      # (3,) array
ag.center_of_geometry()  # unweighted center
ag.radius_of_gyration()  # scalar, Å
ag.total_mass()          # scalar, amu
ag.total_charge()        # total charge (if available)
ag.n_atoms               # int
ag.n_residues            # int
ag.residues              # ResidueGroup
ag.segments              # SegmentGroup

# Distances between two groups
from MDAnalysis.analysis.distances import dist
d = dist(ag1, ag2)       # (3, N) — dist, resid1, resid2
```

---

## Merging Universes / Combining AtomGroups

```python
# Concatenate AtomGroups (preserves order)
combined = ag1 + ag2

# Merge two Universes into one
merged = mda.core.universe.Merge(u1.atoms, u2.atoms)
merged.write('merged.pdb')
```
