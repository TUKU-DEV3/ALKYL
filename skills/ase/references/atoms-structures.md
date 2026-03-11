# ASE — Atoms, Structures & I/O

## The Atoms Object

```python
from ase import Atoms
import numpy as np

# Manual construction: symbols + positions (Å)
atoms = Atoms(
    symbols='H2O',
    positions=[(0, 0, 0), (0.96, 0, 0), (-0.24, 0.93, 0)],
    cell=[10, 10, 10],       # box size in Å
    pbc=True                 # periodic boundary conditions
)

# Key attributes
atoms.positions          # (N, 3) array, Å
atoms.get_positions()    # copy of positions
atoms.set_positions(pos) # update positions
atoms.symbols            # ChemicalSymbols object ('H2O')
atoms.numbers            # atomic numbers array [8, 1, 1]
atoms.cell               # Cell object, 3×3 matrix
atoms.pbc                # [True, True, True] or per-axis
atoms.get_volume()       # cell volume, Å³
len(atoms)               # number of atoms
```

## Building Common Structures

### Molecules (gas phase)

```python
from ase.build import molecule

# Fetch from G2 database
water   = molecule('H2O')
benzene = molecule('C6H6')
co2     = molecule('CO2')
nh3     = molecule('NH3')

# List all available molecules
from ase.collections import g2
print(g2.names)
```

### Bulk Crystals

```python
from ase.build import bulk

# FCC metals
cu = bulk('Cu', 'fcc', a=3.61)
al = bulk('Al', 'fcc', a=4.05)

# BCC
fe = bulk('Fe', 'bcc', a=2.87)

# Diamond / Zincblende
si = bulk('Si', 'diamond', a=5.43)
gaas = bulk('GaAs', 'zincblende', a=5.65)

# Rocksalt, hexagonal
nacl = bulk('NaCl', 'rocksalt', a=5.64)
mg   = bulk('Mg', 'hcp', a=3.21, c=5.21)

# Supercell
from ase.build import make_supercell
P = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
super_cu = make_supercell(cu, P)
```

### Surfaces and Slabs

```python
from ase.build import fcc111, bcc110, hcp0001, surface, add_adsorbate

# FCC(111) — 4 layers, 10 Å vacuum
slab = fcc111('Cu', size=(3, 3, 4), vacuum=10.0)

# Add adsorbate
co = molecule('CO')
add_adsorbate(slab, co, height=2.0, position='fcc')  # fcc hollow site

# Fix bottom 2 layers
from ase.constraints import FixAtoms
mask = [atom.tag > 2 for atom in slab]
slab.set_constraint(FixAtoms(mask=mask))
```

### Nanoparticles

```python
from ase.cluster import Icosahedron, FaceCenteredCubic

# Icosahedral nanoparticle
np_ico = Icosahedron('Au', noshells=4)

# FCC nanoparticle with custom layers
np_fcc = FaceCenteredCubic('Pd', [(1,0,0), (1,1,0), (1,1,1)], [4, 5, 4])
```

## Cell and Periodicity

```python
# Set cell
atoms.set_cell([a, b, c])              # orthorhombic
atoms.set_cell([a, b, c, α, β, γ])    # angles in degrees
atoms.set_cell([[ax, ay, az],          # 3×3 matrix
                [bx, by, bz],
                [cx, cy, cz]])

# Center in cell
atoms.center(vacuum=5.0)              # add 5Å vacuum on all sides
atoms.center(about=[0, 0, 0])

# Periodic boundary conditions
atoms.pbc = True                       # all periodic
atoms.pbc = [True, True, False]        # slab (2D periodic)
atoms.pbc = False                      # isolated molecule
```

## File I/O

```python
from ase import io

# Read any format (auto-detected from extension)
atoms = io.read('structure.xyz')
atoms = io.read('POSCAR')              # VASP
atoms = io.read('geometry.in')        # FHI-aims
atoms = io.read('mol.sdf')            # SDF (requires RDKit or OpenBabel)
atoms = io.read('protein.pdb')        # PDB
atoms = io.read('crystal.cif')        # CIF

# Read multiple frames
traj = io.read('md.traj', index=':')  # all frames
last = io.read('md.traj', index=-1)   # last frame
slice_ = io.read('md.traj', index='::10')  # every 10th

# Write
io.write('output.xyz', atoms)
io.write('relaxed.cif', atoms)
io.write('POSCAR', atoms, vasp5=True)
io.write('traj.traj', atoms)

# Write trajectory (list of Atoms)
io.write('all.xyz', list_of_atoms)

# Supported formats
io.formats.show_all()
```

## Trajectory Files

```python
from ase.io.trajectory import Trajectory

# Write during simulation
traj = Trajectory('md.traj', 'w', atoms)
traj.write()       # snapshot now
traj.close()

# Append mode
with Trajectory('md.traj', 'a', atoms) as traj:
    traj.write()

# Read
traj = Trajectory('md.traj')
print(len(traj))         # number of frames
frame = traj[0]          # first frame (Atoms object)
last  = traj[-1]         # last frame
```

## Manipulating Atoms

```python
# Indexing and slicing
atom = atoms[0]          # single Atom object
sub  = atoms[1:4]        # Atoms subset
sub  = atoms[[0, 2, 5]]  # by index list

# Add / remove atoms
atoms.append(Atom('H', position=[1, 2, 3]))
del atoms[0]             # remove first atom
atoms += other_atoms     # concatenate (use copy() to avoid mutation)

# Translate, rotate
atoms.translate([1, 0, 0])
atoms.rotate(45, 'z', center='COP')  # 45° around z through COP
atoms.rotate(90, [1, 0, 0])

# Distances
d = atoms.get_distance(0, 1)
d_pbc = atoms.get_distance(0, 1, mic=True)   # minimum image convention
dmat = atoms.get_all_distances(mic=True)

# Neighbor list
from ase.neighborlist import NeighborList, natural_cutoffs
cutoffs = natural_cutoffs(atoms)
nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
nl.update(atoms)
indices, offsets = nl.get_neighbors(0)   # neighbors of atom 0
```

## ASE Database

```python
from ase.db import connect

db = connect('results.db')

# Store structure + metadata
db.write(atoms, name='water_opt', energy=atoms.get_potential_energy())

# Query
for row in db.select('name=water_opt'):
    print(row.energy, row.toatoms())

# Update
db.update(row.id, key=value)
```
