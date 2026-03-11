# MDAnalysis — Contacts & Hydrogen Bonds

## Native Contacts (Q value)

Measures fraction of native contacts preserved — canonical folding/binding metric.

```python
from MDAnalysis.analysis import contacts

u   = mda.Universe('system.prmtop', 'traj.nc')
ref = mda.Universe('native.pdb')      # native/crystal structure

# Q1: contacts within one selection
ca   = 'protein and name CA'
q1   = contacts.Contacts(
    u,
    select=ca,
    refgroup=(ref.select_atoms(ca),),
    radius=8.0,
    method='hard_cut'   # or 'soft_cut', 'radius_cut'
)
q1.run()
q_values = q1.results.q   # fraction per frame [0,1]

# Q2: inter-group contacts (e.g., protein–ligand)
prot_sel  = 'protein and not name H*'
lig_sel   = 'resname LIG and not name H*'

q2 = contacts.Contacts(
    u,
    select=(prot_sel, lig_sel),
    refgroup=(ref.select_atoms(prot_sel), ref.select_atoms(lig_sel)),
    radius=4.5,
    method='soft_cut'
)
q2.run()

# Plot
import matplotlib.pyplot as plt
plt.plot(u.trajectory.times, q2.results.q)
plt.xlabel('Time (ps)')
plt.ylabel('Q (fraction of native contacts)')
plt.show()
```

### Contact methods

| Method | Description |
|--------|-------------|
| `'hard_cut'` | 1 if distance ≤ radius, 0 otherwise |
| `'soft_cut'` | Smooth sigmoidal function (default β=5, λ=1.8) |
| `'radius_cut'` | 1 if within `radius * reference_distance` |

---

## Hydrogen Bond Analysis

```python
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA

u = mda.Universe('system.psf', 'traj.dcd')  # topology with bonds preferred

# Full protein H-bonds
hbonds = HBA(
    universe=u,
    # Optional: custom selections
    donors_sel=None,      # auto-detect
    hydrogens_sel=None,   # auto-detect
    acceptors_sel=None,   # auto-detect
    d_h_cutoff=1.2,       # Å — donor-H bond max length
    d_a_cutoff=3.0,       # Å — donor-acceptor cutoff
    d_h_a_angle_cutoff=150.0,  # degrees — min angle
    update_selections=True     # recompute selections each frame
)
hbonds.run()

# Results: array columns = [frame, donor_ix, H_ix, acceptor_ix, dist, angle]
result = hbonds.results.hbonds
```

### Protein–Ligand H-bonds

```python
# Targeted: only protein ↔ ligand H-bonds
prot_donors = u.select_atoms(
    'protein and (name N* O* S*) and not name H*'
)
lig_acceptors = u.select_atoms('resname LIG and (name N* O* S* F*)')
lig_donors    = u.select_atoms('resname LIG and (name N* O*) and not name H*')
prot_acceptors = u.select_atoms('protein and (name N* O* S*) and not name H*')

hbonds = HBA(
    universe=u,
    between=['protein', 'resname LIG'],   # restrict to these two groups
    d_a_cutoff=3.5,
    d_h_a_angle_cutoff=140.0
)
hbonds.hydrogens_sel = hbonds.guess_hydrogens('protein or resname LIG')
hbonds.acceptors_sel = hbonds.guess_acceptors('protein or resname LIG')
hbonds.run()
```

### Summarizing H-bond Results

```python
# H-bonds per timestep
counts = hbonds.count_by_time()
# Returns: [[frame0_time, n_hbonds], [frame1_time, n_hbonds], ...]

# Most persistent H-bonds (sorted by frequency)
table = hbonds.count_by_ids()
# Columns: [donor_ix, H_ix, acceptor_ix, frequency (0–1)]
print(table[:10])   # top 10 most persistent

# By residue type pair
type_counts = hbonds.count_by_type()

# Lifetime / autocorrelation
lifetimes = hbonds.lifetime()

# Plot H-bonds per frame
times = [r[0] for r in counts]
nhb   = [r[1] for r in counts]
plt.plot(times, nhb)
plt.xlabel('Time (ps)'); plt.ylabel('# H-bonds')
plt.show()
```

---

## Distances — Residue/Atom Level

```python
from MDAnalysis.analysis.distances import dist, between
import MDAnalysis.analysis.distances as D
import numpy as np

# Pairwise minimum distances between two groups
d = D.distance_array(
    ag1.positions,
    ag2.positions
)   # shape (n1, n2)

# Minimum distance per frame
min_dists = []
for ts in u.trajectory:
    d = D.distance_array(ag1.positions, ag2.positions)
    min_dists.append(d.min())

# Centre-of-mass distance between protein and ligand
prot = u.select_atoms('protein')
lig  = u.select_atoms('resname LIG')

com_dists = []
for ts in u.trajectory:
    d = np.linalg.norm(prot.center_of_mass() - lig.center_of_mass())
    com_dists.append(d)
```

---

## Binding Site Residues (Dynamic)

```python
# Residues within 5 Å of ligand — re-evaluated each frame
binding_residues_over_time = set()
contact_counts = {}

lig = u.select_atoms('resname LIG')

for ts in u.trajectory:
    # Select protein residues with any atom within 5Å of ligand
    nearby = u.select_atoms(
        'protein and (around 5.0 resname LIG)',
        updating=False   # evaluated once per loop iteration — correct
    )
    for res in nearby.residues:
        key = (res.resname, res.resid)
        contact_counts[key] = contact_counts.get(key, 0) + 1
        binding_residues_over_time.add(key)

# Sort by frequency
n_frames = len(u.trajectory)
persistent = {k: v/n_frames for k, v in contact_counts.items() if v/n_frames > 0.5}
print("Residues in contact > 50% of trajectory:")
for (resname, resid), freq in sorted(persistent.items(), key=lambda x: -x[1]):
    print(f"  {resname}{resid}: {freq:.1%}")
```

---

## Water Bridge Analysis

```python
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import WaterBridgeAnalysis

wb = WaterBridgeAnalysis(
    u,
    selection1='resname LIG',
    selection2='protein',
    water_selection='resname HOH TIP3',
    order=1           # 1 = single water bridge
)
wb.run()
wb.generate_table()
print(wb.table[:5])
```
