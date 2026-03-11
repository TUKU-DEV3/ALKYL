# MDAnalysis — Protein-Ligand Binding Analysis

## Complete Workflow

```
Load traj → Align → RMSD (ligand) → Contacts → H-bonds → Binding site → Export
```

---

## Setup

```python
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align, contacts
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
import MDAnalysis.analysis.distances as D
import numpy as np
import matplotlib.pyplot as plt

u   = mda.Universe('complex.prmtop', 'traj.nc')
ref = mda.Universe('crystal.pdb')   # reference structure

protein = u.select_atoms('protein')
ligand  = u.select_atoms('resname LIG')
```

---

## 1. Ligand Stability (RMSD)

```python
# Align protein backbone, track ligand RMSD
rmsd = rms.RMSD(
    u,
    reference=ref,
    select='backbone',
    groupselections=['resname LIG'],
    ref_frame=0
)
rmsd.run()

times     = rmsd.results.rmsd[:, 1]
prot_rmsd = rmsd.results.rmsd[:, 2]
lig_rmsd  = rmsd.results.rmsd[:, 3]

fig, ax = plt.subplots()
ax.plot(times, prot_rmsd, label='Protein backbone')
ax.plot(times, lig_rmsd,  label='Ligand', linestyle='--')
ax.set(xlabel='Time (ps)', ylabel='RMSD (Å)', title='Complex stability')
ax.legend()
plt.tight_layout()
plt.savefig('rmsd.png', dpi=150)
```

---

## 2. Protein-Ligand Distance (CoM)

```python
prot_heavy = u.select_atoms('protein and not name H*')
lig_heavy  = u.select_atoms('resname LIG and not name H*')

com_dists = []
for ts in u.trajectory:
    d = np.linalg.norm(prot_heavy.center_of_mass() - lig_heavy.center_of_mass())
    com_dists.append(d)

print(f"Mean COM distance: {np.mean(com_dists):.2f} ± {np.std(com_dists):.2f} Å")
```

---

## 3. Per-Residue Contact Frequency

```python
cutoff = 5.0   # Å
contact_freq = {}

for ts in u.trajectory:
    # Atoms within cutoff of ligand
    nearby = u.select_atoms(f'protein and around {cutoff} resname LIG')
    for res in nearby.residues:
        key = (res.resid, res.resname)
        contact_freq[key] = contact_freq.get(key, 0) + 1

n_frames = len(u.trajectory)
df = [(resid, resname, count / n_frames)
      for (resid, resname), count in contact_freq.items()]
df.sort(key=lambda x: -x[2])

print(f"{'Residue':<10} {'Freq':>6}")
print("-" * 18)
for resid, resname, freq in df:
    if freq > 0.2:   # > 20% of frames
        print(f"{resname}{resid:<6}   {freq:.1%}")
```

---

## 4. Contact Map (protein–ligand, per frame)

```python
prot_res = u.select_atoms('protein and not name H*').split('residue')
lig_heavy = u.select_atoms('resname LIG and not name H*')

n_res    = len(prot_res)
n_frames = len(u.trajectory)
contact_map = np.zeros((n_res, n_frames))

for fi, ts in enumerate(u.trajectory):
    for ri, res_atoms in enumerate(prot_res):
        d = D.distance_array(res_atoms.positions, lig_heavy.positions)
        contact_map[ri, fi] = (d.min() < cutoff)

# Heatmap
import seaborn as sns
resids = [r.resids[0] for r in prot_res]
ax = sns.heatmap(contact_map, cmap='Blues', vmin=0, vmax=1,
                 yticklabels=resids, xticklabels=False)
ax.set(xlabel='Frame', ylabel='Residue', title='Ligand contact map')
plt.tight_layout()
plt.savefig('contact_map.png', dpi=150)
```

---

## 5. Hydrogen Bonds (Protein–Ligand)

```python
hbonds = HBA(
    universe=u,
    between=['protein', 'resname LIG'],
    d_a_cutoff=3.5,
    d_h_a_angle_cutoff=140.0,
    update_selections=True
)
hbonds.hydrogens_sel = hbonds.guess_hydrogens('protein or resname LIG')
hbonds.acceptors_sel = hbonds.guess_acceptors('protein or resname LIG')
hbonds.run()

# H-bonds per frame
counts = hbonds.count_by_time()
times  = [c[0] for c in counts]
nhb    = [c[1] for c in counts]

print(f"Mean H-bonds: {np.mean(nhb):.2f} ± {np.std(nhb):.2f}")

# Persistent H-bonds
table = hbonds.count_by_ids()
print("\nPersistent H-bonds (>30% of trajectory):")
for row in table:
    freq = row[3]
    if freq > 0.3:
        d_idx = int(row[0])
        a_idx = int(row[2])
        donor_atom    = u.atoms[d_idx]
        acceptor_atom = u.atoms[a_idx]
        print(f"  {donor_atom.resname}{donor_atom.resid}:{donor_atom.name} → "
              f"{acceptor_atom.resname}{acceptor_atom.resid}:{acceptor_atom.name}  "
              f"freq={freq:.1%}")
```

---

## 6. Binding Site RMSF

```python
from MDAnalysis.analysis import rms, align

# Align on protein core, calculate RMSF of binding site residues
aligner = align.AlignTraj(u, u, select='backbone', in_memory=True).run()

# Identify binding site: residues in contact >50% of frames
binding_resids = [resid for resid, _, freq in df if freq > 0.5]
bs_sel = ' or '.join([f'resid {r}' for r in binding_resids])
bs_atoms = u.select_atoms(f'({bs_sel}) and name CA')

rmsf_bs = rms.RMSF(bs_atoms).run()
print("Binding site RMSF (Å):")
for atom, val in zip(bs_atoms, rmsf_bs.results.rmsf):
    print(f"  {atom.resname}{atom.resid}: {val:.2f}")
```

---

## 7. Pocket Volume (HOLE/Voronoi)

```python
# Via MDAnalysis.analysis.hole2 (requires external HOLE program)
from MDAnalysis.analysis.hole2 import HoleAnalysis

with HoleAnalysis(u, select='protein') as ha:
    ha.run()
    ha.create_vmd_surface(filename='tunnel.vmd')
    print(ha.results.profiles[0])   # profile of frame 0
```

---

## 8. Summary Report

```python
print("=" * 50)
print("PROTEIN-LIGAND BINDING ANALYSIS SUMMARY")
print("=" * 50)
print(f"Trajectory:    {len(u.trajectory)} frames, "
      f"{u.trajectory.totaltime:.0f} ps")
print(f"Protein RMSD:  {np.mean(prot_rmsd):.2f} ± {np.std(prot_rmsd):.2f} Å")
print(f"Ligand RMSD:   {np.mean(lig_rmsd):.2f} ± {np.std(lig_rmsd):.2f} Å")
print(f"COM distance:  {np.mean(com_dists):.2f} ± {np.std(com_dists):.2f} Å")
print(f"Avg H-bonds:   {np.mean(nhb):.2f} ± {np.std(nhb):.2f}")
print(f"Contact residues (>50%): {len([r for r in df if r[2] > 0.5])}")
```

---

## Exporting for Visualization

```python
# Write last frame of trajectory as PDB
u.trajectory[-1]
u.atoms.write('last_frame.pdb')

# Write only protein + ligand
(protein + ligand).write('complex_last.pdb')

# Write trajectory subset (every 10th frame, protein only)
with mda.Writer('protein_thin.xtc', protein.n_atoms) as w:
    for ts in u.trajectory[::10]:
        w.write(protein)
```
