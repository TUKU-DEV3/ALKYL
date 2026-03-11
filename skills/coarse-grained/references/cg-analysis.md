# CG Trajectory Analysis

## Loading CG Trajectories with MDAnalysis

```python
import MDAnalysis as mda
import numpy as np

# MARTINI GROMACS output
u = mda.Universe('md.tpr', 'md_whole.xtc')

# Selection language works identically to AA
protein = u.select_atoms('protein')           # backbone beads
membrane = u.select_atoms('resname POPC CHOL')
water = u.select_atoms('resname W')
upper_leaflet = u.select_atoms('resname POPC and prop z > 6.0')  # rough split
```

## Membrane Thickness

```python
def compute_membrane_thickness(u, lipid_sel='resname POPC',
                                headgroup_name='PO4',
                                n_frames=None):
    """
    Membrane thickness = mean distance between upper and lower PO4 leaflets.
    PO4 is the phosphate bead in POPC/POPE/POPS.
    """
    po4 = u.select_atoms(f'{lipid_sel} and name {headgroup_name}')
    thicknesses = []

    for ts in u.trajectory[-n_frames:] if n_frames else u.trajectory:
        z = po4.positions[:, 2]
        z_com = po4.center_of_geometry()[2]

        upper_z = z[z > z_com].mean()
        lower_z = z[z < z_com].mean()
        thicknesses.append(upper_z - lower_z)

    return {
        'mean': np.mean(thicknesses),
        'std': np.std(thicknesses),
        'values': thicknesses
    }
# Reference: POPC bilayer ~4.0 nm (all-atom ~3.7 nm, CG slightly thicker)
```

## Area Per Lipid (APL)

```python
def compute_apl(u, lipid_sel='resname POPC CHOL',
                n_frames=None):
    """
    Area per lipid = (box_x × box_y) / n_lipids_per_leaflet
    """
    n_lipids = len(u.select_atoms(lipid_sel).residues)
    n_per_leaflet = n_lipids // 2

    apls = []
    for ts in u.trajectory[-n_frames:] if n_frames else u.trajectory:
        box = u.dimensions
        area = box[0] * box[1]   # Å²
        apl = area / n_per_leaflet   # Å² per lipid
        apls.append(apl)

    return {
        'mean_A2': np.mean(apls),
        'std_A2': np.std(apls),
        'mean_nm2': np.mean(apls) / 100,
        'values': apls
    }
# Reference: POPC APL ~65 Å² (AA: 68 Å², CG slightly underestimates)
```

## Lipid Order Parameter (Scd)

```python
def compute_scd_martini(u, lipid_resname='POPC', tail_beads=None):
    """
    Compute MARTINI order parameter analog (P2 order parameter).
    For POPC tail: C1A-C2A-C3A-C4A (sn-1), C1B-C2B-C3B-C4B (sn-2).
    P2 = 0.5 * (3*cos²θ - 1), θ = angle between bond vector and membrane normal.
    """
    if tail_beads is None:
        # POPC MARTINI 3 tail bead names
        tail_beads = {
            'sn1': ['C1A', 'C2A', 'C3A', 'C4A'],
            'sn2': ['C1B', 'C2B', 'C3B', 'C4B']
        }

    membrane_normal = np.array([0, 0, 1])
    results = {tail: [] for tail in tail_beads}

    for ts in u.trajectory:
        for tail_name, beads in tail_beads.items():
            tail_p2 = []
            for i in range(len(beads) - 1):
                b1 = u.select_atoms(f'resname {lipid_resname} and name {beads[i]}')
                b2 = u.select_atoms(f'resname {lipid_resname} and name {beads[i+1]}')
                if len(b1) == 0 or len(b2) == 0:
                    continue

                bond_vectors = b2.positions - b1.positions
                # Normalize
                norms = np.linalg.norm(bond_vectors, axis=1, keepdims=True)
                bond_vectors /= norms + 1e-9

                cos_theta = bond_vectors @ membrane_normal
                p2 = 0.5 * (3 * cos_theta**2 - 1)
                tail_p2.append(p2.mean())

            results[tail_name].append(tail_p2)

    # Average over frames
    for tail_name in results:
        arr = np.array(results[tail_name])
        results[tail_name] = {
            'mean': arr.mean(axis=0),
            'std': arr.std(axis=0)
        }
    return results
```

## Lateral Diffusion Coefficient

```python
from MDAnalysis.analysis import msd as MSD

def compute_lateral_diffusion(u, sel='resname POPC and name PO4',
                               start=None, stop=None, step=1):
    """
    Compute lateral (2D) diffusion coefficient from MSD.
    D_lateral = MSD_xy / (4 * t)  [nm²/ns → ×10⁻⁹ m²/s for SI]
    """
    atoms = u.select_atoms(sel)
    msd_calc = MSD.EinsteinMSD(u, select=sel, msd_type='xy',
                                fft=True)
    msd_calc.run(start=start, stop=stop, step=step)

    msd = msd_calc.results.timeseries   # Å²
    times = msd_calc.results.times      # ps

    # Linear fit in the diffusive regime (skip ballistic at t<1 ps)
    # D = slope / 4 (2D MSD = 4Dt)
    diffusive_mask = (times > 1000) & (times < 0.5 * times[-1])  # ps
    if diffusive_mask.sum() < 5:
        diffusive_mask = np.ones(len(times), dtype=bool)

    from scipy.stats import linregress
    slope, intercept, r, p, se = linregress(times[diffusive_mask],
                                              msd[diffusive_mask])
    D_A2_ps = slope / 4        # Å²/ps
    D_nm2_ns = D_A2_ps / 10   # nm²/ns = 10⁻⁹ m²/s

    return {
        'D_nm2_ns': D_nm2_ns,
        'D_cm2_s': D_nm2_ns * 1e-5,  # cm²/s
        'D_r2': r**2,
        'times_ps': times,
        'msd_A2': msd
    }
# Reference: POPC lateral diffusion ~5×10⁻⁸ cm²/s (exp ~1×10⁻⁸)
# CG overestimates by ×4 (same time-scaling factor)
```

## Membrane Density Profile

```python
def density_profile(u, selections_dict, n_bins=100, axis=2):
    """
    Compute number density profiles along membrane normal (z).
    selections_dict: {'label': 'selection string', ...}
    """
    import pandas as pd

    box_dim = u.dimensions[axis]
    bins = np.linspace(0, box_dim, n_bins + 1)
    bin_centers = (bins[:-1] + bins[1:]) / 2

    profiles = {'z_A': bin_centers}

    for label, sel in selections_dict.items():
        atoms = u.select_atoms(sel)
        counts = np.zeros(n_bins)

        for ts in u.trajectory:
            box = u.dimensions[axis]
            positions = atoms.positions[:, axis] % box
            hist, _ = np.histogram(positions, bins=bins)
            counts += hist

        counts /= len(u.trajectory)   # average over frames
        volume_per_bin = u.dimensions[0] * u.dimensions[1] * (bins[1] - bins[0])
        profiles[label] = counts / volume_per_bin * 1e3   # number/nm³

    return pd.DataFrame(profiles)

# Example usage:
# profiles = density_profile(u, {
#     'PO4': 'resname POPC and name PO4',
#     'Protein_BB': 'protein and name BB',
#     'Water': 'resname W',
# })
```

## Protein CG Analysis

```python
from MDAnalysis.analysis import rms, align

def cg_protein_rmsd(u, ref_universe=None, select='backbone'):
    """
    RMSD of CG protein backbone (BB beads in MARTINI).
    In MARTINI: 'backbone' maps to 'name BB'.
    """
    if 'BB' in u.select_atoms('protein').names:
        cg_sel = 'name BB'   # MARTINI backbone beads
    else:
        cg_sel = 'backbone'

    if ref_universe is None:
        ref_universe = u

    R = rms.RMSD(u, ref_universe, select=cg_sel,
                  groupselections=[cg_sel])
    R.run()
    return R.results.rmsd   # (frame, time, RMSD) array

def cg_protein_rmsf(u, select='name BB', align_first=True):
    """RMSF of CG backbone beads."""
    if align_first:
        align.AlignTraj(u, u, select=select).run()
    from MDAnalysis.analysis import rms
    rmsf = rms.RMSF(u.select_atoms(select))
    rmsf.run()
    return rmsf.results.rmsf
```

## Cholesterol Flip-Flop Rate

```python
def cholesterol_flip_flop_rate(u, chol_sel='resname CHOL',
                                headgroup='ROH', leaflet_cutoff_nm=0.5):
    """
    Count cholesterol flip-flop events (crossing membrane midplane).
    Returns flip-flop rate per lipid per µs.
    """
    chol = u.select_atoms(f'{chol_sel} and name {headgroup}')
    box_z = u.dimensions[2]
    midplane = box_z / 2

    flip_flops = 0
    prev_leaflet = None

    for ts in u.trajectory:
        z = chol.positions[:, 2]
        leaflet = (z > midplane).astype(int)   # 1=upper, 0=lower

        if prev_leaflet is not None:
            changes = (leaflet != prev_leaflet).sum()
            flip_flops += changes

        prev_leaflet = leaflet.copy()

    n_chol = len(chol.residues)
    total_time_ns = (u.trajectory[-1].time - u.trajectory[0].time) / 1000
    rate = flip_flops / n_chol / total_time_ns * 1000   # per lipid per µs

    return {'n_events': flip_flops, 'rate_per_us': rate,
            'total_time_ns': total_time_ns}
```

## Visualization of CG Trajectories

```python
import nglview as nv
import MDAnalysis as mda

def view_cg_membrane(gro_file, xtc_file=None, n_frames=100):
    """
    Visualize MARTINI membrane in NGLview.
    NGLview preferred over py3Dmol for trajectories (see py3Dmol skill).
    """
    u = mda.Universe(gro_file, xtc_file or gro_file)

    if xtc_file:
        # Stride to n_frames
        stride = max(1, len(u.trajectory) // n_frames)
        view = nv.show_mdanalysis(u)
    else:
        view = nv.show_file(gro_file)

    # Protein: cartoon
    view.add_cartoon(selection='protein', color_scheme='residueindex')
    # Lipid headgroups: ball+stick
    view.add_ball_and_stick(selection='resname POPC and name PO4',
                             color='red', sphere_scale=1.5)
    # Lipid tails: line
    view.add_line(selection='resname POPC and not name PO4', color='grey')
    # Cholesterol
    view.add_ball_and_stick(selection='resname CHOL', color='yellow')

    return view
```

## Analysis Summary Report

```python
def cgmd_analysis_report(u, lipid_sel='resname POPC', output='cg_report.txt'):
    """Run all standard membrane analyses and write summary."""
    lines = ['CG-MD Analysis Report', '='*40]

    # Thickness
    thick = compute_membrane_thickness(u, lipid_sel)
    lines.append(f'Membrane thickness: {thick["mean"]:.2f} ± {thick["std"]:.2f} Å')

    # APL
    apl = compute_apl(u, lipid_sel)
    lines.append(f'Area per lipid: {apl["mean_A2"]:.1f} ± {apl["std_A2"]:.1f} Å²')

    # Diffusion (PO4 headgroup)
    diff = compute_lateral_diffusion(u, f'{lipid_sel} and name PO4')
    lines.append(f'Lateral diffusion: {diff["D_cm2_s"]:.2e} cm²/s '
                 f'(R²={diff["D_r2"]:.3f})')
    lines.append(f'  Corrected (÷4): {diff["D_cm2_s"]/4:.2e} cm²/s')

    with open(output, 'w') as f:
        f.write('\n'.join(lines))
    print('\n'.join(lines))
```
