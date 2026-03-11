# MARTINI 3 Protein Setup

## martinize2

The primary tool for converting AA protein structures to MARTINI 3.

```bash
# Install
pip install martinize2
# or: conda install -c conda-forge martinize2

# Basic usage
martinize2 \
  -f protein.pdb \
  -o protein_cg.itp \
  -x protein_cg.pdb \
  -ff martini3001 \          # MARTINI 3.0.0.1
  -dssp dssp \               # secondary structure from DSSP
  -elastic \                 # add elastic network (ElNeDyn)
  -ef 500 \                  # elastic force constant (kJ/mol/nm²)
  -el 0.5 \                  # elastic bond lower cutoff (nm)
  -eu 0.9 \                  # elastic bond upper cutoff (nm)
  -p backbone \              # position restraints on backbone beads
  -pf 1000                   # restraint force constant
```

### martinize2 Python API

```python
import subprocess
from pathlib import Path

def martinize_protein(input_pdb, output_dir='.', ff='martini3001',
                       elastic=True, ef=500, el=0.5, eu=0.9):
    """
    Convert AA protein PDB to MARTINI 3 CG topology.
    Returns paths to generated files.
    """
    out_dir = Path(output_dir)
    out_dir.mkdir(exist_ok=True)
    prefix = Path(input_pdb).stem

    cmd = [
        'martinize2',
        '-f', str(input_pdb),
        '-o', str(out_dir / f'{prefix}.itp'),
        '-x', str(out_dir / f'{prefix}_cg.pdb'),
        '-ff', ff,
        '-dssp', 'dssp',
    ]
    if elastic:
        cmd += ['-elastic', '-ef', str(ef), '-el', str(el), '-eu', str(eu)]

    result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(out_dir))
    if result.returncode != 0:
        raise RuntimeError(f"martinize2 failed:\n{result.stderr}")

    return {
        'itp': out_dir / f'{prefix}.itp',
        'cg_pdb': out_dir / f'{prefix}_cg.pdb',
        'stdout': result.stdout
    }
```

## Elastic Network Models

Elastic networks (ElNeDyn) maintain protein secondary/tertiary structure
by adding harmonic bonds between nearby backbone beads.

```
V_elastic = 0.5 * ef * (r - r0)²
```
Applied to all backbone bead pairs with r₀ between el and eu cutoffs.

```python
def parse_elastic_bonds(itp_file):
    """Parse elastic bonds from martinize2-generated .itp file."""
    bonds = []
    in_bonds = False
    with open(itp_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith('[ bonds ]'):
                in_bonds = True
                continue
            if line.startswith('[') and in_bonds:
                in_bonds = False
            if in_bonds and line and not line.startswith(';'):
                parts = line.split()
                if len(parts) >= 4:
                    bonds.append({
                        'i': int(parts[0]),
                        'j': int(parts[1]),
                        'type': int(parts[2]),
                        'r0': float(parts[3]),
                        'k': float(parts[4]) if len(parts) > 4 else None
                    })
    return bonds
```

## Go-MARTINI (Structure-Based CG)

For intrinsically disordered proteins (IDPs) or proteins where
secondary structure is not well-conserved — uses contact map from native state.

```bash
# Go-MARTINI (Poma et al. 2017, updated 2022)
# Requires: martinize2 + go-martini scripts
# https://github.com/marrink-lab/vermouth-martinize

martinize2 \
  -f protein.pdb \
  -o protein_go.itp \
  -x protein_go.pdb \
  -ff martini3001 \
  -go protein.pdb \          # native structure for Go contacts
  -go-eps 9.414 \            # contact strength (kJ/mol) — default
  -go-res-dist 3 \           # exclude i,i+1,i+2 contacts
  -go-dist-cutoff 0.6        # contact cutoff (nm)
```

## GROMACS CG Protein Simulation

```bash
# 1. Generate CG topology with martinize2
martinize2 -f protein.pdb -o protein.itp -x protein_cg.pdb \
           -ff martini3001 -elastic -ef 500 -el 0.5 -eu 0.9

# 2. Download MARTINI 3 force field files
# https://cgmartini.nl → Download → MARTINI 3
# Unzip to working directory: martini_v3.0.0.itp, water, ions etc.

# 3. Create system.top
cat > system.top << 'EOF'
#include "martini_v3.0.0.itp"
#include "martini_v3.0.0_solvents_v1.itp"
#include "martini_v3.0.0_ions_v1.itp"
#include "protein.itp"

[ system ]
CG Protein in Water

[ molecules ]
Protein    1
W          5000
NA         20
CL         20
EOF

# 4. Solvate with MARTINI water (W beads)
gmx solvate -cp protein_cg.pdb -cs water.gro -p system.top -o solvated.gro

# 5. Add ions
gmx grompp -f ions.mdp -c solvated.gro -p system.top -o ions.tpr
gmx genion -s ions.tpr -o ionized.gro -p system.top -pname NA -nname CL -neutral

# 6. Energy minimization
gmx grompp -f em.mdp -c ionized.gro -p system.top -o em.tpr
gmx mdrun -v -deffnm em

# 7. Equilibration + production
gmx grompp -f eq.mdp -c em.gro -p system.top -o eq.tpr
gmx mdrun -v -deffnm eq

gmx grompp -f md.mdp -c eq.gro -p system.top -o md.tpr
gmx mdrun -v -deffnm md -ntmpi 1 -ntomp 4 -gpu_id 0
```

## MARTINI 3 mdp Parameters for Proteins

```ini
; em.mdp (energy minimization)
integrator = steep
nsteps = 10000
emtol = 100.0
nstxout = 100

; CG non-bonded parameters
cutoff-scheme = Verlet
nstlist = 20
ns_type = grid
rlist = 1.2
rcoulomb = 1.1
rvdw = 1.1
vdw-modifier = Potential-shift
coulombtype = reaction-field
epsilon_r = 15
epsilon_rf = 0

; eq.mdp / md.mdp
integrator = md
dt = 0.020                ; 20 fs timestep
nsteps = 50000000         ; 1 µs (×4 = 4 µs effective)

; Thermostat (v-rescale mandatory for MARTINI)
tcoupl = v-rescale
tc-grps = Protein W_and_ions
tau_t = 1.0 1.0
ref_t = 310 310

; Barostat
pcoupl = Parrinello-Rahman
tau_p = 12.0
compressibility = 3e-4
ref_p = 1.0
```

## OpenMM MARTINI 3

```python
import openmm as mm
import openmm.app as app
from openmm import unit

def setup_martini_openmm(gro_file, top_file, platform='CUDA'):
    """
    Load MARTINI system in OpenMM.
    Note: use GromacsGroFile + GromacsTopFile.
    """
    gro = app.GromacsGroFile(gro_file)
    top = app.GromacsTopFile(
        top_file,
        periodicBoxVectors=gro.getPeriodicBoxVectors(),
        includeDir='/path/to/martini3_itp/'
    )

    system = top.createSystem(
        nonbondedMethod=app.CutoffPeriodic,
        nonbondedCutoff=1.1 * unit.nanometer,
        constraints=None,    # no constraints in CG
        removeCMMotion=True
    )

    # Reaction-field electrostatics (MARTINI standard)
    for force in system.getForces():
        if isinstance(force, mm.NonbondedForce):
            force.setReactionFieldDielectric(15.0)

    # v-rescale thermostat (approximate with Langevin for OpenMM)
    integrator = mm.LangevinMiddleIntegrator(
        310 * unit.kelvin,
        1.0 / unit.picosecond,   # friction ~ v-rescale τ
        0.020 * unit.picoseconds  # 20 fs
    )

    platform_obj = mm.Platform.getPlatformByName(platform)
    sim = app.Simulation(top.topology, system, integrator, platform_obj)
    sim.context.setPositions(gro.positions)
    return sim
```

## Common martinize2 Errors

| Error | Cause | Fix |
|-------|-------|-----|
| "Could not read DSSP" | DSSP not installed/found | `pip install dssp` or `conda install dssp` |
| "Atom not recognized" | Non-standard residue | Use `-ignore` flag or patch .pdb |
| Broken secondary structure | Missing residues in PDB | Fill with modeller/pdbfixer first |
| Elastic bonds too few | el/eu cutoffs too tight | Increase eu to 1.0–1.2 nm |
| Protein unfolds quickly | Insufficient elastic network | Try Go-MARTINI or increase ef |
