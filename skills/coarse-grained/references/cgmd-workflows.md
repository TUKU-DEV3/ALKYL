# CG MD Workflows and Backmapping

## Complete GROMACS MARTINI Workflow

```bash
#!/bin/bash
# Complete CG protein-membrane MD workflow

PREFIX="system"
FF_DIR="/path/to/martini3_itp"
T=310       # temperature (K)
P=1.0       # pressure (bar)

# ── 1. Build membrane with protein ─────────────────────────────────────────
insane \
  -f protein_cg.pdb \
  -l POPC:0.7,CHOL:0.3 \
  -u POPC:0.7,CHOL:0.3 \
  -sol W -salt 0.15 \
  -x 20 -y 20 -z 15 \
  -o ${PREFIX}.gro \
  -p ${PREFIX}.top

# ── 2. Energy minimization ─────────────────────────────────────────────────
gmx grompp -f em.mdp -c ${PREFIX}.gro -p ${PREFIX}.top -o em.tpr -maxwarn 5
gmx mdrun -v -deffnm em -ntmpi 1 -ntomp 4

# ── 3. Equilibration (NPT, restraints) ────────────────────────────────────
gmx grompp -f eq.mdp -c em.gro -p ${PREFIX}.top -o eq.tpr -r em.gro -maxwarn 5
gmx mdrun -v -deffnm eq -ntmpi 1 -ntomp 4

# ── 4. Production MD ───────────────────────────────────────────────────────
gmx grompp -f md.mdp -c eq.gro -p ${PREFIX}.top -o md.tpr -maxwarn 5
gmx mdrun -v -deffnm md -ntmpi 1 -ntomp 4 -gpu_id 0 \
          -nsteps 50000000    # 1 µs CG = 4 µs effective

# ── 5. Analysis ────────────────────────────────────────────────────────────
gmx trjconv -f md.xtc -s md.tpr -o md_whole.xtc -pbc whole -b 0
```

## mdp Files

```ini
; ── em.mdp ─────────────────────────────────────────────────────────────────
integrator = steep
nsteps = 10000
emtol = 100.0
emstep = 0.01

cutoff-scheme = Verlet
nstlist = 20
rlist = 1.2
rcoulomb = 1.1
rvdw = 1.1
vdw-modifier = Potential-shift
coulombtype = reaction-field
epsilon_r = 15
epsilon_rf = 0
```

```ini
; ── eq.mdp ─────────────────────────────────────────────────────────────────
integrator = md
dt = 0.010                  ; start conservative (10 fs)
nsteps = 100000

cutoff-scheme = Verlet
nstlist = 20
rlist = 1.2
rcoulomb = 1.1
rvdw = 1.1
vdw-modifier = Potential-shift
coulombtype = reaction-field
epsilon_r = 15
epsilon_rf = 0

tcoupl = v-rescale
tc-grps = Protein Membrane Solvent
tau_t = 1.0 1.0 1.0
ref_t = 310 310 310

pcoupl = Berendsen           ; use Berendsen for equilibration
tau_p = 3.0
compressibility = 3e-4
ref_p = 1.0
pcoupltype = semiisotropic   ; for membrane: x-y coupled, z independent

; Position restraints on protein backbone
define = -DPOSRES
```

```ini
; ── md.mdp ─────────────────────────────────────────────────────────────────
integrator = md
dt = 0.020                  ; 20 fs (MARTINI standard)
nsteps = 50000000           ; 1 µs CG

cutoff-scheme = Verlet
nstlist = 20
rlist = 1.2
rcoulomb = 1.1
rvdw = 1.1
vdw-modifier = Potential-shift
coulombtype = reaction-field
epsilon_r = 15
epsilon_rf = 0

; Output
nstxout = 0
nstvout = 0
nstfout = 0
nstxout-compressed = 5000   ; every 100 ps
nstenergy = 5000

tcoupl = v-rescale
tc-grps = Protein Membrane Solvent
tau_t = 1.0 1.0 1.0
ref_t = 310 310 310

pcoupl = Parrinello-Rahman   ; production: Parrinello-Rahman
tau_p = 12.0
compressibility = 3e-4
ref_p = 1.0
pcoupltype = semiisotropic
```

## Backmapping: CG → All-Atom

Reconstruct atomic detail from CG trajectory for subsequent AA MD.

### backward.py (Wassenaar et al. 2014)

```bash
# Install
git clone https://github.com/cgmartini/backward.git
cd backward

# Basic backmapping
python initram.sh \
  -f cg_frame.gro \          # CG structure (single frame)
  -o aa_backmapped.gro \     # AA output
  -to charmm36 \             # target AA force field
  -p system.top              # MARTINI topology

# Multiple frames (extract key frame first)
echo "0" | gmx trjconv \
  -f md.xtc -s md.tpr \
  -o frame_1us.gro \
  -b 1000000 -e 1000000    # extract at 1 µs
```

### Python Backmapping Wrapper

```python
import subprocess
from pathlib import Path

def backmap_cg_frame(cg_gro, cg_top, output_aa_gro,
                      target_ff='charmm36', backward_dir='/opt/backward'):
    """
    Backmap single CG frame to all-atom using backward.py.
    """
    cmd = [
        'python', f'{backward_dir}/initram.sh',
        '-f', str(cg_gro),
        '-o', str(output_aa_gro),
        '-to', target_ff,
        '-p', str(cg_top)
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Backmapping failed:\n{result.stderr}")
    return output_aa_gro

def extract_and_backmap(xtc_file, tpr_file, cg_top,
                         time_points_ns, output_dir='backmapped'):
    """Extract CG frames at specified times and backmap all."""
    Path(output_dir).mkdir(exist_ok=True)
    results = []

    for t_ns in time_points_ns:
        t_ps = t_ns * 1000
        gro_out = f'{output_dir}/frame_{t_ns}ns_cg.gro'
        aa_out = f'{output_dir}/frame_{t_ns}ns_aa.gro'

        # Extract frame
        subprocess.run([
            'gmx', 'trjconv',
            '-f', xtc_file, '-s', tpr_file,
            '-o', gro_out,
            '-b', str(t_ps), '-e', str(t_ps),
            '-dump', str(t_ps)
        ], input='0\n', capture_output=True, text=True)

        # Backmap
        try:
            backmap_cg_frame(gro_out, cg_top, aa_out)
            results.append({'time_ns': t_ns, 'aa_gro': aa_out, 'status': 'ok'})
        except Exception as e:
            results.append({'time_ns': t_ns, 'aa_gro': None, 'status': str(e)})

    return results
```

## Post-Backmapping Refinement

After backmapping, always minimize and briefly equilibrate the AA structure:

```python
def refine_backmapped_structure(aa_gro, ff='charmm36m'):
    """
    Minimize and briefly equilibrate backmapped AA structure.
    Uses OpenMM (see force-fields skill).
    """
    from openmm.app import GromacsGroFile, GromacsTopFile, Simulation
    from openmm.app import PME, HBonds
    from openmm import LangevinMiddleIntegrator, MonteCarloBarostat
    from openmm import unit
    import openmm as mm

    gro = GromacsGroFile(aa_gro)
    # Note: need AA topology; use pdbfixer + forcefield assignment
    # See force-fields skill for full setup

    # Staged minimization
    # 1. Minimize with heavy atom restraints (100 kJ/mol/nm²)
    # 2. Minimize without restraints
    # 3. Short NVT 100 ps at 100K → 310K
    # See force-fields skill: openmm-basics.md staged minimization
    pass
```

## Troubleshooting CG Simulations

| Problem | Likely cause | Fix |
|---------|-------------|-----|
| Crash at step 0 | Atom overlap in initial structure | Run energy minimization first; increase emtol |
| LINCS warnings | Bonds too stretched | Check insane output; reduce dt to 5 fs |
| Membrane buckles | Wrong box size or pressure coupling | Use semiisotropic pcoupl; check compressibility |
| Protein unfolds | Insufficient elastic network | Increase ef (500→700) or use Go-MARTINI |
| Protein leaves membrane | Wrong hydrophobic length or orientation | Verify TM helix orientation; adjust z offset in insane |
| Very slow run on CPU | MARTINI needs GPU | Use -gpu_id 0; MARTINI 3 has good GPU performance |
| Lipids flip-flop too fast | Normal in CG | Expected; reduce flip-flop with higher leaflet barriers if needed |

## Small Molecule CG Parameterization (MARTINI 3)

```python
# For drug-like small molecules: use CGSmiles + SMILES-based mapping
# pip install cgsmiles

import cgsmiles
from cgsmiles import MoleculeResolver

# Define CG mapping using CGSmiles notation
# Example: ibuprofen
cg_ibuprofen = "{[#SC5r6]1[#SC5r6][#SC5r6][#SC5r6][#SC5r6][#SC5r6]1[#SC3][#SC3]}.{[#SC3][#C3]=[#Na]}"

resolver = MoleculeResolver.from_string(cg_ibuprofen)
mol_graph = resolver.resolve()

# Or use automated parameterization via Polyply
# pip install polyply
# polyply gen_params -name IBU -smiles 'CC(C)Cc1ccc(cc1)C(C)C(=O)O' -o ibu.itp
```
