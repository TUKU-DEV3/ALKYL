# Homology Modeling — MODELLER 10.x

Comparative homology modeling: build a 3D model from a target sequence using one or more structurally known templates.

---

## Prerequisites and License

```bash
# Free academic license required: https://salilab.org/modeller/registration.html
conda install -c salilab modeller

# Set key (also add to ~/.bashrc)
export KEY_MODELLER="XXXXXXXX"

# Verify
python -c "from modeller import Environ; print('MODELLER OK')"
```

---

## PIR Alignment Format (Critical)

MODELLER requires a strict PIR format. Any deviation causes a segmentation fault or silent wrong output.

```
>P1;TEMPLATE_CODE
structureX:TEMPLATE_CODE:  1 :A: 300:A:protein:organism: 2.0: 0.20
MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVK*

>P1;TARGET_SEQ
sequence:TARGET_SEQ:   :  :   :  :protein:organism:  :
MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVK*
```

**PIR record type fields** (colon-separated after `structureX:`):
1. File/code name (must match the `.pdb` filename without extension)
2. First residue number
3. First chain ID
4. Last residue number
5. Last chain ID
6. Protein type
7. Source organism
8. Resolution
9. R-factor

**Critical rules:**
- Template line uses `structureX:` (known 3D structure), target uses `sequence:`
- Each sequence MUST end with `*`
- Gap character is `-` in aligned region
- No trailing spaces on sequence lines
- Code names are case-sensitive and must exactly match PDB file names on disk

**Common pitfall — segfault from bad PIR:**
```python
# Before running, validate PIR:
from modeller import Environ, Alignment

env = Environ()
aln = Alignment(env)
aln.append(file='alignment.pir', align_codes='ALL')
aln.write(file='alignment_check.pir')  # if this succeeds, PIR is valid
```

---

## Basic Single-Template Workflow

```python
from modeller import Environ
from modeller.automodel import AutoModel

# 1. Initialize environment
env = Environ()
env.io.atom_files_directory = ['.', './templates']
# Suppress HETATM records from template (remove cofactors, ligands)
env.io.hetatm = False     # True = keep HETATM (cofactors, metals)
env.io.water  = False     # True = keep water molecules

# 2. Build models
a = AutoModel(env,
              alnfile  = 'alignment.pir',
              knowns   = '1ABC',          # template PDB code (file: 1ABC.pdb)
              sequence = 'MY_TARGET')     # must match >P1;MY_TARGET in .pir

a.starting_model = 1   # model index start
a.ending_model   = 10  # generate 10 models (more = better sampling)

# Optional: increase optimization cycles for better geometry
a.md_level = None          # no MD (fast)
# a.md_level = refine.slow  # slow MD (better but 5× slower)

a.make()

# 3. Select best model by DOPE score (lower = better)
results = [m for m in a.outputs if m['failure'] is None]
results.sort(key=lambda m: m['DOPE score'])
best = results[0]
print(f"Best model: {best['name']}")
print(f"DOPE score: {best['DOPE score']:.2f}")
# Copy/rename best model
import shutil
shutil.copy(best['name'], 'best_model.pdb')
```

---

## Loop Refinement

Loops (missing/poorly conserved regions) have poor geometry from automodel. Refine them explicitly.

```python
from modeller import Environ
from modeller.automodel import LoopModel, refine

env = Environ()
env.io.atom_files_directory = ['.', './templates']

class MyLoop(LoopModel):
    def select_loop_atoms(self):
        """Define loop regions to refine (residue ranges in target)."""
        from modeller import Selection
        return Selection(
            self.residue_range('45:A', '52:A'),   # loop 1
            self.residue_range('120:A', '131:A'), # loop 2
        )

a = MyLoop(env,
           alnfile      = 'alignment.pir',
           knowns       = '1ABC',
           sequence     = 'MY_TARGET',
           loop_assess_methods = 'DOPE')  # assess loops by DOPE

a.starting_model      = 1
a.ending_model        = 5
a.loop.starting_model = 1
a.loop.ending_model   = 10  # 10 loop conformations per initial model
a.loop.md_level       = refine.fast   # refine.slow for production

a.make()

# Select best loop model
loop_results = [m for m in a.loop.outputs if m['failure'] is None]
loop_results.sort(key=lambda m: m['DOPE score'])
print(f"Best loop model: {loop_results[0]['name']}")
```

---

## Multi-Template Modeling

Use multiple templates to cover different regions or improve accuracy.

```python
from modeller import Environ
from modeller.automodel import AutoModel

env = Environ()
env.io.atom_files_directory = ['.', './templates']

# Multiple templates: list them all in knowns
a = AutoModel(env,
              alnfile  = 'multi_alignment.pir',
              knowns   = ('1ABC', '2XYZ', '3DEF'),   # tuple of template codes
              sequence = 'MY_TARGET')

a.starting_model = 1
a.ending_model   = 10
a.make()
```

**Multi-template PIR format:**
```
>P1;1ABC
structureX:1ABC:  1:A:250:A:::::
MKTAYIAKQRQISFVKSHFSRQLEERLGL---------SRVGDGTQDNLSGAEKAVQVK*

>P1;2XYZ
structureX:2XYZ:  1:A:280:A:::::
-----------QRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK----*

>P1;MY_TARGET
sequence:MY_TARGET::::::::
MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVK*
```

**Template weighting** (downweight a lower-quality template):
```python
from modeller import Environ
from modeller.automodel import AutoModel
from modeller import restraints

class WeightedModel(AutoModel):
    def special_restraints(self, aln):
        rsr = self.restraints
        # Downweight 3DEF template to 50% contribution
        at = aln.positions
        # (Advanced: modify weight matrix — typically not needed for > 30% identity)

# For most cases, MODELLER weights templates automatically based on sequence similarity
```

---

## Custom Restraints

```python
from modeller import Environ
from modeller.automodel import AutoModel
from modeller import physical

class RestrainedModel(AutoModel):
    def special_restraints(self, aln):
        rsr = self.restraints
        at = self.atoms

        # Distance restraint between two atoms (e.g., catalytic dyad)
        # atom1 = 'NE2:57:A:HIS'  (atom_name:resnum:chain:resname)
        rsr.add(restraints.forms.gaussian(
            group    = physical.xy_distance,
            feature  = features.distance(at['NE2:57:A'], at['OD1:102:A']),
            mean     = 3.2,   # target distance in Angstrom
            stdev    = 0.3
        ))

        # Dihedral restraint (e.g., fix a known helix)
        rsr.add(restraints.forms.gaussian(
            group   = physical.phi_dihedral,
            feature = features.dihedral(
                at['C:10:A'], at['N:11:A'], at['CA:11:A'], at['C:11:A']),
            mean    = -60.0,  # degrees
            stdev   = 10.0
        ))
```

---

## DOPE Score Assessment

```python
from modeller import Environ, Selection
from modeller.scripts import complete_pdb

env = Environ()
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')

# Score a single model file
mdl = complete_pdb(env, 'model.B99990001.pdb')
sel = Selection(mdl)
mdl.assess_dope(output='ENERGY', file='model_dope.profile',
                assess_method=mdl.assess.DOPE)

# DOPE per-residue profile → identify poorly modeled regions
# Read profile
import numpy as np
residues, dope_vals = [], []
with open('model_dope.profile') as f:
    for line in f:
        if line.startswith('#') or not line.strip():
            continue
        parts = line.split()
        if len(parts) >= 2:
            residues.append(int(parts[0]))
            dope_vals.append(float(parts[1]))

# Peaks (high DOPE = poor geometry) → loop regions to refine
import matplotlib.pyplot as plt
plt.plot(residues, dope_vals)
plt.axhline(y=-0.1, color='r', linestyle='--', label='Threshold')
plt.xlabel('Residue'); plt.ylabel('DOPE score')
plt.title('Per-residue DOPE profile')
plt.savefig('dope_profile.png', dpi=150)
```

---

## HETATM and Cofactor Issues

```python
# Keep metal cofactor in template (e.g., Zn2+ in active site)
env.io.hetatm = True   # keep all HETATM
# Then manually remove unwanted ligands from template PDB before running

# Alternatively, keep only specific HETATM by residue name:
# Pre-filter the template PDB:
def keep_hetatm(template_pdb, out_pdb, keep_resnames=('ZN', 'MG', 'CA')):
    """Keep only specified HETATM records in template."""
    with open(template_pdb) as f, open(out_pdb, 'w') as out:
        for line in f:
            record = line[:6].strip()
            if record == 'HETATM':
                resname = line[17:20].strip()
                if resname in keep_resnames:
                    out.write(line)
            else:
                out.write(line)

keep_hetatm('1ABC.pdb', '1ABC_clean.pdb', keep_resnames=('ZN',))
```

---

## Common Errors and Fixes

| Error | Cause | Fix |
|-------|-------|-----|
| Segmentation fault | Malformed PIR file | Validate with `aln.write()` before `make()` |
| `IOError: 'XXXX.pdb' not found` | PDB file not in `atom_files_directory` | Add path to `env.io.atom_files_directory` list |
| `ModellerError: HETATM record` | HETATM in template when `io.hetatm=False` | Set `env.io.hetatm = True` or pre-clean template |
| `failure` in output dict | Optimization diverged | Reduce model complexity or increase `a.ending_model` |
| Huge DOPE score (> 0) | Very bad geometry (clashes) | Check alignment for frame-shifts; try loop refinement |
| Missing chain ID in PIR | Chain not specified | Always set chain ID in PIR structureX record |
