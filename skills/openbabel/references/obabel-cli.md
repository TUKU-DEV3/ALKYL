# OpenBabel — obabel Command-Line Tool

## Basic Syntax

```bash
obabel [-i <fmt>] infile [-o <fmt>] [-O outfile] [options]

# Format auto-detected from extension when not specified
obabel mol.sdf -O mol.pdb
obabel mol.smi -O mol.mol2

# Explicit format specification
obabel -isdf input.sdf -opdb -O output.pdb
```

---

## Format Conversion Examples

```bash
# SMILES → SDF (no 3D)
obabel smiles.smi -O out.sdf

# SMILES → SDF with 3D
obabel smiles.smi -O out.sdf --gen3d

# SDF → SMILES (canonical)
obabel library.sdf -O library.smi -ocan

# SDF → MOL2
obabel library.sdf -O library.mol2

# PDB → SDF
obabel protein.pdb -O protein.sdf

# CIF → XYZ (crystal → Cartesian)
obabel structure.cif -O structure.xyz

# InChI → SMILES
echo "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3" | obabel -iinchi -osmi

# SMILES from stdin
echo "CCO" | obabel -ismi -osdf -O ethanol.sdf --gen3d

# Multiple input files → single output
obabel mol1.sdf mol2.sdf mol3.sdf -O merged.sdf

# All SDF files in directory → single SDF
obabel *.sdf -O all.sdf
```

---

## 3D Coordinate Generation

```bash
# --gen3d: SMILES → 3D (default quality: 'med')
obabel smiles.smi -O out.sdf --gen3d

# Quality levels (slowest = best geometry)
obabel smiles.smi -O out.sdf --gen3d fastest   # no cleanup
obabel smiles.smi -O out.sdf --gen3d fast       # 100 FF cycles
obabel smiles.smi -O out.sdf --gen3d med        # + fast rotor search (default)
obabel smiles.smi -O out.sdf --gen3d slow       # better rotor search
obabel smiles.smi -O out.sdf --gen3d best       # slowest, best quality

# Force field selection
obabel smiles.smi -O out.sdf --gen3d --ff mmff94   # default for organics
obabel smiles.smi -O out.sdf --gen3d --ff uff      # for unusual atoms
obabel smiles.smi -O out.sdf --gen3d --ff gaff     # AMBER GAFF
```

---

## Conformer Searching

```bash
# Generate conformers from existing 3D structure
obabel mol.sdf -O confs.sdf --conformer --nconf 50 --score energy

# From SMILES: gen3D first, then conformers
obabel mol.smi -O confs.sdf --gen3d --conformer --nconf 20 --score rmsd

# Parameters:
#   --nconf N          number of conformers to generate
#   --score energy     rank by force field energy (default)
#   --score rmsd       rank by RMSD diversity
#   --weighted         weighted rotor search (slower, better)
#   --ff <name>        force field (mmff94, uff, gaff)
#   --log              print energy log to stderr
```

---

## Hydrogen Management

```bash
# Add all explicit H
obabel mol.sdf -O mol_H.sdf -h

# Remove all H (implicit only)
obabel mol_H.sdf -O mol.sdf -d

# Add H for specific pH (protonation state)
obabel mol.sdf -O mol_ph74.sdf -p 7.4    # physiological pH
obabel mol.sdf -O mol_ph5.sdf  -p 5.0    # acidic pH

# Add only polar H (N, O, S)
obabel mol.sdf -O mol.sdf --addpolarh
```

---

## Filtering Libraries

```bash
# Filter by descriptor values (MW < 500 AND logP < 5)
obabel library.sdf -O ro5.sdf --filter "MW<500 logP<5"

# Multiple conditions (AND implicit)
obabel library.sdf -O filtered.sdf \
  --filter "MW<500 logP>=0 logP<=5 TPSA<140 HBD<=5 HBA<=10"

# OR condition
obabel library.sdf -O out.sdf --filter "MW<200 || logP>4"

# Negate with !
obabel library.sdf -O out.sdf --filter "!logP>5"

# SMARTS substructure filter (keep matching)
obabel library.sdf -O amines.sdf -s '[NH2]'

# SMARTS exclusion filter (remove matching)
obabel library.sdf -O clean.sdf -v '[CH]=O'     # remove aldehydes

# Range selection (molecules 10-50)
obabel library.sdf -O subset.sdf -f 10 -l 50

# Remove duplicates (by InChI)
obabel library.sdf -O dedup.sdf --unique

# Remove duplicates by SMILES
obabel library.sdf -O dedup.sdf --unique cansmi
```

---

## Output Splitting

```bash
# Split multi-molecule SDF into individual files
obabel library.sdf -O mol.sdf -m
# → mol1.sdf, mol2.sdf, mol3.sdf, ...

# Split to MOL files named by title
obabel library.sdf -O %.mol -m   # % = molecule title

# Batch convert all XYZ → PDB individually
obabel *.xyz -opdb -m
```

---

## 2D Coordinate Generation

```bash
# --gen2d: generate 2D layout coordinates
obabel mol.smi -O mol.sdf --gen2d

# Useful for 2D depiction, before writing SVG/PNG
obabel mol.smi -O mol.svg --gen2d
obabel mol.smi -O mol.png --gen2d
```

---

## InChI and InChIKey

```bash
# SMILES → InChI
echo "CCO" | obabel -ismi -oinchi
# InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3

# SMILES → InChIKey
echo "CCO" | obabel -ismi -oinchikey
# LFQSCWFLJHTTHZ-UHFFFAOYSA-N

# InChI → SMILES
echo "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3" | obabel -iinchi -ocan
# CCO
```

---

## Useful Queries

```bash
# List all supported input formats
obabel -L formats read

# List all output formats
obabel -L formats write

# Get format-specific options
obabel -H sdf

# List available descriptors
obabel -L descriptors

# List force fields
obabel -L forcefields

# Verbose output (useful for debugging)
obabel mol.sdf -O mol.pdb -v
```

---

## Calling obabel from Python

```python
import subprocess

def convert(input_file, output_file, extra_args=None):
    cmd = ['obabel', input_file, '-O', output_file]
    if extra_args:
        cmd.extend(extra_args)
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"obabel error: {result.stderr}")
    return result.stdout

# Examples
convert('input.sdf', 'output.pdb')
convert('smiles.smi', 'out.sdf', ['--gen3d', 'slow'])
convert('library.sdf', 'filtered.sdf', ['--filter', 'MW<500 logP<5'])

# SMILES string to SDF file
def smiles_to_sdf(smiles, output_path, gen3d=True):
    cmd = ['obabel', '-ismi', '-osdf', f'-O{output_path}']
    if gen3d:
        cmd += ['--gen3d']
    result = subprocess.run(cmd, input=smiles, capture_output=True, text=True)
    return output_path
```
