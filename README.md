# ALKYL

```
 ░▒▓██████▓▒░░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░
░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░
░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░
░▒▓████████▓▒░▒▓█▓▒░      ░▒▓███████▓▒░ ░▒▓██████▓▒░░▒▓█▓▒░
░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░  ░▒▓█▓▒░   ░▒▓█▓▒░
░▒▓█▓▒░░▒▓█▓▒░▒▓█▓▒░      ░▒▓█▓▒░░▒▓█▓▒░  ░▒▓█▓▒░   ░▒▓█▓▒░
░▒▓█▓▒░░▒▓█▓▒░▒▓████████▓▒░▒▓█▓▒░░▒▓█▓▒░  ░▒▓█▓▒░   ░▒▓████████▓▒░
```

**Computational Chemistry plugin for [Claude Code](https://claude.ai/code).**

Adds molecular modeling, quantum chemistry, and cheminformatics context to every Claude Code session — no separate command, no fork, just `claude`.

## Install

```bash
git clone https://github.com/YOUR_USERNAME/alkyl
cd alkyl
bash install.sh
```

That's it. Open a new `claude` session and ALKYL is active.

## What it adds

- **Identity** — Claude introduces itself as ALKYL, a computational chemistry assistant
- **Chemistry context** — IUPAC nomenclature, SMILES notation, computational cost awareness
- **Skills** — RDKit patterns, cheminformatics workflows, molecular representations

## Stack

`RDKit` · `ASE` · `ORCA` · `Gaussian` · `OpenBabel` · `py3Dmol` · `MDAnalysis` · `DeepChem`

## Uninstall

```bash
bash uninstall.sh
```

## Development

Tests require `rdkit`. Run them from a virtual environment that has it installed:

```bash
# Option A: use the Pepflex venv if available
/home/de/Bureau/Pepflex/venv/bin/python -m pytest tests/ -m "not network" -v

# Option B: create a local venv
python3 -m venv .venv
source .venv/bin/activate
pip install rdkit pytest
python -m pytest tests/ -m "not network" -v
```

## Project structure

```
ALKYL/
├── install.sh          # injects chemistry context into ~/.claude/CLAUDE.md
├── uninstall.sh        # removes the injected block
├── config/
│   └── CLAUDE.md       # chemistry identity and behavior
└── skills/
    ├── rdkit.md        # RDKit patterns and usage
    └── cheminformatics.md  # molecular representations & workflows
```
