# ALKYL Scripts — Design Document

**Date:** 2026-03-01
**Status:** Approved
**Scope:** Pre-built chemistry scripts invoked by Claude (ALKYL) via Bash tool

---

## Context

ALKYL currently provides: skills (rdkit, deepchem, synkit, torchdrug, nextflow, pepflex),
MCP servers (pubmed, chembl, arxiv), and install/uninstall infrastructure.

The missing layer: **executable scripts** that Claude can call directly during a session
to perform chemistry tasks without writing ad-hoc code each time.

---

## Design Decisions

| Decision | Choice | Rationale |
|---|---|---|
| Invocation mode | Claude calls via Bash tool | Scripts as tools for ALKYL, not CLI for users |
| Interface | Independent scripts with argparse | Simple, no dependencies between scripts |
| Output | JSON on stdout | Claude parses directly; structured, unambiguous |
| File output | `--out FILE` flag | Optional, when user needs a file artifact |
| Molecule input | `--smiles`, `--sdf`, `--inchi` | Flexible; SMILES is the default text form |
| Errors | stderr + non-zero exit code | Claude detects failure and reports clearly |
| Core dependency | RDKit (required) | ASE, cclib imported on demand |

---

## Architecture

```
ALKYL/
  scripts/
    chem_convert.py    # inter-format I/O
    chem_props.py      # physico-chemical descriptors + fingerprints
    chem_3d.py         # conformer generation, RMSD, alignment
    chem_qm.py         # QM input/output (ORCA, Gaussian)
    chem_fetch.py      # external data fetch (PubChem, ChEMBL)
```

Scripts are named `chem_*.py` — 5 names, one per functional axis. Each stays under ~300 lines.

---

## Script Specifications

### `chem_convert.py` — Format I/O

```
--smiles STR | --sdf FILE | --inchi STR | --mol2 FILE | --pdb FILE
--to [smiles|sdf|inchi|mol2|pdb|svg|png]
--sanitize / --no-sanitize
--out FILE    (default: stdout)
```

Output: converted molecule in requested format, or JSON metadata for text formats.

---

### `chem_props.py` — Descriptors & Drug-likeness

```
--smiles STR | --sdf FILE
--descriptors [all | mw,logp,hbd,hba,tpsa,rotbonds,rings,fsp3,...]
--lipinski           → pass/fail + violation list
--pains              → structural alert flags
--fingerprint [morgan|maccs|rdkit]  --radius 2  --nbits 2048
--similarity FILE.sdf --ref SMILES  → Tanimoto scores
```

Output JSON example:
```json
{"mw": 180.16, "logp": 1.19, "hbd": 1, "hba": 3, "tpsa": 63.6,
 "lipinski": {"pass": true, "violations": []}}
```

---

### `chem_3d.py` — Conformers & Geometry

```
--smiles STR | --sdf FILE
--conformers N  --method [etkdg|etdg]   → multi-conf SDF
--minimize [mmff94|uff]
--rmsd --ref conf.sdf                    → RMSD matrix
--align --ref REF.sdf                    → rigid alignment
```

Output: SDF file (via `--out`) or JSON with RMSD matrix / conformer energies.

---

### `chem_qm.py` — QM Interface (ORCA / Gaussian)

```
--smiles STR | --xyz FILE | --sdf FILE
--engine [orca|gaussian]
--task [sp|opt|freq|scan]
--method STR   (e.g. "B3LYP")
--basis STR    (e.g. "6-31G*")
--charge INT   --mult INT
--parse OUTPUT_FILE   → JSON: energies, frequencies, optimized geometry
--out FILE            → write input file
```

Output (generation): writes input file.
Output (parsing): JSON with energies, thermochemistry, geometry.

---

### `chem_fetch.py` — External Data

```
--source [pubchem|chembl]
--name STR | --smiles STR | --cid INT | --chembl-id STR
--properties [mw,logp,inchi,...]    → subset of fields
--out FILE.sdf / FILE.json
--batch FILE.txt                    → list of names/IDs, one per line
```

Output: SDF or JSON depending on `--out` extension; stdout JSON if no `--out`.

---

## CLAUDE.md Integration

A `## Scripts` section is added to `config/CLAUDE.md`:

```markdown
## Scripts disponibles

Scripts dans `$ALKYL_SCRIPTS/` (résolu à l'install).
Appelle via : `python $ALKYL_SCRIPTS/<script>.py`

| Tâche | Script |
|---|---|
| Convertir format moléculaire | chem_convert.py |
| Calculer MW, LogP, TPSA, fingerprints | chem_props.py |
| Vérifier Lipinski / PAINS | chem_props.py --lipinski --pains |
| Générer conformères 3D | chem_3d.py --conformers N |
| Préparer input ORCA/Gaussian | chem_qm.py --engine orca |
| Parser output QM | chem_qm.py --parse output.log |
| Récupérer molécule PubChem/ChEMBL | chem_fetch.py |

Règles :
- Parser le JSON stdout avant de répondre
- Si RDKit absent : signaler, ne pas halluciner les valeurs
- `--help` disponible sur chaque script
```

`install.sh` resolves the absolute path of `scripts/` and exports `$ALKYL_SCRIPTS`
in the injected block — Claude always knows the correct path on any machine.

---

## Priority Order

1. `chem_convert.py` + `chem_props.py` — cheminformatics de base
2. `chem_fetch.py` — données externes (dépend souvent des autres)
3. `chem_3d.py` — conformères
4. `chem_qm.py` — QM (dépendances optionnelles: ASE, cclib)

---

## Out of Scope

- GUI or web interface
- Script-to-script pipelines (Claude orchestrates, not scripts)
- Batch processing beyond `--batch` flag in chem_fetch
- MD simulation setup (future phase)
