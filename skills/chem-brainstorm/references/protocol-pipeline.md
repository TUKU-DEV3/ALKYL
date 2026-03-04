# Protocol — Pipeline Construction (Rigid)

Use when building a reproducible, multi-step computational workflow over a library or dataset.

---

## Stage 1 — Define Pipeline Goal

Answer:
- **Input**: what format? (SMILES list, SDF, CSV, ChEMBL query)
- **Output**: what do we want? (filtered set, ranked list, trained model, QM data)
- **Scale**: how many molecules? (< 1K / 1K–100K / > 100K)
- **Reproducibility**: one-off or recurring?

## Stage 2 — Pipeline Skeleton

Sketch the stages as a DAG:
```
fetch → standardize → filter → analyze → score → select → (ML or QM)
```

Map each stage to an ALKYL script:
```
fetch        → chem_fetch.py (PubChem/ChEMBL) or ChEMBL MCP
standardize  → chem_standardize.py
filter       → chem_filter.py (Lipinski/PAINS/Veber)
analyze      → chem_batch.py (bulk descriptors)
search       → chem_search.py (substructure/similarity)
diversity    → chem_diversity.py (MaxMin subset)
3D prep      → chem_3d.py (conformers)
score        → chem_metabolism.py / chem_pka.py / ML model
QM           → chem_qm.py (ORCA/Gaussian inputs)
```

## Stage 3 — Fetch Input Data ⚡

```
ChEMBL.compound_search(query="[target or scaffold]", limit=500)
→ starting set with known activity

ChEMBL.get_bioactivity(target_chembl_id="CHEMBLXXX", standard_type="IC50")
→ quantitative SAR dataset for ML

ZINC database (scientific-skills:zinc-database)
→ purchasable, 3D-ready compounds for docking
```

## Stage 4 — Standardize + Filter ⚡

```bash
# For each SMILES in input:
python chem_standardize.py --smiles "..." → canonical form
python chem_batch.py --input library.smi --descriptors all --lipinski --pains \
  --skip-invalid --out results.json
```

Checkpoints:
- Log n_valid vs n_total
- Flag PAINS hits separately (don't discard, annotate)
- Ghose filter requires ≥ 20 heavy atoms — skip for fragments

## Stage 5 — Diversity Selection ⚡ (if library too large)

```bash
python chem_diversity.py --input filtered.smi --n 1000 --fingerprint morgan
# MaxMin selection: representative diverse subset
```

## Stage 6 — Enrichment / Scoring ⚡–⚡⚡

Choose enrichment strategy:
- **Rule-based**: `chem_filter.py` + `chem_metabolism.py`
- **Similarity-based**: `chem_search.py --mode similarity --threshold 0.6`
- **ML-based**: `deepchem` skill (QSAR model on ChEMBL data)
- **Docking**: `scientific-skills:diffdock` or prepare with `chem_3d.py`

## Stage 7 — Automation (if recurring) ⚡⚡

For HPC or cloud deployment, load `nextflow` skill:
```
nextflow skill → references/pipeline-patterns.md
→ DSL2 workflow, containers, resume, parallel execution
```

For Python orchestration:
```python
# chem_batch.py handles bulk; chain with subprocess or snakemake
```

## Stage 8 — Output & Reporting

Outputs to always produce:
1. **Full results JSON** (`--out results.json`)
2. **Summary statistics** (n_input, n_passed, n_flagged, property distributions)
3. **Top-N list** (ranked by QED or activity prediction)
4. **Reproducibility info** (script versions, fingerprint type, thresholds used)

## Scale Guidance

| Library size | Approach |
|-------------|---------|
| < 100 | run scripts individually, inspect manually |
| 100–10K | `chem_batch.py`, in-process Python |
| 10K–1M | `chem_batch.py` + parallel chunks, or Nextflow |
| > 1M | Nextflow + HPC, or cloud (Modal/DnaNexus) |
