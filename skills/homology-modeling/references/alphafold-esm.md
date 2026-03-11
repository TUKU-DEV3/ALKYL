# Homology Modeling — AlphaFold2 / ColabFold / ESMFold

AI-based structure prediction: no template required, captures remote homologs and de novo folds.

---

## ColabFold — Local CLI (Recommended)

ColabFold wraps AlphaFold2 with faster MMseqs2 MSA generation (100× faster than genetic databases).

### Installation

```bash
# Via pip (GPU required for inference)
pip install "colabfold[alphafold]"

# Via conda (more reproducible)
conda install -c conda-forge -c bioconda colabfold

# Download weights once (~3 GB)
python -m colabfold.download
```

### Basic prediction

```bash
# Single sequence (FASTA file)
colabfold_batch target.fasta output_dir/ \
    --num-models 5 \
    --num-recycle 3 \
    --amber \
    --use-gpu-relax

# Multiple sequences (one per line or multi-FASTA)
colabfold_batch sequences.fasta output_dir/ \
    --num-models 3 \
    --num-recycle 6 \
    --amber

# Complex / multimer
colabfold_batch complex.fasta output_dir/ \
    --num-models 5 \
    --num-recycle 3 \
    --model-type alphafold2_multimer_v3 \
    --amber
```

### Key options

| Flag | Default | Notes |
|------|---------|-------|
| `--num-models` | 5 | Number of model parameters to use (1–5) |
| `--num-recycle` | 3 | Recycling iterations; 6–12 for harder targets |
| `--amber` | False | AMBER relaxation after prediction (strongly recommended) |
| `--use-gpu-relax` | False | GPU-accelerated AMBER relax |
| `--msa-mode` | `mmseqs2_uniref_env` | Use `single_sequence` for ESMFold-like (no MSA) |
| `--model-type` | `auto` | `alphafold2_ptm` (monomer) or `alphafold2_multimer_v3` |
| `--templates` | True | Use PDB templates; `--no-templates` for ab initio |
| `--max-msa` | 512:1024 | Reduce for faster/less accurate runs |

### Output files

```
output_dir/
  target_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_000.pdb
  target_relaxed_rank_001_alphafold2_ptm_model_1_seed_000.pdb  ← USE THIS
  target_scores_rank_001_alphafold2_ptm_model_1_seed_000.json
  target_coverage.png    ← MSA depth per residue
  target_PAE.png         ← Predicted Aligned Error matrix
  target_plddt.png       ← pLDDT per residue
```

---

## Parsing ColabFold / AlphaFold2 Output

### Load best model

```python
from pathlib import Path
import json

def load_best_af2_model(output_dir: str) -> tuple:
    """Find and load rank_001 relaxed PDB + scores."""
    out = Path(output_dir)

    # Find rank_001 relaxed PDB
    pdb_files = sorted(out.glob("*relaxed_rank_001*.pdb"))
    if not pdb_files:
        pdb_files = sorted(out.glob("*unrelaxed_rank_001*.pdb"))
    if not pdb_files:
        raise FileNotFoundError(f"No rank_001 PDB in {output_dir}")

    pdb_path = pdb_files[0]

    # Load matching scores JSON
    stem = pdb_path.stem.replace("relaxed_", "").replace("unrelaxed_", "")
    score_file = out / f"{stem.replace('rank_001_', 'scores_rank_001_')}.json"
    # Fallback: find any scores_rank_001 json
    if not score_file.exists():
        score_files = sorted(out.glob("*scores_rank_001*.json"))
        score_file = score_files[0] if score_files else None

    scores = {}
    if score_file and score_file.exists():
        with open(score_file) as f:
            scores = json.load(f)

    return str(pdb_path), scores

pdb_path, scores = load_best_af2_model("output_dir/")
print(f"Best model: {pdb_path}")
print(f"Mean pLDDT: {sum(scores.get('plddt', [0])) / max(len(scores.get('plddt', [1])), 1):.1f}")
```

### Extract pLDDT per residue (from B-factor column)

AlphaFold2 stores pLDDT in the B-factor column of the PDB.

```python
from Bio.PDB import PDBParser
import numpy as np
import matplotlib.pyplot as plt

def extract_plddt(pdb_file: str) -> tuple:
    """Extract per-residue pLDDT from AF2 PDB (stored in B-factor)."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("model", pdb_file)

    residue_nums = []
    plddt_vals   = []

    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] != ' ':
                    continue  # skip HETATM
                # Take CA atom bfactor as pLDDT
                if 'CA' in residue:
                    plddt = residue['CA'].get_bfactor()
                    residue_nums.append(residue.id[1])
                    plddt_vals.append(plddt)
        break  # first model only

    return np.array(residue_nums), np.array(plddt_vals)

resids, plddt = extract_plddt("best_model_relaxed.pdb")

# Plot pLDDT
plt.figure(figsize=(12, 4))
plt.plot(resids, plddt, linewidth=1)
plt.fill_between(resids, plddt, alpha=0.3)
plt.axhline(90, color='green',  linestyle='--', label='High (>90)')
plt.axhline(70, color='orange', linestyle='--', label='Medium (>70)')
plt.axhline(50, color='red',    linestyle='--', label='Low (<50)')
plt.xlabel('Residue'); plt.ylabel('pLDDT')
plt.title('Per-residue confidence (pLDDT)')
plt.legend(); plt.tight_layout()
plt.savefig('plddt_profile.png', dpi=150)

# Flag low-confidence regions
low_conf = resids[plddt < 70]
print(f"Low-confidence residues (pLDDT < 70): {len(low_conf)}")
print(f"Disordered residues (pLDDT < 50): {(plddt < 50).sum()}")
```

### pLDDT interpretation

| pLDDT range | Confidence | Meaning |
|-------------|-----------|---------|
| > 90 | Very high | Confident, use for docking |
| 70–90 | High | Generally reliable backbone |
| 50–70 | Low | Flexible/uncertain, check carefully |
| < 50 | Very low | Likely disordered — DO NOT dock here |

**For docking:** require pLDDT > 80 in the binding pocket residues.

---

## PAE (Predicted Aligned Error) Matrix

PAE quantifies inter-residue / inter-domain position confidence.

```python
import json
import numpy as np
import matplotlib.pyplot as plt

def plot_pae(scores_json: str, out_png: str = "pae.png") -> np.ndarray:
    """Load and plot PAE matrix from ColabFold scores JSON."""
    with open(scores_json) as f:
        scores = json.load(f)

    # ColabFold key: 'pae' (N x N matrix of floats)
    pae = np.array(scores['pae'])

    fig, ax = plt.subplots(figsize=(8, 7))
    im = ax.imshow(pae, vmin=0, vmax=30, cmap='Greens_r')
    plt.colorbar(im, ax=ax, label='PAE (Å)')
    ax.set_xlabel('Scored residue')
    ax.set_ylabel('Aligned residue')
    ax.set_title('Predicted Aligned Error')
    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    return pae

pae = plot_pae("scores_rank_001.json", "pae.png")

# Inter-domain confidence: check off-diagonal blocks
# Low PAE (dark) between domains = reliable relative orientation
# High PAE (light) = uncertain relative placement (e.g., flexible linker)
n = pae.shape[0]
domain_boundary = 150   # example: two domains split at residue 150
d1_d2_pae = pae[:domain_boundary, domain_boundary:].mean()
print(f"Domain1–Domain2 mean PAE: {d1_d2_pae:.1f} Å")
print(f"  < 5 Å: reliable interface  |  > 15 Å: uncertain orientation")
```

---

## ESMFold — Single-Sequence Prediction

ESMFold uses an ESM-2 language model instead of MSA. Faster, but less accurate for remote homologs.

### Python API

```python
import torch
import esm

# Load model (downloads ~2.5 GB on first run, cached in ~/.cache/torch/hub)
model = esm.pretrained.esmfold_v1()
model = model.eval()

# Move to GPU if available
device = 'cuda' if torch.cuda.is_available() else 'cpu'
model = model.to(device)

# Single sequence prediction
sequence = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVK"

with torch.no_grad():
    output = model.infer_pdb(sequence)

with open("esmfold_model.pdb", "w") as f:
    f.write(output)
print("ESMFold model saved.")
```

### Batch prediction

```python
import torch, esm
from pathlib import Path

def esmfold_batch(sequences: dict, out_dir: str) -> None:
    """
    Predict structures for multiple sequences.
    sequences: {name: amino_acid_sequence}
    """
    Path(out_dir).mkdir(exist_ok=True)

    model = esm.pretrained.esmfold_v1()
    model = model.eval().to('cuda' if torch.cuda.is_available() else 'cpu')

    for name, seq in sequences.items():
        # ESMFold has a practical limit ~1000 residues before OOM on 16 GB GPU
        if len(seq) > 1000:
            print(f"Warning: {name} has {len(seq)} residues — may OOM")
        with torch.no_grad():
            pdb_str = model.infer_pdb(seq)
        out_path = Path(out_dir) / f"{name}.pdb"
        out_path.write_text(pdb_str)
        print(f"Saved: {out_path}")

seqs = {
    "protein_A": "MKTAYIAKQRQISFVKSHFSRQ...",
    "protein_B": "MSKLFRQIAVKHLRNLEEQARG...",
}
esmfold_batch(seqs, "esmfold_output/")
```

### Extract pLDDT from ESMFold output

ESMFold also stores pLDDT in B-factor column — same extraction code as AlphaFold2 applies.

---

## AF2 vs ESMFold — When to Use Each

| Criterion | AlphaFold2 / ColabFold | ESMFold |
|-----------|------------------------|---------|
| Accuracy (general) | ★★★★★ | ★★★ |
| Accuracy (orphan proteins, no homologs) | ★★★★ | ★★★ |
| Speed | 10–60 min / protein | < 1 min / protein |
| MSA required | Yes (MMseqs2, automatic) | No |
| GPU VRAM | 10 GB (monomer) | 15 GB (full model) |
| Multimer support | Yes (multimer_v3) | No |
| Batch throughput | Moderate | High |
| Best use case | Production modeling | Rapid screening, mutant scan |

**Rule of thumb:** use ESMFold for initial screening or when MSA fails; use ColabFold for final model.

---

## LocalColabFold / OpenFold for Reproducibility

```bash
# LocalColabFold: installs AF2 weights locally, no internet required during inference
# Install: https://github.com/YoshitakaMo/localcolabfold
./install_colabfold.sh /opt/localcolabfold

# Run (same CLI as colabfold_batch)
/opt/localcolabfold/colabfold_batch/bin/colabfold_batch \
    target.fasta output_dir/ \
    --num-models 5 --num-recycle 3 --amber

# OpenFold: JAX-free PyTorch reimplementation, reproducible training
# pip install openfold  (requires separate weight download)
```

---

## Common Pitfalls

| Problem | Symptom | Fix |
|---------|---------|-----|
| pLDDT < 50 in binding pocket | Unreliable docking results | Use MODELLER with template for that region |
| OOM on ESMFold | CUDA out of memory | Reduce sequence length or use CPU (slow) |
| ColabFold slow MSA | Hang at MSA step | Add `--msa-mode single_sequence` for testing |
| Multimer chains wrong order | Chain A/B swapped | Specify chains in FASTA: `>A\nSEQ_A\n>B\nSEQ_B` |
| relaxed PDB missing | `--amber` not set or AMBER failed | Check AMBER install: `conda install -c conda-forge openmm` |
| Wrong residue numbering | AF2 renumbers from 1 | Use `--pdb-id` option or renumber after with Biopython |
