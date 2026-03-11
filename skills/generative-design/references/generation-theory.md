# Generation Theory — Molecular Space, Representations, and Evaluation

## Molecular Space

- **Drug-like chemical space**: ~10^60 molecules (Lipinski-compliant)
- **ZINC20**: ~1.4 billion purchasable, drug-like compounds
- **Practical explored space**: ~10^8–10^9 molecules from HTS history
- **Goal of generative models**: efficiently navigate this space guided by property oracles

## Molecular Representations for Generation

### SMILES (Simplified Molecular Input Line Entry System)

```
CC(=O)Oc1ccccc1C(=O)O    # aspirin
```

- Linear string over alphabet ~65 tokens
- Multiple valid SMILES per molecule → augmentation trick (randomized SMILES)
- **Problem**: ~10-40% of LM-generated SMILES are chemically invalid
- **Validity check**: `Chem.MolFromSmiles(smi) is not None`

### SELFIES (Self-Referencing Embedded Strings)

```python
import selfies as sf

smiles = "CC(=O)Oc1ccccc1C(=O)O"
enc = sf.encoder(smiles)     # "[C][C][=Branch1][C][=O][O][C:1][=C][C][=C][C][=C][C:1][C][=Branch1][C][=O][O]"
dec = sf.decoder(enc)        # back to valid SMILES
```

- **Alphabet**: semantic derivation rules — every sequence decodes to a valid molecule
- **100% validity** by construction; critical for LM generation
- Version 2.x alphabet: `sf.get_semantic_robust_alphabet()` (~800 tokens)
- Limitation: some complex stereochemistry not fully representable; SMILES more expressive

### Molecular Graphs

- **Nodes**: atoms (atom type, charge, H count); **Edges**: bonds (single/double/triple/aromatic)
- Used natively by GNN-based generators (GCPN, GraphAF, JT-VAE)
- **Junction Tree (JT)**: decompose molecule into ring/functional group vocabulary, then assemble tree
- Most expressive representation; no string validity issues; slower to decode

### 3D Point Clouds / Voxels

- Atom type + (x, y, z) coordinates
- Used by diffusion models (DiffSBDD, TargetDiff)
- Requires conformer generation; conditioned on protein pocket grid

## Benchmark Suites

### MOSES (Molecular Sets)

```bash
pip install molsets
```

| Metric | Definition | Typical target |
|--------|-----------|----------------|
| Validity | % valid SMILES | > 90% |
| Uniqueness | % unique (if valid) | > 99% |
| Novelty | % not in train set | > 99% |
| Filters | % passing structural filters | > 99% |
| FCD | Fréchet ChemNet Distance | < 0.5 (excellent) |
| IntDiv | Mean 1 - Tanimoto(gen_i, gen_j) | > 0.85 |
| SNN | Mean Tanimoto to nearest train molecule | < 0.6 |
| Scaffold | % unique Murcko scaffolds / N_valid | Higher = better |

```python
from moses import get_all_metrics

# train_smiles: list of training SMILES (reference distribution)
# gen_smiles: list of generated SMILES
metrics = get_all_metrics(gen=gen_smiles, train=train_smiles, test=test_smiles)
print(metrics)  # dict with all MOSES metrics
```

### GuacaMol

```bash
pip install guacamol
```

Two evaluation modes:

**Distributional learning benchmark** (is the output distribution drug-like?):
```python
from guacamol.distribution_learning_benchmark import (
    KLDivBenchmark, FrechetBenchmark
)
from guacamol.assess_distribution_learning import assess_distribution_learning
```

**Goal-directed benchmark** (17 tasks: rediscovery, similarity, MPO, scaffold decoration):
```python
from guacamol.assess_goal_directed_generation import assess_goal_directed_generation
from guacamol.goal_directed_generator import GoalDirectedGenerator

class MyGenerator(GoalDirectedGenerator):
    def generate_optimized_molecules(self, scoring_function, number_molecules,
                                      starting_population=None):
        # scoring_function.score(smiles) → float ∈ [0, 1]
        # return top-N SMILES sorted by score
        ...

results = assess_goal_directed_generation(MyGenerator(), json_output_file="results.json")
```

Key goal-directed tasks:
- `rediscovery_*`: recover exact molecule
- `similarity_*`: maximize Tanimoto to target
- `MPO_*`: multi-property optimization (celecoxib MPO, fexofenadine MPO, ...)
- `scaffold_deco_*`: fill SMILES scaffold with substituents
- `median_*`: find molecule mediating two targets

## Fréchet ChemNet Distance (FCD)

- **Architecture**: ChemNet (2-layer LSTM trained on 72M ChEMBL SMILES) — extracts 512-dim embedding
- Compute FCD between two distributions using Gaussian approximation (like FID for images)
- **Formula**: FCD = ||μ₁ - μ₂||² + Tr(Σ₁ + Σ₂ - 2√(Σ₁Σ₂))

```python
from fcd_torch import FCD  # pip install fcd_torch

fcd = FCD(device='cuda', n_jobs=8)
score = fcd(gen_smiles, ref_smiles)  # lower = more similar distributions
# Typically: < 0.5 = excellent; 1-3 = good; > 5 = poor
```

## Property Scoring Oracles

These are used as reward signals for RL or as evaluation filters:

```python
from rdkit import Chem
from rdkit.Chem import QED, Descriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

# QED: drug-likeness composite (0-1, higher = better)
def qed_score(smi):
    mol = Chem.MolFromSmiles(smi)
    return QED.qed(mol) if mol else 0.0

# SA score: synthetic accessibility (1=easy, 10=hard, target ≤ 4)
# Requires sascorer (from RDKit Contrib or pip install sascorer)
from sascorer import calculateScore  # pip install sascorer
def sa_score(smi):
    mol = Chem.MolFromSmiles(smi)
    return calculateScore(mol) if mol else 10.0

# LogP
def logp_score(smi):
    mol = Chem.MolFromSmiles(smi)
    return Descriptors.MolLogP(mol) if mol else None

# Penalized LogP (common RL benchmark objective)
def penalized_logp(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None: return -100.0
    log_p = Descriptors.MolLogP(mol)
    sa = calculateScore(mol)
    # Penalize for large rings
    from rdkit.Chem import rdMolDescriptors
    cycle_list = mol.GetRingInfo().AtomRings()
    large_ring_penalty = max(0, max((len(j) - 6) for j in cycle_list) if cycle_list else 0)
    return log_p - (sa - 1) / 9 - large_ring_penalty
```

## When to Use Which Generation Paradigm

| Situation | Recommended approach |
|-----------|---------------------|
| Fast screening, no structure | SELFIES + fine-tuned GPT/LSTM |
| Property optimization from scratch | REINVENT 4 (RL) |
| Smooth interpolation, BO in latent | JT-VAE / MolVAE |
| Protein pocket available | DiffSBDD / TargetDiff |
| Fragment growing / linker | DiffLinker / DeLinker |
| Benchmarking comparison | MOSES + GuacaMol (both) |
| Graph-native, custom architecture | TorchDrug (GCPN, GraphAF) |

## Scaffold Diversity Analysis

```python
from rdkit.Chem.Scaffolds import MurckoScaffold

def get_scaffold(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None: return None
    return MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)

def scaffold_diversity(smiles_list):
    scaffolds = [get_scaffold(s) for s in smiles_list]
    scaffolds = [s for s in scaffolds if s]
    return len(set(scaffolds)) / len(scaffolds)  # 1.0 = fully diverse
```

## Key Numbers to Remember

| Chemical space fact | Value |
|--------------------|----|
| Drug-like space (Lipinski) | ~10^60 |
| ZINC20 purchasable | ~1.4 × 10^9 |
| ChEMBL 34 compounds | ~2.4 × 10^6 |
| SMILES LM typical validity | 85-99% |
| SELFIES validity | 100% |
| Good FCD (GuacaMol) | < 0.5 |
| MOSES IntDiv target | > 0.85 |
