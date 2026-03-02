# TorchDrug — Molecular Generation

## Overview

TorchDrug supports unconditional generation (explore chemical space) and conditional generation (optimize for target properties). Two main paradigms: autoregressive (GCPN with RL) and normalizing flows (GraphAutoregressiveFlow).

---

## Tasks

### AutoregressiveGeneration

Generates molecules atom-by-atom via sequential actions (add atom, add bond, terminate).

```python
from torchdrug import tasks, models

model = models.GIN(input_dim=..., hidden_dims=[256, 256])

task = tasks.AutoregressiveGeneration(
    model,
    atom_types=[6, 7, 8, 9, 16],   # C, N, O, F, S
    max_node=38,                     # max atoms per molecule
    criterion="nll"
)
```

**Generation strategies:**
- `beam_size=10` — beam search (quality)
- `num_sample=100` — stochastic sampling (diversity)

### GCPNGeneration (policy-based RL)

Reinforcement learning to directly optimize non-differentiable objectives (QED, SA score, binding affinity).

```python
from torchdrug import tasks

task = tasks.GCPNGeneration(
    model,
    atom_types=[6, 7, 8, 9, 16],
    max_edge_unroll=12,
    max_node=38,
    criterion="ppo",               # proximal policy optimization
    reward_temperature=1.0
)
```

---

## Generative Models

### GraphAutoregressiveFlow

Normalizing flow — exact likelihood, stable training, invertible transformations.

```python
from torchdrug import models

model = models.GraphAutoregressiveFlow(
    input_dim=...,
    hidden_dims=[256, 256, 256],
    num_flow=32,
    use_edge_feat=True
)
```

**Advantages over VAE/GAN:**
- Exact log-likelihood (no ELBO approximation)
- No mode collapse or training instability
- Invertible: can map between molecule and latent space

---

## Datasets for Generation

| Dataset | Size | Purpose |
|---------|------|---------|
| `ZINC250k` | 250k | Drug-like molecules, standard benchmark |
| `ZINC2M` | 2M | Large-scale drug-like pre-training |
| `QM9` | 133k | Small organics (≤ 9 heavy atoms) |
| `ChEMBL` | millions | Target-annotated bioactive compounds |

```python
from torchdrug import datasets

dataset = datasets.ZINC250k("~/datasets/")
```

---

## Unconditional Generation Workflow

```python
import torch
from torch.utils.data import DataLoader
from torchdrug import datasets, models, tasks

dataset = datasets.ZINC250k("~/datasets/")
train_set, _, _ = dataset.split()

model = models.GraphAutoregressiveFlow(
    input_dim=dataset.node_feature_dim,
    hidden_dims=[256, 256, 256],
    num_flow=32
)

task = tasks.AutoregressiveGeneration(model, atom_types=[6,7,8,9,16], max_node=38)
optimizer = torch.optim.Adam(task.parameters(), lr=1e-3)

for epoch in range(20):
    task.train()
    for batch in DataLoader(train_set, batch_size=32, shuffle=True):
        loss = task(batch)
        optimizer.zero_grad(); loss.backward(); optimizer.step()

# Sample new molecules
task.train(mode=False)
with torch.no_grad():
    mols = task.generate(num_sample=100)
# mols: list of torchdrug Molecule objects

# Validate with RDKit
from rdkit import Chem
valid_mols = [m for m in mols if Chem.MolFromSmiles(m.to_smiles()) is not None]
print(f"Validity: {len(valid_mols)/len(mols):.2%}")
```

---

## Conditional Generation (Property Optimization)

### Single-objective with GCPN

```python
from torchdrug import tasks
from rdkit.Chem import QED

def reward_fn(mol):
    """Example: maximize QED (drug-likeness)."""
    try:
        rdkit_mol = mol.to_molecule()
        return float(QED.qed(rdkit_mol))
    except Exception:
        return 0.0

task = tasks.GCPNGeneration(
    model,
    atom_types=[6, 7, 8, 9, 16],
    max_node=38,
    criterion="ppo",
    reward_function=reward_fn
)
```

### Multi-objective reward

```python
from rdkit.Chem import QED, Descriptors
from torchdrug.utils import sascore   # synthetic accessibility

def multi_reward(mol):
    try:
        rdkit_mol = mol.to_molecule()
        qed  = QED.qed(rdkit_mol)
        sa   = 1 - sascore.calculateScore(rdkit_mol) / 10.0   # normalize 0-1
        logp = Descriptors.MolLogP(rdkit_mol)
        logp_score = 1.0 if 0 <= logp <= 5 else 0.0           # Lipinski

        reward = 0.5 * qed + 0.3 * sa + 0.2 * logp_score
        return reward
    except Exception:
        return 0.0

task = tasks.GCPNGeneration(model, atom_types=[6,7,8,9,16], max_node=38,
                             criterion="ppo", reward_function=multi_reward)
```

---

## Validation and Filtering

After generation, always validate:

```python
from rdkit import Chem
from rdkit.Chem import QED, Descriptors

def filter_candidates(molecules):
    results = []
    for mol in molecules:
        try:
            rdkit_mol = mol.to_molecule()
            if rdkit_mol is None:
                continue

            mw   = Descriptors.ExactMolWt(rdkit_mol)
            logp = Descriptors.MolLogP(rdkit_mol)
            hbd  = Descriptors.NumHDonors(rdkit_mol)
            hba  = Descriptors.NumHAcceptors(rdkit_mol)

            # Lipinski Ro5
            if mw <= 500 and logp <= 5 and hbd <= 5 and hba <= 10:
                results.append({
                    "smiles": Chem.MolToSmiles(rdkit_mol),
                    "qed":    QED.qed(rdkit_mol),
                    "mw":     mw,
                    "logp":   logp
                })
        except Exception:
            continue
    return sorted(results, key=lambda x: x["qed"], reverse=True)
```

---

## Evaluation Metrics

| Metric | Description | Tool |
|--------|-------------|------|
| **Validity** | % chemically valid molecules | RDKit `MolFromSmiles` |
| **Uniqueness** | % unique SMILES among valid | set() deduplication |
| **Novelty** | % not in training set | fingerprint comparison |
| **Diversity** | Internal Tanimoto dissimilarity | Morgan FP + BulkSimilarity |
| **FCD** | Fréchet ChemNet Distance | `fcd` library |
| **KL divergence** | Property distribution match | scipy |

---

## Scaffold-Based Generation

Generate analogs around a fixed core:

```python
from torchdrug import data

# Define scaffold as SMILES with attachment points (*)
scaffold_smiles = "c1ccc(NC(=O)*)cc1"   # aniline scaffold
scaffold = data.Molecule.from_smiles(scaffold_smiles)

# Condition generation on scaffold
task = tasks.AutoregressiveGeneration(
    model,
    scaffold=scaffold,
    atom_types=[6, 7, 8, 9, 16],
    max_node=38
)
```

---

## Best Practices

- **Start unconditional** → validate quality → add property optimization
- **GCPN reward**: always include validity + drug-likeness as baseline components
- **Validity**: expect 85–99% for flow models; 60–90% for RL-based GCPN
- **ZINC250k** for drug-like; **QM9** for small organics; **ChEMBL** for target-specific
- **SMILES augmentation** during training improves generalization
- **Post-filter** every generated molecule with RDKit before ranking
- **Synthesizability**: use SA score < 4.0 as filter; run retrosynthesis on top candidates
- **Diversity**: use MaxMin picker (RDKit) to select representative subset for wet lab
