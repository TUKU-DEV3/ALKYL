# RL-Based Optimization — REINVENT 4

## REINVENT 4 Overview

REINVENT 4 (AstraZeneca, open-source) is the production-grade RL framework for de novo molecular design. Uses an LSTM prior (frozen) + LSTM agent (trainable) with REINFORCE updates guided by a multi-component scoring function.

```bash
pip install reinvent
# Requires: rdkit, torch, numpy, pandas, toml
```

**Architecture**:
```
Prior (frozen pre-trained LSTM) ──┐
                                  ├──→ REINFORCE update ──→ Agent (updated each step)
Agent (copy of prior, trainable) ─┘
                                       ↑
                         Scoring function (oracle calls)
```

**Key insight**: Prior acts as KL-regularizer: `Loss = -Score × (log P_agent - log P_prior)`. Without prior constraint, agent diverges to invalid/unrealistic molecules.

## TOML Configuration

REINVENT 4 is fully config-driven via TOML files.

### Minimal RL run

```toml
[parameters]
run_type = "reinforcement_learning"
device = "cuda"
tb_logdir = "tb_logs"
json_out_config = "_run_config.json"

[parameters.prior]
filename = "models/reinvent.prior"

[parameters.agent]
filename = "models/reinvent.prior"   # agent starts from prior

[parameters.reinforcement_learning]
batch_size = 128
n_steps = 1000
sigma = 128                          # RL temperature: higher = more exploitation
learning_rate = 0.0001
starting_kl_coefficient = 0.0       # adaptive KL penalty (0 = off)

[scoring]
type = "custom_product"             # geometric mean of components
parallel = false

[[scoring.component]]
type = "qed"
name = "QED"
weight = 1.0

[[scoring.component]]
type = "sa_score"
name = "SA Score"
weight = 1.0

[scoring.component.transform]
type = "reverse_sigmoid"
low = 1.0
high = 6.0
k = 0.5

[diversity_filter]
type = "IdenticalMurckoScaffold"    # penalize repeated scaffolds
minscore = 0.4                      # minimum score to enter scaffold memory
bucket_size = 25                    # max molecules per scaffold bucket

[inception]
memory_size = 20                    # top molecules kept in replay buffer
sample_size = 5                     # molecules sampled from buffer each step
min_score = 0.4
```

```bash
reinvent -l rl_run.log rl_config.toml
```

## Scoring Component Types

### Built-in components

| Type | Description | Range |
|------|-------------|-------|
| `qed` | Quantitative drug-likeness | 0-1 |
| `sa_score` | Synthetic accessibility (inverted) | 0-1 (1=easy) |
| `num_rotatable_bonds` | Rotatable bond count | raw int |
| `tpsa` | Topological polar surface area | raw float |
| `molecular_weight` | MW in Da | raw float |
| `alogp` | Calculated LogP | raw float |
| `num_hbd` / `num_hba` | H-bond donors/acceptors | raw int |
| `matching_substructure` | SMARTS substructure match | 0 or 1 |
| `custom_alerts` | SMARTS filter (PAINS, etc.) | 0 or 1 |

### Docking oracle (Vina/Gnina integration)

```toml
[[scoring.component]]
type = "dockstream"                 # DockStream plugin for Vina/Glide
name = "Vina_docking"
weight = 1.0

[scoring.component.params]
configuration_path = "dockstream_config.json"
docker_script_path = "path/to/dockstream/docker/AZdock/AZdock.py"
environment_path = "path/to/conda/env"
receptor_path = "receptor.pdbqt"
grid_center = [10.5, 20.3, -5.1]
grid_size = [20.0, 20.0, 20.0]
```

### Custom Python scoring component

```python
# custom_component.py — implement score() method

from rdkit import Chem
from rdkit.Chem import Descriptors

class CustomScoringComponent:
    """Custom scoring component for REINVENT 4."""

    def __init__(self, params: dict):
        self.target_mw = params.get("target_mw", 400.0)

    def score(self, smiles: list[str]) -> list[float]:
        scores = []
        for smi in smiles:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                scores.append(0.0)
            else:
                mw = Descriptors.ExactMolWt(mol)
                # Gaussian reward centered on target MW
                score = float(np.exp(-0.5 * ((mw - self.target_mw) / 50) ** 2))
                scores.append(score)
        return scores
```

```toml
[[scoring.component]]
type = "custom"
name = "MW_target"
weight = 1.0

[scoring.component.params]
module = "custom_component"
class = "CustomScoringComponent"
target_mw = 400.0
```

## Scoring Aggregation

REINVENT supports three aggregation modes:

| Type | Formula | When to use |
|------|---------|-------------|
| `custom_product` | ∏ score_i^(w_i) | All components must be satisfied |
| `custom_sum` | Σ w_i × score_i | Partial satisfaction acceptable |
| `arithmetic_mean` | mean(scores) | Equal importance |

**Recommendation**: Use `custom_product` — it forces all constraints to be satisfied simultaneously (score collapses to ~0 if any component fails).

## Score Transforms

Raw oracle outputs → normalized [0, 1] via transforms:

```python
# sigmoid: score_norm = 1 / (1 + exp(-k*(x - inflection)))
# reverse_sigmoid: 1 - sigmoid (for SA: low SA = good)
# step: threshold function (hard cutoff)
# double_sigmoid: window function (MW between 300-500)
# gaussian: Gaussian around target value
```

```toml
# MW window 300-500 Da
[scoring.component.transform]
type = "double_sigmoid"
low = 300.0
high = 500.0
coef_div = 500.0
coef_si = 20.0
coef_se = 20.0
```

## Diversity Filters

Critical to prevent mode collapse (RL converges to 1 scaffold).

| Filter | Behavior | Use when |
|--------|----------|----------|
| `IdenticalMurckoScaffold` | Penalize exact Murcko scaffold repeats | Default; most common |
| `IdenticalTopologicalScaffold` | Stricter: penalize topological scaffold | When Murcko too loose |
| `ScaffoldSimilarity` | Penalize similar (Tanimoto > threshold) scaffolds | Broad scaffold diversity |
| `NoFilter` | No diversity enforcement | Testing only |

```toml
[diversity_filter]
type = "IdenticalMurckoScaffold"
minscore = 0.4          # only high-scoring mols enter scaffold memory
bucket_size = 25        # after 25 mols with same scaffold, score → 0
```

## Staged Learning (Curriculum)

Multi-stage training: start broad → progressively tighten constraints.

```toml
[parameters.staged_learning]
stages = [
    # Stage 1: broad QED + SA filter (500 steps)
    {type = "reinforcement_learning", n_steps = 500,
     scoring = {components = [{type = "qed", weight = 1.0}]}},
    # Stage 2: add MW constraint (500 steps)
    {type = "reinforcement_learning", n_steps = 500,
     scoring = {components = [
         {type = "qed", weight = 1.0},
         {type = "molecular_weight", transform = {type = "double_sigmoid", low = 300, high = 500}}
     ]}},
]
```

## Transfer Learning Mode (prior fine-tuning)

Before RL, fine-tune prior on a seed library:

```toml
[parameters]
run_type = "transfer_learning"

[parameters.transfer_learning]
filename = "models/reinvent.prior"
input_smiles_path = "seed_library.smi"
num_epochs = 10
batch_size = 64
learning_rate = 0.0001
starting_epoch = 0
save_every_n_epochs = 1
sample_after_epoch = true
```

## Sampling from Trained Agent

```python
import torch
from reinvent.models import RNN  # or use reinvent CLI

# CLI sampling (easiest):
# reinvent --sample rl_checkpoint.chkpt --n 10000 --output generated.csv

# Python API:
from reinvent.chemistry.library_design.reaction_filters.reaction_filter_factory import ReactionFilterFactory
from reinvent.models.reinvent.models.model import Model

model = Model.load_from_file("rl_checkpoint.chkpt", mode="inference")
smiles, nlls = model.sample_smiles(n=1000)  # nlls = negative log-likelihoods (diversity proxy)
```

## REINFORCE Algorithm (Manual Implementation)

For custom RL loops not using REINVENT infrastructure:

```python
import torch
from rdkit import Chem

def reinforce_step(agent, prior, smiles_batch, scores, sigma=128, device='cpu'):
    """Single REINFORCE update step."""
    # Tokenize and get log-likelihoods from agent and prior
    agent_nlls = agent.likelihood(smiles_batch)    # negative log-likelihoods
    prior_nlls = prior.likelihood(smiles_batch)

    # Augmented NLL (penalize deviation from prior)
    augmented_nll = prior_nlls + sigma * scores    # lower score = higher augmented NLL

    # REINFORCE loss
    loss = torch.mean((augmented_nll - agent_nlls) ** 2)

    return loss

# Training loop
optimizer = torch.optim.Adam(agent.parameters(), lr=1e-4)
for step in range(n_steps):
    smiles, agent_nlls = agent.sample(batch_size=128)
    scores = scoring_function(smiles)   # oracle: array of floats [0, 1]
    loss = reinforce_step(agent, prior, smiles, torch.tensor(scores))
    optimizer.zero_grad(); loss.backward(); optimizer.step()
```

## Multi-Objective: Pareto Front

To maintain diverse multi-objective solutions (not just top-1):

```python
import numpy as np

def pareto_front(scores_matrix):
    """Find Pareto-optimal indices. scores_matrix: (N, n_objectives), higher=better."""
    N = len(scores_matrix)
    dominated = np.zeros(N, dtype=bool)
    for i in range(N):
        for j in range(N):
            if i == j: continue
            # j dominates i if j is >= on all objectives and > on at least one
            if np.all(scores_matrix[j] >= scores_matrix[i]) and \
               np.any(scores_matrix[j] > scores_matrix[i]):
                dominated[i] = True
                break
    return np.where(~dominated)[0]

# Usage: scores_matrix shape (batch_size, n_objectives)
# objectives: [qed_score, -sa_score, tanimoto_to_target, -docking_score]
pareto_idx = pareto_front(np.column_stack([qed, sa, tani, dock]))
```

## Key Pitfalls

- **No diversity filter → mode collapse**: RL always collapses to 1-3 scaffolds; always add `IdenticalMurckoScaffold`
- **Reward hacking**: very high `sigma` → agent ignores prior → generates chemically unusual structures
- **SA score is essential**: without it, RL generates high-scoring but unsynthesizable molecules (SA > 7)
- **Docking oracle must be fast**: Vina ~10 s/mol → batch 128 takes 20 min/step; use Gnina rescoring on top-20% only
- **`custom_product` with 0.0 component score = 0 total**: set minimum score floor (e.g., 0.05) to avoid dead-end penalties

## Related Tools

- **REINVENT 4 GitHub**: `MolecularAI/REINVENT` (Apache 2.0)
- **DockStream**: Vina/Glide/AutoDock integration for REINVENT oracle
- **GuacaMol**: standard benchmarking for goal-directed generation
- **molscore**: `pip install molscore` — modular scoring system compatible with REINVENT
