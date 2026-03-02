# TorchDrug — Model Reference Catalog

## GNN Models (Molecular)

### GCN — Graph Convolutional Network

**Paper:** Kipf & Welling, ICLR 2017

```python
from torchdrug import models

model = models.GCN(
    input_dim=dataset.node_feature_dim,
    hidden_dims=[256, 256, 256],
    edge_input_dim=dataset.edge_feature_dim,  # optional
    batch_norm=True,
    activation="relu",
    dropout=0.2,
    readout="mean"
)
```

| Aspect | Detail |
|--------|--------|
| Aggregation | Normalized adjacency (symmetric) |
| Best for | Baselines, simple graphs |
| Weakness | Homophily assumption; no edge types |
| Speed | Fast |

---

### GAT — Graph Attention Network

**Paper:** Veličković et al., ICLR 2018

```python
model = models.GAT(
    input_dim=dataset.node_feature_dim,
    hidden_dims=[256, 256, 256],
    num_heads=8,
    negative_slope=0.2,   # LeakyReLU slope
    dropout=0.2,
    readout="mean"
)
```

| Aspect | Detail |
|--------|--------|
| Aggregation | Learned attention weights per neighbor |
| Best for | Interpretability; heterogeneous neighborhoods |
| Output | Multi-head attention → concatenate or average |
| Speed | Medium |

---

### GIN — Graph Isomorphism Network ★ (default choice)

**Paper:** Xu et al., ICLR 2019

```python
model = models.GIN(
    input_dim=dataset.node_feature_dim,
    hidden_dims=[256, 256, 256, 256],
    edge_input_dim=dataset.edge_feature_dim,
    batch_norm=True,
    readout="sum",   # sum is theoretically motivated for GIN
    eps=0            # learnable epsilon
)
```

| Aspect | Detail |
|--------|--------|
| Aggregation | Injective (maximally expressive) |
| Best for | Molecular property prediction — state-of-the-art |
| Key insight | Distinguishes graph structures GCN/GAT cannot |
| Speed | Medium |

---

### RGCN — Relational GCN

**Paper:** Schlichtkrull et al., ESWC 2018

```python
model = models.RGCN(
    input_dim=dataset.node_feature_dim,
    num_relation=dataset.num_bond_type,
    hidden_dims=[256, 256, 256],
    num_bases=30,     # basis decomposition (reduces parameters)
    dropout=0.2
)
```

| Aspect | Detail |
|--------|--------|
| Aggregation | Relation-specific weight matrices |
| Best for | Knowledge graphs; retrosynthesis (multi-bond types) |
| Key param | `num_bases` — lower = fewer params, more sharing |
| Speed | Medium-slow (scales with num_relation) |

---

### MPNN — Message Passing Neural Network

**Paper:** Gilmer et al., ICML 2017

```python
model = models.MPNN(
    input_dim=dataset.node_feature_dim,
    hidden_dim=256,
    edge_input_dim=dataset.edge_feature_dim,
    num_layer=3,
    num_mlp_layer=2,
    readout="set2set"   # Set2Set readout (more expressive)
)
```

| Aspect | Detail |
|--------|--------|
| Aggregation | Custom message + GRU update |
| Best for | Quantum chemistry (QM9); edge-feature-rich tasks |
| Key feature | Set2Set readout for graph-level tasks |
| Speed | Medium |

---

### SchNet — Continuous-Filter CNN

**Paper:** Schütt et al., NeurIPS 2017

```python
model = models.SchNet(
    input_dim=dataset.node_feature_dim,
    hidden_dims=[256, 256, 256],
    num_gaussian=25,   # RBF basis for inter-atomic distances
    cutoff=10.0        # interaction cutoff in Å
)
```

| Aspect | Detail |
|--------|--------|
| Requires | 3D coordinates |
| Invariance | Rotation + translation invariant |
| Best for | QM9, molecular dynamics, protein-ligand (3D) |
| Speed | Slow (3D graph construction) |

---

### NFP — Neural Fingerprint

**Paper:** Duvenaud et al., NeurIPS 2015

```python
model = models.NFP(
    input_dim=dataset.node_feature_dim,
    output_dim=256,
    hidden_dims=[256, 256],
    num_layer=3
)
```

| Aspect | Detail |
|--------|--------|
| Output | Differentiable molecular fingerprint |
| Best for | When interpretability or ECFP drop-in needed |
| Speed | Fast |

---

### ChebNet — Chebyshev Spectral GNN

**Paper:** Defferrard et al., NeurIPS 2016

```python
model = models.ChebNet(
    input_dim=dataset.node_feature_dim,
    hidden_dims=[256, 256],
    num_cheb=3   # Chebyshev polynomial order
)
```

| Aspect | Detail |
|--------|--------|
| Aggregation | Spectral (graph Laplacian) |
| Best for | Tasks requiring global structure information |
| Speed | Fast |

---

## Protein-Specific Models

### GearNet — Geometry-Aware Relational Graph Network ★

**Paper:** Zhang et al., ICLR 2023 (doi:10.48550/arXiv.2203.06125)

```python
model = models.GearNet(
    input_dim=21,          # 20 AA + 1 virtual
    hidden_dims=[512, 512, 512, 512, 512, 512],
    num_relation=7,        # seq-forward, seq-backward, radius, KNN, ...
    edge_input_dim=59,     # geometric features (distances, angles, dihedrals)
    batch_norm=True,
    readout="sum",
    dropout=0.2
)
```

| Aspect | Detail |
|--------|--------|
| Edge types | Sequential + spatial (radius/KNN) + contact |
| Geometric features | Distances, angles, dihedrals in edge_input |
| Best for | Structure-based protein tasks — SOTA |
| Pre-training | MultiviewContrast on PDB structures |

---

### ESM — Evolutionary Scale Modeling ★ (sequence tasks)

**Paper:** Rives et al., PNAS 2021 (doi:10.1073/pnas.2016239118)

```python
model = models.ESM(path="esm1b_t33_650M_UR50S.pt")

# ESM-2 variants by size:
# esm2_t6_8M_UR50D.pt       (8M params, fast)
# esm2_t12_35M_UR50D.pt     (35M)
# esm2_t30_150M_UR50D.pt    (150M)
# esm2_t33_650M_UR50D.pt    (650M, recommended)
# esm2_t36_3B_UR50D.pt      (3B, large)
```

| Aspect | Detail |
|--------|--------|
| Architecture | Transformer (BERT-style) |
| Pre-training | 250M+ UniRef50 sequences |
| Best for | Any sequence-only protein task |
| Fine-tuning lr | 1e-5 on backbone; 1e-3 on head |

---

### ProteinBERT

```python
model = models.ProteinBERT(path="proteinbert.pt")
```

Lighter than ESM; masked language model pre-training on UniProt.

---

### ProteinLSTM / ProteinCNN / ProteinResNet

```python
# Fast baselines for sequence tasks
model = models.ProteinLSTM(input_dim=21, hidden_dim=512, num_layers=3)
model = models.ProteinCNN(input_dim=21, hidden_dims=[512]*5, kernel_size=3)
model = models.ProteinResNet(input_dim=21, hidden_dims=[512]*6, kernel_size=3)
```

| Model | Speed | Long-range | Use when |
|-------|-------|-----------|----------|
| ProteinLSTM | Medium | Good | Sequential tasks |
| ProteinCNN | Fast | Poor | Local motifs |
| ProteinResNet | Fast | Medium | Deep CNN with residuals |
| ESM | Slow | Excellent | Default sequence choice |
| GearNet | Slow | Excellent | Structure available |

---

## Knowledge Graph Embedding Models

### Quick comparison

| Model | Relation patterns | `embedding_dim` | Notes |
|-------|------------------|-----------------|-------|
| `TransE` | 1-to-1 | any | Simplest, fastest |
| `DistMult` | Symmetric | any | Cannot model antisymmetry |
| `ComplEx` | Sym + antisym | any | Good general baseline |
| `RotatE` | All patterns | even only | Best overall ★ |
| `SimplE` | All patterns | any | Slightly more params than ComplEx |

```python
# All share the same constructor pattern:
model = models.RotatE(
    num_entity=dataset.num_entity,
    num_relation=dataset.num_relation,
    embedding_dim=256    # must be even for RotatE
)
```

---

## Generative Models

### GraphAutoregressiveFlow

```python
model = models.GraphAutoregressiveFlow(
    input_dim=dataset.node_feature_dim,
    hidden_dims=[256, 256, 256],
    num_flow=32,
    use_edge_feat=True
)
```

Normalizing flow — exact likelihood, invertible, stable training.

---

## Pre-training / Self-supervised Models

### InfoGraph (molecular)

```python
from torchdrug import models as td_models
pretrain = td_models.InfoGraph(model)
```

Maximizes mutual information between graph-level and node-level representations. Use to pre-train molecule encoders on unlabeled ZINC/ChEMBL.

### MultiviewContrast (protein)

```python
from torchdrug import tasks

pretrain_task = tasks.MultiviewContrast(
    model,
    views=["sequence", "structure"],
    temperature=0.07
)
```

Contrasts sequence and structure views of the same protein. Best for pre-training GearNet on PDB.

---

## Model Selection by Task

| Task | Dataset type | Recommended |
|------|-------------|-------------|
| Molecular property (class/reg) | 2D mol graph | **GIN** |
| Molecular property (3D available) | 3D structure | **SchNet** |
| Protein function (no structure) | Sequence | **ESM** |
| Protein function (structure) | 3D structure | **GearNet** |
| Knowledge graph completion | Entity-relation | **RotatE** |
| Retrosynthesis center | Reaction graph | **RGCN** |
| Retrosynthesis synthon | Product graph | **GIN** |
| Molecular generation (flow) | ZINC | **GraphAutoregressiveFlow** |
| Molecular generation (RL) | ZINC | **GCPNGeneration + GIN** |

---

## Model Selection by Dataset Size

| Dataset size | Strategy |
|-------------|---------|
| < 500 | Pre-trained model, frozen backbone |
| 500–5k | Pre-trained + fine-tune, strong regularization |
| 5k–100k | GIN/GAT from scratch, scaffold split |
| > 100k | Any model; deeper architectures viable |

---

## Custom Model Template

```python
from torchdrug import core, layers
import torch.nn as nn

@core.register("models.MyModel")
class MyModel(nn.Module, core.Configurable):
    def __init__(self, input_dim, hidden_dim, output_dim):
        super().__init__()
        self.input_dim  = input_dim
        self.output_dim = output_dim
        self.conv = layers.GraphConv(input_dim, hidden_dim)
        self.pool = nn.Linear(hidden_dim, output_dim)

    def forward(self, graph, input, all_loss=None, metric=None):
        node_feat = self.conv(graph, input)
        graph_feat = graph.node2graph(node_feat, reduction="mean")
        return {"node_feature": node_feat, "graph_feature": graph_feat}
```

Key rules:
1. Set `self.input_dim` and `self.output_dim`
2. Forward returns a dict with `"node_feature"` and/or `"graph_feature"`
3. Inherit `core.Configurable` for serialization
4. Decorate with `@core.register(...)` for string-based lookup
