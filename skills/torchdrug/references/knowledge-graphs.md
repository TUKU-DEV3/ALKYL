# TorchDrug — Knowledge Graph Reasoning

## Overview

Knowledge graph (KG) completion predicts missing (head, relation, tail) triples. In biomedical KGs this enables drug repurposing, disease mechanism discovery, and gene-disease link prediction.

---

## Dataset Catalog

### General KGs

| Dataset | Entities | Relations | Triples | Notes |
|---------|----------|-----------|---------|-------|
| `FB15k` | 14,951 | 1,345 | 592k | Freebase subset |
| `FB15k-237` | 14,541 | 237 | 310k | FB15k without inverse leakage |
| `WN18` | 40,943 | 18 | 151k | WordNet |
| `WN18RR` | 40,943 | 11 | 93k | WN18 without inverse leakage |

### Biomedical KGs

| Dataset | Entities | Relations | Triples | Content |
|---------|----------|-----------|---------|---------|
| `Hetionet` | 45,158 | 11 | 2.25M | Drugs, diseases, genes, pathways, anatomy |

**Hetionet** covers compounds, diseases, genes, pathways, anatomy, GO terms with 11 relation types including `Compound-treats-Disease`.

---

## Loading Datasets

```python
from torchdrug import datasets

dataset = datasets.FB15k237("~/kg-datasets/")
# or
dataset = datasets.Hetionet("~/kg-datasets/")

print(dataset.num_entity)     # total entity count
print(dataset.num_relation)   # total relation types
print(dataset[0])             # (head_idx, relation_idx, tail_idx)
```

---

## Task: KnowledgeGraphCompletion

```python
from torchdrug import datasets, models, tasks
import torch

dataset = datasets.FB15k237("~/kg-datasets/")
train_set, valid_set, test_set = dataset.split()

model = models.RotatE(
    num_entity=dataset.num_entity,
    num_relation=dataset.num_relation,
    embedding_dim=256,       # must be even for RotatE
    max_score=9.0
)

task = tasks.KnowledgeGraphCompletion(
    model,
    num_negative=256,
    adversarial_temperature=1.0,
    strict_negative=True
)

optimizer = torch.optim.Adam(task.parameters(), lr=2e-4)
```

### Full training loop

```python
from torch.utils.data import DataLoader

train_loader = DataLoader(train_set, batch_size=1024, shuffle=True)
valid_loader = DataLoader(valid_set, batch_size=256)

for epoch in range(1000):
    task.train()
    for batch in train_loader:
        loss = task(batch)
        optimizer.zero_grad(); loss.backward(); optimizer.step()

    if epoch % 50 == 0:
        task.train(mode=False)
        with torch.no_grad():
            preds   = torch.cat([task.predict(b) for b in valid_loader])
            targets = torch.cat([task.target(b)  for b in valid_loader])
        metrics = task.evaluate(preds, targets)
        print(f"Epoch {epoch}: MRR={metrics['mrr']:.4f}, Hits@10={metrics['hits@10']:.4f}")
```

---

## Embedding Models

### TransE
```python
# h + r ≈ t — simple, fast, 1-to-1 relations
model = models.TransE(num_entity=..., num_relation=..., embedding_dim=200, p_norm=1)
```

### RotatE (recommended)
```python
# Relations as rotations in complex space — all relation patterns
model = models.RotatE(num_entity=..., num_relation=..., embedding_dim=256, max_score=9.0)
```

### DistMult
```python
# Bilinear; fast; symmetric relations only
model = models.DistMult(num_entity=..., num_relation=..., embedding_dim=256)
```

### ComplEx
```python
# Complex-valued; asymmetric + symmetric; good balance
model = models.ComplEx(num_entity=..., num_relation=..., embedding_dim=256)
```

### SimplE
```python
# Two embeddings per entity; fully expressive
model = models.SimplE(num_entity=..., num_relation=..., embedding_dim=256)
```

### RGCN (GNN-based)
```python
# Relational graph convolution — uses graph structure
model = models.RGCN(
    num_entity=dataset.num_entity,
    num_relation=dataset.num_relation,
    hidden_dims=[256, 256],
    num_bases=30
)
```

---

## Model Selection Guide

| Relation pattern | Best model |
|-----------------|------------|
| General (unknown) | RotatE |
| Symmetric only | DistMult |
| Asymmetric + symmetric | ComplEx |
| Large graph, memory limit | TransE |
| GNN-based, multi-relational | RGCN |

---

## Evaluation Metrics

| Metric | Description | Target |
|--------|-------------|--------|
| `mr` | Mean rank of correct entity | lower |
| `mrr` | Mean reciprocal rank | higher |
| `hits@1` | Correct in top-1 | higher |
| `hits@3` | Correct in top-3 | higher |
| `hits@10` | Correct in top-10 | higher |

```python
metrics = task.evaluate(preds, targets)
# {"mr": ..., "mrr": ..., "hits@1": ..., "hits@3": ..., "hits@10": ...}
```

---

## Drug Repurposing with Hetionet

```python
dataset = datasets.Hetionet("~/kg-datasets/")

# Find relation index for "Compound-treats-Disease"
rel2idx = {rel: i for i, rel in enumerate(dataset.relation2id)}
treats_idx = rel2idx["Compound-treats-Disease"]

# After training, rank all compounds against a target disease
task.train(mode=False)
disease_idx = dataset.entity2id["Alzheimer's disease"]
# Build (compound_idx, treats_idx, disease_idx) triples for all compounds
# and use model.score() to rank them
```

---

## Negative Sampling

```python
task = tasks.KnowledgeGraphCompletion(
    model,
    num_negative=256,            # more = harder training signal
    adversarial_temperature=1.0, # self-adversarial negative sampling
    strict_negative=True,        # exclude known positives
)
```

---

## Hyperparameter Recommendations

| Setting | FB15k-237 | Hetionet |
|---------|-----------|---------|
| Model | RotatE | RotatE |
| `embedding_dim` | 256 | 200 |
| `num_negative` | 256 | 128 |
| `adversarial_temperature` | 1.0 | 0.5 |
| `lr` | 2e-4 | 1e-4 |
| `batch_size` | 1024 | 512 |
| Epochs | 1000 | 500 |

---

## Best Practices

- Use **RotatE** by default — handles symmetric, antisymmetric, inverse, composition patterns
- Use **FB15k-237** and **WN18RR** — originals (FB15k / WN18) have inverse relation leakage
- **Self-adversarial sampling** (`adversarial_temperature=1.0`) consistently improves MRR
- **Hetionet**: filter repurposing predictions by known mechanism and pathway coverage
- Log **MRR** (primary) + **Hits@10**; mean rank (MR) is sensitive to outliers
- Embedding dim ≥ 200 for competitive performance; 500 for state-of-the-art
