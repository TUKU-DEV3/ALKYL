# TorchDrug — Molecular Property Prediction

## Dataset Catalog

### Classification (drug discovery)

| Dataset | Molecules | Task | Targets |
|---------|-----------|------|---------|
| `BBBP` | 2,039 | Binary | Blood-brain barrier penetration |
| `BACE` | 1,513 | Binary | β-secretase (BACE-1) inhibition |
| `HIV` | 41,127 | Binary | HIV replication inhibition |
| `ClinTox` | 1,478 | Binary | Clinical trial toxicity (2 tasks) |
| `Tox21` | 7,831 | Multi-label | 12 toxicology targets |
| `ToxCast` | 8,576 | Multi-label | 617 assays |
| `SIDER` | 1,427 | Multi-label | 27 drug side-effect organ classes |
| `MUV` | 93,087 | Multi-label | 17 targets (maximum unbiased validation) |

### Regression

| Dataset | Molecules | Target | Unit |
|---------|-----------|--------|------|
| `ESOL` | 1,128 | Aqueous solubility | log mol/L |
| `FreeSolv` | 642 | Hydration free energy | kcal/mol |
| `Lipophilicity` | 4,200 | LogD (pH 7.4) | — |
| `SAMPL` | 643 | Solvation free energy | kcal/mol |

### Quantum chemistry

| Dataset | Molecules | Properties |
|---------|-----------|------------|
| `QM7` | 7,165 | Atomization energy, polarizability |
| `QM8` | 21,786 | Electronic spectra, excited states |
| `QM9` | 133,885 | 12 geometric/energetic/thermodynamic |
| `PCQM4M` | 3.8M | HOMO-LUMO gap (OGB benchmark) |

### Large-scale

| Dataset | Molecules | Purpose |
|---------|-----------|---------|
| `ZINC250k` | 250k | Drug-like, generative modeling |
| `ZINC2M` | 2M | Drug-like, large-scale pre-training |

---

## Tasks

### PropertyPrediction

General purpose — classification or regression at graph level.

```python
from torchdrug import datasets, models, tasks

# Classification example: BBBP
dataset = datasets.BBBP("~/datasets/")
train_set, valid_set, test_set = dataset.split()

model = models.GIN(
    input_dim=dataset.node_feature_dim,
    hidden_dims=[256, 256, 256, 256],
    edge_input_dim=dataset.edge_feature_dim,
    batch_norm=True, readout="mean"
)

task = tasks.PropertyPrediction(
    model,
    task=dataset.tasks,          # ["p_np"]
    criterion="bce",
    metric=("auprc", "auroc"),
    num_mlp_layer=2
)
```

### MultipleBinaryClassification

For multi-label datasets with missing labels (Tox21, SIDER, ToxCast):

```python
dataset = datasets.Tox21("~/datasets/")

task = tasks.MultipleBinaryClassification(
    model,
    task=list(range(len(dataset.tasks))),
    criterion="bce",
    metric=["auroc", "auprc"]
)
```

### Regression example (ESOL)

```python
dataset = datasets.ESOL("~/datasets/")

task = tasks.PropertyPrediction(
    model,
    task=dataset.tasks,
    criterion="mse",
    metric=["mae", "rmse"]
)
```

---

## Data Splitting

```python
# Scaffold split (default, recommended for drug discovery)
train, valid, test = dataset.split()

# Random split
train, valid, test = dataset.split(ratios=[0.8, 0.1, 0.1], random=True)

# Stratified (classification, balance labels)
train, valid, test = dataset.split(ratios=[0.8, 0.1, 0.1], stratified=True)
```

**Use scaffold split** for realistic drug discovery evaluation — prevents molecules with the same scaffold appearing in both train and test (which inflates performance with random splits).

---

## Model Selection Guide

| Scenario | Recommended model | Notes |
|----------|-------------------|-------|
| Small dataset (< 1k) | `GIN` + pre-training | Or pre-trained InfoGraph |
| General classification | `GIN` | Best balance of accuracy/speed |
| Interpretable predictions | `GAT` | Attention weights per neighbor |
| Quantum properties (2D) | `MPNN` | Edge features essential |
| Quantum properties (3D) | `SchNet` | Requires 3D coordinates |
| Large-scale (> 100k) | `GIN` deeper | Train from scratch |

---

## Full Training Pipeline

```python
import torch
from torch.utils.data import DataLoader
from torchdrug import datasets, models, tasks, transforms

# Optional: add virtual node for better graph readout
transform = transforms.VirtualNode()
dataset = datasets.BBBP("~/datasets/", transform=transform)
train_set, valid_set, test_set = dataset.split()

model = models.GIN(
    input_dim=dataset.node_feature_dim,
    hidden_dims=[256, 256, 256],
    edge_input_dim=dataset.edge_feature_dim,
    batch_norm=True, readout="mean"
)

task = tasks.PropertyPrediction(
    model, task=dataset.tasks,
    criterion="bce", metric=["auroc", "auprc"]
)

optimizer = torch.optim.Adam(task.parameters(), lr=1e-3, weight_decay=1e-5)
scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, patience=10)

train_loader = DataLoader(train_set, batch_size=32, shuffle=True,  num_workers=4)
valid_loader = DataLoader(valid_set, batch_size=64, shuffle=False, num_workers=4)

best_auroc, patience = 0, 0

for epoch in range(200):
    task.train()
    for batch in train_loader:
        loss = task(batch)
        optimizer.zero_grad(); loss.backward(); optimizer.step()

    task.train(mode=False)
    with torch.no_grad():
        preds   = torch.cat([task.predict(b) for b in valid_loader])
        targets = torch.cat([task.target(b)  for b in valid_loader])
    metrics = task.evaluate(preds, targets)
    scheduler.step(metrics["auroc"])

    if metrics["auroc"] > best_auroc:
        best_auroc = metrics["auroc"]
        torch.save(task.state_dict(), "best_model.pth")
        patience = 0
    else:
        patience += 1
        if patience > 20:
            break

print(f"Best validation AUROC: {best_auroc:.4f}")
```

---

## Self-Supervised Pre-training

Pre-train on unlabeled molecules before fine-tuning on small datasets:

```python
from torchdrug import tasks

# Attribute masking (mask atom features, predict them)
pretrain_task = tasks.AttributeMasking(model, mask_rate=0.15)

# Edge prediction (contrastive)
pretrain_task = tasks.EdgePrediction(model)

# Context prediction (subgraph/graph contrast)
pretrain_task = tasks.ContextPrediction(model, context_size=3)

# Or use InfoGraph (mutual information maximization)
from torchdrug import models as td_models
pretrain_task = td_models.InfoGraph(model)
```

Then fine-tune on downstream task:

```python
# Load pre-trained weights
model.load_state_dict(torch.load("pretrained.pth"), strict=False)

# Fine-tune with smaller learning rate
optimizer = torch.optim.Adam([
    {"params": model.parameters(), "lr": 1e-4},   # pre-trained
    {"params": task.mlp.parameters(), "lr": 1e-3}  # new head
])
```

---

## Custom Feature Engineering

```python
from torchdrug import data, transforms

# Add virtual node + edges (recommended for global features)
transform = transforms.Compose([
    transforms.VirtualNode(),
    transforms.VirtualEdge(),
])

# Build dataset with custom node features via subclassing
class CustomMolDataset(data.MoleculeDataset):
    def get_item(self, index):
        mol = super().get_item(index)
        # Add your custom features here
        return mol
```

---

## Best Practices

- **Scaffold split** over random split for realistic evaluation
- **`VirtualNode` transform** improves global information flow
- **Multi-task** (Tox21, SIDER) → use `MultipleBinaryClassification`, not `PropertyPrediction`
- **AUROC** for imbalanced binary; **RMSE + MAE** for regression
- **Pre-train** on ZINC or ChEMBL before fine-tuning on < 1k molecule datasets
- **Early stopping** with validation AUROC patience = 20 epochs
- **Ensemble** 5 models with different seeds for production
