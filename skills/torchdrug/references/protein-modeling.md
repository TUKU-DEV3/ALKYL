# TorchDrug — Protein Modeling

## Dataset Catalog

### Function prediction

| Dataset | Proteins | Task | Classes |
|---------|----------|------|---------|
| `EnzymeCommission` | 17,562 | Multi-label EC number | 538 |
| `GeneOntology` | 46,796 | Multi-label GO terms | BP / MF / CC |
| `BetaLactamase` | 5,864 | Regression (activity) | — |
| `Fluorescence` | 54,025 | Regression (GFP intensity) | — |
| `Stability` | 53,614 | Regression (thermostability) | — |
| `Solubility` | 62,478 | Binary | membrane / soluble |
| `BinaryLocalization` | 22,168 | Binary | membrane / soluble |
| `SubcellularLocalization` | 8,943 | 10-class | localization |

### Structure

| Dataset | Proteins | Task |
|---------|----------|------|
| `Fold` | 16,712 | Fold classification (1,195 classes) |
| `SecondaryStructure` | 8,678 | 3-state or 8-state SS per residue |

### Interaction

| Dataset | Entries | Task |
|---------|---------|------|
| `HumanPPI` | 6,584 pairs | PPI binary classification |
| `YeastPPI` | 6,451 pairs | PPI binary classification |
| `PPIAffinity` | 2,156 pairs | Binding affinity regression |
| `PDBBind` | 20,000+ | Protein-ligand affinity (3D) |
| `BindingDB` | ~1.5M | Protein-ligand affinity |

---

## Protein Representation

### Sequence-based models

```python
from torchdrug import models

# ESM (recommended for sequence-only tasks)
model = models.ESM(path="esm1b_t33_650M_UR50S.pt")   # ESM-1b
# or ESM-2 variants: esm2_t6_8M, esm2_t33_650M, esm2_t36_3B ...

# ProteinBERT (lighter pre-trained)
model = models.ProteinBERT(path="proteinbert.pt")

# ProteinLSTM (fast baseline)
model = models.ProteinLSTM(input_dim=21, hidden_dim=512, num_layers=2)

# ProteinCNN (local pattern detection)
model = models.ProteinCNN(input_dim=21, hidden_dims=[512, 512], kernel_size=3)

# ProteinResNet (deep CNN with residual blocks)
model = models.ProteinResNet(input_dim=21, hidden_dims=[512]*6, kernel_size=3)
```

### Structure-based models

```python
from torchdrug import models, transforms

# GearNet — best for structure tasks
model = models.GearNet(
    input_dim=21,
    hidden_dims=[512, 512, 512],
    num_relation=7,           # edge type count (seq + radius + knn)
    edge_input_dim=59,        # geometric features
    batch_norm=True,
    readout="sum"
)

# SchNet — 3D continuous-filter convolutions
model = models.SchNet(
    input_dim=21,
    hidden_dims=[256, 256, 256],
    num_gaussian=25,
    cutoff=10.0
)
```

---

## Graph Construction for Proteins

```python
from torchdrug import data, transforms

# Load structure
protein = data.Protein.from_pdb("1a3x.pdb")

# Residue-level graph with multiple edge types
graph = protein.residue_graph(
    node_position="ca",
    edge_types=["sequential", "radius", "knn"],
    radius_cutoff=10.0,    # Å
    knn=10                 # K nearest neighbors
)
```

Use `TruncateProtein` transform to cap sequence length for large proteins:

```python
from torchdrug import transforms, datasets

transform = transforms.TruncateProtein(max_length=500)
dataset = datasets.EnzymeCommission("~/datasets/", transform=transform)
```

---

## Task Types

### PropertyPrediction (protein-level)

```python
from torchdrug import tasks, datasets, models

dataset = datasets.EnzymeCommission("~/datasets/")
train_set, valid_set, test_set = dataset.split()

model = models.ESM(path="esm1b_t33_650M_UR50S.pt")

task = tasks.PropertyPrediction(
    model,
    task=dataset.tasks,
    criterion="bce",
    metric=["auroc", "auprc"]
)
```

### NodePropertyPrediction (residue-level)

For secondary structure, contact maps, binding sites:

```python
dataset = datasets.SecondaryStructure("~/datasets/")

task = tasks.NodePropertyPrediction(
    model,
    task=dataset.tasks,
    criterion="ce",          # 3-class or 8-class SS
    metric=["acc"]
)
```

### InteractionPrediction (protein pairs or protein-ligand)

```python
dataset = datasets.HumanPPI("~/datasets/")

task = tasks.InteractionPrediction(
    model,
    task="interaction",
    criterion="bce",
    metric=["auroc"]
)
```

---

## Pre-training Strategies

### Load and fine-tune ESM

```python
import torch
from torchdrug import models, tasks, datasets

# Load pre-trained ESM
model = models.ESM(path="esm1b_t33_650M_UR50S.pt")

# Fine-tune on stability prediction
dataset = datasets.Stability("~/datasets/")
task = tasks.PropertyPrediction(
    model, task=dataset.tasks,
    criterion="mse", metric=["mae", "rmse"]
)

optimizer = torch.optim.Adam([
    {"params": model.parameters(),   "lr": 1e-5},   # frozen-ish
    {"params": task.mlp.parameters(), "lr": 1e-3},   # fresh head
])
```

### GearNet with MultiviewContrast pre-training

```python
from torchdrug import tasks as td_tasks

# Pre-train GearNet on large structure databases
pretrain_task = td_tasks.MultiviewContrast(
    model,
    views=["sequence", "structure"],
    temperature=0.07
)

# Then fine-tune on downstream task
```

---

## Full Workflow: Enzyme Function Prediction

```python
import torch
from torch.utils.data import DataLoader
from torchdrug import datasets, models, tasks, transforms

transform = transforms.TruncateProtein(max_length=500)
dataset = datasets.EnzymeCommission("~/datasets/", transform=transform)
train_set, valid_set, test_set = dataset.split()

# Use GearNet (structure available in EnzymeCommission)
model = models.GearNet(
    input_dim=21,
    hidden_dims=[512, 512, 512],
    num_relation=7,
    edge_input_dim=59,
    batch_norm=True,
    readout="sum"
)

task = tasks.PropertyPrediction(
    model, task=dataset.tasks,
    criterion="bce", metric=["auroc", "auprc"]
)

optimizer = torch.optim.Adam(task.parameters(), lr=1e-3)

for epoch in range(50):
    task.train()
    for batch in DataLoader(train_set, batch_size=4, shuffle=True):
        loss = task(batch)
        optimizer.zero_grad(); loss.backward(); optimizer.step()

    task.train(mode=False)
    with torch.no_grad():
        preds   = torch.cat([task.predict(b) for b in DataLoader(valid_set, batch_size=4)])
        targets = torch.cat([task.target(b)  for b in DataLoader(valid_set, batch_size=4)])
    print(f"Epoch {epoch}: {task.evaluate(preds, targets)}")
```

---

## Integration with External Tools

### AlphaFold predicted structures

```python
from torchdrug import data

# Load AF2 predicted PDB (same API as experimental)
protein = data.Protein.from_pdb("AF-P12345-F1-model_v4.pdb")

# Use pLDDT b-factor as confidence filter
# (access via protein.residue_number or custom PDB parsing)
```

### ESMFold (sequence → structure → TorchDrug)

```python
# Generate structure with ESMFold
import esm
model_esm, alphabet = esm.pretrained.esmfold_v1()
structure_pdb = model_esm.infer_pdb("MKTIIALSYIFCLVFA")

# Load in TorchDrug
from torchdrug import data
protein = data.Protein.from_pdb(structure_pdb)
```

---

## Model Selection Guide

| Task | Sequence? | Structure? | Recommended |
|------|-----------|------------|-------------|
| Function (EC, GO) | ✓ only | — | ESM-1b / ESM-2 |
| Function (EC, GO) | ✓ | ✓ | GearNet |
| Stability / fluorescence | ✓ | — | ESM fine-tune |
| Secondary structure | ✓ | — | ProteinResNet or ESM |
| PPI | ✓ | — | ProteinBERT + interaction head |
| Protein-ligand binding | ✓ | ✓ | SchNet on complex |
| Small dataset (< 500 proteins) | ✓ | — | Freeze ESM, train head only |

---

## Best Practices

- **Sequence tasks**: start with pre-trained ESM; fine-tune with lr 1e-5 on backbone, 1e-3 on head
- **Structure tasks**: GearNet with `["sequential", "radius", "knn"]` edges
- **Small datasets**: freeze pre-trained backbone; only train prediction head
- **Batch size**: 4–8 for proteins (large graphs); accumulate gradients if needed
- **Truncate**: always apply `TruncateProtein(500)` to avoid OOM on full sequences
- **Regression** (Fluorescence, Stability): MSE loss + MAE + Spearman correlation
- **Multi-label** (EC, GO): BCE loss + AUROC per ontology branch
