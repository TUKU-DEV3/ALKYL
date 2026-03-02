# TorchDrug — Core Data Structures & Architecture

## Data Structures at a Glance

| Class | Inherits | Represents | Key extra attributes |
|-------|----------|------------|----------------------|
| `Graph` | — | Any graph | `num_node`, `num_edge`, `node_feature`, `edge_feature`, `edge_list` |
| `Molecule` | `Graph` | Small molecule | `atom_type`, `bond_type`, `formal_charge`, `explicit_hs` |
| `Protein` | `Graph` | Protein | `residue_type`, `atom_name`, `residue_number`, `chain_id` |
| `PackedGraph` | — | Batched graphs | `num_nodes`, `num_edges`, `graph_ind` |

---

## Molecule

```python
from torchdrug import data

# From SMILES
mol = data.Molecule.from_smiles("CCO")
print(mol.atom_type)   # tensor([6, 6, 8])   C, C, O
print(mol.bond_type)   # tensor([1, 1])        single bonds

# From / to RDKit
from rdkit import Chem
rdkit_mol = Chem.MolFromSmiles("c1ccccc1")
mol = data.Molecule.from_molecule(rdkit_mol)
back = mol.to_molecule()   # → RDKit Mol
smiles = mol.to_smiles()
```

### Node features (auto-extracted)
- atom type, formal charge, explicit/implicit H, hybridization, aromaticity, chirality

### Edge features (auto-extracted)
- bond type (single/double/triple/aromatic), stereochemistry, conjugation, ring membership

---

## Protein

```python
from torchdrug import data

# Load from PDB
protein = data.Protein.from_pdb("1a3x.pdb")

# From sequence only (no structure)
protein = data.Protein.from_sequence("MKTIIALSYIFCLVFA")

# Build residue-level graph (nodes = residues, not atoms)
graph = protein.residue_graph(
    node_position="ca",                        # Cα atoms
    edge_types=["sequential", "radius"],       # backbone + spatial
    radius_cutoff=10.0                         # Å
)

protein.to_pdb("output.pdb")
```

### Edge types for proteins
| Type | Description |
|------|-------------|
| `sequential` | Adjacent residues in sequence |
| `radius` | Cα within distance cutoff |
| `knn` | K nearest neighbors in 3D |
| `contact` | Heavy-atom distance < 8 Å |

---

## Graph (generic)

```python
from torchdrug import data
import torch

# Construct manually
edge_list = torch.tensor([[0,1],[1,2],[2,0]])
node_feat = torch.randn(3, 10)
graph = data.Graph(edge_list=edge_list, node_feature=node_feat, num_node=3)

# Manipulate
g_undir = graph.undirected()
sub = graph.node_mask(torch.tensor([True, True, False]))
```

---

## PackedGraph (batching)

DataLoader handles batching automatically — graphs of different sizes are packed into one disconnected graph:

```python
from torch.utils.data import DataLoader
from torchdrug import datasets

dataset = datasets.BBBP("~/datasets/")
loader = DataLoader(dataset, batch_size=32, shuffle=True)

for batch in loader:
    # batch["graph"] is a PackedGraph
    graph = batch["graph"]
    print(graph.num_node)   # total nodes across 32 molecules
    print(graph.graph_ind)  # which graph each node belongs to
```

---

## Model Interface

All TorchDrug models follow this signature:

```python
def forward(self, graph, input, all_loss=None, metric=None):
    """
    graph  : PackedGraph (batched)
    input  : node features [num_node, input_dim]
    returns: dict with "node_feature" and/or "graph_feature"
    """
    ...
    return {"node_feature": node_out, "graph_feature": graph_out}
```

Every model must expose:
- `model.input_dim` — expected node feature dimension
- `model.output_dim` — output embedding dimension

```python
from torchdrug import models, datasets

dataset = datasets.BBBP("~/datasets/")
model = models.GIN(
    input_dim=dataset.node_feature_dim,    # must match!
    hidden_dims=[256, 256, 256],
    edge_input_dim=dataset.edge_feature_dim,
    batch_norm=True,
    readout="mean"     # "sum" | "mean" | "max"
)
print(model.output_dim)  # 256 (last hidden_dim)
```

---

## Task Interface

```python
from torchdrug import tasks, models

model = models.GIN(input_dim=..., hidden_dims=[256, 256])
task  = tasks.PropertyPrediction(
    model,
    task=dataset.tasks,          # list of task names
    criterion="bce",             # "bce" | "ce" | "mse" | "mae"
    metric=["auroc", "auprc"],   # evaluation metrics
    num_mlp_layer=2
)
```

Tasks implement:

| Method | Purpose |
|--------|---------|
| `task(batch)` | Forward → returns scalar loss |
| `task.predict(batch)` | → predictions tensor |
| `task.target(batch)` | → ground truth tensor |
| `task.evaluate(pred, target)` | → metrics dict |
| `task.preprocess(train, valid, test)` | optional dataset prep |

---

## Standard Training Loop

```python
import torch
from torch.utils.data import DataLoader
from torchdrug import datasets, models, tasks

dataset = datasets.BBBP("~/datasets/")
train_set, valid_set, test_set = dataset.split()  # scaffold split by default

model = models.GIN(input_dim=dataset.node_feature_dim,
                   hidden_dims=[256, 256, 256],
                   edge_input_dim=dataset.edge_feature_dim,
                   batch_norm=True, readout="mean")

task = tasks.PropertyPrediction(model, task=dataset.tasks,
                                 criterion="bce", metric=["auroc", "auprc"])

optimizer = torch.optim.Adam(task.parameters(), lr=1e-3)
train_loader = DataLoader(train_set, batch_size=32, shuffle=True)
valid_loader = DataLoader(valid_set, batch_size=64)

for epoch in range(100):
    task.train()
    for batch in train_loader:
        loss = task(batch)
        optimizer.zero_grad(); loss.backward(); optimizer.step()

    task.train(mode=False)   # switch to evaluation/inference mode
    preds, targets = [], []
    with torch.no_grad():
        for batch in valid_loader:
            preds.append(task.predict(batch))
            targets.append(task.target(batch))
    metrics = task.evaluate(torch.cat(preds), torch.cat(targets))
    print(f"Epoch {epoch}: {metrics}")
```

---

## Serialization (core.Configurable)

```python
from torchdrug import core, models
import torch

model = models.GIN(input_dim=10, hidden_dims=[256, 256])

# Save / load configuration (architecture only)
config = model.config_dict()
model2 = core.Configurable.load_config_dict(config)

# Save / load weights
torch.save(task.state_dict(), "task.pth")
task.load_state_dict(torch.load("task.pth"))
```

---

## Transforms

```python
from torchdrug import transforms, datasets

# Add virtual node connected to all atoms (improves global readout)
transform = transforms.Compose([
    transforms.VirtualNode(),
    transforms.VirtualEdge(),
])
dataset = datasets.BBBP("~/datasets/", transform=transform)

# Protein: truncate long sequences
transform = transforms.TruncateProtein(max_length=500)
```

---

## PyTorch Lightning Integration

```python
import pytorch_lightning as pl
import torch

class LightningTask(pl.LightningModule):
    def __init__(self, task):
        super().__init__()
        self.task = task

    def training_step(self, batch, batch_idx):
        return self.task(batch)

    def validation_step(self, batch, batch_idx):
        pred   = self.task.predict(batch)
        target = self.task.target(batch)
        return {"pred": pred, "target": target}

    def validation_epoch_end(self, outputs):
        preds   = torch.cat([o["pred"]   for o in outputs])
        targets = torch.cat([o["target"] for o in outputs])
        self.log_dict(self.task.evaluate(preds, targets))

    def configure_optimizers(self):
        return torch.optim.Adam(self.parameters(), lr=1e-3)
```

---

## Loss Functions & Metrics Reference

| Criterion string | Use for |
|-----------------|---------|
| `"bce"` | Binary classification (single or multi-label) |
| `"ce"` | Multi-class classification |
| `"mse"` | Regression |
| `"mae"` | Regression (robust to outliers) |

| Metric string | Domain |
|---------------|--------|
| `"auroc"` | Classification (imbalanced) |
| `"auprc"` | Classification (very imbalanced) |
| `"acc"` | Classification (balanced) |
| `"mae"`, `"rmse"`, `"r2"` | Regression |
| `"mr"`, `"mrr"`, `"hits@K"` | Knowledge graph ranking |

---

## Common Pitfalls

| Problem | Fix |
|---------|-----|
| Dimension mismatch | Check `dataset.node_feature_dim` == `model.input_dim` |
| Data leakage | Use scaffold split, not random, for drug discovery |
| Poor performance | Add `VirtualNode` transform; increase `hidden_dims` depth |
| Memory error | Reduce `batch_size`; use gradient accumulation |
| Generated molecules invalid | Post-validate with RDKit `Chem.MolFromSmiles` |
| Forgetting `input_dim` / `output_dim` | Always define these in custom models |
