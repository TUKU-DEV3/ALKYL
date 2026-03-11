# VAE and Graph-Based Generation — JT-VAE, Latent BO, TorchDrug

## Variational Autoencoders for Molecules

VAEs learn a continuous latent space Z where: encode(mol) → z ~ N(μ, σ²), decode(z) → mol.

**Advantage over LMs**: smooth latent space enables gradient-based optimization and Bayesian optimization (BO) in latent space.

```
Molecule → Encoder → (μ, σ) → z = μ + ε·σ → Decoder → Molecule
             ↑                                    ↑
         Graph/SMILES                       Tree + Graph decode
         encoder                            (JT-VAE)
```

## Junction Tree VAE (JT-VAE)

JT-VAE (Jin et al., 2018) represents molecules as a junction tree of ring/functional group vocabulary, then reconstructs at the graph level. Achieves ~100% valid reconstructions (vs. ~80% for SMILES VAE).

### Junction tree decomposition

```python
# Install: pip install rdkit torchdrug
# JT-VAE original: https://github.com/wengong-jin/icml18-jtnn
# TorchDrug implements JT-VAE natively

from torchdrug import data as td_data, models, tasks
from torchdrug.data import MoleculeDataset

# 1. Load dataset
dataset = MoleculeDataset()
dataset.load_smiles(train_smiles, targets={}, verbose=True)

# 2. JT-VAE model
model = models.JTVAE(
    hidden_dim=450,
    latent_dim=56,           # latent space dimension
    depth=3,                 # message passing depth
)

# 3. Generation task
task = tasks.VaeGeneration(
    model,
    atom_feature="default",
    bond_feature="default",
)

# 4. Sample from latent space
with torch.no_grad():
    smiles_list = task.generate(num_sample=100)
```

### Manual JT-VAE encoding (for latent BO)

```python
import torch

# Encode a molecule to latent vector
def encode_molecule(smiles, model, task):
    mol = td_data.Molecule.from_smiles(smiles)
    batch = td_data.PackedMolecule.pack([mol])
    with torch.no_grad():
        z_mean, z_log_var = task.model.encode(batch)
    return z_mean.squeeze(0)  # latent vector

# Decode a latent vector to SMILES
def decode_latent(z, task, n_trials=10):
    with torch.no_grad():
        # Decode attempts (VAE decoder is stochastic)
        smiles_list = []
        for _ in range(n_trials):
            smi = task.model.decode(z.unsqueeze(0))
            if smi: smiles_list.append(smi)
    return smiles_list
```

## Latent Space Bayesian Optimization (LSBO)

The key application of VAEs in drug design: optimize a property f(z) in continuous latent space, then decode the optimal z back to a molecule.

```python
# pip install botorch gpytorch

import torch
from botorch.models import SingleTaskGP
from botorch.fit import fit_gpytorch_mll
from botorch.acquisition import ExpectedImprovement
from botorch.optim import optimize_acqf
from gpytorch.mlls import ExactMarginalLogLikelihood

def latent_bo_loop(
    encode_fn,      # smiles → z tensor (latent)
    decode_fn,      # z tensor → smiles (or None)
    oracle_fn,      # smiles → float score (property to maximize)
    seed_smiles,    # list of initial molecules
    n_iterations=20,
    n_candidates=10,
    latent_dim=56,
):
    # 1. Initialize with seed molecules
    Z_train = torch.stack([encode_fn(s) for s in seed_smiles])
    Y_train = torch.tensor([[oracle_fn(s)] for s in seed_smiles], dtype=torch.float64)

    best_score = Y_train.max().item()
    history = []

    for iteration in range(n_iterations):
        # 2. Fit GP surrogate on (Z, Y)
        gp = SingleTaskGP(Z_train.double(), Y_train)
        mll = ExactMarginalLogLikelihood(gp.likelihood, gp)
        fit_gpytorch_mll(mll)

        # 3. Optimize acquisition function
        ei = ExpectedImprovement(model=gp, best_f=Y_train.max())
        bounds = torch.stack([Z_train.min(0).values, Z_train.max(0).values])

        z_candidate, _ = optimize_acqf(
            acq_function=ei,
            bounds=bounds.double(),
            q=n_candidates,
            num_restarts=10,
            raw_samples=512,
        )

        # 4. Decode, score, and add to training set
        new_smiles = [decode_fn(z) for z in z_candidate]
        valid = [(s, z) for s, z in zip(new_smiles, z_candidate) if s is not None]

        for smi, z in valid:
            score = oracle_fn(smi)
            Z_train = torch.cat([Z_train, z.unsqueeze(0)], dim=0)
            Y_train = torch.cat([Y_train, torch.tensor([[score]])], dim=0)
            if score > best_score:
                best_score = score
                history.append((iteration, smi, score))
                print(f"Iter {iteration}: {smi} → score={score:.3f}")

    return history
```

## SMILES VAE (simpler alternative to JT-VAE)

When JT-VAE is overkill (no graph-level constraints needed):

```python
import torch
import torch.nn as nn

class MolVAE(nn.Module):
    """SMILES-based VAE with GRU encoder/decoder."""

    def __init__(self, vocab_size, embed_dim=64, hidden_dim=256, latent_dim=128, max_len=120):
        super().__init__()
        self.latent_dim = latent_dim
        self.max_len = max_len

        # Encoder: bidirectional GRU
        self.embed = nn.Embedding(vocab_size, embed_dim)
        self.encoder_rnn = nn.GRU(embed_dim, hidden_dim, batch_first=True, bidirectional=True)
        self.fc_mu = nn.Linear(hidden_dim * 2, latent_dim)
        self.fc_logvar = nn.Linear(hidden_dim * 2, latent_dim)

        # Decoder: unidirectional GRU
        self.z2h = nn.Linear(latent_dim, hidden_dim)
        self.decoder_rnn = nn.GRU(embed_dim + latent_dim, hidden_dim, batch_first=True)
        self.head = nn.Linear(hidden_dim, vocab_size)

    def encode(self, x):
        emb = self.embed(x)
        _, h = self.encoder_rnn(emb)              # h: (2, B, H)
        h = torch.cat([h[0], h[1]], dim=-1)       # concat fwd+bwd: (B, 2H)
        return self.fc_mu(h), self.fc_logvar(h)

    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std

    def decode(self, z, target=None, teacher_forcing=True):
        B = z.size(0)
        h = self.z2h(z).unsqueeze(0)              # (1, B, H)
        inp = torch.zeros(B, 1, dtype=torch.long, device=z.device)  # [BOS] token
        outputs = []
        for t in range(self.max_len):
            emb = torch.cat([self.embed(inp), z.unsqueeze(1)], dim=-1)
            out, h = self.decoder_rnn(emb, h)
            logits = self.head(out)
            outputs.append(logits)
            inp = target[:, t:t+1] if (teacher_forcing and target is not None) \
                  else logits.argmax(-1)
        return torch.cat(outputs, dim=1)

    def forward(self, x):
        mu, logvar = self.encode(x)
        z = self.reparameterize(mu, logvar)
        logits = self.decode(z, target=x)
        return logits, mu, logvar

def vae_loss(logits, target, mu, logvar, beta=1.0):
    """ELBO loss with β-VAE weighting."""
    # Reconstruction loss
    B, T, V = logits.shape
    recon = nn.functional.cross_entropy(logits.reshape(-1, V), target.reshape(-1), ignore_index=0)
    # KL divergence: -0.5 * Σ(1 + logvar - μ² - exp(logvar))
    kl = -0.5 * torch.mean(1 + logvar - mu.pow(2) - logvar.exp())
    return recon + beta * kl
```

### β-VAE scheduling (critical for training)

```python
# Warm up β from 0 → 1 over first N steps to prevent posterior collapse
def get_beta(step, warmup_steps=10000):
    return min(1.0, step / warmup_steps)

# Training step
for step, batch in enumerate(loader):
    beta = get_beta(step)
    logits, mu, logvar = model(batch)
    loss = vae_loss(logits, batch, mu, logvar, beta=beta)
    ...
```

## TorchDrug Graph Generative Models

For graph-native generation (no string intermediates). See `torchdrug` skill for full implementation.

| Model | Architecture | Notes |
|-------|-------------|-------|
| **GCPN** | Graph Convolutional Policy Network (RL + GNN) | Property-guided; atom/bond actions |
| **GraphAF** | Autoregressive normalizing flow on graphs | Exact likelihood; fast sampling |
| **GraphDF** | Discrete flow | Better validity than GraphAF |
| **JT-VAE** | Junction tree VAE | Best validity; slowest |

```python
# Quick GCPN setup (TorchDrug)
from torchdrug import models, tasks

model = models.RGCN(
    input_dim=td_data.Molecule.node_feature_dim,
    hidden_dims=[64, 64, 64],
    num_relation=td_data.Molecule.edge_feature_dim,
)
task = tasks.GCPNGeneration(
    model,
    atom_types=["C", "N", "O", "S", "F", "Cl"],
    max_edge_unroll=12,
    max_node=38,
    criterion="nll",
)
# → see torchdrug/references/molecular-generation.md for full training loop
```

## Latent Space Interpolation

Smooth interpolation between two molecules via VAE latent space:

```python
def interpolate_molecules(smi_a, smi_b, encode_fn, decode_fn, n_steps=10):
    """Generate molecules interpolating between A and B in latent space."""
    z_a = encode_fn(smi_a)
    z_b = encode_fn(smi_b)
    results = []
    for alpha in torch.linspace(0, 1, n_steps):
        z_interp = (1 - alpha) * z_a + alpha * z_b
        smi = decode_fn(z_interp)
        results.append((float(alpha), smi))
    return results

# Applications:
# - Scaffold hopping between two known actives
# - Property gradient exploration (LogP, MW along interpolation path)
# - Fragment merging (two fragments → merged core)
```

## When to Use VAE vs. RL vs. LM

| Situation | Recommendation |
|-----------|---------------|
| Bayesian optimization, few oracle calls | VAE + LSBO (JT-VAE preferred) |
| Multi-property optimization, many oracle calls | REINVENT 4 (RL) |
| Fast focused library generation | SELFIES + fine-tuned GPT |
| Graph-native, custom atom constraints | TorchDrug GCPN/GraphAF |
| Scaffold interpolation, analogue design | MolVAE (SMILES) |

## Key Pitfalls

- **Posterior collapse in VAE**: KL term → 0, encoder ignores input; fix with β-VAE warmup
- **JT-VAE tree vocabulary** must cover training set; out-of-vocabulary fragments = decode failure
- **BO in high-dimensional latent**: latent_dim > 64 → BO struggles (curse of dimensionality); use PCA projection or lower latent_dim
- **Decoder stochasticity**: same z → different SMILES on each call; average over multiple decode attempts
- **TorchDrug GCPN convergence**: sensitive to reward scaling; normalize all rewards to [0, 1] before training
