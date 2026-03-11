# SELFIES + Language Models for Molecular Generation

## SELFIES Grammar (Always-Valid Generation)

SELFIES (Self-Referencing Embedded Strings) encodes molecules as token sequences where every possible sequence decodes to a chemically valid molecule.

```python
import selfies as sf

# Core conversion
smiles = "c1ccc(cc1)C(=O)O"  # benzoic acid
selfies_str = sf.encoder(smiles)
decoded = sf.decoder(selfies_str)

# Verify roundtrip
from rdkit import Chem
assert Chem.MolFromSmiles(decoded) is not None

# Get alphabet (all valid tokens, ~800)
alphabet = list(sf.get_semantic_robust_alphabet())
stoi = {tok: i for i, tok in enumerate(["[nop]"] + alphabet)}  # [nop] = pad
itos = {i: tok for tok, i in stoi.items()}

# Tokenize a SELFIES string
tokens = list(sf.split_selfies(selfies_str))  # e.g. ["[C]", "[=O]", ...]
indices = [stoi[t] for t in tokens]

# Decode token indices back to molecule
tokens_back = [itos[i] for i in indices]
smiles_back = sf.decoder("".join(tokens_back))  # always valid
```

### SELFIES Constraints (advanced)

```python
# Customize allowed atomic valences
constraints = sf.get_preset_constraints("default")
constraints["N"] = 3      # cap nitrogen valence at 3
sf.set_semantic_constraints(constraints)

# Random SELFIES generation (always valid by construction)
import random
def random_selfies(n_tokens=20):
    tokens = random.choices(alphabet, k=n_tokens)
    return sf.decoder("".join(tokens))  # guaranteed valid SMILES
```

## SMILES Data Augmentation (Randomized SMILES)

Randomized SMILES: same molecule, different atom traversal order → more training diversity.

```python
from rdkit import Chem
import random

def randomize_smiles(smiles, n=5):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None: return []
    results = set()
    for _ in range(n * 10):
        new_atom_order = list(range(mol.GetNumAtoms()))
        random.shuffle(new_atom_order)
        new_mol = Chem.RenumberAtoms(mol, new_atom_order)
        smi = Chem.MolToSmiles(new_mol, canonical=False)
        results.add(smi)
        if len(results) == n: break
    return list(results)

# Use 10x augmented SMILES for LM training — improves validity and novelty
```

## LSTM-based SMILES/SELFIES Generator

Classic approach (still competitive for focused libraries):

```python
import torch
import torch.nn as nn

class MolLSTM(nn.Module):
    def __init__(self, vocab_size, embed_dim=128, hidden_dim=512, n_layers=3, dropout=0.1):
        super().__init__()
        self.embed = nn.Embedding(vocab_size, embed_dim, padding_idx=0)
        self.lstm = nn.LSTM(embed_dim, hidden_dim, n_layers, batch_first=True, dropout=dropout)
        self.head = nn.Linear(hidden_dim, vocab_size)
        self.dropout = nn.Dropout(dropout)

    def forward(self, x, hidden=None):
        # x: (batch, seq_len) token indices
        emb = self.dropout(self.embed(x))
        out, hidden = self.lstm(emb, hidden)
        logits = self.head(self.dropout(out))
        return logits, hidden

    @torch.no_grad()
    def sample(self, start_token, end_token, max_len=100, temperature=1.0, device='cpu'):
        self.eval()
        tokens = [start_token]
        hidden = None
        x = torch.tensor([[start_token]], device=device)
        for _ in range(max_len):
            logits, hidden = self(x, hidden)
            probs = (logits[0, -1] / temperature).softmax(-1)
            next_tok = torch.multinomial(probs, 1).item()
            tokens.append(next_tok)
            if next_tok == end_token: break
            x = torch.tensor([[next_tok]], device=device)
        return tokens

def train_epoch(model, loader, optimizer, criterion, device):
    model.train()
    total_loss = 0
    for batch in loader:
        x = batch[:, :-1].to(device)   # input: all tokens except last
        y = batch[:, 1:].to(device)    # target: all tokens except first
        logits, _ = model(x)
        loss = criterion(logits.reshape(-1, logits.size(-1)), y.reshape(-1))
        optimizer.zero_grad(); loss.backward(); optimizer.step()
        total_loss += loss.item()
    return total_loss / len(loader)
```

## GPT-style (Transformer) Molecular LM with HuggingFace

### Fine-tuning a pre-trained model on a focused library

```python
from transformers import (
    GPT2Config, GPT2LMHeadModel, GPT2TokenizerFast,
    DataCollatorForLanguageModeling, Trainer, TrainingArguments
)
from datasets import Dataset
import selfies as sf

# 1. Prepare SELFIES corpus
def smiles_to_selfies_tokens(smiles_list):
    all_tokens = []
    for smi in smiles_list:
        enc = sf.encoder(smi)
        if enc:
            tokens = list(sf.split_selfies(enc))
            all_tokens.append(" ".join(tokens))
    return all_tokens

# Build custom tokenizer from SELFIES alphabet
alphabet = list(sf.get_semantic_robust_alphabet())
special_tokens = ["[PAD]", "[BOS]", "[EOS]", "[UNK]"]
vocab = {tok: i for i, tok in enumerate(special_tokens + alphabet)}

# 2. GPT2 config (small for molecular LM)
config = GPT2Config(
    vocab_size=len(vocab),
    n_embd=256,
    n_layer=4,
    n_head=8,
    n_positions=128,
    bos_token_id=vocab["[BOS]"],
    eos_token_id=vocab["[EOS]"],
    pad_token_id=vocab["[PAD]"],
)
model = GPT2LMHeadModel(config)

# 3. Dataset prep
tokenized = smiles_to_selfies_tokens(train_smiles)
dataset = Dataset.from_dict({"text": tokenized})

# 4. Training
training_args = TrainingArguments(
    output_dir="mol_gpt", num_train_epochs=20,
    per_device_train_batch_size=64, learning_rate=3e-4,
    warmup_steps=500, weight_decay=0.01,
    logging_steps=100, save_steps=500,
)
trainer = Trainer(
    model=model, args=training_args,
    train_dataset=dataset,
    data_collator=DataCollatorForLanguageModeling(tokenizer=..., mlm=False),
)
trainer.train()
```

### Using ChemGPT (pre-trained, HuggingFace Hub)

```python
from transformers import AutoTokenizer, AutoModelForCausalLM
import torch

# ncfrey/ChemGPT-1.2B trained on ~1.1B SMILES (PubChem)
tokenizer = AutoTokenizer.from_pretrained("ncfrey/ChemGPT-1.2B")
model = AutoModelForCausalLM.from_pretrained("ncfrey/ChemGPT-1.2B")

prompt = "CC(=O)"  # acetyl group — generate amides/ketones
inputs = tokenizer(prompt, return_tensors="pt")
outputs = model.generate(
    **inputs,
    max_new_tokens=50,
    do_sample=True,
    temperature=1.0,
    top_p=0.95,
    num_return_sequences=100,
)
gen_smiles = [tokenizer.decode(o, skip_special_tokens=True) for o in outputs]

# Validity filter
from rdkit import Chem
valid = [s for s in gen_smiles if Chem.MolFromSmiles(s) is not None]
print(f"Validity: {len(valid)/len(gen_smiles):.1%}")
```

## Sampling Strategies

### Temperature sampling

```
Lower T (0.5-0.8) → conservative, reproduce training distribution
Higher T (1.0-1.3) → diverse, risk invalid/unusual structures
Optimal for molecules: 0.7 - 1.2
```

### Nucleus (Top-p) sampling

```
top_p=0.9: sample from smallest set covering 90% probability mass
Prevents degenerate tokens while keeping diversity
Recommended: temperature=0.9 + top_p=0.95
```

### Beam search (NOT recommended for molecular generation)

```
Beam search finds high-likelihood → repetitive outputs
Use temperature sampling for structural diversity
```

### Constrained generation (prefix forcing)

```python
# Fix a SMARTS-derived prefix then sample completion
# e.g., force Boc group "CC(C)(C)OC(=O)" as prefix
# Model completes with N-containing fragments

from transformers import LogitsProcessor

class PrefixForceProcessor(LogitsProcessor):
    """Force specific prefix tokens at each step."""
    def __init__(self, prefix_ids, n_forced):
        self.prefix_ids = prefix_ids
        self.n_forced = n_forced

    def __call__(self, input_ids, scores):
        step = input_ids.shape[-1]
        if step < self.n_forced:
            forced_id = self.prefix_ids[step]
            scores[:, :] = -float('inf')
            scores[:, forced_id] = 0
        return scores
```

## Transfer Learning on Focused Library

**Strategy**:
1. **Pre-train** on large diverse set (ZINC/ChEMBL, >100k molecules) → learn molecular grammar
2. **Fine-tune** on focused library (target-active compounds, 100-10k molecules) → adapt distribution
3. **Monitor FCD** against held-out set during fine-tuning (stop before memorization)

```python
# Fine-tune only top transformer blocks (last N layers)
# Freeze embedding + first 2 layers for stability
for name, param in model.named_parameters():
    if "h.0." in name or "h.1." in name or "wte" in name or "wpe" in name:
        param.requires_grad = False
# Fine-tuning time: 5-20 epochs on focused library (100-10k mols)
# Learning rate: 1e-4 to 5e-5 (lower than from-scratch)
```

## Key Pitfalls

- **SMILES LMs produce 10-30% invalid**: always filter with `Chem.MolFromSmiles`; use SELFIES for 100% validity
- **SELFIES alphabet version matters**: use `selfies >= 2.1`; old alphabets differ
- **Batch encoding**: `sf.encoder` is slow for large datasets — use `multiprocessing.Pool`
- **Overfitting on small focused libs**: add dropout + weight decay; track FCD on held-out set
- **Temperature too low → mode collapse** (same scaffold); too high → invalid SMILES
- **Randomized SMILES augmentation** can 5-10× improve novelty without extra data
