# Template Search and Sequence Alignment

Finding the right template is the most critical step in comparative modeling. A better template (higher identity, more coverage) directly yields a better model.

---

## Template Reliability Rules of Thumb

| Sequence identity | Expected model quality |
|------------------|----------------------|
| > 50% | Excellent — backbone reliable, sidechains mostly correct |
| 30–50% | Good — backbone ~1 Å RMSD, some loop errors |
| 25–30% | Twilight zone — backbone ~2–3 Å, loops unreliable |
| < 25% | Threading/fold recognition only — use AlphaFold2 instead |
| > 90% | Near-identical — model ≈ template (use for conformational variants) |

---

## HHblits / HHpred (Best for Remote Homologs)

HHblits builds an MSA via profile–profile comparison — far more sensitive than BLAST at <30% identity.

```bash
# Install: conda install -c bioconda hhsuite
# Database: download UniRef30 from https://wwwuser.gwdg.de/~compbiol/

# Step 1: Build MSA for query
hhblits -i query.fasta \
        -d /data/UniRef30_2023_02/UniRef30_2023_02 \
        -o query.hhr \
        -oa3m query.a3m \
        -n 3 \           # iterations
        -e 1e-3 \        # E-value cutoff
        -cpu 8

# Step 2: Search PDB (for template identification)
hhblits -i query.a3m \
        -d /data/pdb70/pdb70 \
        -o query_vs_pdb.hhr \
        -n 1 \
        -e 1e-3 \
        -cpu 8

# Parse HHR output for top templates
python - <<'EOF'
import re

def parse_hhr(hhr_file, top_n=10):
    templates = []
    with open(hhr_file) as f:
        in_hits = False
        for line in f:
            if line.startswith(' No '):
                in_hits = True
                continue
            if in_hits and re.match(r'\s+\d+\s+', line):
                parts = line.split()
                if len(parts) >= 9:
                    templates.append({
                        'rank': int(parts[0]),
                        'pdb_id': parts[1][:4].lower(),
                        'chain': parts[1][4] if len(parts[1]) > 4 else 'A',
                        'prob': float(parts[2]),    # probability (%)
                        'e_value': parts[3],
                        'identity': float(parts[7].replace('%', '')),
                        'coverage': parts[8],
                    })
            if in_hits and line.strip() == '':
                break

    return templates[:top_n]

templates = parse_hhr('query_vs_pdb.hhr')
for t in templates:
    print(f"Rank {t['rank']}: {t['pdb_id']}{t['chain']} "
          f"prob={t['prob']:.0f}% id={t['identity']:.0f}%")
EOF
```

---

## Biopython BLAST (Quick Template Search)

```python
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO

def blast_against_pdb(fasta_file, max_hits=10, e_threshold=0.001):
    """Run BLAST against PDB database and return top hits."""
    with open(fasta_file) as f:
        query = SeqIO.read(f, 'fasta')

    print(f"BLASTing {query.id} ({len(query.seq)} aa) against PDB...")
    result_handle = NCBIWWW.qblast(
        'blastp', 'pdb', str(query.seq),
        hitlist_size=max_hits,
        expect=e_threshold,
    )

    blast_records = list(NCBIXML.parse(result_handle))
    hits = []
    for record in blast_records:
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                pdb_id = alignment.title.split('|')[3][:4].lower()
                chain = alignment.title.split('|')[3][4] if len(alignment.title.split('|')[3]) > 4 else 'A'
                identity = hsp.identities / hsp.align_length * 100
                coverage = hsp.align_length / len(query.seq) * 100
                hits.append({
                    'pdb_id': pdb_id,
                    'chain': chain,
                    'identity': identity,
                    'coverage': coverage,
                    'e_value': hsp.expect,
                    'score': hsp.score,
                })

    return sorted(hits, key=lambda x: x['e_value'])

hits = blast_against_pdb('query.fasta')
for h in hits:
    print(f"{h['pdb_id']}{h['chain']}: id={h['identity']:.0f}% "
          f"cov={h['coverage']:.0f}% E={h['e_value']:.2e}")
```

---

## PDB Download (Programmatic)

```python
from Bio.PDB import PDBList
import os

def download_templates(pdb_ids, output_dir='templates/'):
    """Download PDB structures for templates."""
    os.makedirs(output_dir, exist_ok=True)
    pdbl = PDBList()

    downloaded = []
    for pdb_id in pdb_ids:
        path = pdbl.retrieve_pdb_file(
            pdb_id.lower(),
            file_type='pdb',
            pdir=output_dir,
        )
        # Rename from ent format to .pdb
        new_path = os.path.join(output_dir, f'{pdb_id.lower()}.pdb')
        if os.path.exists(path) and path != new_path:
            os.rename(path, new_path)
        downloaded.append(new_path)
        print(f"Downloaded: {pdb_id} → {new_path}")

    return downloaded

template_files = download_templates(['1TQN', '2V00', '3EML'])
```

---

## Sequence Alignment (Biopython)

```python
from Bio import Align, SeqIO
from Bio.Seq import Seq

def pairwise_align(query_seq, template_seq, mode='global'):
    """Pairwise alignment with BLOSUM62 scoring."""
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = Align.substitution_matrices.load('BLOSUM62')
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.mode = mode  # 'global' or 'local'

    alignments = aligner.align(query_seq, template_seq)
    best = alignments[0]

    # Compute statistics
    n_id = sum(a == b and a != '-' for a, b in zip(*best))
    n_aligned = sum(1 for a, b in zip(*best) if a != '-' and b != '-')
    identity = n_id / n_aligned * 100 if n_aligned > 0 else 0

    print(f"Score: {best.score:.1f}")
    print(f"Identity: {identity:.1f}% ({n_id}/{n_aligned})")
    print(best)
    return best

# From FASTA files
with open('query.fasta') as f:
    query = str(SeqIO.read(f, 'fasta').seq)
with open('template.fasta') as f:
    template = str(SeqIO.read(f, 'fasta').seq)

alignment = pairwise_align(query, template)
```

---

## PIR Alignment Format (MODELLER)

MODELLER requires a specific PIR-format alignment file. Errors here are the most common source of MODELLER crashes.

```python
def write_pir_alignment(query_id, query_seq,
                         template_id, template_seq,
                         template_pdb, template_chain,
                         start_res, end_res,
                         outfile='alignment.ali'):
    """
    Write PIR alignment file for MODELLER.

    PIR format rules:
    - >P1;name: sequence identifier line
    - structureX:pdb:start:chain:end:chain:org:name:resolution:rfactor
    - sequence ends with '*'
    - gaps represented as '-'
    """
    with open(outfile, 'w') as f:
        # Template (structure)
        f.write(f">P1;{template_id}\n")
        f.write(f"structureX:{template_pdb}:"
                f"{start_res}:{template_chain}:"
                f"{end_res}:{template_chain}:undefined:undefined:-1.00:-1.00\n")
        # Wrap sequence at 75 chars
        for i in range(0, len(template_seq), 75):
            f.write(template_seq[i:i+75] + '\n')
        if not template_seq.endswith('*'):
            f.write('*\n')
        f.write('\n')

        # Query (sequence to model)
        f.write(f">P1;{query_id}\n")
        f.write(f"sequence:{query_id}:::::::0.00: 0.00\n")
        for i in range(0, len(query_seq), 75):
            f.write(query_seq[i:i+75] + '\n')
        if not query_seq.endswith('*'):
            f.write('*\n')

    print(f"PIR alignment written: {outfile}")
    return outfile

# Example
write_pir_alignment(
    query_id='target',
    query_seq='MKVLWAALLVTFLAGCQ---AQIRDEFF*',
    template_id='1abc_A',
    template_seq='MKVLWAALLVTFLAGCQSTAQIRDEFF*',
    template_pdb='1abc',
    template_chain='A',
    start_res=1,
    end_res=28,
)
```

---

## Structural Alignment (Biopython Superimposer)

```python
from Bio.PDB import PDBParser, Superimposer, PDBIO

def structural_align(ref_pdb, mobile_pdb, output='aligned.pdb'):
    """Align mobile structure to reference by Cα atoms."""
    parser = PDBParser(QUIET=True)
    ref = parser.get_structure('ref', ref_pdb)
    mobile = parser.get_structure('mob', mobile_pdb)

    # Get Cα atoms (matching residues)
    ref_cas = [r['CA'] for r in ref.get_residues()
               if 'CA' in r and r.id[0] == ' ']
    mob_cas = [r['CA'] for r in mobile.get_residues()
               if 'CA' in r and r.id[0] == ' ']

    n = min(len(ref_cas), len(mob_cas))
    ref_cas, mob_cas = ref_cas[:n], mob_cas[:n]

    sup = Superimposer()
    sup.set_atoms(ref_cas, mob_cas)
    sup.apply(list(mobile.get_atoms()))

    print(f"Cα RMSD: {sup.rms:.3f} Å ({n} residues)")

    io = PDBIO()
    io.set_structure(mobile)
    io.save(output)
    return sup.rms

rmsd = structural_align('template.pdb', 'model.pdb')
```

---

## Multi-Template Selection Strategy

```python
def select_templates(hits, max_templates=5,
                      min_identity=25.0, min_coverage=60.0):
    """
    Select complementary templates:
    - High identity templates for core regions
    - Additional templates for poorly covered loops
    """
    # Filter by minimum thresholds
    viable = [h for h in hits
              if h['identity'] >= min_identity
              and h['coverage'] >= min_coverage]

    if not viable:
        print("WARNING: No templates above thresholds — consider AlphaFold2")
        return []

    # Rank: primary = highest identity × coverage
    viable.sort(key=lambda h: -h['identity'] * h['coverage'] / 100)

    # Select up to max_templates (diverse coverage)
    selected = viable[:max_templates]

    print(f"Selected {len(selected)} templates:")
    for i, t in enumerate(selected):
        print(f"  {i+1}. {t['pdb_id']}{t['chain']} "
              f"id={t['identity']:.0f}% cov={t['coverage']:.0f}%")

    return selected
```
