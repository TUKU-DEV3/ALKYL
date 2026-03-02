# Nextflow Core Concepts

## Processes

A process wraps a script or command. It is the basic unit of work.

```nextflow
process MY_PROCESS {
    // Directives (optional)
    label 'big_mem'
    cpus  4
    memory '8 GB'

    input:
    val  sample_id      // plain value
    path reads          // staged file

    output:
    path "*.bam",     emit: bam
    path "*.log",     emit: log

    script:
    """
    aligner --threads ${task.cpus} -i ${reads} -o ${sample_id}.bam 2> ${sample_id}.log
    """
}
```

### Script types

| Section | Use case |
|---------|----------|
| `script:` | Default — Bash (or other via shebang) |
| `exec:` | Groovy/Java code — no subprocess, no container |
| `shell:` | *(deprecated 24.11)* Use `!{var}` for NF vars, `$BASH_VAR` for Bash |
| `stub:` | Dry-run skeleton — activated with `-stub` flag |

**Multi-language example:**
```nextflow
process PYTHON_CALC {
    input:
    path sdf_file

    output:
    path "descriptors.csv"

    script:
    """
    #!/usr/bin/env python3
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    import csv

    suppl = Chem.SDMolSupplier('${sdf_file}')
    with open('descriptors.csv', 'w') as f:
        w = csv.writer(f)
        w.writerow(['smiles', 'mw', 'logp'])
        for mol in suppl:
            if mol:
                w.writerow([Chem.MolToSmiles(mol),
                            Descriptors.MolWt(mol),
                            Descriptors.MolLogP(mol)])
    """
}
```

**Variable escaping — critical:**
```nextflow
process BLAST {
    script:
    """
    # \$DB  → Bash variable (escape with backslash)
    # $query → Nextflow variable (no escape)
    blastp -db \$BLASTDB -query ${query} -out results.txt
    """
}
```

### Conditional scripts

```nextflow
process ALIGN {
    input:
    path seqs
    val  aligner   // 'mafft' or 'clustalo'

    script:
    if (aligner == 'mafft')
        """
        mafft --auto $seqs > aligned.fa
        """
    else if (aligner == 'clustalo')
        """
        clustalo -i $seqs -o aligned.fa
        """
    else
        error "Unknown aligner: ${aligner}"
}
```

---

## Input Qualifiers

| Qualifier | Receives | Use for |
|-----------|----------|---------|
| `val x` | Any value | Strings, numbers, maps |
| `path f` | File/dir | Staged into task workdir |
| `env VAR` | String | Set env variable in script |
| `stdin` | String | Piped to process stdin |
| `tuple val(id), path(f)` | Multiple | Group related inputs (sample metadata + file) |
| `each x` | List/channel | Run process once per element |

**Tuple pattern (most common in bioinformatics):**
```nextflow
process PROCESS_SAMPLE {
    input:
    tuple val(meta), path(reads)

    script:
    """
    tool -s ${meta.id} -i ${reads} -o ${meta.id}.out
    """
}
```

**Fix file name on staging:**
```nextflow
input:
path query, name: 'input.fa'   // always staged as input.fa
```

**`each` — cartesian product:**
```nextflow
process GRID_SEARCH {
    input:
    path mol
    each cutoff   // runs once per cutoff value

    script:
    """
    dock --mol $mol --cutoff $cutoff > ${cutoff}.out
    """
}
workflow {
    mols     = channel.fromPath('*.sdf')
    cutoffs  = [5.0, 7.5, 10.0]
    GRID_SEARCH(mols, cutoffs)
}
```

---

## Output Qualifiers

| Qualifier | Captures | Notes |
|-----------|----------|-------|
| `path "*.txt"` | File glob | Most common |
| `path "dir/"` | Directory | Whole directory |
| `val x` | Computed value | From `exec:` or `when:` |
| `env VAR` | Env variable | After script runs |
| `stdout` | Process stdout | Capture as string |
| `eval 'cmd'` | Shell command output | Run after script |
| `tuple val(x), path(f)` | Multiple | Pair metadata with file |

**Naming outputs (emit):**
```nextflow
output:
path "*.bam",  emit: bam
path "*.bai",  emit: index
path "*.log",  emit: log
```

Access in workflow: `MY_PROCESS.out.bam`

---

## Process Directives (key ones)

```nextflow
process MY_PROCESS {
    // Resources
    cpus   8
    memory '32 GB'
    time   '4h'
    disk   '100 GB'          // scratch space

    // Scheduling
    label  'gpu'             // match config withLabel
    queue  'bigmem'          // HPC queue name
    clusterOptions '--gres=gpu:1'   // raw scheduler flags

    // Containers
    container 'docker://nfcore/rnaseq:3.14'

    // Output management
    publishDir 'results/bam', mode: 'copy'
    publishDir 'results/qc',  pattern: '*.log', mode: 'copy'

    // Retry logic
    maxRetries 3
    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }
    memory { 8.GB * task.attempt }   // scale memory on retry

    // Cache
    cache 'lenient'    // ignore timestamps — good for NFS

    // Misc
    tag   "${meta.id}"          // label in log/report
    storeDir '/data/cache'      // persistent cache across runs
}
```

---

## Groovy / Script Language Basics

Nextflow scripts use Groovy syntax.

```nextflow
// Variables — always use def in closures
def name = 'aspirin'
def mw   = 180.16

// String interpolation
def msg = "Molecule: ${name}, MW: ${mw}"

// Lists
def smiles_list = ['CC(=O)O', 'c1ccccc1', 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C']
println smiles_list.size()          // 3
println smiles_list[0]              // CC(=O)O
println smiles_list.collect { it.length() }  // [8, 6, 37]
println smiles_list.findAll { it.length() > 10 }

// Maps (metadata pattern)
def meta = [id: 'SAMP001', group: 'treatment', replicate: 1]
println meta.id
println meta['group']

// Closures
def double_it = { x -> x * 2 }
println double_it(21)   // 42

// File paths
def p = file('/data/molecules/aspirin.sdf')
println p.name          // aspirin.sdf
println p.baseName      // aspirin
println p.extension     // sdf
println p.parent        // /data/molecules
```

---

## task Special Variable

Inside a script, `task` exposes runtime metadata:

```nextflow
script:
"""
echo "Running with ${task.cpus} CPUs and ${task.memory} RAM"
echo "Attempt number: ${task.attempt}"
"""
```

| Property | Value |
|----------|-------|
| `task.cpus` | Allocated CPUs |
| `task.memory` | Allocated memory (MemoryUnit) |
| `task.time` | Allocated time |
| `task.attempt` | Retry attempt number (1-based) |
| `task.index` | Task index in channel |
| `task.name` | Full task name |
| `task.workDir` | Absolute path to task work directory |
| `task.ext.*` | Custom properties from config `ext` scope |
