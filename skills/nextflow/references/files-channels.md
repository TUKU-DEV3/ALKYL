# Nextflow Files, Channels & Operators

## Channel Factories

```nextflow
// Single value
ch = channel.of('aspirin')
ch = channel.value(file('/data/ref.fa'))   // Value — reused forever

// Multiple values
ch = channel.of('CC(=O)O', 'c1ccccc1', 'CCO')

// From a list
smiles = ['CC(=O)O', 'c1ccccc1']
ch = channel.fromList(smiles)

// Files matching glob
ch = channel.fromPath('data/*.sdf')
ch = channel.fromPath('data/**/*.mol2', hidden: false)

// File pairs (common in NGS/MD)
ch = channel.fromFilePairs('data/*_{1,2}.fastq')
// → ['sample1', [file('sample1_1.fastq'), file('sample1_2.fastq')]], ...

// From a CSV/TSV samplesheet
ch = channel.fromPath('samplesheet.csv')
    | splitCsv(header: true)
    | map { row -> tuple(row.id, file(row.path), row.condition) }

// From a value file (one value per line)
ch = channel.fromPath('ids.txt')
    | splitText(strip: true)

// Timer
ch = channel.interval('1s')   // emits 0, 1, 2, … every second

// Watching filesystem (live)
ch = channel.watchPath('data/*.sdf', 'create,modify')
```

---

## Essential Operators

### Transform

```nextflow
ch.map { v -> v * 2 }
ch.flatMap { v -> [v, v.toUpperCase()] }   // one-to-many
ch.buffer(size: 10, remainder: true)        // emit lists of N items
ch.collate(3)                               // alias for buffer
ch.splitFasta()                             // one sequence per item
ch.splitText(by: 100, file: true)          // N lines per chunk
ch.splitCsv(header: true)                  // CSV → map per row
ch.splitJson()                             // JSON array → items
```

### Filter / Select

```nextflow
ch.filter { it > 5 }
ch.filter(String)                          // type filter
ch.filter ~/.*\.sdf/                       // regex filter
ch.unique()                                // deduplicate
ch.distinct()                              // consecutive dedup
ch.first()                                 // Value with first item
ch.last()                                  // Value with last item
ch.take(10)                                // first N items
ch.until { v -> v == 'STOP' }
```

### Collect / Reduce

```nextflow
ch.collect()                   // all items → list (Value)
ch.toList()                    // same as collect
ch.toSortedList()
ch.count()                     // Value with count
ch.sum()
ch.min()
ch.max()
ch.reduce(0) { acc, v -> acc + v }
```

### Combine / Join

```nextflow
// Join on first element (key)
bam_ch.join(bai_ch)
// → [key, bam, bai]

// Join on a specific index
ch1.join(ch2, by: [0, 1])   // match on first two elements

// Cartesian product
mols_ch.combine(params_ch)
// → every mol paired with every param set

// Concat — preserve order
ch1.concat(ch2, ch3)   // ch2 only starts after ch1 closes

// Mix — no order guarantee
ch1.mix(ch2, ch3)

// Cross — like combine but based on key
ch1.cross(ch2)   // first element must be key
```

### Group

```nextflow
// Group [key, val] pairs by key
channel.of(
    ['A', 1], ['B', 2], ['A', 3]
).groupTuple()
// → ['A', [1, 3]], ['B', [2]]

// Group with size — emit when N items collected
ch.groupTuple(size: 3, remainder: true)

// Buffer into fixed-size windows
ch.buffer(size: 4)
```

### Branch (conditional routing)

```nextflow
ch.branch {
    pass: it > 5
    fail: it <= 5
}
// Access: ch.pass, ch.fail

// With map transform per branch
ch.branch {
    large: it.size() > 1e6
        return it
    small: true
        return it
}
```

---

## Working with Files

### Path attributes

```nextflow
def p = file('/data/molecules/aspirin.sdf.gz')

p.name           // 'aspirin.sdf.gz'
p.baseName       // 'aspirin.sdf'        (without last extension)
p.extension      // 'gz'
p.parent         // /data/molecules
p.simpleName     // 'aspirin'            (without all extensions)
p.nameWithoutExtension // same as baseName

p.exists()
p.isFile()
p.isDirectory()
p.size()         // bytes
p.text           // read entire file as String
p.readLines()    // List<String>
p.bytes          // byte[]
```

### File operations in scripts

```nextflow
// Read
content = file('molecules.txt').text
lines   = file('molecules.txt').readLines()

// Write
outFile = file('results.txt')
outFile.text = 'MW,LogP\n'
outFile << '180.16,1.19\n'   // append

// Line by line (memory-efficient)
file('large.sdf').eachLine { line ->
    if (line =~ /^M  END/) { /* ... */ }
}

// Copy / move
file('source.sdf').copyTo('dest.sdf')
file('source.sdf').moveTo('archive/')

// Path building
def base = file('/data')
def mol  = base / 'molecules' / 'aspirin.sdf'   // /data/molecules/aspirin.sdf
```

### Remote files

```nextflow
// HTTP/FTP
pdb = file('https://files.rcsb.org/download/1ABC.pdb')

// S3
mol  = file('s3://my-bucket/data/mol.sdf')
ref  = file('s3://my-bucket/refs/genome.fa')

// Azure Blob
f = file('az://my-container/data/file.csv')

// Google Cloud Storage
g = file('gs://my-bucket/data/file.csv')
```

---

## Samplesheet Pattern (nf-core style)

```nextflow
// samplesheet.csv:
// id,smiles,activity
// cpd001,CC(=O)O,1
// cpd002,c1ccccc1,0

workflow {
    ch_input = channel.fromPath(params.samplesheet)
        | splitCsv(header: true)
        | map { row ->
            def meta = [id: row.id, activity: row.activity.toInteger()]
            [meta, row.smiles]
          }

    PROCESS_MOL(ch_input)
}

process PROCESS_MOL {
    tag "${meta.id}"

    input:
    tuple val(meta), val(smiles)

    output:
    tuple val(meta), path("${meta.id}.csv"), emit: results

    script:
    """
    python calc.py --smiles '${smiles}' --id ${meta.id} --out ${meta.id}.csv
    """
}
```

---

## publishDir Patterns

```nextflow
process MY_PROCESS {
    // Copy all outputs to results/
    publishDir "${params.outdir}", mode: 'copy'

    // Different dirs per file type
    publishDir "${params.outdir}/logs",   pattern: '*.log',  mode: 'copy'
    publishDir "${params.outdir}/output", pattern: '*.csv',  mode: 'copy'

    // Dynamic path based on metadata
    publishDir { "${params.outdir}/${meta.id}" }, mode: 'copy'

    // Save space — symlink instead of copy
    publishDir "${params.outdir}", mode: 'symlink'

    // Compress on publish
    publishDir "${params.outdir}", saveAs: { fn -> fn.endsWith('.log') ? null : fn }
}
```

| Mode | Behavior |
|------|----------|
| `copy` | Full copy — safest |
| `move` | Move file out of workDir — breaks resume |
| `link` | Hard link (same filesystem only) |
| `symlink` | Symbolic link (default) |
| `copyNoFollow` | Copy without following symlinks |

---

## Splitting Large Inputs

```nextflow
// Split a large FASTA into chunks of 1000 sequences
channel.fromPath('proteome.fa')
    | splitFasta(by: 1000, file: true)
    | BLAST_SEARCH

// Split SDF into batches of 500 molecules
channel.fromPath('library.sdf')
    | splitText(by: 500 * 5, file: true)   // ~500 molecules (5 lines/mol approx)
    | VIRTUAL_SCREEN
```
