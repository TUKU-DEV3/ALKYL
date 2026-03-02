# Nextflow Pipeline Patterns

## Workflow Structure

```nextflow
// Entry workflow — one per script, no name
workflow {
    // Legacy params (DSL2 classic)
    params.input  = 'data/*.fastq'
    params.outdir = 'results'

    reads = channel.fromFilePairs(params.input)
    QC(reads)
    ALIGN(QC.out.trimmed, params.genome)
}

// Named workflow — reusable subworkflow
workflow VARIANT_CALLING {
    take:
    bam_ch   // Channel<tuple(meta, bam, bai)>

    main:
    CALL_VARIANTS(bam_ch)
    FILTER_VARIANTS(CALL_VARIANTS.out.vcf)

    emit:
    vcf = FILTER_VARIANTS.out.vcf
    stats = CALL_VARIANTS.out.stats
}
```

### Typed parameters (DSL2 strict, NXF ≥ 25.10)

```bash
export NXF_SYNTAX_PARSER=v2
```

```nextflow
params {
    input:  Path          // required — no default
    outdir: Path = 'results'
    max_cpus: Integer = 8
    save_bams: Boolean    // defaults to false if no value
}
```

---

## Dataflow: Channels and Values

### Channel vs Value

| Type | Behavior | Created by |
|------|----------|------------|
| `Channel` | Async stream, consumed once | `channel.of()`, `channel.fromPath()`, … |
| `Value` | Async singleton, reused forever | `channel.value()`, process output with one item |

```nextflow
// Value — shared across all tasks
ref = channel.value(file('/data/genome.fa'))

// Channel — one task per item
reads = channel.fromFilePairs('data/*_{1,2}.fastq')

ALIGN(reads, ref)   // ref reused for each reads tuple
```

### Key channel operators

```nextflow
ch = channel.of(1, 2, 3, 4, 5)

ch.map    { v -> v * 2 }          // transform each element
ch.filter { v -> v > 2 }          // keep matching elements
ch.flatMap { v -> [v, v*10] }     // one-to-many transform
ch.collect()                       // gather all into a list (Value)
ch.toList()                        // same as collect
ch.view   { v -> "item: $v" }     // print and pass through

// Joining channels on a key
ch1 = channel.of(['s1', file('a.bam')], ['s2', file('b.bam')])
ch2 = channel.of(['s1', file('a.bai')], ['s2', file('b.bai')])
ch1.join(ch2)   // → [s1, a.bam, a.bai], [s2, b.bam, b.bai]

// Group by key
ch.groupTuple()   // group [key, val] pairs by key

// Mix multiple channels
ch1.mix(ch2, ch3)   // merge, order not guaranteed

// Combine (cartesian product)
mols.combine(params_ch)   // all mol × param combinations
```

---

## Pipe `|` and And `&` Operators

```nextflow
// Pipe — chain process/operator calls
workflow {
    channel.fromPath('*.fa')
        | TRIM
        | ALIGN
        | map { bam -> tuple(bam.baseName, bam) }
        | SORT
        | view
}

// And — fan out to multiple processes with same input
workflow {
    channel.fromPath('*.sdf')
        | ( COMPUTE_DESCRIPTORS & DOCK & FILTER )
        | mix
        | view
}
```

---

## Modules

### Including definitions

```nextflow
// main.nf
include { TRIM_READS }                from './modules/trim'
include { ALIGN; ALIGN as ALIGN_V2 } from './modules/align'   // alias
include { QC; STATS }                from './modules/utils'

workflow {
    reads = channel.fromFilePairs(params.input)
    TRIM_READS(reads)
    ALIGN(TRIM_READS.out.trimmed)
}
```

### Module directory layout

```
pipeline/
├── main.nf
├── nextflow.config
└── modules/
    ├── trim/
    │   ├── main.nf        ← process TRIM_READS
    │   └── templates/
    │       └── trim.sh
    ├── align/
    │   └── main.nf
    └── utils.nf           ← multiple processes in one file
```

**Module file (modules/trim/main.nf):**
```nextflow
process TRIM_READS {
    tag "${meta.id}"
    label 'process_medium'
    container 'quay.io/biocontainers/trim-galore:0.6.7'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_trimmed.fq.gz"), emit: trimmed
    path "*.log",                              emit: log

    script:
    """
    trim_galore --cores ${task.cpus} ${reads} 2> ${meta.id}.log
    """
}
```

### Using the same process twice (aliases)

```nextflow
include { ALIGN as ALIGN_DNA } from './modules/align'
include { ALIGN as ALIGN_RNA } from './modules/align'

workflow {
    ALIGN_DNA(dna_reads, dna_ref)
    ALIGN_RNA(rna_reads, rna_ref)
}
```

---

## Workflow Composition

```nextflow
// workflows/qc.nf
workflow QC_WORKFLOW {
    take:
    reads

    main:
    FASTQC(reads)
    MULTIQC(FASTQC.out.zip.collect())

    emit:
    report = MULTIQC.out.report
}

// main.nf
include { QC_WORKFLOW }       from './workflows/qc'
include { VARIANT_CALLING }   from './workflows/variant_calling'

workflow {
    reads = channel.fromFilePairs(params.input)
    QC_WORKFLOW(reads)
    VARIANT_CALLING(reads)

    QC_WORKFLOW.out.report.view { "QC done: $it" }
}
```

---

## Process Invocation Rules

- A process/workflow can only be called **once per workflow scope**
- Use aliases to call the same process multiple times in the same workflow
- Calling order does NOT define execution order — dataflow does
- Named outputs accessed as: `PROCESS.out.name` or `PROCESS.out[0]`

---

## Workflow Outputs (NXF ≥ 24.04, stable in 25.10)

Replaces `publishDir` with a declarative output block:

```nextflow
workflow {
    main:
    ch_bam = ALIGN(reads)

    publish:
    bam_files = ch_bam
}

output {
    bam_files {
        path 'alignments'
        mode 'copy'
    }
}
```

**Dynamic path + index file:**
```nextflow
output {
    samples {
        path { meta -> "samples/${meta.id}/" }
        index {
            path   'samples.csv'
            header true
        }
    }
}
```

---

## Recursion (preview feature)

```nextflow
// Enable in config: nextflow.preview.recursion = true
process ITERATE {
    input:  val n
    output: val result
    exec:   result = n - 1
}

workflow {
    ITERATE
        .recurse(10)
        .until { v -> v <= 0 }
        .view()
}
```
