# Nextflow for Chemistry & Bioinformatics

## Typical Computational Chemistry Pipeline Structure

```
main.nf
nextflow.config
modules/
├── rdkit/main.nf         ← molecular preprocessing
├── dock/main.nf          ← docking (AutoDock, Vina, Glide)
├── orca/main.nf          ← quantum chemistry (ORCA, Gaussian)
├── md/main.nf            ← molecular dynamics (GROMACS, AMBER)
└── analysis/main.nf      ← postprocessing / ML scoring
data/
├── library.sdf           ← compound library
└── receptor.pdb          ← protein target
conf/
├── base.config
├── hpc.config
└── containers.config
```

---

## Pattern 1 — Virtual Screening Pipeline

```nextflow
// main.nf
params.library  = 'data/library.sdf'
params.receptor = 'data/receptor.pdb'
params.outdir   = 'results'
params.batch_size = 500

workflow {
    library  = channel.fromPath(params.library)
    receptor = channel.value(file(params.receptor))

    // Split library into batches
    batches = library | splitText(by: params.batch_size * 5, file: true)

    // Prepare molecules
    PREPARE_LIGANDS(batches)

    // Dock all batches against receptor
    VINA_DOCK(PREPARE_LIGANDS.out.prepared, receptor)

    // Collect and rank
    RANK_RESULTS(VINA_DOCK.out.poses.collect())
}
```

```nextflow
// modules/dock/main.nf
process VINA_DOCK {
    tag "${ligand_batch.baseName}"
    label 'process_medium'
    container 'quay.io/biocontainers/autodock-vina:1.2.3'
    publishDir "${params.outdir}/docking", mode: 'copy'

    input:
    path ligand_batch
    path receptor

    output:
    path "poses_${ligand_batch.baseName}.pdbqt", emit: poses
    path "log_${ligand_batch.baseName}.txt",     emit: log

    script:
    """
    vina \\
        --receptor ${receptor} \\
        --ligand   ${ligand_batch} \\
        --out      poses_${ligand_batch.baseName}.pdbqt \\
        --log      log_${ligand_batch.baseName}.txt \\
        --cpu      ${task.cpus} \\
        --exhaustiveness 8
    """
}
```

---

## Pattern 2 — ORCA Quantum Chemistry (HPC)

```nextflow
// modules/orca/main.nf
process ORCA_OPT {
    tag "${meta.id}"
    label 'process_high'       // 32 CPUs, 256 GB, 72h in config
    container 'docker://orca-mpi:5.0.4'

    input:
    tuple val(meta), path(xyz)

    output:
    tuple val(meta), path("${meta.id}.out"),   emit: output
    tuple val(meta), path("${meta.id}_opt.xyz"), emit: optimized
    tuple val(meta), path("${meta.id}.engrad"), emit: gradient, optional: true

    script:
    def inp = """
    ! ${meta.functional} ${meta.basis} Opt
    %pal nprocs ${task.cpus} end
    %maxcore ${(task.memory.toMega() / task.cpus).intValue()}
    * xyzfile ${meta.charge} ${meta.multiplicity} ${xyz}
    """.stripIndent()

    """
    cat > ${meta.id}.inp <<EOF
${inp}
EOF
    orca ${meta.id}.inp > ${meta.id}.out
    """
}
```

```nextflow
// Usage in workflow
workflow {
    molecules = channel.fromPath('structures/*.xyz')
        | map { xyz ->
            def meta = [
                id:           xyz.baseName,
                functional:   'PBE0',
                basis:        'def2-TZVP',
                charge:       0,
                multiplicity: 1
            ]
            tuple(meta, xyz)
          }

    ORCA_OPT(molecules)
    EXTRACT_ENERGIES(ORCA_OPT.out.output)
}
```

---

## Pattern 3 — GROMACS MD Pipeline

```nextflow
process GROMACS_MINIMIZE {
    tag "${meta.id}"
    label 'process_high'
    container 'quay.io/biocontainers/gromacs:2023.2'

    input:
    tuple val(meta), path(gro), path(top), path(mdp)

    output:
    tuple val(meta), path("em.gro"),  emit: structure
    tuple val(meta), path("em.edr"),  emit: energy
    tuple val(meta), path("em.log"),  emit: log

    script:
    """
    gmx grompp -f ${mdp} -c ${gro} -p ${top} -o em.tpr -maxwarn 5
    gmx mdrun -v -deffnm em -ntmpi 1 -ntomp ${task.cpus}
    """
}

process GROMACS_NVT {
    tag "${meta.id}"
    label 'gpu'

    input:
    tuple val(meta), path(gro), path(top), path(mdp)

    output:
    tuple val(meta), path("nvt.gro"), emit: structure
    tuple val(meta), path("nvt.xtc"), emit: trajectory
    tuple val(meta), path("nvt.edr"), emit: energy

    script:
    """
    gmx grompp -f ${mdp} -c ${gro} -p ${top} -o nvt.tpr
    gmx mdrun -v -deffnm nvt -ntmpi 1 -ntomp ${task.cpus} -gpu_id 0
    """
}
```

---

## Pattern 4 — RDKit Molecular Processing

```nextflow
process CALC_DESCRIPTORS {
    tag "${meta.id}"
    container 'quay.io/biocontainers/rdkit:2023.09.4'
    publishDir "${params.outdir}/descriptors", mode: 'copy'

    input:
    tuple val(meta), val(smiles)

    output:
    tuple val(meta), path("${meta.id}_desc.csv"), emit: descriptors

    script:
    """
    python3 - <<'PYEOF'
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors, QED, Crippen
    import csv, sys

    mol = Chem.MolFromSmiles('${smiles}')
    if mol is None:
        print(f"Invalid SMILES: ${smiles}", file=sys.stderr)
        sys.exit(1)

    with open('${meta.id}_desc.csv', 'w') as f:
        w = csv.writer(f)
        w.writerow(['id','smiles','mw','logp','hbd','hba','tpsa','qed','rings'])
        w.writerow([
            '${meta.id}',
            '${smiles}',
            round(Descriptors.MolWt(mol), 3),
            round(Crippen.MolLogP(mol), 3),
            rdMolDescriptors.CalcNumHBD(mol),
            rdMolDescriptors.CalcNumHBA(mol),
            round(Descriptors.TPSA(mol), 2),
            round(QED.qed(mol), 4),
            rdMolDescriptors.CalcNumRings(mol),
        ])
    PYEOF
    """
}
```

---

## Pattern 5 — BLAST for Enzyme Discovery

```nextflow
process BLAST_SEARCH {
    tag "${query.baseName}"
    label 'process_medium'
    container 'quay.io/biocontainers/blast:2.14.0'

    input:
    path query      // protein FASTA
    path db_dir     // BLAST database directory (Value)

    output:
    path "${query.baseName}_hits.txt", emit: hits
    path "${query.baseName}.blast6",   emit: table

    script:
    """
    blastp \\
        -query ${query} \\
        -db ${db_dir}/nr \\
        -out ${query.baseName}.blast6 \\
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \\
        -num_threads ${task.cpus} \\
        -evalue 1e-10 \\
        -max_target_seqs 500

    awk '\$3 >= 30' ${query.baseName}.blast6 | cut -f2 | sort -u > ${query.baseName}_hits.txt
    """
}
```

---

## Pattern 6 — Parallel DFT Grid Search

Run all molecules × all method combinations:

```nextflow
workflow {
    molecules  = channel.fromPath('structures/*.xyz')
    functionals = channel.of('PBE', 'B3LYP', 'PBE0', 'wB97X-D3')
    bases       = channel.of('def2-SVP', 'def2-TZVP')

    // Cartesian product: each molecule × functional × basis
    molecules
        .combine(functionals)
        .combine(bases)
        .map { xyz, func, basis ->
            def meta = [
                id:         "${xyz.baseName}_${func}_${basis}",
                functional: func,
                basis:      basis,
                charge:     0,
                multiplicity: 1
            ]
            tuple(meta, xyz)
          }
        | ORCA_SP
        | COLLECT_ENERGIES
}
```

---

## nf-core Chemistry Pipelines

| Pipeline | Purpose | URL |
|----------|---------|-----|
| `nf-core/proteinfold` | AlphaFold2, ESMFold, RoseTTAFold | nf-co.re/proteinfold |
| `nf-core/drugresponseeval` | Drug response prediction | nf-co.re/drugresponseeval |
| `nf-core/ampliseq` | Amplicon sequencing | nf-co.re/ampliseq |
| `nf-core/rnaseq` | RNA-seq reference | nf-co.re/rnaseq |
| `nf-core/sarek` | Variant calling | nf-co.re/sarek |

**Use nf-core modules as building blocks:**
```bash
nf-core modules install blast/blastp
nf-core modules install rdkit/mol_to_fingerprint
```

---

## Config for Chemistry HPC (SLURM example)

```groovy
// conf/hpc.config — chemistry cluster
process {
    withLabel: process_low    { cpus = 4;  memory = '16 GB'; time = '4h'  }
    withLabel: process_medium { cpus = 16; memory = '64 GB'; time = '24h' }
    withLabel: process_high   { cpus = 32; memory = '256 GB'; time = '72h'; queue = 'highmem' }
    withLabel: gpu {
        cpus          = 8
        memory        = '32 GB'
        time          = '24h'
        queue         = 'gpu'
        clusterOptions = '--gres=gpu:v100:1'
    }
    withName: 'ORCA_.*' {
        // MPI jobs — use Intel MPI or OpenMPI
        clusterOptions = '--ntasks-per-node=32 --nodes=1'
        module = 'openmpi/4.1 orca/5.0.4'   // environment modules
    }
}

singularity {
    enabled    = true
    autoMounts = true
    cacheDir   = '/shared/containers'
}

params {
    max_cpus = 128
    max_mem  = '1.TB'
    max_time = '168.h'
}
```

---

## Useful One-liners

```bash
# Run with SLURM and Singularity, resuming
nextflow run main.nf -profile slurm,singularity -resume

# Dry run — check pipeline logic without running
nextflow run main.nf -stub

# Override params on CLI
nextflow run main.nf --library data/lib2.sdf --outdir results_lib2

# Show execution report
nextflow run main.nf -with-report reports/report.html -with-timeline reports/timeline.html

# List cached runs
nextflow log

# Clean a specific run from cache
nextflow clean -f run_name

# Watch running jobs
watch nextflow log -f
```
