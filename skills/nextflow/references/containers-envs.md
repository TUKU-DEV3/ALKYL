# Nextflow Containers & Software Environments

## Container Runtimes Supported

| Runtime | ID | HPC-friendly | Root-free | Notes |
|---------|----|:---:|:---:|-------|
| Docker | `docker` | No | No | Development; not allowed on most HPC |
| Singularity | `singularity` | Yes | Yes | Widely deployed on HPC clusters |
| Apptainer | `apptainer` | Yes | Yes | Open-source fork of Singularity (preferred) |
| Podman | `podman` | Yes | Yes | Rootless Docker-compatible |
| Charliecloud | `charliecloud` | Yes | Yes | User-namespace based |
| Shifter | `shifter` | Yes | Yes | NERSC/Cray environments |
| Conda | — | Yes | Yes | Software envs, no isolation |
| Spack | — | Yes | Yes | HPC package manager |

---

## Docker

```groovy
// nextflow.config
docker {
    enabled    = true
    runOptions = '-u $(id -u):$(id -g)'   // run as current user
    temp       = 'auto'                   // mount /tmp
}
```

Per-process container:
```nextflow
process RUN_RDKIT {
    container 'docker.io/continuumio/miniconda3'
    // or just the image name — docker:// prefix optional
    container 'nfcore/rnaseq:3.14'
    ...
}
```

Global container:
```groovy
process.container = 'docker://my-chem-image:1.0'
docker.enabled    = true
```

---

## Singularity / Apptainer (HPC)

**Apptainer is the recommended runtime for HPC** — runs without root, supports autofs.

```groovy
// nextflow.config
apptainer {
    enabled    = true
    autoMounts = true     // auto-mount host paths (requires user bind control)
    cacheDir   = '/shared/containers/nxf-cache'   // shared across nodes
}
```

Pull from Docker Hub automatically:
```groovy
process.container = 'docker://quay.io/biocontainers/blast:2.14.0'
apptainer.enabled = true
// Nextflow will pull and convert to .sif automatically
```

**Force local file only:**
```groovy
process.container = 'file:///path/to/my_image.sif'
```

**Pull from different registries:**
```groovy
process.container = 'shub://my-singularity-hub/image'   // SingularityHub
process.container = 'oras://registry.example.com/img'   // OCI registry
```

**Library vs cache dir:**
- `libraryDir` — read-only repository of pre-downloaded images
- `cacheDir` — writable, where new pulls are stored

```groovy
apptainer {
    libraryDir = '/shared/containers/library'  // checked first (read-only)
    cacheDir   = '/scratch/${USER}/sif-cache'  // fallback (writable)
}
```

---

## Conda

```nextflow
process RDKIT_CALC {
    conda 'rdkit=2023.09.4 pandas=2.0'   // package spec
    // or
    conda 'environment.yml'               // environment file
    // or
    conda 'bioconda::blast=2.14.0'        // channel::package

    script:
    """
    python calc_properties.py ${input}
    """
}
```

```groovy
// nextflow.config — global conda settings
conda {
    enabled     = true
    useMamba    = true     // faster solver
    cacheDir    = '/shared/conda-envs'
    createTimeout = '30 min'
}
```

---

## Per-Process Containers (Multiple Containers)

```groovy
process {
    withName: 'BLAST_SEARCH' {
        container = 'quay.io/biocontainers/blast:2.14.0'
    }
    withName: 'RDKIT_FILTER' {
        container = 'quay.io/biocontainers/rdkit:2023.09.4'
    }
    withName: 'ORCA_.*' {
        container = 'docker://orca-mpi:5.0.4'
        clusterOptions = '--gres=gpu:1'
    }
}
apptainer.enabled = true
```

---

## Wave (Seqera — on-demand containers)

Wave builds and provisions containers on the fly from a `Dockerfile` or Conda spec — no manual image builds.

```groovy
// nextflow.config
wave {
    enabled = true
}
tower {
    accessToken = secrets.TOWER_TOKEN
}
```

```nextflow
process CUSTOM_TOOL {
    conda 'rdkit=2023 openbabel=3.1'   // Wave builds the container automatically

    script:
    """
    python process.py ${input}
    """
}
```

---

## BioContainers — Chemistry/Bioinformatics Images

Recommended image sources for chemistry pipelines:

| Tool | Image |
|------|-------|
| RDKit | `quay.io/biocontainers/rdkit:2023.09.4` |
| OpenBabel | `quay.io/biocontainers/openbabel:3.1.1` |
| BLAST | `quay.io/biocontainers/blast:2.14.0` |
| AutoDock Vina | `quay.io/biocontainers/autodock-vina:1.2.3` |
| GROMACS | `quay.io/biocontainers/gromacs:2023.2` |
| AmberTools | `quay.io/biocontainers/ambertools:22.0` |
| Multiqc | `quay.io/biocontainers/multiqc:1.21` |
| FastQC | `quay.io/biocontainers/fastqc:0.12.1` |
| DeepChem | `docker.io/deepchemio/deepchem:2.7.1` |

Search: `https://quay.io/repository/biocontainers/<tool>`

---

## Container Best Practices

```groovy
// ✓ Pin exact tags — never use :latest
container = 'quay.io/biocontainers/rdkit:2023.09.4--py311hc857e9e_1'

// ✓ Use a shared cache for Apptainer/Singularity on HPC
apptainer.cacheDir = '/shared/nxf-apptainer-cache'

// ✓ Check image has /bin/bash and ps (required by Nextflow)

// ✓ Automate container selection per process via config, not hardcoding
process.withName: 'MY_PROC' { container = 'image:tag' }

// ✗ Avoid mounting /home or / into containers
// ✗ Avoid running as root inside containers on shared HPC
```

---

## Minimal Dockerfile for Chemistry Process

```dockerfile
FROM condaforge/miniforge3:latest

RUN conda install -y -c conda-forge -c rdkit \
    rdkit=2023.09.4 \
    openbabel=3.1.1 \
    pandas=2.0 \
    && conda clean -afy

# Nextflow requirements
RUN apt-get update && apt-get install -y procps && rm -rf /var/lib/apt/lists/*

ENTRYPOINT ["/bin/bash"]
```
