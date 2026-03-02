# Nextflow Execution & Configuration

## Configuration File Hierarchy

Files loaded in order (later overrides earlier):

1. `$HOME/.nextflow/config` — global user defaults
2. `nextflow.config` in project dir
3. `nextflow.config` in launch dir
4. `-c extra.config` on CLI

**Force a specific config only:** `nextflow run main.nf -C prod.config`

---

## Configuration Syntax

```groovy
// nextflow.config

// ── Params ──────────────────────────────────────────────────────────────────
params {
    input    = 'data/*.sdf'
    outdir   = 'results'
    max_cpus = 16
    max_mem  = '128 GB'
    max_time = '24h'
}

// ── Process defaults ────────────────────────────────────────────────────────
process {
    executor     = 'slurm'
    queue        = 'normal'
    cpus         = 1
    memory       = '4 GB'
    time         = '1h'
    errorStrategy = 'retry'
    maxRetries   = 2

    // Label-based resource tiers
    withLabel: process_low {
        cpus   = 2
        memory = '8 GB'
        time   = '2h'
    }
    withLabel: process_medium {
        cpus   = 8
        memory = '32 GB'
        time   = '8h'
    }
    withLabel: process_high {
        cpus   = 16
        memory = '64 GB'
        time   = '24h'
    }
    withLabel: gpu {
        cpus          = 8
        memory        = '32 GB'
        clusterOptions = '--gres=gpu:1'
        queue         = 'gpu'
    }

    // Process-specific override
    withName: 'ORCA_CALC' {
        cpus   = 32
        memory = '256 GB'
        time   = '72h'
        queue  = 'highmem'
    }
    // Regex selector — all processes starting with ALIGN
    withName: 'ALIGN.*' {
        memory = '16 GB'
    }
}

// ── Executor ─────────────────────────────────────────────────────────────────
executor {
    name      = 'slurm'
    queueSize = 200       // max concurrent jobs
    pollInterval = '30s'
    retry.maxAttempt = 3
}

// ── Work directory ───────────────────────────────────────────────────────────
workDir = '/scratch/${USER}/nextflow-work'

// ── Reports ──────────────────────────────────────────────────────────────────
report {
    enabled = true
    file    = 'reports/report.html'
}
timeline {
    enabled = true
    file    = 'reports/timeline.html'
}
trace {
    enabled = true
    file    = 'reports/trace.txt'
    fields  = 'task_id,name,status,cpus,memory,time,realtime,%cpu,%mem'
}
dag {
    enabled = true
    file    = 'reports/dag.html'
}
```

---

## Profiles

Define alternate execution environments selected with `-profile`:

```groovy
profiles {
    local {
        process.executor = 'local'
        process.cpus     = 4
    }

    slurm {
        process.executor = 'slurm'
        process.queue    = 'normal'
        docker.enabled   = false
        singularity.enabled = true
        singularity.autoMounts = true
    }

    docker {
        docker.enabled = true
        docker.runOptions = '-u $(id -u):$(id -g)'
    }

    test {
        params.input = "${projectDir}/test/data/*.sdf"
        params.outdir = 'test-results'
        process.cpus   = 1
        process.memory = '2 GB'
    }

    cloud_aws {
        process.executor  = 'awsbatch'
        process.container = 'public.ecr.aws/my-org/my-image:latest'
        aws.region        = 'eu-west-1'
        workDir           = 's3://my-bucket/nextflow-work'
    }
}
```

Activate multiple: `nextflow run main.nf -profile slurm,docker`

---

## Executors Reference

| Executor | ID | Notes |
|----------|----|-------|
| Local | `local` | Default — your machine |
| SLURM | `slurm` | HPC — most common |
| LSF | `lsf` | IBM Platform LSF |
| SGE/UGE | `sge` | Sun/Univa Grid Engine |
| PBS/Torque | `pbs` | PBS Works / Torque |
| HTCondor | `condor` | Shared filesystem required |
| HyperQueue | `hq` | Lightweight HPC (≥ 0.20) |
| Flux | `flux` | Flux Framework |
| AWS Batch | `awsbatch` | Requires Docker, S3 workDir |
| Azure Batch | `azurebatch` | Requires Docker, Azure Blob workDir |
| Google Batch | `google-batch` | Requires Docker, GCS workDir |
| Kubernetes | `k8s` | PVC-backed shared storage required |

**SLURM example:**
```groovy
process {
    executor     = 'slurm'
    queue        = 'compute'
    clusterOptions = '--account=chem-lab --mail-user=user@inst.edu'
}
executor {
    queueSize    = 100
    submitRateLimit = '10/1min'
}
```

---

## Cache and Resume

### Enable resume
```bash
nextflow run main.nf -resume
```

### What gets hashed (task cache key)
- Session ID
- Process name + calling workflow name
- Container image
- Conda/Spack environment
- All input files (path + mtime + size by default)
- Process script content
- Global variables referenced in the script

### Cache modes

| Mode | Hash includes | Use when |
|------|---------------|----------|
| `standard` (default) | path, mtime, size | Local filesystem |
| `lenient` | path, size (no mtime) | NFS / shared FS with unreliable timestamps |
| `deep` | file content (MD5) | Maximum correctness, slow |

```groovy
process {
    cache = 'lenient'   // global
}
// or per-process
process MY_PROCESS {
    cache 'deep'
}
```

### Troubleshooting resume failures
```bash
# View what changed between runs
nextflow log            # list all runs
nextflow log run-name   # tasks in a run
nextflow log run-name -f 'name,hash,status'

# Force re-run a specific process
nextflow run main.nf -resume --force_rerun MY_PROCESS
```

**Common causes of cache miss:**
- Input file modified (mtime changed, even if content unchanged → use `lenient`)
- Process script changed (even whitespace)
- Global variable race condition (use `def x` in closures)
- Different session ID (new run without `-resume`)
- Container image tag updated (`:latest` is a bad practice)

---

## Dynamic Resources (retry scaling)

```groovy
process {
    errorStrategy = { task.exitStatus in [104, 134, 137, 139, 143, 247] ? 'retry' : 'finish' }
    maxRetries    = 3
    memory        = { check_max( 8.GB  * task.attempt, 'memory') }
    time          = { check_max( 4.h   * task.attempt, 'time') }
    cpus          = { check_max( 2     * task.attempt, 'cpus') }
}

// Helper function to cap resources
def check_max(obj, type) {
    if (type == 'memory') {
        return obj.compareTo(params.max_mem as nextflow.util.MemoryUnit) == 1
            ? params.max_mem as nextflow.util.MemoryUnit : obj
    }
    if (type == 'time') {
        return obj.compareTo(params.max_time as nextflow.util.Duration) == 1
            ? params.max_time as nextflow.util.Duration : obj
    }
    if (type == 'cpus') {
        return Math.min(obj as int, params.max_cpus as int)
    }
}
```

---

## Environment Variables (key NXF vars)

| Variable | Effect |
|----------|--------|
| `NXF_WORK` | Override default work directory |
| `NXF_TEMP` | Temp directory for tasks |
| `NXF_EXECUTOR` | Default executor |
| `NXF_SYNTAX_PARSER=v2` | Enable strict DSL2 parser |
| `NXF_CLOUDCACHE_PATH` | Use S3/GCS as cache store |
| `NXF_APPTAINER_CACHEDIR` | Where Apptainer images are cached |
| `NXF_SINGULARITY_CACHEDIR` | Where Singularity images are cached |
| `NXF_ANSI_LOG=false` | Disable ANSI colors |
| `NXF_DEBUG=1` | Enable debug log |

---

## Config Includes and Secrets

```groovy
// nextflow.config
includeConfig 'conf/base.config'
includeConfig 'conf/containers.config'

// Conditional include
if (params.use_hpc) {
    includeConfig 'conf/hpc.config'
}

// Secrets (nextflow secrets set MY_KEY value)
process.environment = [
    API_KEY: secrets.MY_KEY
]
```

---

## Workflow-level Config

```groovy
workflow {
    output {
        mode      = 'copy'       // 'copy', 'move', 'link', 'symlink'
        overwrite = true
    }
}

// Output dir
outputDir = 'my-results'   // default: 'results'
```
