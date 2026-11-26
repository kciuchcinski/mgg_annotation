# MGG Annotation Pipeline

A Nextflow pipeline for annotating phage genomes using PHROGs-enriched profile HMM searches across multiple databases (PHROGs, PFAM, ECOD) for sensitive and precise functional annotation.

## Quick Start (macOS, Apple Silicon or Intel)

You need [Homebrew](https://brew.sh), [Docker Desktop](https://www.docker.com/products/docker-desktop/) and ~30 GB of free disk space.

```bash
# 1. Install Nextflow (Homebrew installs a suitable Java automatically)
brew install nextflow

# 2. Start Docker Desktop, then build the pipeline image (native on Apple Silicon, ~2 min)
git clone https://github.com/kciuchcinski/mgg_annotation.git
cd mgg_annotation
docker build -t phage_annotation:latest .

# 3. Download the databases into ./databases (~6.5 GB download, ~25 GB unpacked)
./setup_databases.sh databases

# 4. Run on the bundled test genomes (~10 min on a MacBook), results in test/output/
nextflow run main.nf -params-file test/config.yml
```

To annotate your own genomes, edit `PHAGES_DIR` and `OUTPUT_DIR` in `config.yml` and run `nextflow run main.nf -params-file config.yml`. No `-profile` is needed: Docker is the default.

Linux workstations use the same steps (install Nextflow from [nextflow.io](https://www.nextflow.io/docs/latest/install.html) and Docker Engine). For HPC clusters with Apptainer see [Option B](#3-set-up-the-execution-environment).

See [Installation](#installation), [Database setup](#database-setup), and [Configuration](#configuration) for full details.

## Overview

The pipeline takes phage genome assemblies (FASTA) as input and produces per-ORF functional annotations by:

1. **Preprocessing** — validates input genomes and filters by minimum length. Files are skipped (and listed in `dropped_inputs.tsv`) if they hold more than one sequence or are shorter than `PHAGE_MIN_LENGTH`. Contig IDs (the first word of the FASTA header) must be unique across all input files; if two files share one, the run stops and lists the clashing files.
2. **ORF calling** — predicts open reading frames using both Prodigal-gv (meta mode) and Glimmer3, then reconciles predictions.
3. **Protein clustering** — clusters translated ORFs with MMseqs2 to reduce redundancy and generate multiple sequence alignments (MSAs).
4. **MSA enrichment** — enriches MSAs against the PHROGs database using HHblits to improve search sensitivity.
5. **Database searches** — searches enriched profiles against three HHsuite-formatted databases: PHROGs, PFAM, and ECOD.
6. **Annotation & export** — collects and filters hits, builds a per-ORF annotation table, and exports GenBank files.

## Requirements

- **OS**: macOS (Apple Silicon or Intel) or Linux
- **Nextflow** ≥ 25.04 (requires Java 17+; `brew install nextflow` takes care of both)
- **Container runtime** (one of):
  - [Docker](https://www.docker.com/) — default; recommended for macOS and Linux workstations
  - [Apptainer](https://apptainer.org/) ≥ 1.1 — recommended for HPC
  - [Conda](https://docs.conda.io/) / [Micromamba](https://mamba.readthedocs.io/) — Linux x86_64 only (Glimmer has no macOS / ARM conda build)
- **Disk**: ~30 GB for the databases, plus space for `work/` during a run

## Installation

### 1. Install Nextflow

On macOS:

```bash
brew install nextflow
```

On Linux (needs Java 17 or newer):

```bash
curl -s https://get.nextflow.io | bash
chmod +x nextflow
mkdir -p ~/.local/bin && mv nextflow ~/.local/bin/   # make sure ~/.local/bin is on your $PATH
```

See the [official Nextflow docs](https://www.nextflow.io/docs/latest/install.html) for other installation methods.

### 2. Clone the repository

```bash
git clone https://github.com/kciuchcinski/mgg_annotation.git
cd mgg_annotation
```

### 3. Set up the execution environment

Choose **one** of the three options below.

<details open>
<summary><strong>Option A — Docker (default; macOS and Linux workstations)</strong></summary>

Start Docker Desktop (macOS) or the Docker daemon (Linux), then build the image:

```bash
docker build -t phage_annotation:latest .
```

The image builds natively for the machine it is built on (`arm64` on Apple Silicon, `amd64` on Intel/AMD), so no emulation is involved.

**Memory:** the pipeline uses at most `MAX_MEMORY` (default `7 GB`) and `MAX_CPUS` (default `8`) at once. Docker Desktop must have at least that much memory (Settings → Resources → Memory; the default on most Macs is 8 GB). On a larger machine raise both in `config.yml` together with Docker Desktop's limits, e.g. `MAX_CPUS: 12` and `MAX_MEMORY: '24 GB'`.

</details>

<details>
<summary><strong>Option B — Apptainer (HPC)</strong></summary>

Build the Apptainer image from the provided definition file:

```bash
apptainer build phage_annotation.sif phage_annotation.def
```

The resulting `.sif` file must be in the pipeline root directory (or update the `container` path in `nextflow.config`). Host paths (including `DB_ROOT`) are mounted automatically via `autoMounts`. Run with `-profile apptainer`.

</details>

<details>
<summary><strong>Option C — Conda / Micromamba (Linux x86_64 only)</strong></summary>

Create the environment from the provided YAML file:

```bash
# Using Micromamba (faster, recommended)
micromamba create -f environment.yml
micromamba activate phage_annotation

# OR using Conda
conda env create -f environment.yml
conda activate phage_annotation
```

After activation, update the `conda` path in `nextflow.config` to point to your environment:

```groovy
conda {
    process {
        withName: '.*' {
            conda = '/path/to/your/envs/phage_annotation'
        }
    }
}
```

Replace `/path/to/your/envs/phage_annotation` with the actual path (e.g. the output of `micromamba env list` or `conda env list`). If using Micromamba, ensure `useMicromamba = true` is set in `nextflow.config`. Run with `-profile conda`.

</details>

## Database setup

The pipeline requires three HHsuite-formatted profile databases, a Glimmer training file, and one metadata table. All paths are resolved relative to `DB_ROOT` (set in `config.yml`; a relative `DB_ROOT` is taken from the directory you launch Nextflow in), but individual paths can be overridden. The pipeline checks that all databases exist before it starts.

> **Disk space:** the downloads are PHROGs ~0.16 GB · Pfam-A ~2.25 GB · ECOD ~4.2 GB (~6.5 GB), but unpacked they take **~25 GB** (Pfam ~5.5 GB, ECOD ~12.6 GB of `ffdata` alone). The setup script checks for ~30 GB of free space before it starts.

### Quick setup (recommended)

Run the provided setup script to download and prepare all databases:

```bash
bash setup_databases.sh /path/to/databases
```

The script downloads each database, extracts it, and renames directories to the version-agnostic paths expected by the pipeline defaults. It works with the tools that ship with macOS (`curl`) and Linux (`curl` or `wget`). Existing databases are skipped and interrupted downloads resume, so the script is safe to re-run.

After completion, set `DB_ROOT` in `config.yml` to the path you provided.

### Manual setup

If you prefer to set up databases manually, download each one into `DB_ROOT`:

| Database | Download | Target directory |
|---|---|---|
| PHROGs v4 | [phrogs_hhsuite_db.tar.gz](https://phrogs.lmge.uca.fr/downloads_from_website/phrogs_hhsuite_db.tar.gz) | `phrogs/` |
| PHROGs metadata | [phrogs_annot_v4.tsv](https://phrogs.lmge.uca.fr/downloads_from_website/phrog_annot_v4.tsv) | `tables/phrog_annot_v4.tsv` |
| Pfam-A | [pfamA_35.0.tar.gz](https://wwwuser.gwdguser.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/pfamA_35.0.tar.gz) | `pfam/` |
| ECOD | [ecod.v294.F40.hhm_db.tar.gz](http://prodata.swmed.edu/ecod/distributions/ecod.v294.F40.hhm_db.tar.gz) | `ecod/` |
| Glimmer ICM | bundled in `resources/` (copied by setup script) | `training-file_refseq.icm` |

After extraction, rename the versioned directories to their version-agnostic names (e.g. `pfamA_35.0/` to `pfam/`). The HHsuite database prefix inside each directory must also match the directory name (e.g. `pfam/pfam_a3m.ffdata`).

Alternatively, keep the original directory names and override the paths in `config.yml`:

```yaml
HHSUITE:
  PFAM: "/path/to/databases/pfamA_35.0/pfam"
  ECOD: "/path/to/databases/ECOD_F70_v294/ecod"
```

### Expected directory layout

After running `setup_databases.sh` (or manual setup with default names):

```
DB_ROOT/
├── training-file_refseq.icm
├── phrogs/
│   └── phrogs{_a3m.ffdata, _a3m.ffindex, _hhm.ffdata, _hhm.ffindex, ...}
├── pfam/
│   └── pfam{_a3m.ffdata, _a3m.ffindex, _hhm.ffdata, _hhm.ffindex, ...}
├── ecod/
│   └── ecod{_a3m.ffdata, _a3m.ffindex, ...}
└── tables/
    └── phrog_annot_v4.tsv
```

## Configuration

All user-configurable parameters are set in `config.yml`. Edit it directly before running.

### Minimal configuration

Every parameter has a default (see `nextflow.config`), so you only need to set three paths:

```yaml
# Input directory containing phage .fasta files
PHAGES_DIR: /path/to/your/phage/fastas

# Output directory
OUTPUT_DIR: /path/to/results

# Root directory containing all databases (see "Database setup")
DB_ROOT: /path/to/databases
```

### Full parameter reference

#### Basic parameters

| Parameter | Default | Description |
|---|---|---|
| `PHAGES_DIR` | `test/data` | Directory containing input FASTA files (one genome per file) |
| `INPUT_EXTENSION` | `fasta` | File extension to look for in `PHAGES_DIR` |
| `OUTPUT_DIR` | `results` | Directory for pipeline outputs |
| `DB_ROOT` | `databases` (in the pipeline directory) | Root directory for all databases |
| `PHAGE_MIN_LENGTH` | `2000` | Minimum genome length in bp to process |
| `WRITE_SEARCH_TABLE` | `true` | Write raw HHsuite hits to `search.tsv` |
| `WRITE_ANNOTATION_TABLE` | `true` | Write per-ORF annotation to `annotation.tsv` |
| `N_FUNCTIONS_PER_DB` | `2` | Max number of functions reported per protein per database |
| `PUBLISH_INTERMEDIATES` | `false` | Also copy intermediate files (ORF tables, proteins, clusters, MSAs, HHR files) to `OUTPUT_DIR/intermediates/` |

#### Annotation filtering thresholds

These thresholds control which HHsuite hits are retained in the final annotation. Separate thresholds are applied to PHROGs and to the domain databases (PFAM, ECOD).

| Parameter | Default | Description |
|---|---|---|
| `ANNOTATION_FILTERS.PHROGS.MIN_PROB` | `0.7` | Minimum HHsearch probability for PHROGs hits |
| `ANNOTATION_FILTERS.PHROGS.MIN_QCOV` | `0.1` | Minimum query coverage for PHROGs hits |
| `ANNOTATION_FILTERS.PHROGS.MIN_TCOV` | `0.1` | Minimum template coverage for PHROGs hits |
| `ANNOTATION_FILTERS.DOMAIN_DBS.MIN_PROB` | `0.7` | Minimum HHsearch probability for domain DB hits |
| `ANNOTATION_FILTERS.DOMAIN_DBS.MIN_QCOV` | `0.1` | Minimum query coverage for domain DB hits |
| `ANNOTATION_FILTERS.DOMAIN_DBS.MIN_TCOV` | `0.1` | Minimum template coverage for domain DB hits |

#### Clustering parameters

Controls MMseqs2 protein clustering behavior.

| Parameter | Default | Description |
|---|---|---|
| `CLUSTERING.IDENTITY` | `0.5` | Minimum sequence identity |
| `CLUSTERING.COVERAGE` | `0.8` | Minimum alignment coverage |
| `CLUSTERING.EVAL` | `0.001` | Maximum e-value |
| `CLUSTERING.SENSITIVITY` | `7.5` | MMseqs2 sensitivity (1–7.5) |

#### Execution parameters

| Parameter | Default | Description |
|---|---|---|
| `BATCH_SIZE` | `20` | Number of protein clusters per HHsuite search batch |
| `THREADS_PER_BATCH` | `4` | CPU threads allocated per HHsuite search task |
| `SEARCH_TOOL` | `hhblits_omp` | HHsuite search binary (`hhblits_omp` or `hhblits`) |
| `ENRICH_CPUS` | `4` | CPU threads for the PHROGs enrichment step |
| `PARSING_CPUS` | `4` | CPU threads for clustering, MSA preparation, hit collection and annotation building |
| `MAX_CPUS` | `8` | Docker profile: total CPUs the pipeline may use at once (also caps every task) |
| `MAX_MEMORY` | `7 GB` | Docker profile: total memory the pipeline may use at once (also caps every task) |

> **Note:** The per-process parameters control CPU allocation per task. `MAX_CPUS` / `MAX_MEMORY` set the global pool for the default Docker profile; for other profiles see the section below.

#### Resource limits

Nextflow uses the **local executor**, which runs all tasks on the current machine. With the default Docker profile the pool is set by `MAX_CPUS` and `MAX_MEMORY`. For the `apptainer` profile it is set to 64 CPUs / 256 GB in `nextflow.config`; change it by editing the `executor` block there:

```groovy
executor {
    cpus = 16        // max total CPUs across all concurrent tasks
    memory = '64 GB' // max total memory across all concurrent tasks
}
```

Nextflow will schedule tasks up to these limits based on each process's CPU and memory requests. For example, with `executor.cpus = 16` and `THREADS_PER_BATCH = 4`, up to four HHsuite search batches can run in parallel.

On an **HPC cluster** (SLURM, PBS, SGE, etc.), add the appropriate executor settings to `nextflow.config`:

```groovy
process {
    executor = 'slurm'
    queue = 'your-partition'
}
```

See the [Nextflow executors documentation](https://www.nextflow.io/docs/latest/executor.html) for the full list of supported schedulers and their configuration options.

## Running the pipeline

```bash
# Docker (default; macOS / Linux workstations)
nextflow run main.nf -params-file config.yml

# Apptainer (HPC)
nextflow run main.nf -profile apptainer -params-file config.yml

# Conda (Linux x86_64)
nextflow run main.nf -profile conda -params-file config.yml
```

### Resuming a failed run

If a run fails or is interrupted, relaunch with `-resume` to restart from the last successful step:

```bash
nextflow run main.nf -params-file config.yml -resume
```

Nextflow caches intermediate results in the `work/` directory. Do not delete it if you intend to resume.

## Resource requirements

Resource usage is dominated by two steps: HHsuite searches (CPU-bound, many parallel batches) and annotation building (RAM-bound, single process). The figures below were measured on an HPC cluster using the default `BATCH_SIZE=20` and `THREADS_PER_BATCH=4`.

### HHsuite (search and enrichment)

Each HHsuite batch runs independently and uses a fixed amount of memory regardless of dataset size:

| Step | RAM per instance |
|---|---|
| `ENRICH_MSA_PHROGS` | ~2 GB |
| `HHSUITE_SEARCH_BATCH` | ~4 GB |

Total HHsuite compute scales with the number of protein clusters produced by MMseqs2:

| Input genomes | Total HHsuite CPU-hours |
|---|---|
| 10 | ~0.8 h |
| 100 | ~12 h |
| 1 000 | ~67 h |

Wall time depends on how many batches run in parallel. With `THREADS_PER_BATCH=4`, each batch slot needs 4 CPUs and ~4 GB RAM. Allocate CPUs and memory in multiples of that ratio.

### Annotation building

`BUILD_ANNOTATION` loads all hits and ORF sequences into a single in-memory table. Peak RAM scales with protein cluster diversity (not raw genome count) — highly redundant datasets benefit significantly from MMseqs2 clustering:

| Input genomes | Peak RAM |
|---|---|
| 10 | ~10 GB |
| 100 | ~35 GB |
| 1 000 | ~130 GB |

> These peaks were measured when `BUILD_ANNOTATION` and `COLLECT_HITS` always started 160 worker processes. They now start `PARSING_CPUS` workers, which needs much less memory on small machines.

For large real-world phage collections (e.g. ~30 000 genomes from INPHARED), peak RAM was well under 400 GB due to high sequence redundancy converging the cluster count.

### Suggested job allocation

| Input size | CPUs | RAM |
|---|---|---|
| up to ~100 genomes | 16 | 64 GB |
| up to ~1 000 genomes | 32 | 256 GB |
| up to ~10 000 genomes | 64 | 400 GB |

Set matching limits in `nextflow.config`:

```groovy
executor {
    cpus   = 32
    memory = '256 GB'
}
```

## Output

Results are written to the directory specified by `OUTPUT_DIR`.

| File | Description |
|---|---|
| `annotation.tsv` | Per-ORF annotation table with top hits from each database |
| `report.tsv` | Filtered hits passing the probability and coverage thresholds |
| `dropped_inputs.tsv` | Input files that were skipped, with the reason (multi-record FASTA, too short, empty) |
| `failed_hhr.tsv` | HH-suite result files that could not be parsed (database, PC, error); their hits are missing. Header only when all files parsed |
| `search.tsv` | Full unfiltered hits table (if `WRITE_SEARCH_TABLE` is enabled) |
| `PCs2proteins.tsv` | Mapping of protein clusters (PCs) to individual protein IDs |
| `genbanks/<GENOME_ID>.gb` | GenBank-formatted files with annotated features for each input genome |

### Intermediate files

Nextflow stores intermediate files in the `work/` directory. This directory can grow large and may be safely deleted after a successful run. Set `PUBLISH_INTERMEDIATES: true` to also copy the main intermediate files to `OUTPUT_DIR/intermediates/`, using the same folder names as the original Snakemake pipeline.

## Known issues

- **ECOD annotations with the database from `setup_databases.sh`.** The script downloads ECOD v294 F40, whose profile names look like `3520895|A0A077YYV3_F1_nD5|1075.4.1.9`. The ECOD annotation code expects the HHpred ECOD_F70 name format (`... | X: ..., H: ..., T: ..., F: ...`). With v294 F40 every hit gets X-group `-` and its raw ID as the "function", so only one ECOD hit is kept per cluster and the label is not human-readable. Use the HHpred ECOD_F70 database (set `HHSUITE.ECOD`) if you need ECOD functions.
- **Pfam version.** `setup_databases.sh` downloads Pfam-A 35.0; the original MGG pipeline used Pfam-A 32.0.

## Software versions

The containerized environment includes:

| Tool | Version | Purpose |
|---|---|---|
| Python | 3.11 | Pipeline scripts |
| Biopython | 1.79 | Sequence I/O and GenBank export |
| HH-suite | 3.3.0 | Profile-profile searches and MSA enrichment |
| MMseqs2 | 16.747c6 | Protein clustering |
| Prodigal-gv | 2.11.0 | ORF prediction (viral/metagenomic) |
| Glimmer | 3.02 | ORF prediction (trained model); compiled from source in the Docker image |
| csb | 1.2.5 | HHR parsing |
| pandas / NumPy | 3.0.2 / 2.4.4 | Tables |

## Citation

If you use this pipeline, please cite:

> *Citation will be added upon publication.*

## License

This project is licensed under the MIT License.

```
MIT License

Copyright (c) 2024 Mostowy Lab

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```
