# MGG Annotation Pipeline

A Nextflow pipeline for annotating phage genomes using PHROGs-enriched profile HMM searches across multiple databases for sensitive and precise functional annotation.

## Overview

The pipeline takes phage genome assemblies (FASTA) as input and produces per-ORF functional annotations by:

1. **Preprocessing** — validates input genomes and filters by minimum length.
2. **ORF calling** — predicts open reading frames using both Prodigal-gv (meta mode) and Glimmer3, then reconciles predictions.
3. **Protein clustering** — clusters translated ORFs with MMseqs2 to reduce redundancy and generate multiple sequence alignments (MSAs).
4. **MSA enrichment** — enriches MSAs against the PHROGs database using HHblits to improve search sensitivity.
5. **Database searches** — searches enriched profiles against four HHsuite-formatted databases: PHROGs, PFAM, ECOD, and AlanDB.
6. **Annotation & export** — collects and filters hits, builds a per-ORF annotation table, and exports GenBank files.

## Requirements

- **OS**: Linux (recommended) or macOS
- **Nextflow** ≥ 22.10 (requires Java 11+)
- **Container runtime** (one of):
  - [Apptainer](https://apptainer.org/) ≥ 1.1 — recommended for HPC / Linux
  - [Docker](https://www.docker.com/) — recommended for macOS / local workstations
  - [Conda](https://docs.conda.io/) / [Micromamba](https://mamba.readthedocs.io/) — fallback if containers are not available

## Installation

### 1. Install Nextflow

```bash
curl -s https://get.nextflow.io | bash
chmod +x nextflow

# Optionally 
mv nextflow ~/.local/bin/   # or any directory on your $PATH
```

See the [official Nextflow docs](https://www.nextflow.io/docs/latest/install.html) for alternative installation methods and Java requirements.

### 2. Clone the repository

```bash
git clone https://github.com/kciuchcinski/mgg_annotation.git
cd mgg_annotation
```

### 3. Set up the execution environment

Choose **one** of the three options below.

<details>
<summary><strong>Option A — Apptainer (recommended for HPC / Linux)</strong></summary>

Build the Apptainer image from the provided definition file:

```bash
apptainer build phage_annotation.sif phage_annotation.def
```

The resulting `.sif` file must be in the pipeline root directory (or update the `container` path in `nextflow.config`). Host paths (including `DB_ROOT`) are mounted automatically via `autoMounts`.

</details>

<details>
<summary><strong>Option B — Docker (recommended for macOS)</strong></summary>

Build the Docker image:

```bash
docker build --platform=linux/amd64 -t phage_annotation:latest .
```

> **Note:** On Apple Silicon Macs, the `--platform=linux/amd64` flag is required because several bioinformatics tools lack native ARM builds.

Make sure Docker Desktop is running before launching the pipeline.

</details>

<details>
<summary><strong>Option C — Conda / Micromamba (fallback)</strong></summary>

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

Replace `/path/to/your/envs/phage_annotation` with the actual path (e.g. the output of `micromamba env list` or `conda env list`). If using Micromamba, ensure `useMicromamba = true` is set in `nextflow.config`.

</details>

## Database setup

The pipeline requires the following databases and metadata files to be available locally. All database paths are resolved relative to a single root directory (`DB_ROOT` in `config.yml`), but individual paths can be overridden.

### Required databases

| Database | Expected path (relative to `DB_ROOT`) | Description |
|---|---|---|
| Glimmer training file | `training-file_refseq.icm` | Pre-trained ICM for Glimmer3 ORF prediction |
| PHROGs (HHsuite) | `PHROGS_v4/phrogs` | Prokaryotic Virus Remote Homologous Groups v4, HHsuite-formatted |
| PFAM (HHsuite) | `pfamA_32/pfam` | Pfam-A v32, HHsuite-formatted |
| ECOD (HHsuite) | `ECOD_F70_20230309/ECOD_F70_20230309` | ECOD domain database (F70, 2023-03-09), HHsuite-formatted |
| AlanDB (HHsuite) | `AlanDavidson/profile-db/all_proteins` | Alan Davidson phage protein database, HHsuite-formatted |

### Required metadata tables

| File | Expected path (relative to `DB_ROOT`) | Description |
|---|---|---|
| PHROGs annotation table | `tables/phrog_annot_v4.tsv` | Functional category mappings for PHROGs |
| AlanDB annotation table | `tables/alan_annot.tsv` | Functional annotations for AlanDB entries |

### Expected directory layout

```
DB_ROOT/
├── training-file_refseq.icm
├── PHROGS_v4/
│   └── phrogs{_a3m.ffdata, _a3m.ffindex, _hhm.ffdata, _hhm.ffindex, ...}
├── pfamA_32/
│   └── pfam{_a3m.ffdata, _a3m.ffindex, _hhm.ffdata, _hhm.ffindex, ...}
├── ECOD_F70_20230309/
│   └── ECOD_F70_20230309{_a3m.ffdata, _a3m.ffindex, ...}
├── AlanDavidson/
│   └── profile-db/
│       └── all_proteins{_a3m.ffdata, _a3m.ffindex, ...}
└── tables/
    ├── phrog_annot_v4.tsv
    └── alan_annot.tsv
```

<!-- TODO: Add download links / setup script for each database before publication -->

## Configuration

All user-configurable parameters are set in `config.yml`. Copy and edit the example before running:

```bash
cp config.yml my_config.yml
```

### Minimal configuration

You only need to set three paths:

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
| `PHAGES_DIR` | `test/data/` | Directory containing input FASTA files |
| `INPUT_EXTENSION` | `fasta` | File extension to look for in `PHAGES_DIR` |
| `OUTPUT_DIR` | `test/output` | Directory for pipeline outputs |
| `DB_ROOT` | `databases/` | Root directory for all databases |
| `PHAGE_MIN_LENGTH` | `2000` | Minimum genome length in bp to process |
| `WRITE_SEARCH_TABLE` | `true` | Write raw HHsuite hits to `search.tsv` |
| `WRITE_ANNOTATION_TABLE` | `true` | Write per-ORF annotation to `annotation.tsv` |
| `N_FUNCTIONS_PER_DB` | `2` | Max number of functions reported per protein per database |

#### Annotation filtering thresholds

These thresholds control which HHsuite hits are retained in the final annotation. Separate thresholds are applied to PHROGs and to the domain databases (PFAM, ECOD, AlanDB).

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
| `PARSING_CPUS` | `16` | CPU threads for hit collection and annotation parsing |

> **Note:** These parameters control CPU allocation per process, not the total resources available to the pipeline. See the section below on how to set the global resource pool.

#### Resource limits

By default, Nextflow uses the **local executor**, which runs all tasks on the current machine. You can cap the total resources it may use by adding an `executor` block to `nextflow.config`:

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

Select the profile matching your execution environment:

```bash
# Apptainer (HPC / Linux)
nextflow run main.nf -profile apptainer -params-file config.yml

# Docker (macOS / local)
nextflow run main.nf -profile docker -params-file config.yml

# Conda
nextflow run main.nf -profile conda -params-file config.yml
```

### Resuming a failed run

If a run fails or is interrupted, relaunch with `-resume` to restart from the last successful step:

```bash
nextflow run main.nf -profile apptainer -params-file config.yml -resume
```

Nextflow caches intermediate results in the `work/` directory. Do not delete it if you intend to resume.

## Output

Results are written to the directory specified by `OUTPUT_DIR`.

| File | Description |
|---|---|
| `annotation.tsv` | Per-ORF annotation table with top hits from each database (if `WRITE_ANNOTATION_TABLE` is enabled) |
| `report.tsv` | Filtered hits passing the probability and coverage thresholds |
| `search.tsv` | Full unfiltered hits table (if `WRITE_SEARCH_TABLE` is enabled) |
| `PCs2proteins.tsv` | Mapping of protein clusters (PCs) to individual protein IDs |
| `genbanks/<GENOME_ID>.gb` | GenBank-formatted files with annotated features for each input genome |

### Intermediate files

Nextflow stores intermediate files in the `work/` directory. This directory can grow large and may be safely deleted after a successful run.

## Pipeline DAG

```
FASTA files
    │
    ▼
CHECK_INPUT ──────────────────────────────────────┐
    │                                              │
    ├──▶ PRODIGAL ──▶ CONCAT_PRODIGAL ────────┐   │
    │                                          │   │
    └──▶ GLIMMER3 ──▶ CONCAT_GLIMMER          │   │
              │                                │   │
              ▼                                │   │
         FILTER_GLIMMER                        │   │
              │                                │   │
              ▼                                ▼   │
         PROCESSING ◀─────────────────────────────┘
              │
              ▼
    EXTRACT_TRANSLATE_ORFS
         │          │
         ▼          ▼
  REORGANIZE    CONCAT_PROTEINS
    TABLE            │
      │              ▼
      │         CLUSTERING
      │          │       │
      │          ▼       ▼
      │    PCS2PROTEINS  PREPARE_A3M
      │          │           │
      │          │           ▼
      │          │     BUILD_FFINDEX
      │          │           │
      |          │           ▼
      │          │    ENRICH_MSA_PHROGS
      │          │           │
      │          │           ▼
      │          │    HHSUITE_SEARCH_BATCH (×4 DBs)
      │          │           │
      │          │           ▼
      │          │       UNPACK_HHR
      │          │           │
      │          │           ▼
      │          │       GATHER_HHR
      │          │           │
      │          ▼           ▼
      │        COLLECT_HITS
      │              │
      │              ▼
      │         FILTER_HITS
      │              │
      │              ▼
      └───▶  BUILD_ANNOTATION
                     │
                     ▼
                  GENBANK
```

## Software versions

The containerized environment includes:

| Tool | Version | Purpose |
|---|---|---|
| Python | 3.11 | Pipeline scripts |
| Biopython | 1.79 | Sequence I/O and GenBank export |
| HH-suite | (latest via bioconda) | Profile-profile searches and MSA enrichment |
| MMseqs2 | 13.45111 | Protein clustering |
| Prodigal-gv | 2.11.0 | ORF prediction (viral/metagenomic) |
| Glimmer | 3.02 | ORF prediction (trained model) |
| Clustal Omega | 1.2.4 | Multiple sequence alignment |
| BEDTools | 2.30.0 | Genomic interval operations |
| FFindex | (latest via bioconda) | Flat-file index operations for HHsuite |

## Citation

<!-- TODO: Add citation when published -->

If you use this pipeline, please cite:

> *Citation pending.*

## License

<!-- TODO: Add license information -->
