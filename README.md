# MGG_ANNOTATION_PIPELINE
A Nextflow pipeline for annotating phage genomes using PHROGs-based enrichment and HMM-based searches for precise and sensitive annotation.

## Setup
You will need a POSIX-compatible system (Linux or macOS) with the following installed:

1. Nextflow

```bash
curl -s https://get.nextflow.io | bash
chmod +x nextflow
mv nextflow /usr/local/bin/  # Optional: move to path
```

Please note that nextflow requires Java to run. Check out [official nextflow docs](https://www.nextflow.io/docs/latest/install.html#install-page) for recommended install methods.


2. Apptainer (recommended)

Linux: Follow official installation instructions [here](https://github.com/apptainer/apptainer/blob/release-1.4/INSTALL.md)

Before you run the pipeline, build the container with:

```bash
apptainer build phage_annotation.sif phage_annotation.def
```

2. Conda (alternative)

If, for some reason, you cannot use apptainer, you can fall back to using Conda. We recommend using Micromamba or Mamba for faster installation, but standard Conda works too (see [this](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html) for installation instructions)

```bash
# Using Micromamba
micromamba create -f environment.yml
micromamba activate phage_annotation

# OR using Conda
conda env create -f environment.yml
conda activate phage_annotation
```

Then, set the `conda` param in `nextflow.config` to point to your env:
```
  conda {
    process {
      withName: '.*' {
        conda = '/path/to/micromamba/envs/phage_annotation'
      }
    }
    ...
```

Note: if using micromamba, make sure to set `useMicromamba: true` in `nextflow.config`

After the environment is created, you can skip setting up the container (step 4). When running the pipeline, replace the `-profile apptainer` flag with `-profile conda` and point to the environment file.

```bash
./nextflow run main.nf \
  -profile conda \
  -params-file config.yml
```

2. Docker (recommended for macOS)

On macOS (especially Apple Silicon), using Docker is often simpler than Apptainer or full Conda setups. This pipeline provides a Docker image that mirrors the Apptainer environment.

Start by installing docker - [official docs](https://docs.docker.com/desktop/setup/install/mac-install/)

Then, build the Docker image:

```bash
docker build --platform=linux/amd64 -t phage_annotation:latest .
```

Finally, run the pipeline with Docker:

```bash
./nextflow run main.nf \
  -profile docker \
  -params-file config.yml
```

On macOS, make sure Docker Desktop is installed and running before starting the pipeline.

3. Clone the repository:
```bash
git clone <your-repo-url>
cd <your-repo-name>
```


## Configuration
CRITICAL: You must edit the configuration files to match your local environment before running the pipeline. The default configuration contains hardcoded paths specific to a cluster environment.

1. Update config.yml

Open config.yml and modify the following paths:

PHAGES_DIR: Path to your directory containing input .fasta files.

OUTPUT_DIR: Desired output location.

DB_ROOT: (Important) Update this to the absolute path where you have stored the required MGG databases.

2. Update nextflow.config

Open nextflow.config and locate the runOptions line inside the singularity scope.

Action: Change the bind path (-B) to match your local storage, or remove it if not needed.

Example: Change -B /net/storage:/net/storage to -B /path/to/my/data:/data.


## Running the Pipeline
Once configured, run the pipeline using the command below. This command mounts the container and loads your custom parameters.

bash
./nextflow run main.nf \
  -profile apptainer \
  -params-file config.yml

If a run fails, you can try the `-resume` flag. It will restart the pipeline from the last successful step instead of beginning from scratch.

## Outputs
Results will be saved to the directory specified in config.yml (Default: output_test).

annotation.tsv: Detailed per-ORF annotation table.

report.tsv: All hits passing the filtering step

search.tsv: Full raw hits table (if enabled).

genbanks/GENOME_ID.gb: GenBank files with per-phage annotation results

Intermediate files: Located in the work/ directory (can be deleted after successful completion).