# another phage annotation tool

__WARNING:__ Use sample fasta headers, like __accession numbers__ or __simple names__ e.g., `PHAGE0001`, `PHAGE0002`

## Installation and Execution

### 1. Install apptainer, snakemake and dependencies

The simplest way is to use `conda`:
```bash
conda create -n snakemake -c conda-forge -c bioconda snakemake pandas biopython=1.79 mamba
conda activate snakemake
```

You can also install the dependencies with python in a virtual environment:
``bash
python3 -m venv annotation_env
source annotation_env/bin/activate
pip3 install snakemake pandas BioPython
```

For apptainer, refer to the installation instruction available in [Apptainer Docs](https://apptainer.org/docs/admin/main/installation.html#)

### 2. Build the Singularity container

Build the container from the definition file. This will create a `phage_annotation.sif` file containing all the necessary dependencies.

```bash
singularity build phage_annotation.sif phage_annotation.def
```

### 3. Configure the pipeline

The pipeline requires access to several databases. You need to provide the paths to these databases in a configuration file.

First, copy the template configuration file:

```bash
cp config.template.yml config.yml
```

Then, open `config.yml` in a text editor and replace the placeholder paths (e.g., `/path/to/your/databases/...`) with the actual paths to databases and input/output directories

### 4. Run the pipeline

You can now run the pipeline using the `snakemake` command.

```bash
snakemake --cores all --snakefile ANNOTATION --configfile config.yml --sdm apptainer
```