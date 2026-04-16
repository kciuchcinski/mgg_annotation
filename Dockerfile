FROM mambaorg/micromamba:latest

# Become root to install OS packages
USER root

# System packages
RUN apt-get update \
 && apt-get install -y --no-install-recommends procps \
 && rm -rf /var/lib/apt/lists/*

# Ensure conda 'base' environment is active during RUN instructions
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# Back to the micromamba user (defined in the base image)
USER $MAMBA_USER

# Install conda and pip packages into the base environment
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        python=3.11 \
        biopython=1.79 \
        pandas=3.0.2 \
        numpy=2.4.4 \
        seaborn \
        matplotlib=3.10.8 \
        hhsuite=3.3.0 \
        mmseqs2=16.747c6 \
        clustalo=1.2.4 \
        bedtools=2.30.0 \
        glimmer=3.02 \
        prodigal-gv=2.11.0 \
        ffindex=0.98 \
        pip && \
    micromamba run -n base pip install \
        phylotreelib \
        csb==1.2.5 && \
    micromamba clean --all --yes

# Match the %environment section
ENV PATH="/opt/conda/bin:${PATH}"