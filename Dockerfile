# Multi-arch image (linux/amd64 and linux/arm64, e.g. Apple Silicon Macs).
# Build with:  docker build -t phage_annotation:latest .
ARG MICROMAMBA_TAG=2.9.0-debian13

# ---------------------------------------------------------------------------
# Stage 1: Glimmer 3.02 from source (bioconda has no linux-aarch64 build).
# Same tarball and checksum as the bioconda recipe. CXXDEFS is cleared because
# the default (-D__cplusplus) breaks modern libstdc++ headers.
# ---------------------------------------------------------------------------
FROM mambaorg/micromamba:${MICROMAMBA_TAG} AS glimmer
USER root
RUN apt-get update \
 && apt-get install -y --no-install-recommends ca-certificates curl g++ make \
 && rm -rf /var/lib/apt/lists/*
WORKDIR /build
RUN curl -fsSL -o glimmer302b.tar.gz https://ccb.jhu.edu/software/glimmer/glimmer302b.tar.gz \
 && echo "344d012ae12596de905866fe9eb7f16c  glimmer302b.tar.gz" | md5sum -c - \
 && tar xzf glimmer302b.tar.gz \
 && make -C glimmer3.02/src all CXXDEFS= \
 && mkdir -p /opt/glimmer/bin \
 && cp glimmer3.02/bin/* glimmer3.02/scripts/* /opt/glimmer/bin/

# ---------------------------------------------------------------------------
# Stage 2: runtime image
# ---------------------------------------------------------------------------
FROM mambaorg/micromamba:${MICROMAMBA_TAG}

# procps provides `ps`, which Nextflow uses for task metrics
USER root
RUN apt-get update \
 && apt-get install -y --no-install-recommends procps \
 && rm -rf /var/lib/apt/lists/*

COPY --from=glimmer /opt/glimmer/bin/ /usr/local/bin/

# Ensure conda 'base' environment is active during RUN instructions
ARG MAMBA_DOCKERFILE_ACTIVATE=1
USER $MAMBA_USER

# Keep in sync with environment.yml and phage_annotation.def (glimmer comes from stage 1)
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        python=3.11 \
        biopython=1.79 \
        pandas=3.0.2 \
        numpy=2.4.4 \
        hhsuite=3.3.0 \
        mmseqs2=16.747c6 \
        prodigal-gv=2.11.0 \
        pip \
 && micromamba run -n base pip install --no-cache-dir csb==1.2.5 \
 && micromamba clean --all --yes

ENV PATH="/opt/conda/bin:${PATH}"
