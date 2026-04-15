# ─────────────────────────────────────────────────────────────────────────────
# LycianWay — Dockerfile
#
# Multi-stage build:
#   Stage 1 (rust-builder): Compiles Autocycler from source
#   Stage 2 (main):         Installs all conda envs + the FastAPI web app
#
# Key design decisions:
#   - CONDA_BASE=/opt/conda  (micromamba default)
#   - fastp + kraken2 installed via apt (system-wide, no conda env needed)
#   - Each bioinformatics conda env is a separate RUN layer for Docker cache
#   - Databases are NEVER embedded — they must be mounted as volumes
#   - Bandage installed from bioconda inside the 'analiz' env
#   - Autocycler compiled from source in the rust-builder stage
# ─────────────────────────────────────────────────────────────────────────────

# ── Stage 1: Build Autocycler ─────────────────────────────────────────────────
FROM rust:1.78-slim-bookworm AS rust-builder

RUN apt-get update -qq && apt-get install -y --no-install-recommends \
      git ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN git clone --depth 1 https://github.com/rrwick/Autocycler.git /opt/autocycler-src \
    && cd /opt/autocycler-src \
    && cargo build --release


# ── Stage 2: Main application image ──────────────────────────────────────────
FROM mambaorg/micromamba:1.5.10-noble

LABEL maintainer="gultekinnunal@gmail.com"
LABEL description="LycianWay — Automated bacterial genomic analysis platform"
LABEL version="0.1.0"

USER root

# ── System packages ───────────────────────────────────────────────────────────
# fastp and kraken2 are installed system-wide (same as the host setup)
RUN apt-get update -qq && apt-get install -y --no-install-recommends \
      fastp \
      kraken2 \
      curl \
      wget \
      procps \
      ca-certificates \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# ── Autocycler binary from builder stage ─────────────────────────────────────
COPY --from=rust-builder /opt/autocycler-src/target/release/autocycler /usr/local/bin/autocycler

# ── Conda environments ────────────────────────────────────────────────────────
# Each env is a separate layer so Docker can cache independently.
# Channels: conda-forge (priority), bioconda, defaults

# analiz: Flye (long-read assembly), AMRFinderPlus, Abricate, NanoPlot, Bandage
RUN micromamba create -n analiz -y \
      -c conda-forge -c bioconda -c defaults \
      "flye>=2.9" \
      "ncbi-amrfinderplus>=3.12" \
      abricate \
      nanoplot \
      bandage \
    && micromamba clean --all --yes

# shovill: Illumina short-read assembly (wraps SPAdes)
RUN micromamba create -n shovill -y \
      -c conda-forge -c bioconda -c defaults \
      shovill \
    && micromamba clean --all --yes

# mlst: MLST sequence typing
RUN micromamba create -n mlst -y \
      -c conda-forge -c bioconda -c defaults \
      mlst \
    && micromamba clean --all --yes

# mobsuite: MOB-Suite plasmid typing and reconstruction
RUN micromamba create -n mobsuite -y \
      -c conda-forge -c bioconda -c defaults \
      mob_suite \
    && micromamba clean --all --yes

# quast5: QUAST assembly quality assessment
RUN micromamba create -n quast5 -y \
      -c conda-forge -c bioconda -c defaults \
      "quast>=5.2" \
    && micromamba clean --all --yes

# checkM: CheckM2 genome completeness/contamination
RUN micromamba create -n checkM -y \
      -c conda-forge -c bioconda -c defaults \
      checkm2 \
    && micromamba clean --all --yes

# checkv: CheckV viral genome quality (optional, requires DB volume)
RUN micromamba create -n checkv -y \
      -c conda-forge -c bioconda -c defaults \
      checkv \
    && micromamba clean --all --yes

# bakta: Bakta rapid bacterial genome annotation (optional)
RUN micromamba create -n bakta -y \
      -c conda-forge -c bioconda -c defaults \
      bakta \
    && micromamba clean --all --yes

# ── Application user ──────────────────────────────────────────────────────────
RUN useradd -m -u 1000 -s /bin/bash appuser

# ── Python web application ────────────────────────────────────────────────────
WORKDIR /app
COPY requirements.txt /app/

# Install into micromamba's base Python environment
RUN /opt/conda/bin/pip install --no-cache-dir -r requirements.txt

COPY backend/   /app/backend/
COPY frontend/  /app/frontend/

RUN mkdir -p /app/data/{uploads,results,logs} \
    && chown -R appuser:appuser /app

# ── Runtime environment variables ─────────────────────────────────────────────
# These point to tools/paths inside the container.
# All database paths must be provided via docker-compose environment or .env
ENV CONDA_BASE=/opt/conda
ENV BANDAGE_BIN=/opt/conda/envs/analiz/bin/Bandage
ENV AUTOCYCLER_BIN=/usr/local/bin/autocycler
ENV PYTHONUNBUFFERED=1

USER appuser

EXPOSE 8000

HEALTHCHECK --interval=30s --timeout=10s --start-period=60s --retries=3 \
  CMD curl -f http://localhost:8000/health || exit 1

WORKDIR /app/backend

CMD ["/opt/conda/bin/python", "-m", "uvicorn", "main:app", \
     "--host", "0.0.0.0", "--port", "8000", "--log-level", "info"]
