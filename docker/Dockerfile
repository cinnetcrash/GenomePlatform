# ─────────────────────────────────────────────────────────────────────────────
# GenomePlatform — Dockerfile
# Strateji: Python + temel bioinformatik araçları Conda ile kurulur.
# Ağır araçlar (Flye, Shovill) için host'un Conda ortamları bind-mount edilir.
# ─────────────────────────────────────────────────────────────────────────────
FROM mambaorg/micromamba:1.5.8-jammy

LABEL maintainer="gultekinnunal@gmail.com"
LABEL description="GenomePlatform - Automated genomic analysis platform"
LABEL version="0.1.0"

# Root olarak sistem bağımlılıklarını kur
USER root
RUN apt-get update -qq && apt-get install -y --no-install-recommends \
      primer3 \
      mafft \
      fastp \
      mlst \
      curl \
      procps \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Uygulama kullanıcısı oluştur (root olarak çalışma — güvenlik)
RUN useradd -m -u 1000 -s /bin/bash appuser

# Micromamba ortamı oluştur
COPY docker/environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml \
    && micromamba clean --all --yes

# Uygulama dosyalarını kopyala
WORKDIR /app
COPY backend/   /app/backend/
COPY frontend/  /app/frontend/
COPY requirements.txt /app/

# Python bağımlılıklarını kur
RUN /opt/conda/bin/pip install --no-cache-dir -r requirements.txt

# Veri dizinlerini oluştur
RUN mkdir -p /app/data/{uploads,results,logs} \
    && chown -R appuser:appuser /app/data

# Güvenli kullanıcıya geç
USER appuser

# Sağlık kontrolü
HEALTHCHECK --interval=30s --timeout=10s --start-period=30s --retries=3 \
  CMD curl -f http://localhost:8000/health || exit 1

EXPOSE 8000

# Başlatma komutu
CMD ["/opt/conda/bin/python", "-m", "uvicorn", "main:app", \
     "--host", "0.0.0.0", "--port", "8000", \
     "--log-level", "info"]

WORKDIR /app/backend
