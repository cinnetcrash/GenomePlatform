"""
GenomePlatform Configuration
Tüm ayarlar tek yerden yönetilir.
"""
import os
from pathlib import Path

# === Dizinler ===
BASE_DIR       = Path(__file__).parent.parent
DATA_DIR       = BASE_DIR / "data"
UPLOAD_DIR     = DATA_DIR / "uploads"
RESULTS_DIR    = DATA_DIR / "results"
LOGS_DIR       = DATA_DIR / "logs"
SCRIPTS_DIR    = Path("/home/analysis/Desktop/Scripts")

for d in [UPLOAD_DIR, RESULTS_DIR, LOGS_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# === Güvenlik ===
MAX_FILE_SIZE_MB   = 2000          # 2 GB (FASTQ.gz için)
ALLOWED_EXTENSIONS = {".fastq", ".fastq.gz", ".fq", ".fq.gz"}
ALLOWED_MIME_TYPES = {"application/gzip", "text/plain", "application/octet-stream"}
MAX_JOBS_PER_IP    = 3             # Aynı IP'den aynı anda en fazla 3 iş
RATE_LIMIT         = "5/minute"    # Upload rate limiti

# FASTQ magic bytes (gzip: 1f 8b, plain: '@' veya '#')
GZIP_MAGIC = b"\x1f\x8b"

# === İş Yaşam Döngüsü ===
JOB_EXPIRY_HOURS = 24              # 24 saat sonra otomatik silinir

# === Conda Ortamları (mevcut kuruluma göre doğrulandı) ===
CONDA_BASE = Path("/home/analysis/miniconda3")
CONDA_ENVS = {
    "qc":       "analiz",    # NanoPlot, Flye
    "assembly": "shovill",   # Shovill, SPAdes
    "bakta":    "bakta",     # Bakta annotasyonu
    "mlst":     "mlst",      # CGE MLST
}

# === Araç Yolları (sistem PATH ve conda env'lerinden) ===
TOOLS = {
    "fastp":     "/usr/bin/fastp",           # sistem PATH'te
    "nanoplot":  "NanoPlot",                 # analiz env'inde
    "flye":      "flye",                     # analiz env'inde
    "shovill":   "shovill",                  # shovill env'inde
    "mlst":      "mlst",                     # mlst env'inde
    "amrfinder": "amrfinder",
    "bakta":     "bakta",                    # bakta env'inde
    "mafft":     str(CONDA_BASE / "bin" / "mafft"),  # base'de
    "primer3":   "/usr/bin/primer3_core",    # sistem PATH'te
}

# === Claude API ===
ANTHROPIC_API_KEY = os.getenv("ANTHROPIC_API_KEY", "")
CLAUDE_MODEL      = "claude-sonnet-4-6"

# === Sunucu ===
HOST = "0.0.0.0"
PORT = 8000
