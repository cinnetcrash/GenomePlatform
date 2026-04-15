"""
LycianWay — Central Configuration
"""
import os
import math
import multiprocessing
from pathlib import Path

# === Directories ===
BASE_DIR    = Path(__file__).parent.parent
DATA_DIR    = BASE_DIR / "data"
UPLOAD_DIR  = DATA_DIR / "uploads"
RESULTS_DIR = DATA_DIR / "results"
LOGS_DIR    = DATA_DIR / "logs"
SCRIPTS_DIR = Path("/home/analysis/Desktop/Scripts")

for d in [UPLOAD_DIR, RESULTS_DIR, LOGS_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# === CPU Budget (80% of available cores) ===
_total_cpus  = multiprocessing.cpu_count()
MAX_THREADS  = max(1, math.floor(_total_cpus * 0.80))   # e.g. 16 cores → 12 threads
print(f"[config] CPU budget: {MAX_THREADS}/{_total_cpus} threads (80%)")

# === Security ===
MAX_FILE_SIZE_MB   = 2000
ALLOWED_EXTENSIONS = {".fastq", ".fastq.gz", ".fq", ".fq.gz"}
ALLOWED_MIME_TYPES = {"application/gzip", "text/plain", "application/octet-stream"}
MAX_JOBS_PER_IP    = 3
RATE_LIMIT         = "5/minute"

# FASTQ magic bytes (gzip: 1f 8b)
GZIP_MAGIC = b"\x1f\x8b"

# === Job Lifecycle ===
JOB_EXPIRY_HOURS = 24   # Files auto-deleted after 24 hours

# === Conda Environments (validated against this machine) ===
CONDA_BASE = Path("/home/analysis/miniconda3")
CONDA_ENVS = {
    "qc":       "analiz",    # NanoPlot, Flye
    "assembly": "shovill",   # Shovill, SPAdes
    "bakta":    "bakta",     # Bakta annotation
    "mlst":     "mlst",      # CGE MLST
}

# === Tool Paths ===
TOOLS = {
    "fastp":     "/usr/bin/fastp",
    "nanoplot":  "NanoPlot",
    "flye":      "flye",
    "shovill":   "shovill",
    "mlst":      "mlst",
    "amrfinder": "amrfinder",
    "bakta":     "bakta",
    "mafft":     str(CONDA_BASE / "bin" / "mafft"),
    "primer3":   "/usr/bin/primer3_core",
}

# === Claude API ===
ANTHROPIC_API_KEY = os.getenv("ANTHROPIC_API_KEY", "")
CLAUDE_MODEL      = "claude-sonnet-4-6"

# === Server ===
HOST = "0.0.0.0"
PORT = 8000
