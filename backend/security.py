"""
Security layer:
- File validation (magic bytes, extension, size)
- Command injection prevention
- Path traversal protection
- Rate limiting helpers
"""
import re
import uuid
import hashlib
import logging
from pathlib import Path

from config import (
    ALLOWED_EXTENSIONS, MAX_FILE_SIZE_MB,
    GZIP_MAGIC, UPLOAD_DIR
)

logger = logging.getLogger("security")


# ─── File Validation ──────────────────────────────────────────────────────────

def validate_filename(filename: str) -> str:
    """
    Returns a sanitised filename.
    - Strips path traversal characters
    - Allows only permitted extensions
    - Rejects empty or suspicious names
    """
    if not filename or len(filename) > 255:
        raise ValueError("Invalid filename.")

    safe_name = re.sub(r"[^\w.\-]", "_", Path(filename).name)

    lower = safe_name.lower()
    allowed = any(lower.endswith(ext) for ext in ALLOWED_EXTENSIONS)
    if not allowed:
        raise ValueError(
            f"Unsupported file type. Allowed extensions: {ALLOWED_EXTENSIONS}"
        )

    return safe_name


def validate_file_size(size_bytes: int) -> None:
    """Checks file size against the configured limit."""
    max_bytes = MAX_FILE_SIZE_MB * 1024 * 1024
    if size_bytes > max_bytes:
        raise ValueError(
            f"File too large ({size_bytes / 1e6:.0f} MB). "
            f"Maximum allowed: {MAX_FILE_SIZE_MB} MB."
        )


def validate_magic_bytes(data: bytes, filename: str) -> None:
    """
    Validates the first bytes of the file.
    .fastq.gz → must start with gzip magic bytes (1f 8b)
    .fastq    → must start with '@' (FASTQ format)
    """
    lower = filename.lower()
    if lower.endswith(".gz"):
        if not data.startswith(GZIP_MAGIC):
            raise ValueError(
                "File has a .gz extension but is not a valid gzip file. "
                "Please upload a real FASTQ.gz file."
            )
    else:
        text_start = data[:4].decode("ascii", errors="ignore").strip()
        if not text_start.startswith("@"):
            raise ValueError(
                "File does not appear to be in FASTQ format. "
                "The first line must start with '@'."
            )


# ─── Job ID Generation ────────────────────────────────────────────────────────

def generate_job_id() -> str:
    """Generates an unpredictable, URL-safe job ID."""
    return uuid.uuid4().hex


def job_upload_dir(job_id: str) -> Path:
    """Returns an isolated upload directory for a job — path traversal impossible."""
    if not re.fullmatch(r"[0-9a-f]{32}", job_id):
        raise ValueError("Invalid job ID.")
    return UPLOAD_DIR / job_id


# ─── Command Safety ───────────────────────────────────────────────────────────

def sanitize_path(path: Path) -> str:
    """
    Validates a path before use in subprocess commands.
    Not strictly necessary when using list args (shell=False),
    but provides an extra layer of defence.
    """
    resolved = path.resolve()
    allowed_roots = [UPLOAD_DIR.resolve(), Path("/home/analysis/Desktop").resolve()]
    if not any(str(resolved).startswith(str(root)) for root in allowed_roots):
        raise ValueError(f"Access to path not allowed: {resolved}")
    return str(resolved)


def safe_sample_name(name: str) -> str:
    """Returns a shell-safe sample name (alphanumeric, hyphens, underscores only)."""
    clean = re.sub(r"[^\w\-]", "_", name)
    return clean[:64]


# ─── Checksum ─────────────────────────────────────────────────────────────────

def compute_md5(path: Path, chunk_size: int = 1 << 20) -> str:
    """Computes the MD5 hash of a file for integrity checking."""
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()
