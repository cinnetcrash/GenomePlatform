"""
Güvenlik katmanı:
- Dosya validasyonu (magic bytes, uzantı, boyut)
- Komut enjeksiyonu önleme
- Path traversal koruması
- Rate limiting yardımcıları
"""
import re
import uuid
import hashlib
import logging
from pathlib import Path
from typing import Optional

from config import (
    ALLOWED_EXTENSIONS, MAX_FILE_SIZE_MB,
    GZIP_MAGIC, UPLOAD_DIR
)

logger = logging.getLogger("security")


# ─── Dosya Validasyonu ────────────────────────────────────────────────────────

def validate_filename(filename: str) -> str:
    """
    Güvenli dosya adı döndürür.
    - Path traversal karakterlerini temizler
    - Sadece izin verilen uzantılara izin verir
    - Boş veya şüpheli isimleri reddeder
    """
    if not filename or len(filename) > 255:
        raise ValueError("Geçersiz dosya adı.")

    # Sadece alfanümerik, nokta, tire, alt çizgiye izin ver
    safe_name = re.sub(r"[^\w.\-]", "_", Path(filename).name)

    # Uzantı kontrolü (.fastq.gz gibi çift uzantıları da yakala)
    lower = safe_name.lower()
    allowed = any(lower.endswith(ext) for ext in ALLOWED_EXTENSIONS)
    if not allowed:
        raise ValueError(
            f"Desteklenmeyen dosya türü. İzin verilenler: {ALLOWED_EXTENSIONS}"
        )

    return safe_name


def validate_file_size(size_bytes: int) -> None:
    """Dosya boyutu limitini kontrol eder."""
    max_bytes = MAX_FILE_SIZE_MB * 1024 * 1024
    if size_bytes > max_bytes:
        raise ValueError(
            f"Dosya çok büyük ({size_bytes / 1e6:.0f} MB). "
            f"Maksimum: {MAX_FILE_SIZE_MB} MB."
        )


def validate_magic_bytes(data: bytes, filename: str) -> None:
    """
    Dosyanın ilk byte'larını kontrol eder.
    .fastq.gz → gzip magic (1f 8b)
    .fastq    → '@' ile başlamalı (FASTQ format)
    """
    lower = filename.lower()
    if lower.endswith(".gz"):
        if not data.startswith(GZIP_MAGIC):
            raise ValueError(
                "Dosya .gz uzantılı ama gzip formatında değil. "
                "Lütfen gerçek bir FASTQ.gz dosyası yükleyin."
            )
    else:
        # Plain FASTQ: ilk karakter '@' veya boşluk/yorum satırı olabilir
        text_start = data[:4].decode("ascii", errors="ignore").strip()
        if not text_start.startswith("@"):
            raise ValueError(
                "Dosya FASTQ formatında görünmüyor. "
                "İlk satır '@' ile başlamalıdır."
            )


# ─── Job ID Üretimi ───────────────────────────────────────────────────────────

def generate_job_id() -> str:
    """Tahmin edilemez, URL-güvenli iş ID'si üretir."""
    return uuid.uuid4().hex


def job_upload_dir(job_id: str) -> Path:
    """Her iş için izole bir dizin döndürür — path traversal imkânsız."""
    # job_id yalnızca hex karakterlerinden oluşur (uuid4().hex)
    if not re.fullmatch(r"[0-9a-f]{32}", job_id):
        raise ValueError("Geçersiz job ID.")
    return UPLOAD_DIR / job_id


# ─── Komut Güvenliği ─────────────────────────────────────────────────────────

def sanitize_path(path: Path) -> str:
    """
    Shell komutlarında kullanılacak yolu güvenli hale getirir.
    subprocess list argümanı kullanıldığında gerekli değil,
    ama ek güvenlik için yine de kontrol ederiz.
    """
    resolved = path.resolve()
    # Sadece beklenen dizinler altında olmalı
    allowed_roots = [UPLOAD_DIR.resolve(), Path("/home/analysis/Desktop").resolve()]
    if not any(str(resolved).startswith(str(root)) for root in allowed_roots):
        raise ValueError(f"İzin verilmeyen dizin erişimi: {resolved}")
    return str(resolved)


def safe_sample_name(name: str) -> str:
    """
    Örnek adını shell-güvenli hale getirir.
    Sadece alfanümerik + tire + alt çizgiye izin verir.
    """
    clean = re.sub(r"[^\w\-]", "_", name)
    return clean[:64]  # Maksimum 64 karakter


# ─── Checksum ─────────────────────────────────────────────────────────────────

def compute_md5(path: Path, chunk_size: int = 1 << 20) -> str:
    """Dosyanın MD5 hash'ini hesaplar (bütünlük kontrolü için)."""
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()
