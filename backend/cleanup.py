"""
24 saatlik otomatik silme sistemi.
Arka planda çalışan bir thread her 10 dakikada bir
süresi dolmuş işleri ve dosyalarını siler.
"""
import shutil
import logging
import threading
import time
from pathlib import Path

import database as db
from config import UPLOAD_DIR, RESULTS_DIR

logger = logging.getLogger("cleanup")

SWEEP_INTERVAL_SEC = 600   # Her 10 dakikada bir tara


def delete_job_files(job_id: str) -> None:
    """Bir işe ait tüm dosya ve dizinleri güvenli biçimde siler."""
    deleted = []

    for base_dir in [UPLOAD_DIR, RESULTS_DIR]:
        job_dir = base_dir / job_id
        if job_dir.exists() and job_dir.is_dir():
            # Sembolik link değil mi kontrol et (path traversal güvenliği)
            if not job_dir.is_symlink():
                shutil.rmtree(job_dir, ignore_errors=True)
                deleted.append(str(job_dir))

    if deleted:
        logger.info("İş %s silindi: %s", job_id, deleted)
    else:
        logger.debug("İş %s için silinecek dosya bulunamadı.", job_id)


def sweep() -> None:
    """Süresi dolmuş tüm işleri temizler."""
    expired = db.get_expired_jobs()
    if not expired:
        return

    logger.info("%d adet süresi dolmuş iş bulundu.", len(expired))
    for job in expired:
        try:
            delete_job_files(job["id"])
            db.mark_deleted(job["id"])
        except Exception as exc:
            logger.error("İş %s silinemedi: %s", job["id"], exc)


def _sweep_loop() -> None:
    """Arka plan thread döngüsü."""
    logger.info("Temizleme servisi başlatıldı (her %ds).", SWEEP_INTERVAL_SEC)
    while True:
        try:
            sweep()
        except Exception as exc:
            logger.error("Sweep hatası: %s", exc)
        time.sleep(SWEEP_INTERVAL_SEC)


def start_cleanup_thread() -> threading.Thread:
    """Arka plan temizleme thread'ini başlatır ve döndürür."""
    t = threading.Thread(target=_sweep_loop, daemon=True, name="cleanup-sweep")
    t.start()
    return t
