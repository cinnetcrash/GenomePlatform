"""
24-hour automatic file deletion system.
A background thread sweeps every 10 minutes to delete
expired jobs and their associated files.
"""
import shutil
import logging
import threading
import time
from pathlib import Path

import database as db
from config import UPLOAD_DIR, RESULTS_DIR

logger = logging.getLogger("cleanup")

SWEEP_INTERVAL_SEC = 600   # Sweep every 10 minutes


def delete_job_files(job_id: str) -> None:
    """Safely deletes all files and directories associated with a job."""
    deleted = []

    for base_dir in [UPLOAD_DIR, RESULTS_DIR]:
        job_dir = base_dir / job_id
        if job_dir.exists() and job_dir.is_dir():
            # Guard against symlink-based path traversal
            if not job_dir.is_symlink():
                shutil.rmtree(job_dir, ignore_errors=True)
                deleted.append(str(job_dir))

    if deleted:
        logger.info("Job %s deleted: %s", job_id, deleted)
    else:
        logger.debug("No files found to delete for job %s.", job_id)


def sweep() -> None:
    """Cleans up all expired jobs."""
    expired = db.get_expired_jobs()
    if not expired:
        return

    logger.info("Found %d expired job(s).", len(expired))
    for job in expired:
        try:
            delete_job_files(job["id"])
            db.mark_deleted(job["id"])
        except Exception as exc:
            logger.error("Failed to delete job %s: %s", job["id"], exc)


def _sweep_loop() -> None:
    """Background thread loop."""
    logger.info("Cleanup service started (every %ds).", SWEEP_INTERVAL_SEC)
    while True:
        try:
            sweep()
        except Exception as exc:
            logger.error("Sweep error: %s", exc)
        time.sleep(SWEEP_INTERVAL_SEC)


def start_cleanup_thread() -> threading.Thread:
    """Starts and returns the background cleanup thread."""
    t = threading.Thread(target=_sweep_loop, daemon=True, name="cleanup-sweep")
    t.start()
    return t
