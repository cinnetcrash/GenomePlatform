"""
Automatic file cleanup system.

Two sweep strategies run together:
  1. Expiry sweep  — deletes all files + DB row for jobs older than 24 h
  2. Overflow sweep — for jobs beyond the most recent MAX_STORED_REPORTS,
                      deletes upload + result working files but keeps the
                      report HTML and the DB row (history stays visible)
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
MAX_STORED_REPORTS = 10    # Keep working files for the most recent N finished jobs


def _safe_rmtree(path: Path) -> bool:
    """Removes a directory tree if it exists and is not a symlink. Returns True if removed."""
    if path.exists() and path.is_dir() and not path.is_symlink():
        shutil.rmtree(path, ignore_errors=True)
        return True
    return False


def delete_job_files(job_id: str) -> None:
    """Deletes all upload + result files for a job (report HTML is inside results dir)."""
    deleted = []
    for base_dir in [UPLOAD_DIR, RESULTS_DIR]:
        if _safe_rmtree(base_dir / job_id):
            deleted.append(str(base_dir / job_id))
    if deleted:
        logger.info("Job %s: deleted files %s", job_id, deleted)


def prune_job_workfiles(job_id: str, report_path: str | None) -> None:
    """
    Deletes working files (uploads + intermediate results) for an overflow job,
    but preserves the report HTML so history remains viewable.
    """
    # Remove uploads entirely
    _safe_rmtree(UPLOAD_DIR / job_id)

    # Remove everything in results dir except the report HTML
    results_dir = RESULTS_DIR / job_id
    if results_dir.exists() and not results_dir.is_symlink():
        report = Path(report_path) if report_path else None
        for child in list(results_dir.iterdir()):
            if report and child.resolve() == report.resolve():
                continue   # keep the report
            if child.is_dir():
                shutil.rmtree(child, ignore_errors=True)
            else:
                child.unlink(missing_ok=True)

    logger.info("Job %s: pruned working files (report kept)", job_id)


def sweep() -> None:
    """Runs both expiry and overflow sweeps."""
    # 1. Expired jobs (24 h) — full delete including report
    expired = db.get_expired_jobs()
    if expired:
        logger.info("Expiry sweep: %d job(s).", len(expired))
        for job in expired:
            try:
                delete_job_files(job["id"])
                db.mark_deleted(job["id"])
            except Exception as exc:
                logger.error("Failed to delete job %s: %s", job["id"], exc)

    # 2. Overflow jobs — keep files for only the last MAX_STORED_REPORTS
    overflow = db.get_overflow_jobs(keep=MAX_STORED_REPORTS)
    if overflow:
        logger.info("Overflow sweep: pruning working files for %d job(s).", len(overflow))
        for job in overflow:
            try:
                prune_job_workfiles(job["id"], job.get("report_path"))
                db.mark_files_deleted(job["id"])
            except Exception as exc:
                logger.error("Failed to prune job %s: %s", job["id"], exc)


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
