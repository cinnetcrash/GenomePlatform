# backend/scheduler.py
import threading
from config import MAX_THREADS

_lock: threading.Lock = threading.Lock()
_active_jobs: set[str] = set()


def job_start(job_id: str) -> None:
    with _lock:
        _active_jobs.add(job_id)


def job_done(job_id: str) -> None:
    with _lock:
        _active_jobs.discard(job_id)


def threads_for_stage() -> int:
    with _lock:
        n = len(_active_jobs)
    return max(1, MAX_THREADS // max(1, n))


def active_count() -> int:
    with _lock:
        return len(_active_jobs)
