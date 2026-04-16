"""
SQLite database — job status and timestamps.
"""
import sqlite3
import json
import logging
from datetime import datetime, timezone
from pathlib import Path

DB_PATH = Path(__file__).parent.parent / "data" / "jobs.db"
logger = logging.getLogger("database")


def get_conn() -> sqlite3.Connection:
    conn = sqlite3.connect(str(DB_PATH), check_same_thread=False)
    conn.row_factory = sqlite3.Row
    return conn


def init_db() -> None:
    """Creates database tables if they do not exist."""
    with get_conn() as conn:
        conn.executescript("""
            CREATE TABLE IF NOT EXISTS jobs (
                id            TEXT PRIMARY KEY,
                status        TEXT NOT NULL DEFAULT 'queued',
                created_at    TEXT NOT NULL,
                expires_at    TEXT NOT NULL,
                filename      TEXT NOT NULL,
                file_md5      TEXT,
                ip_hash       TEXT NOT NULL,
                read_type     TEXT,
                stages        TEXT DEFAULT '{}',
                error         TEXT,
                report_path   TEXT,
                files_deleted INTEGER NOT NULL DEFAULT 0
            );

            CREATE INDEX IF NOT EXISTS idx_expires ON jobs(expires_at);
            CREATE INDEX IF NOT EXISTS idx_ip      ON jobs(ip_hash);
        """)
        # Migrate existing databases — safe no-op if column already present
        try:
            conn.execute(
                "ALTER TABLE jobs ADD COLUMN files_deleted INTEGER NOT NULL DEFAULT 0"
            )
        except Exception:
            pass
    logger.info("Database ready: %s", DB_PATH)


def create_job(job_id: str, filename: str, expires_at: str,
               ip_hash: str, file_md5: str = "") -> None:
    now = datetime.now(timezone.utc).isoformat()
    with get_conn() as conn:
        conn.execute(
            """INSERT INTO jobs
               (id, status, created_at, expires_at, filename, file_md5, ip_hash)
               VALUES (?, 'queued', ?, ?, ?, ?, ?)""",
            (job_id, now, expires_at, filename, file_md5, ip_hash)
        )


def update_job_status(job_id: str, status: str,
                      error: str = None,
                      read_type: str = None,
                      report_path: str = None) -> None:
    fields, vals = [], []
    fields.append("status = ?");     vals.append(status)
    if error is not None:
        fields.append("error = ?");  vals.append(error)
    if read_type is not None:
        fields.append("read_type = ?"); vals.append(read_type)
    if report_path is not None:
        fields.append("report_path = ?"); vals.append(report_path)
    vals.append(job_id)
    with get_conn() as conn:
        conn.execute(
            f"UPDATE jobs SET {', '.join(fields)} WHERE id = ?", vals
        )


def update_stage(job_id: str, stage: str, status: str,
                 detail: str = "") -> None:
    """Updates a pipeline stage entry, recording a timestamp on each change."""
    now = datetime.now(timezone.utc).isoformat()
    with get_conn() as conn:
        row = conn.execute(
            "SELECT stages FROM jobs WHERE id = ?", (job_id,)
        ).fetchone()
        if not row:
            return
        stages = json.loads(row["stages"])
        prev = stages.get(stage, {})
        stages[stage] = {
            "status":     status,
            "detail":     detail,
            "started_at": prev.get("started_at", now) if status == "running" else prev.get("started_at"),
            "updated_at": now,
        }
        conn.execute(
            "UPDATE jobs SET stages = ? WHERE id = ?",
            (json.dumps(stages), job_id)
        )


def get_job(job_id: str) -> dict | None:
    with get_conn() as conn:
        row = conn.execute(
            "SELECT * FROM jobs WHERE id = ?", (job_id,)
        ).fetchone()
    return dict(row) if row else None


def get_expired_jobs() -> list[dict]:
    """Returns all jobs whose expiry time has passed."""
    now = datetime.now(timezone.utc).isoformat()
    with get_conn() as conn:
        rows = conn.execute(
            "SELECT * FROM jobs WHERE expires_at < ? AND status != 'deleted'",
            (now,)
        ).fetchall()
    return [dict(r) for r in rows]


def get_recent_jobs(limit: int = 20) -> list[dict]:
    """Returns the most recent jobs ordered by creation time, excluding deleted."""
    with get_conn() as conn:
        rows = conn.execute(
            """SELECT id, status, filename, read_type, created_at, expires_at,
                      error, report_path, stages
               FROM jobs
               WHERE status != 'deleted'
               ORDER BY created_at DESC
               LIMIT ?""",
            (limit,)
        ).fetchall()
    return [dict(r) for r in rows]


def get_overflow_jobs(keep: int = 10) -> list[dict]:
    """
    Returns finished jobs (completed/failed) beyond the most recent `keep`.
    These are candidates for pruning — upload/result files are deleted but
    the DB row and report HTML are preserved.
    """
    with get_conn() as conn:
        rows = conn.execute(
            """SELECT * FROM jobs
               WHERE status IN ('completed', 'failed')
                 AND files_deleted = 0
               ORDER BY created_at DESC
               LIMIT -1 OFFSET ?""",
            (keep,)
        ).fetchall()
    return [dict(r) for r in rows]


def cancel_stale_jobs() -> int:
    """
    Called once at server startup.  Any job still in an in-progress state
    (queued / running / ai_pending) was interrupted by a crash or restart and
    will never complete.  Mark them cancelled so they don't block new uploads.
    Returns the number of rows updated.
    """
    with get_conn() as conn:
        cur = conn.execute(
            """UPDATE jobs
               SET status = 'cancelled',
                   error  = 'Server restarted while job was in progress.'
               WHERE status IN ('queued', 'running', 'ai_pending')"""
        )
        n = cur.rowcount
    if n:
        logger.warning("Cancelled %d stale job(s) left over from previous run.", n)
    return n


def mark_files_deleted(job_id: str) -> None:
    """Marks a job's working files as deleted (report is kept)."""
    with get_conn() as conn:
        conn.execute(
            "UPDATE jobs SET files_deleted = 1 WHERE id = ?", (job_id,)
        )


def mark_deleted(job_id: str) -> None:
    with get_conn() as conn:
        conn.execute(
            "UPDATE jobs SET status = 'deleted' WHERE id = ?", (job_id,)
        )


def cancel_job(job_id: str) -> bool:
    """
    Marks a job as cancelled.  Returns True if the job was active, False otherwise.
    """
    with get_conn() as conn:
        cur = conn.execute(
            """UPDATE jobs SET status = 'cancelled', error = 'Cancelled by user.'
               WHERE id = ? AND status IN ('queued', 'running', 'ai_pending')""",
            (job_id,)
        )
        return cur.rowcount > 0


def get_job_status(job_id: str) -> str | None:
    """Lightweight single-field fetch — avoids loading stages JSON for cancel checks."""
    with get_conn() as conn:
        row = conn.execute(
            "SELECT status FROM jobs WHERE id = ?", (job_id,)
        ).fetchone()
    return row["status"] if row else None


def count_active_jobs_for_ip(ip_hash: str) -> int:
    """Returns the number of active jobs for a given IP hash."""
    with get_conn() as conn:
        row = conn.execute(
            """SELECT COUNT(*) as n FROM jobs
               WHERE ip_hash = ? AND status NOT IN ('completed','failed','deleted','cancelled')""",
            (ip_hash,)
        ).fetchone()
    return row["n"] if row else 0
