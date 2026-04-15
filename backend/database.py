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
                id          TEXT PRIMARY KEY,
                status      TEXT NOT NULL DEFAULT 'queued',
                created_at  TEXT NOT NULL,
                expires_at  TEXT NOT NULL,
                filename    TEXT NOT NULL,
                file_md5    TEXT,
                ip_hash     TEXT NOT NULL,
                read_type   TEXT,
                stages      TEXT DEFAULT '{}',
                error       TEXT,
                report_path TEXT
            );

            CREATE INDEX IF NOT EXISTS idx_expires ON jobs(expires_at);
            CREATE INDEX IF NOT EXISTS idx_ip      ON jobs(ip_hash);
        """)
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
    """Updates a pipeline stage entry."""
    with get_conn() as conn:
        row = conn.execute(
            "SELECT stages FROM jobs WHERE id = ?", (job_id,)
        ).fetchone()
        if not row:
            return
        stages = json.loads(row["stages"])
        stages[stage] = {"status": status, "detail": detail}
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


def mark_deleted(job_id: str) -> None:
    with get_conn() as conn:
        conn.execute(
            "UPDATE jobs SET status = 'deleted' WHERE id = ?", (job_id,)
        )


def count_active_jobs_for_ip(ip_hash: str) -> int:
    """Returns the number of active jobs for a given IP hash."""
    with get_conn() as conn:
        row = conn.execute(
            """SELECT COUNT(*) as n FROM jobs
               WHERE ip_hash = ? AND status NOT IN ('completed','failed','deleted')""",
            (ip_hash,)
        ).fetchone()
    return row["n"] if row else 0
