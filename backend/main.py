"""
LycianWay — Main FastAPI Application
Secure FASTQ upload, pipeline triggering, status polling, report download.
"""
import hashlib
import logging
import sys
import threading
from datetime import datetime, timedelta, timezone
from pathlib import Path

import aiofiles
from fastapi import FastAPI, File, Form, HTTPException, Request, UploadFile
from typing import Optional
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from slowapi import Limiter
from slowapi.errors import RateLimitExceeded
from slowapi.util import get_remote_address

sys.path.insert(0, str(Path(__file__).parent))
import cleanup
import database as db
import system_check as sc
from ai_interpreter import interpret
from config import (
    JOB_EXPIRY_HOURS, KRAKEN2_DEFAULT_DB, MAX_FILE_SIZE_MB,
    RATE_LIMIT, RESULTS_DIR, UPLOAD_DIR,
)
from pipeline import run_pipeline, get_cancel_event, PipelineCancelledError
from comparison_pipeline import run_comparison
from report_generator import generate_html_report, save_report
from security import (
    compute_md5, generate_job_id, job_upload_dir,
    validate_file_size, validate_filename, validate_magic_bytes,
)

# ─── Logging ──────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(name)s] %(levelname)s: %(message)s",
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler(str(Path(__file__).parent.parent / "data" / "logs" / "app.log")),
    ],
)
logger = logging.getLogger("main")

# ─── FastAPI Setup ────────────────────────────────────────────────────────────
limiter = Limiter(key_func=get_remote_address)

app = FastAPI(
    title="LycianWay",
    description="Automated bacterial genomic analysis platform",
    version="0.1.0",
    docs_url=None,
    redoc_url=None,
)

app.state.limiter = limiter

@app.exception_handler(RateLimitExceeded)
async def rate_limit_handler(request: Request, exc: RateLimitExceeded):
    return JSONResponse(
        status_code=429,
        content={"error": "Too many requests. Please wait before trying again."},
    )

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["GET", "POST"],
    allow_headers=["Content-Type"],
)

FRONTEND_DIR = Path(__file__).parent.parent / "frontend"
app.mount("/static", StaticFiles(directory=str(FRONTEND_DIR / "static")), name="static")
templates = Jinja2Templates(directory=str(FRONTEND_DIR / "templates"))

# ─── Startup ──────────────────────────────────────────────────────────────────
@app.on_event("startup")
async def startup():
    db.init_db()
    stale = db.cancel_stale_jobs()   # clear jobs orphaned by a crash / restart
    cleanup.start_cleanup_thread()
    logger.info("LycianWay started. (stale jobs cancelled: %d)", stale)


# ─── Helpers ──────────────────────────────────────────────────────────────────

def _ip_hash(request: Request) -> str:
    """Hash the client IP — raw IP is never stored (GDPR)."""
    ip = get_remote_address(request)
    return hashlib.sha256(ip.encode()).hexdigest()[:16]


def _resolve_kraken2_db(raw: str) -> Path | None:
    """
    Validates and returns the Kraken2 database path.
    Accepts user-supplied path or falls back to the configured default.
    Returns None if no valid DB is found.
    """
    candidate = Path(raw.strip()) if raw and raw.strip() else KRAKEN2_DEFAULT_DB
    if candidate.is_dir() and (candidate / "hash.k2d").exists():
        return candidate
    # Try default
    if KRAKEN2_DEFAULT_DB.is_dir() and (KRAKEN2_DEFAULT_DB / "hash.k2d").exists():
        return KRAKEN2_DEFAULT_DB
    return None


def _run_full_pipeline(job_id: str, fastq_path: Path,
                       fastq_r2: Path | None = None,
                       kraken2_db: Path | None = None,
                       sample_type: str = "auto") -> None:
    """
    Runs the full pipeline + AI interpretation in a background thread.
    Even on failure a partial report is generated so the user can see
    what completed before the error.
    """
    logger.info("[%s] Pipeline starting: %s%s", job_id, fastq_path,
                f" + {fastq_r2.name}" if fastq_r2 else "")
    out_dir = RESULTS_DIR / job_id
    out_dir.mkdir(parents=True, exist_ok=True)

    pipeline_results: dict = {}
    ai_results:       dict = {}
    job_error:        str | None = None

    # ── 1. Genomic pipeline ───────────────────────────────────────────────────
    try:
        pipeline_results = run_pipeline(job_id, fastq_path,
                                        fastq_r2=fastq_r2,
                                        kraken2_db=kraken2_db,
                                        sample_type=sample_type)
        if pipeline_results.get("error"):
            job_error = pipeline_results["error"]
        # If pipeline was cancelled, skip AI and report generation
        if db.get_job_status(job_id) == "cancelled":
            logger.info("[%s] Job cancelled — skipping AI and report.", job_id)
            return
    except PipelineCancelledError:
        logger.info("[%s] Pipeline cancelled.", job_id)
        return
    except Exception as exc:
        logger.error("[%s] Pipeline error: %s", job_id, exc, exc_info=True)
        job_error = str(exc)

    # ── 2. AI interpretation (skip only if pipeline produced nothing) ─────────
    if pipeline_results and not pipeline_results.get("error"):
        try:
            db.update_stage(job_id, "ai", "running")
            ai_results = interpret(pipeline_results)
            db.update_stage(job_id, "ai", "done")
        except Exception as exc:
            logger.error("[%s] AI error: %s", job_id, exc, exc_info=True)
            ai_results = {
                "summary": f"AI interpretation failed: {exc}",
                "risk_level": "UNKNOWN",
            }
            db.update_stage(job_id, "ai", "failed", str(exc))
            # Not a fatal error — continue to report

    # ── 3. Generate report (always, even on partial/failed pipelines) ─────────
    try:
        db.update_stage(job_id, "report", "running")
        html = generate_html_report(job_id, pipeline_results, ai_results,
                                    job_error=job_error)
        report_path = save_report(job_id, html, out_dir)
        db.update_stage(job_id, "report", "done")

        final_status = "failed" if job_error else "completed"
        db.update_job_status(job_id, final_status,
                             report_path=str(report_path),
                             error=job_error)
        logger.info("[%s] %s. Report: %s", job_id,
                    "Failed (partial report)" if job_error else "Completed",
                    report_path)

    except Exception as exc:
        logger.error("[%s] Report generation error: %s", job_id, exc, exc_info=True)
        db.update_job_status(job_id, "failed", error=str(exc))


# ─── Routes ───────────────────────────────────────────────────────────────────

@app.get("/", response_class=HTMLResponse)
async def index(request: Request):
    return templates.TemplateResponse("index.html", {"request": request})


@app.get("/kraken2-databases")
async def list_kraken2_databases():
    """
    Scans common locations for valid Kraken2 databases (contain hash.k2d).
    Returns a list of {label, path} objects for the UI dropdown.
    """
    search_roots = [
        Path("/home/analysis/databases_all"),
        Path("/home/analysis/databases"),
        Path("/databases"),        # Docker volume mount point
        Path("/data/kraken2"),
        Path("/opt/databases"),
    ]
    found = []
    seen = set()
    for root in search_roots:
        if not root.exists():
            continue
        for p in root.rglob("hash.k2d"):
            db_dir = p.parent
            key = str(db_dir.resolve())
            if key in seen:
                continue
            seen.add(key)
            found.append({"label": db_dir.name, "path": str(db_dir)})

    # Always include the configured default if not already found
    default_key = str(KRAKEN2_DEFAULT_DB.resolve())
    if default_key not in seen and KRAKEN2_DEFAULT_DB.is_dir():
        found.insert(0, {"label": KRAKEN2_DEFAULT_DB.name + " (default)",
                         "path": str(KRAKEN2_DEFAULT_DB)})

    return JSONResponse({"databases": found})


async def _save_upload(upload: UploadFile, dest: Path) -> int:
    """Streams an uploaded file to disk; returns total bytes written."""
    total = 0
    first_chunk = await upload.read(8192)
    validate_magic_bytes(first_chunk, dest.name)
    total += len(first_chunk)
    async with aiofiles.open(dest, "wb") as out:
        await out.write(first_chunk)
        while chunk := await upload.read(1 << 20):
            total += len(chunk)
            if total > MAX_FILE_SIZE_MB * 1024 * 1024:
                dest.unlink(missing_ok=True)
                raise HTTPException(
                    status_code=413,
                    detail=f"File exceeds the {MAX_FILE_SIZE_MB} MB limit."
                )
            await out.write(chunk)
    return total


@app.post("/upload")
@limiter.limit(RATE_LIMIT)
async def upload_fastq(
    request: Request,
    file:        UploadFile = File(...),
    file_r2:     Optional[UploadFile] = File(default=None),
    kraken2_db:  str = Form(default=""),
    sample_type: str = Form(default="auto"),
):
    """
    Securely uploads FASTQ/FASTQ.gz file(s) and starts the analysis pipeline.
    Accepts a single file (MinION / SE Illumina) or two files (R1 + R2, paired-end).

    Security checks:
    - Filename sanitisation
    - Magic byte validation
    - Size limit (2 GB per file)
    - Concurrent job limit per IP
    """
    ip_hash = _ip_hash(request)
    if db.count_active_jobs_for_ip(ip_hash) >= 3:
        raise HTTPException(
            status_code=429,
            detail="You can have at most 3 active analyses running at the same time."
        )

    # Validate and generate job
    try:
        safe_name = validate_filename(file.filename or "upload.fastq.gz")
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))

    job_id = generate_job_id()
    up_dir = job_upload_dir(job_id)
    up_dir.mkdir(parents=True, exist_ok=True)
    dest = up_dir / safe_name

    # Save R1
    try:
        total_bytes = await _save_upload(file, dest)
    except HTTPException:
        raise
    except Exception as e:
        dest.unlink(missing_ok=True)
        logger.error("File save error (R1): %s", e)
        raise HTTPException(status_code=500, detail="Failed to save file.")

    # Save R2 (optional)
    dest_r2: Path | None = None
    if file_r2 and file_r2.filename:
        try:
            safe_name_r2 = validate_filename(file_r2.filename)
        except ValueError as e:
            dest.unlink(missing_ok=True)
            raise HTTPException(status_code=400, detail=f"R2: {e}")
        dest_r2 = up_dir / safe_name_r2
        try:
            r2_bytes = await _save_upload(file_r2, dest_r2)
            total_bytes += r2_bytes
        except HTTPException:
            dest.unlink(missing_ok=True)
            raise
        except Exception as e:
            dest.unlink(missing_ok=True)
            dest_r2.unlink(missing_ok=True)
            logger.error("File save error (R2): %s", e)
            raise HTTPException(status_code=500, detail="Failed to save R2 file.")

    md5 = compute_md5(dest)
    expires_at = (
        datetime.now(timezone.utc) + timedelta(hours=JOB_EXPIRY_HOURS)
    ).isoformat()
    db.create_job(job_id, safe_name, expires_at, ip_hash, md5)

    resolved_db = _resolve_kraken2_db(kraken2_db)
    logger.info("New job: %s | %s%s | %d bytes | MD5: %s | Kraken2 DB: %s",
                job_id, safe_name,
                f" + {dest_r2.name}" if dest_r2 else "",
                total_bytes, md5,
                str(resolved_db) if resolved_db else "none")

    t = threading.Thread(
        target=_run_full_pipeline,
        args=(job_id, dest, dest_r2, resolved_db, sample_type),
        daemon=True,
        name=f"pipeline-{job_id[:8]}",
    )
    t.start()

    return JSONResponse({
        "job_id":       job_id,
        "filename":     safe_name,
        "paired_end":   dest_r2 is not None,
        "size_mb":      round(total_bytes / 1e6, 2),
        "md5":          md5,
        "expires_at":   expires_at,
        "kraken2_db":   str(resolved_db) if resolved_db else None,
        "message":      "Analysis started. Poll /status/{job_id} for progress.",
    })


@app.get("/status/{job_id}")
async def job_status(request: Request, job_id: str):
    """Returns the current status and stage details for a job."""
    import re
    if not re.fullmatch(r"[0-9a-f]{32}", job_id):
        raise HTTPException(status_code=400, detail="Invalid job ID.")

    job = db.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found.")

    import json as _json
    stages = _json.loads(job["stages"]) if job["stages"] else {}

    return JSONResponse({
        "job_id":     job["id"],
        "status":     job["status"],
        "filename":   job["filename"],
        "read_type":  job["read_type"],
        "created_at": job["created_at"],
        "expires_at": job["expires_at"],
        "stages":     stages,
        "error":      job["error"],
        "has_report": bool(job["report_path"]),
    })


@app.post("/cancel/{job_id}")
async def cancel_job_endpoint(request: Request, job_id: str):
    """Cancels a running analysis. Interrupts the active subprocess within ~5 s."""
    import re
    if not re.fullmatch(r"[0-9a-f]{32}", job_id):
        raise HTTPException(status_code=400, detail="Invalid job ID.")

    job = db.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found.")

    if job["status"] not in ("running", "queued", "ai_pending"):
        raise HTTPException(
            status_code=409,
            detail=f"Job is not active (status: {job['status']})."
        )

    # Signal the pipeline thread to stop (interrupts current subprocess within 5 s)
    ev = get_cancel_event(job_id)
    if ev:
        ev.set()

    # Mark DB immediately so UI sees "Cancelled" on next poll
    db.cancel_job(job_id)
    logger.info("[%s] Cancelled by user.", job_id)
    return JSONResponse({"status": "cancelled", "job_id": job_id})


def _get_report_path(job_id: str) -> Path:
    """Shared validation logic for report endpoints."""
    import re
    if not re.fullmatch(r"[0-9a-f]{32}", job_id):
        raise HTTPException(status_code=400, detail="Invalid job ID.")
    job = db.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found.")
    if job["status"] not in ("completed", "failed", "cancelled"):
        raise HTTPException(
            status_code=202,
            detail=f"Analysis not yet complete (status: {job['status']})."
        )
    if not job["report_path"]:
        raise HTTPException(status_code=404, detail="Report not found.")
    report_path = Path(job["report_path"])
    if not report_path.exists():
        raise HTTPException(status_code=410, detail="Report has been deleted.")
    return report_path


@app.get("/view/{job_id}", response_class=HTMLResponse)
async def view_report(request: Request, job_id: str):
    """Opens the completed HTML report inline in the browser."""
    report_path = _get_report_path(job_id)
    return HTMLResponse(content=report_path.read_text(encoding="utf-8"))


@app.get("/report/{job_id}")
async def download_report(request: Request, job_id: str):
    """Downloads the completed HTML report as a file."""
    report_path = _get_report_path(job_id)
    return FileResponse(
        path=str(report_path),
        media_type="text/html",
        filename=f"LycianWay_{job_id[:8]}_report.html",
    )


@app.get("/assembly-graph/{job_id}")
async def assembly_graph(request: Request, job_id: str):
    """Serves the Bandage assembly graph PNG for a completed job."""
    import re
    if not re.fullmatch(r"[0-9a-f]{32}", job_id):
        raise HTTPException(status_code=400, detail="Invalid job ID.")
    job = db.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found.")

    import json as _json
    stages = _json.loads(job["stages"]) if job["stages"] else {}
    bandage = stages.get("bandage", {})
    detail = _json.loads(bandage.get("detail", "{}")) if bandage.get("detail") else {}
    image_path = detail.get("image_path")

    if not image_path:
        raise HTTPException(status_code=404, detail="Assembly graph image not available.")
    img = Path(image_path)
    if not img.exists():
        raise HTTPException(status_code=410, detail="Assembly graph image has been deleted.")
    return FileResponse(path=str(img), media_type="image/png",
                        filename=f"assembly_graph_{job_id[:8]}.png")


@app.get("/history")
async def analysis_history():
    """Returns the most recent analyses for the history panel."""
    import json as _json
    jobs = db.get_recent_jobs(limit=20)
    result = []
    for j in jobs:
        stages = _json.loads(j["stages"]) if j["stages"] else {}
        # Calculate elapsed time from first running stage to last updated
        started_at = None
        last_updated = None
        for s in stages.values():
            sa = s.get("started_at")
            ua = s.get("updated_at")
            if sa and (started_at is None or sa < started_at):
                started_at = sa
            if ua and (last_updated is None or ua > last_updated):
                last_updated = ua

        elapsed_sec = None
        if started_at and last_updated:
            try:
                from datetime import datetime, timezone
                fmt = "%Y-%m-%dT%H:%M:%S.%f%z"
                def _parse(s):
                    try:
                        return datetime.fromisoformat(s)
                    except Exception:
                        return None
                t0 = _parse(started_at)
                t1 = _parse(last_updated)
                if t0 and t1:
                    elapsed_sec = int((t1 - t0).total_seconds())
            except Exception:
                pass

        result.append({
            "job_id":       j["id"],
            "status":       j["status"],
            "filename":     j["filename"],
            "read_type":    j["read_type"],
            "created_at":   j["created_at"],
            "elapsed_sec":  elapsed_sec,
            "has_report":   bool(j["report_path"]) and (j["files_deleted"] == 0 or
                             (j["report_path"] and Path(j["report_path"]).exists())),
            "error":        j["error"],
        })
    return JSONResponse({"analyses": result})


@app.get("/system-stats")
async def system_stats():
    """Returns CPU and disk usage for the server status widget."""
    import psutil
    cpu_percent = psutil.cpu_percent(interval=0.2)
    cpu_count   = psutil.cpu_count(logical=True)
    disk = psutil.disk_usage("/")
    return JSONResponse({
        "cpu": {
            "percent":     cpu_percent,
            "total_cores": cpu_count,
            "budget":      __import__("config").MAX_THREADS,
        },
        "disk": {
            "total_gb": round(disk.total / 1e9, 1),
            "used_gb":  round(disk.used  / 1e9, 1),
            "free_gb":  round(disk.free  / 1e9, 1),
            "percent":  disk.percent,
        },
    })


@app.get("/system-check")
async def system_check():
    """Scans installed tools and databases; reports what is present or missing."""
    return JSONResponse(sc.run_system_check())


@app.get("/health")
async def health():
    return {"status": "ok"}


# ─── Multi-sample comparison ───────────────────────────────────────────────────

@app.post("/compare")
@limiter.limit("10/minute")
async def start_comparison(request: Request):
    """
    Accepts a JSON body: {"job_ids": ["id1", "id2", ...]}
    Validates all jobs, then launches the comparison in a background thread.
    """
    body = await request.json()
    job_ids = body.get("job_ids", [])

    if not isinstance(job_ids, list) or len(job_ids) < 2:
        raise HTTPException(status_code=400, detail="Provide at least 2 job IDs.")
    if len(job_ids) > 20:
        raise HTTPException(status_code=400, detail="Maximum 20 samples per comparison.")

    # Validate all jobs exist and are completed
    valid_ids = []
    for jid in job_ids:
        job = db.get_job(jid)
        if not job:
            raise HTTPException(status_code=404, detail=f"Job {jid} not found.")
        if job["status"] != "completed":
            raise HTTPException(
                status_code=400,
                detail=f"Job {jid} is not completed (status: {job['status']})."
            )
        valid_ids.append(jid)

    ip_hash = _ip_hash(request)
    comp_id = generate_job_id()
    db.create_comparison(comp_id, valid_ids, ip_hash)

    def _bg():
        run_comparison(comp_id, valid_ids)

    threading.Thread(target=_bg, daemon=True).start()
    logger.info("Comparison %s started for %d jobs", comp_id, len(valid_ids))
    return {"comparison_id": comp_id}


@app.get("/compare/{comp_id}/status")
async def comparison_status(comp_id: str):
    """Poll comparison status."""
    comp = db.get_comparison(comp_id)
    if not comp:
        raise HTTPException(status_code=404, detail="Comparison not found.")
    return {
        "comparison_id": comp_id,
        "status":        comp["status"],
        "error":         comp.get("error"),
        "has_report":    bool(comp.get("report_path") and
                              Path(comp["report_path"]).exists()),
    }


@app.get("/compare/{comp_id}/report", response_class=HTMLResponse)
async def comparison_report(comp_id: str):
    """View the comparison HTML report."""
    comp = db.get_comparison(comp_id)
    if not comp:
        raise HTTPException(status_code=404, detail="Comparison not found.")
    if comp["status"] != "completed":
        raise HTTPException(status_code=400,
                            detail=f"Comparison not ready (status: {comp['status']}).")
    rp = Path(comp["report_path"])
    if not rp.exists():
        raise HTTPException(status_code=404, detail="Report file not found.")
    return HTMLResponse(content=rp.read_text(encoding="utf-8"))
