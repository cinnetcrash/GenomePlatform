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
from fastapi import FastAPI, File, HTTPException, Request, UploadFile
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
from ai_interpreter import interpret
from config import (
    JOB_EXPIRY_HOURS, MAX_FILE_SIZE_MB,
    RATE_LIMIT, RESULTS_DIR, UPLOAD_DIR,
)
from pipeline import run_pipeline
from primer_designer import design_all_primers
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
    description="Automated bacterial genomic analysis and PCR primer design platform",
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
    cleanup.start_cleanup_thread()
    logger.info("LycianWay started.")


# ─── Helpers ──────────────────────────────────────────────────────────────────

def _ip_hash(request: Request) -> str:
    """Hash the client IP — raw IP is never stored (GDPR)."""
    ip = get_remote_address(request)
    return hashlib.sha256(ip.encode()).hexdigest()[:16]


def _run_full_pipeline(job_id: str, fastq_path: Path) -> None:
    """Runs the full pipeline + AI interpretation + primer design in a background thread."""
    try:
        logger.info("[%s] Pipeline starting: %s", job_id, fastq_path)
        out_dir = RESULTS_DIR / job_id
        out_dir.mkdir(parents=True, exist_ok=True)

        # 1. Genomic pipeline
        pipeline_results = run_pipeline(job_id, fastq_path)
        if pipeline_results.get("error"):
            return

        # 2. AI interpretation
        db.update_stage(job_id, "ai", "running")
        ai_results = interpret(pipeline_results)
        db.update_stage(job_id, "ai", "done")

        # 3. Primer design — skipped if assembly failed or no AMR/AI targets found
        assembly_fasta = pipeline_results.get("assembly_fasta")
        fasta_path = Path(assembly_fasta) if assembly_fasta else None
        amr_count = pipeline_results.get("amr", {}).get("count", 0)
        ai_targets = ai_results.get("pcr_targets", [])
        has_targets = amr_count > 0 or len(ai_targets) > 0
        assembly_ok = fasta_path and fasta_path.exists()

        if not assembly_ok:
            db.update_stage(job_id, "primers", "skipped",
                            "Assembly did not produce a FASTA — primer design skipped.")
            primer_results = []
            logger.warning("[%s] Primer design skipped: no assembly FASTA.", job_id)
        elif not has_targets:
            db.update_stage(job_id, "primers", "skipped",
                            "No AMR genes or AI PCR targets found — primer design skipped.")
            primer_results = []
            logger.info("[%s] Primer design skipped: no targets.", job_id)
        else:
            db.update_stage(job_id, "primers", "running")
            primer_results = design_all_primers(
                pipeline_results.get("amr", {}),
                ai_targets,
                fasta_path,
            )
            db.update_stage(job_id, "primers", "done", f"{len(primer_results)} genes")

        # 4. Generate report
        db.update_stage(job_id, "report", "running")
        html = generate_html_report(job_id, pipeline_results, ai_results, primer_results)
        report_path = save_report(job_id, html, out_dir)
        db.update_stage(job_id, "report", "done")

        db.update_job_status(job_id, "completed", report_path=str(report_path))
        logger.info("[%s] Completed. Report: %s", job_id, report_path)

    except Exception as exc:
        logger.error("[%s] Background error: %s", job_id, exc, exc_info=True)
        db.update_job_status(job_id, "failed", error=str(exc))


# ─── Routes ───────────────────────────────────────────────────────────────────

@app.get("/", response_class=HTMLResponse)
async def index(request: Request):
    return templates.TemplateResponse("index.html", {"request": request})


@app.post("/upload")
@limiter.limit(RATE_LIMIT)
async def upload_fastq(request: Request, file: UploadFile = File(...)):
    """
    Securely uploads a FASTQ/FASTQ.gz file and starts the analysis pipeline.

    Security checks:
    - Filename sanitisation
    - Magic byte validation
    - Size limit (2 GB)
    - Concurrent job limit per IP
    """
    ip_hash = _ip_hash(request)

    active = db.count_active_jobs_for_ip(ip_hash)
    if active >= 3:
        raise HTTPException(
            status_code=429,
            detail="You can have at most 3 active analyses running at the same time."
        )

    try:
        safe_name = validate_filename(file.filename or "upload.fastq.gz")
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))

    first_chunk = await file.read(8192)
    try:
        validate_magic_bytes(first_chunk, safe_name)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))

    job_id  = generate_job_id()
    up_dir  = job_upload_dir(job_id)
    up_dir.mkdir(parents=True, exist_ok=True)
    dest    = up_dir / safe_name

    total_bytes = len(first_chunk)
    try:
        async with aiofiles.open(dest, "wb") as out:
            await out.write(first_chunk)
            while chunk := await file.read(1 << 20):
                total_bytes += len(chunk)
                if total_bytes > MAX_FILE_SIZE_MB * 1024 * 1024:
                    dest.unlink(missing_ok=True)
                    raise HTTPException(
                        status_code=413,
                        detail=f"File exceeds the {MAX_FILE_SIZE_MB} MB limit."
                    )
                await out.write(chunk)
    except HTTPException:
        raise
    except Exception as e:
        dest.unlink(missing_ok=True)
        logger.error("File save error: %s", e)
        raise HTTPException(status_code=500, detail="Failed to save file.")

    md5 = compute_md5(dest)

    expires_at = (
        datetime.now(timezone.utc) + timedelta(hours=JOB_EXPIRY_HOURS)
    ).isoformat()
    db.create_job(job_id, safe_name, expires_at, ip_hash, md5)

    logger.info("New job: %s | File: %s | %d bytes | MD5: %s",
                job_id, safe_name, total_bytes, md5)

    t = threading.Thread(
        target=_run_full_pipeline,
        args=(job_id, dest),
        daemon=True,
        name=f"pipeline-{job_id[:8]}",
    )
    t.start()

    return JSONResponse({
        "job_id":     job_id,
        "filename":   safe_name,
        "size_mb":    round(total_bytes / 1e6, 2),
        "md5":        md5,
        "expires_at": expires_at,
        "message":    "Analysis started. Poll /status/{job_id} for progress.",
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


def _get_report_path(job_id: str) -> Path:
    """Shared validation logic for report endpoints."""
    import re
    if not re.fullmatch(r"[0-9a-f]{32}", job_id):
        raise HTTPException(status_code=400, detail="Invalid job ID.")
    job = db.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found.")
    if job["status"] != "completed":
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


@app.get("/system-stats")
async def system_stats():
    """Returns CPU and disk usage for the server status widget."""
    import psutil
    from config import MAX_THREADS, RESULTS_DIR, UPLOAD_DIR

    cpu_percent  = psutil.cpu_percent(interval=0.2)
    cpu_count    = psutil.cpu_count(logical=True)

    disk = psutil.disk_usage("/")
    data_dir = RESULTS_DIR.parent
    try:
        data_disk = psutil.disk_usage(str(data_dir))
        data_used_gb  = round(data_disk.used / 1e9, 1)
        data_total_gb = round(data_disk.total / 1e9, 1)
        data_percent  = data_disk.percent
    except Exception:
        data_used_gb = data_total_gb = data_percent = None

    return JSONResponse({
        "cpu": {
            "percent":     cpu_percent,
            "total_cores": cpu_count,
            "budget":      MAX_THREADS,
        },
        "disk": {
            "total_gb":  round(disk.total / 1e9, 1),
            "used_gb":   round(disk.used / 1e9, 1),
            "free_gb":   round(disk.free / 1e9, 1),
            "percent":   disk.percent,
        },
    })


@app.get("/health")
async def health():
    return {"status": "ok"}
