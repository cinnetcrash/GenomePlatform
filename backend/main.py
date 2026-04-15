"""
GenomePlatform — Ana FastAPI Uygulaması
Güvenli FASTQ yükleme, pipeline tetikleme, durum sorgulama, rapor indirme.
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

# Yerel modüller
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

# ─── Logging ─────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(name)s] %(levelname)s: %(message)s",
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler(str(Path(__file__).parent.parent / "data" / "logs" / "app.log")),
    ],
)
logger = logging.getLogger("main")

# ─── FastAPI Kurulumu ─────────────────────────────────────────────────────────
limiter = Limiter(key_func=get_remote_address)

app = FastAPI(
    title="GenomePlatform",
    description="Otomatik genomik analiz ve PCR primer tasarım platformu",
    version="0.1.0",
    docs_url=None,      # Prodüksiyon'da Swagger UI'yi kapat
    redoc_url=None,
)

app.state.limiter = limiter

@app.exception_handler(RateLimitExceeded)
async def rate_limit_handler(request: Request, exc: RateLimitExceeded):
    return JSONResponse(
        status_code=429,
        content={"error": "Çok fazla istek gönderdiniz. Lütfen bekleyin."},
    )

# CORS — sadece aynı origin (prodüksiyonda domain'e göre güncelle)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],   # Prototip; prodüksiyonda kısıtla
    allow_methods=["GET", "POST"],
    allow_headers=["Content-Type"],
)

# Statik dosyalar ve şablonlar
FRONTEND_DIR = Path(__file__).parent.parent / "frontend"
app.mount("/static", StaticFiles(directory=str(FRONTEND_DIR / "static")), name="static")
templates = Jinja2Templates(directory=str(FRONTEND_DIR / "templates"))

# ─── Başlatma ────────────────────────────────────────────────────────────────
@app.on_event("startup")
async def startup():
    db.init_db()
    cleanup.start_cleanup_thread()
    logger.info("GenomePlatform başlatıldı.")


# ─── Yardımcılar ─────────────────────────────────────────────────────────────

def _ip_hash(request: Request) -> str:
    """IP'yi hash'le — ham IP'yi kaydetme (GDPR)."""
    ip = get_remote_address(request)
    return hashlib.sha256(ip.encode()).hexdigest()[:16]


def _run_full_pipeline(job_id: str, fastq_path: Path) -> None:
    """Arka planda tam pipeline + AI + Primer çalıştırır."""
    try:
        logger.info("[%s] Pipeline başlıyor: %s", job_id, fastq_path)
        out_dir = RESULTS_DIR / job_id
        out_dir.mkdir(parents=True, exist_ok=True)

        # 1. Genomik pipeline
        pipeline_results = run_pipeline(job_id, fastq_path)
        if pipeline_results.get("error"):
            return  # Hata zaten DB'ye yazıldı

        # 2. AI yorumu
        db.update_stage(job_id, "ai", "running")
        ai_results = interpret(pipeline_results)
        db.update_stage(job_id, "ai", "done")

        # 3. Primer tasarımı
        db.update_stage(job_id, "primers", "running")
        assembly_fasta = pipeline_results.get("assembly_fasta")
        fasta_path = Path(assembly_fasta) if assembly_fasta else None
        primer_results = design_all_primers(
            pipeline_results.get("amr", {}),
            ai_results.get("pcr_targets", []),
            fasta_path,
        )
        db.update_stage(job_id, "primers", "done",
                        f"{len(primer_results)} gen")

        # 4. Rapor üret
        db.update_stage(job_id, "report", "running")
        html = generate_html_report(
            job_id, pipeline_results, ai_results, primer_results
        )
        report_path = save_report(job_id, html, out_dir)
        db.update_stage(job_id, "report", "done")

        db.update_job_status(job_id, "completed",
                             report_path=str(report_path))
        logger.info("[%s] Tamamlandı. Rapor: %s", job_id, report_path)

    except Exception as exc:
        logger.error("[%s] Arka plan hatası: %s", job_id, exc, exc_info=True)
        db.update_job_status(job_id, "failed", error=str(exc))


# ─── Rotalar ─────────────────────────────────────────────────────────────────

@app.get("/", response_class=HTMLResponse)
async def index(request: Request):
    return templates.TemplateResponse("index.html", {"request": request})


@app.post("/upload")
@limiter.limit(RATE_LIMIT)
async def upload_fastq(request: Request,
                       file: UploadFile = File(...)):
    """
    FASTQ/FASTQ.gz dosyasını güvenli biçimde yükler ve analizi başlatır.

    Güvenlik kontrolleri:
    - Dosya adı temizleme
    - Magic byte doğrulama
    - Boyut limiti (2 GB)
    - Aynı IP'den eş zamanlı iş limiti
    """
    ip_hash = _ip_hash(request)

    # IP başına aktif iş kontrolü
    active = db.count_active_jobs_for_ip(ip_hash)
    if active >= 3:
        raise HTTPException(
            status_code=429,
            detail="Aynı anda en fazla 3 aktif analiziniz olabilir."
        )

    # Dosya adı validasyonu
    try:
        safe_name = validate_filename(file.filename or "upload.fastq.gz")
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))

    # İlk bloğu oku (magic byte + boyut kontrolü için)
    first_chunk = await file.read(8192)
    try:
        validate_magic_bytes(first_chunk, safe_name)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))

    # Job oluştur
    job_id  = generate_job_id()
    up_dir  = job_upload_dir(job_id)
    up_dir.mkdir(parents=True, exist_ok=True)
    dest    = up_dir / safe_name

    # Dosyayı kaydet
    total_bytes = len(first_chunk)
    try:
        async with aiofiles.open(dest, "wb") as out:
            await out.write(first_chunk)
            while chunk := await file.read(1 << 20):   # 1 MB parçalar
                total_bytes += len(chunk)
                if total_bytes > MAX_FILE_SIZE_MB * 1024 * 1024:
                    dest.unlink(missing_ok=True)
                    raise HTTPException(
                        status_code=413,
                        detail=f"Dosya {MAX_FILE_SIZE_MB} MB limitini aşıyor."
                    )
                await out.write(chunk)
    except HTTPException:
        raise
    except Exception as e:
        dest.unlink(missing_ok=True)
        logger.error("Dosya kayıt hatası: %s", e)
        raise HTTPException(status_code=500, detail="Dosya kaydedilemedi.")

    # MD5 hesapla
    md5 = compute_md5(dest)

    # DB'ye kaydet
    expires_at = (
        datetime.now(timezone.utc) + timedelta(hours=JOB_EXPIRY_HOURS)
    ).isoformat()
    db.create_job(job_id, safe_name, expires_at, ip_hash, md5)

    logger.info("Yeni iş: %s | Dosya: %s | %d bayt | MD5: %s",
                job_id, safe_name, total_bytes, md5)

    # Pipeline'ı arka planda başlat
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
        "message":    "Analiz başlatıldı. /status/{job_id} ile takip edin.",
    })


@app.get("/status/{job_id}")
async def job_status(request: Request, job_id: str):
    """İş durumunu ve aşama bilgilerini döndürür."""
    import re
    if not re.fullmatch(r"[0-9a-f]{32}", job_id):
        raise HTTPException(status_code=400, detail="Geçersiz iş ID.")

    job = db.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="İş bulunamadı.")

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


@app.get("/report/{job_id}")
async def download_report(request: Request, job_id: str):
    """Tamamlanmış HTML raporunu indir."""
    import re
    if not re.fullmatch(r"[0-9a-f]{32}", job_id):
        raise HTTPException(status_code=400, detail="Geçersiz iş ID.")

    job = db.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="İş bulunamadı.")
    if job["status"] != "completed":
        raise HTTPException(status_code=202,
                            detail=f"Analiz henüz tamamlanmadı ({job['status']}).")
    if not job["report_path"]:
        raise HTTPException(status_code=404, detail="Rapor bulunamadı.")

    report_path = Path(job["report_path"])
    if not report_path.exists():
        raise HTTPException(status_code=410, detail="Rapor silinmiş.")

    return FileResponse(
        path=str(report_path),
        media_type="text/html",
        filename=f"GenomePlatform_{job_id[:8]}_report.html",
    )


@app.get("/health")
async def health():
    return {"status": "ok"}
