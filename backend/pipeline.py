"""
Pipeline motoru:
1. FASTQ tipini tespit et (MinION / Illumina / Paired-end)
2. Aşamaları sırayla çalıştır
3. Her aşamayı veritabanına kaydet
4. Sonuçları ayrıştır ve döndür
"""
import gzip
import json
import logging
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Any

import database as db
from config import CONDA_BASE, RESULTS_DIR, SCRIPTS_DIR
from security import safe_sample_name, sanitize_path

logger = logging.getLogger("pipeline")

# ─── Yardımcılar ─────────────────────────────────────────────────────────────

def _conda_run(env_name: str, cmd: list[str], cwd: Path,
               timeout: int = 7200) -> subprocess.CompletedProcess:
    """
    Belirtilen conda ortamında komutu çalıştırır.
    shell=False kullanarak komut enjeksiyonu engellenir.
    """
    conda_bin = CONDA_BASE / "envs" / env_name / "bin"
    env_python = conda_bin / "python"

    # PATH'i conda env bin dizini ile genişlet
    import os
    env = os.environ.copy()
    env["PATH"] = str(conda_bin) + ":" + env["PATH"]
    env["CONDA_PREFIX"] = str(CONDA_BASE / "envs" / env_name)

    return subprocess.run(
        cmd,
        cwd=str(cwd),
        env=env,
        capture_output=True,
        text=True,
        timeout=timeout,
    )


def _run(cmd: list[str], cwd: Path, timeout: int = 3600,
         env_name: str = None) -> tuple[int, str, str]:
    """Komutu çalıştır; (returncode, stdout, stderr) döndür."""
    try:
        if env_name:
            result = _conda_run(env_name, cmd, cwd, timeout)
        else:
            import os
            result = subprocess.run(
                cmd, cwd=str(cwd), capture_output=True,
                text=True, timeout=timeout, env=os.environ.copy()
            )
        return result.returncode, result.stdout, result.stderr
    except subprocess.TimeoutExpired:
        return -1, "", f"Zaman aşımı ({timeout}s)"
    except FileNotFoundError as e:
        return -1, "", f"Araç bulunamadı: {e}"


# ─── FASTQ Tipi Tespiti ───────────────────────────────────────────────────────

def detect_read_type(fastq_path: Path) -> str:
    """
    İlk 2000 okumanın uzunluk dağılımına bakarak okuma tipini belirler.
    - Medyan ≥ 1000 bp → MinION (long-read)
    - Medyan < 400 bp  → Illumina (short-read)
    """
    lengths = []
    open_fn = gzip.open if str(fastq_path).endswith(".gz") else open

    try:
        with open_fn(fastq_path, "rt", errors="ignore") as fh:
            while len(lengths) < 2000:
                header = fh.readline()
                if not header:
                    break
                seq = fh.readline().strip()
                fh.readline()   # '+'
                fh.readline()   # qual
                if seq:
                    lengths.append(len(seq))
    except Exception as e:
        logger.warning("Okuma tipi tespit hatası: %s", e)
        return "unknown"

    if not lengths:
        return "unknown"

    lengths.sort()
    median = lengths[len(lengths) // 2]
    logger.info("Okuma tipi tespiti: medyan=%d bp, n=%d", median, len(lengths))

    if median >= 1000:
        return "minion"
    elif median < 400:
        return "illumina"
    else:
        return "unknown"


# ─── Pipeline Aşamaları ───────────────────────────────────────────────────────

def stage_qc(job_id: str, fastq: Path, out_dir: Path,
             read_type: str) -> dict[str, Any]:
    """QC aşaması: MinION → NanoPlot, Illumina → FastP"""
    db.update_stage(job_id, "qc", "running")
    qc_dir = out_dir / "qc"
    qc_dir.mkdir(exist_ok=True)
    result: dict[str, Any] = {}

    if read_type == "minion":
        cmd = [
            "NanoPlot", "--fastq", sanitize_path(fastq),
            "--outdir", str(qc_dir),
            "--N50", "--loglength", "--no_static",
            "--threads", "8"
        ]
        rc, out, err = _run(cmd, qc_dir, env_name="analiz")

        # NanoStats.txt ayrıştır
        nano_stats = qc_dir / "NanoStats.txt"
        if nano_stats.exists():
            for line in nano_stats.read_text().splitlines():
                if "Mean read length" in line:
                    result["mean_read_length"] = line.split(":")[-1].strip()
                elif "Mean read quality" in line:
                    result["mean_quality"] = line.split(":")[-1].strip()
                elif "Total bases" in line:
                    result["total_bases"] = line.split(":")[-1].strip()
                elif "N50" in line:
                    result["n50"] = line.split(":")[-1].strip()

    else:  # illumina
        json_out = qc_dir / "fastp.json"
        cmd = [
            "fastp",
            "-i", sanitize_path(fastq),
            "-o", str(qc_dir / "filtered.fastq.gz"),
            "--json", str(json_out),
            "--html", str(qc_dir / "fastp.html"),
            "--thread", "8",
            "--detect_adapter_for_pe",
            "-g", "-x",
            "--length_required", "50",
            "-q", "20"
        ]
        rc, out, err = _run(cmd, qc_dir)

        if json_out.exists():
            data = json.loads(json_out.read_text())
            summary = data.get("summary", {})
            af = summary.get("after_filtering", {})
            result = {
                "total_reads":       af.get("total_reads", 0),
                "total_bases":       af.get("total_bases", 0),
                "q30_rate":          round(af.get("q30_rate", 0) * 100, 1),
                "mean_read_length":  af.get("read1_mean_length", 0),
            }

    db.update_stage(job_id, "qc", "done", json.dumps(result))
    return result


def stage_assembly(job_id: str, fastq: Path, out_dir: Path,
                   read_type: str) -> Path | None:
    """Assembly aşaması: MinION → Flye, Illumina → Shovill"""
    db.update_stage(job_id, "assembly", "running")
    asm_dir = out_dir / "assembly"
    asm_dir.mkdir(exist_ok=True)

    if read_type == "minion":
        cmd = [
            "flye",
            "--nano-raw", sanitize_path(fastq),
            "--out-dir", str(asm_dir),
            "--threads", "16",
            "--genome-size", "5m",
        ]
        env = "analiz"
    else:
        cmd = [
            "shovill",
            "--R1", sanitize_path(fastq),
            "--outdir", str(asm_dir),
            "--cpus", "16",
            "--force",
        ]
        env = "shovill"

    rc, out, err = _run(cmd, asm_dir, timeout=10800, env_name=env)

    # Assembly FASTA bul
    candidates = list(asm_dir.glob("assembly.fasta")) + \
                 list(asm_dir.glob("contigs.fa"))
    if candidates:
        fasta = candidates[0]
        db.update_stage(job_id, "assembly", "done", str(fasta))
        return fasta

    db.update_stage(job_id, "assembly", "failed", err[-500:])
    return None


def stage_mlst(job_id: str, fasta: Path, out_dir: Path) -> dict[str, Any]:
    """MLST tiplemesi"""
    db.update_stage(job_id, "mlst", "running")
    mlst_dir = out_dir / "mlst"
    mlst_dir.mkdir(exist_ok=True)

    cmd = ["mlst", "--quiet", sanitize_path(fasta)]
    rc, out, err = _run(cmd, mlst_dir)

    result: dict[str, Any] = {"raw": out.strip()}
    # Çıktı formatı: dosya\tşema\tST\tallel1\tallel2...
    if out.strip():
        parts = out.strip().split("\t")
        if len(parts) >= 3:
            result["scheme"] = parts[1]
            result["st"]     = parts[2]
            result["alleles"] = parts[3:]

    db.update_stage(job_id, "mlst", "done", json.dumps(result))
    return result


def stage_amr(job_id: str, fasta: Path, out_dir: Path) -> dict[str, Any]:
    """AMR gen tespiti: AMRFinderPlus"""
    db.update_stage(job_id, "amr", "running")
    amr_dir = out_dir / "amr"
    amr_dir.mkdir(exist_ok=True)
    amr_out = amr_dir / "amrfinder.tsv"

    cmd = [
        "amrfinder",
        "--nucleotide", sanitize_path(fasta),
        "--output", str(amr_out),
        "--threads", "8",
        "--plus",
    ]
    rc, out, err = _run(cmd, amr_dir)

    genes = []
    if amr_out.exists():
        lines = amr_out.read_text().splitlines()
        headers = lines[0].split("\t") if lines else []
        for line in lines[1:]:
            parts = line.split("\t")
            if len(parts) >= len(headers):
                row = dict(zip(headers, parts))
                genes.append({
                    "gene":     row.get("Gene symbol", ""),
                    "class":    row.get("Class", ""),
                    "subclass": row.get("Subclass", ""),
                    "identity": row.get("% Identity to reference sequence", ""),
                    "method":   row.get("Method", ""),
                })

    result = {"genes": genes, "count": len(genes)}
    db.update_stage(job_id, "amr", "done", json.dumps(result))
    return result


def stage_annotation(job_id: str, fasta: Path, out_dir: Path,
                      sample_name: str) -> Path | None:
    """Bakta ile genom annotasyonu"""
    db.update_stage(job_id, "annotation", "running")
    ann_dir = out_dir / "annotation"
    ann_dir.mkdir(exist_ok=True)

    cmd = [
        "bakta",
        "--db", str(CONDA_BASE / "envs" / "bakta" / "db"),
        "--output", str(ann_dir),
        "--prefix", safe_sample_name(sample_name),
        "--threads", "8",
        "--force",
        sanitize_path(fasta),
    ]
    rc, out, err = _run(cmd, ann_dir, timeout=7200, env_name="bakta")

    gbk_files = list(ann_dir.glob("*.gbff"))
    if gbk_files:
        db.update_stage(job_id, "annotation", "done", str(gbk_files[0]))
        return gbk_files[0]

    db.update_stage(job_id, "annotation", "skipped",
                    "Bakta tamamlanamadı; rapor AMR bulgularıyla devam eder.")
    return None


# ─── Ana Pipeline ─────────────────────────────────────────────────────────────

def run_pipeline(job_id: str, fastq_path: Path) -> dict[str, Any]:
    """
    Tam pipeline'ı çalıştırır. Tüm sonuçları içeren sözlük döndürür.
    """
    out_dir = RESULTS_DIR / job_id
    out_dir.mkdir(parents=True, exist_ok=True)

    sample_name = fastq_path.stem.replace(".fastq", "").replace(".fq", "")
    sample_name = safe_sample_name(sample_name)

    results: dict[str, Any] = {
        "job_id":      job_id,
        "sample_name": sample_name,
    }

    try:
        db.update_job_status(job_id, "running")

        # 1. Okuma tipi tespiti
        db.update_stage(job_id, "detect", "running")
        read_type = detect_read_type(fastq_path)
        db.update_stage(job_id, "detect", "done", read_type)
        db.update_job_status(job_id, "running", read_type=read_type)
        results["read_type"] = read_type
        logger.info("[%s] Okuma tipi: %s", job_id, read_type)

        # 2. QC
        results["qc"] = stage_qc(job_id, fastq_path, out_dir, read_type)
        logger.info("[%s] QC tamamlandı.", job_id)

        # 3. Assembly
        fasta = stage_assembly(job_id, fastq_path, out_dir, read_type)
        results["assembly_fasta"] = str(fasta) if fasta else None
        if not fasta:
            raise RuntimeError("Assembly başarısız oldu.")
        logger.info("[%s] Assembly: %s", job_id, fasta)

        # 4. MLST
        results["mlst"] = stage_mlst(job_id, fasta, out_dir)
        logger.info("[%s] MLST: %s", job_id, results["mlst"].get("st"))

        # 5. AMR
        results["amr"] = stage_amr(job_id, fasta, out_dir)
        logger.info("[%s] AMR: %d gen bulundu.", job_id,
                    results["amr"].get("count", 0))

        # 6. Annotasyon (opsiyonel — başarısız olsa da devam et)
        try:
            results["annotation"] = str(
                stage_annotation(job_id, fasta, out_dir, sample_name)
            )
        except Exception as ann_err:
            logger.warning("[%s] Annotasyon atlandı: %s", job_id, ann_err)
            results["annotation"] = None

        db.update_job_status(job_id, "ai_pending")
        return results

    except Exception as exc:
        logger.error("[%s] Pipeline hatası: %s", job_id, exc, exc_info=True)
        db.update_job_status(job_id, "failed", error=str(exc))
        results["error"] = str(exc)
        return results
