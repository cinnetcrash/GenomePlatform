"""
Pipeline engine:
1. Detect FASTQ type (MinION / Illumina)
2. Run stages sequentially
3. Record each stage in the database
4. Return consolidated results dict
"""
import gzip
import json
import logging
import subprocess
import sys
from pathlib import Path
from typing import Any

import database as db
from config import CONDA_BASE, MAX_THREADS, RESULTS_DIR, SCRIPTS_DIR
from security import safe_sample_name, sanitize_path

logger = logging.getLogger("pipeline")


# ─── Helpers ──────────────────────────────────────────────────────────────────

def _conda_run(env_name: str, cmd: list[str], cwd: Path,
               timeout: int = 7200) -> subprocess.CompletedProcess:
    """
    Runs a command inside the specified conda environment.
    Uses shell=False to prevent command injection.
    """
    import os
    conda_bin = CONDA_BASE / "envs" / env_name / "bin"
    env = os.environ.copy()
    env["PATH"] = str(conda_bin) + ":" + env["PATH"]
    env["CONDA_PREFIX"] = str(CONDA_BASE / "envs" / env_name)

    return subprocess.run(
        cmd, cwd=str(cwd), env=env,
        capture_output=True, text=True, timeout=timeout,
    )


def _run(cmd: list[str], cwd: Path, timeout: int = 3600,
         env_name: str = None) -> tuple[int, str, str]:
    """Runs a command; returns (returncode, stdout, stderr)."""
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
        return -1, "", f"Tool timed out after {timeout}s"
    except FileNotFoundError as e:
        return -1, "", f"Tool not found: {e}"


# ─── Read Type Detection ──────────────────────────────────────────────────────

def detect_read_type(fastq_path: Path) -> str:
    """
    Determines read type from the length distribution of the first 2000 reads.
    - Median ≥ 1000 bp → MinION (long-read)
    - Median < 400 bp  → Illumina (short-read)
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
        logger.warning("Read type detection failed: %s", e)
        return "unknown"

    if not lengths:
        return "unknown"

    lengths.sort()
    median = lengths[len(lengths) // 2]
    logger.info("Read type: median=%d bp, n=%d", median, len(lengths))

    if median >= 1000:
        return "minion"
    elif median < 400:
        return "illumina"
    else:
        return "unknown"


# ─── Pipeline Stages ──────────────────────────────────────────────────────────

def stage_qc(job_id: str, fastq: Path, out_dir: Path,
             read_type: str) -> dict[str, Any]:
    """QC stage: MinION → NanoPlot, Illumina → FastP"""
    db.update_stage(job_id, "qc", "running")
    qc_dir = out_dir / "qc"
    qc_dir.mkdir(exist_ok=True)
    result: dict[str, Any] = {}

    if read_type == "minion":
        cmd = [
            "NanoPlot", "--fastq", sanitize_path(fastq),
            "--outdir", str(qc_dir),
            "--N50", "--loglength", "--no_static",
            "--threads", str(MAX_THREADS),
        ]
        rc, out, err = _run(cmd, qc_dir, env_name="analiz")

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
            "--thread", str(MAX_THREADS),
            "--detect_adapter_for_pe",
            "-g", "-x",
            "--length_required", "50",
            "-q", "20",
        ]
        rc, out, err = _run(cmd, qc_dir)

        if json_out.exists():
            data = json.loads(json_out.read_text())
            summary = data.get("summary", {})
            af = summary.get("after_filtering", {})
            result = {
                "total_reads":      af.get("total_reads", 0),
                "total_bases":      af.get("total_bases", 0),
                "q30_rate":         round(af.get("q30_rate", 0) * 100, 1),
                "mean_read_length": af.get("read1_mean_length", 0),
            }

    db.update_stage(job_id, "qc", "done", json.dumps(result))
    return result


def _assembly_stats(fasta: Path) -> dict[str, Any]:
    """
    Computes basic assembly statistics from a FASTA file:
    contig count, total length, N50, largest contig, GC%.
    """
    lengths: list[int] = []
    gc_count = 0
    total_bases = 0

    with open(fasta) as fh:
        seq_chunks: list[str] = []
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if seq_chunks:
                    s = "".join(seq_chunks)
                    lengths.append(len(s))
                    s_up = s.upper()
                    gc_count += s_up.count("G") + s_up.count("C")
                    total_bases += len(s)
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if seq_chunks:
            s = "".join(seq_chunks)
            lengths.append(len(s))
            s_up = s.upper()
            gc_count += s_up.count("G") + s_up.count("C")
            total_bases += len(s)

    lengths.sort(reverse=True)

    # N50
    n50 = 0
    cumsum = 0
    half = total_bases / 2
    for ln in lengths:
        cumsum += ln
        if cumsum >= half:
            n50 = ln
            break

    return {
        "total_contigs":    len(lengths),
        "total_length_bp":  total_bases,
        "n50_bp":           n50,
        "largest_contig_bp": lengths[0] if lengths else 0,
        "gc_percent":       round(gc_count / total_bases * 100, 1) if total_bases else 0,
    }


def stage_assembly(job_id: str, fastq: Path, out_dir: Path,
                   read_type: str) -> Path | None:
    """Assembly stage: MinION → Flye, Illumina → Shovill"""
    db.update_stage(job_id, "assembly", "running")
    asm_dir = out_dir / "assembly"
    asm_dir.mkdir(exist_ok=True)

    if read_type == "minion":
        cmd = [
            "flye",
            "--nano-raw", sanitize_path(fastq),
            "--out-dir", str(asm_dir),
            "--threads", str(MAX_THREADS),
            "--genome-size", "5m",
        ]
        env = "analiz"
    else:
        cmd = [
            "shovill",
            "--R1", sanitize_path(fastq),
            "--outdir", str(asm_dir),
            "--cpus", str(MAX_THREADS),
            "--force",
        ]
        env = "shovill"

    rc, out, err = _run(cmd, asm_dir, timeout=10800, env_name=env)

    candidates = (list(asm_dir.glob("assembly.fasta")) +
                  list(asm_dir.glob("contigs.fa")))
    if candidates:
        fasta = candidates[0]
        try:
            stats = _assembly_stats(fasta)
        except Exception as e:
            logger.warning("Assembly stats failed: %s", e)
            stats = {}
        stats["fasta_path"] = str(fasta)
        db.update_stage(job_id, "assembly", "done", json.dumps(stats))
        return fasta

    db.update_stage(job_id, "assembly", "failed", err[-500:])
    return None


def stage_mlst(job_id: str, fasta: Path, out_dir: Path) -> dict[str, Any]:
    """MLST typing"""
    db.update_stage(job_id, "mlst", "running")
    mlst_dir = out_dir / "mlst"
    mlst_dir.mkdir(exist_ok=True)

    cmd = ["mlst", "--quiet", sanitize_path(fasta)]
    rc, out, err = _run(cmd, mlst_dir, env_name="mlst")

    result: dict[str, Any] = {"raw": out.strip()}
    if out.strip():
        parts = out.strip().split("\t")
        if len(parts) >= 3:
            result["scheme"]  = parts[1]
            result["st"]      = parts[2]
            result["alleles"] = parts[3:]

    db.update_stage(job_id, "mlst", "done", json.dumps(result))
    return result


def stage_amr(job_id: str, fasta: Path, out_dir: Path) -> dict[str, Any]:
    """AMR gene detection via AMRFinderPlus"""
    db.update_stage(job_id, "amr", "running")
    amr_dir = out_dir / "amr"
    amr_dir.mkdir(exist_ok=True)
    amr_out = amr_dir / "amrfinder.tsv"

    cmd = [
        "amrfinder",
        "--nucleotide", sanitize_path(fasta),
        "--output", str(amr_out),
        "--threads", str(MAX_THREADS),
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
    """Genome annotation with Bakta"""
    db.update_stage(job_id, "annotation", "running")
    ann_dir = out_dir / "annotation"
    ann_dir.mkdir(exist_ok=True)

    cmd = [
        "bakta",
        "--db", str(CONDA_BASE / "envs" / "bakta" / "db"),
        "--output", str(ann_dir),
        "--prefix", safe_sample_name(sample_name),
        "--threads", str(MAX_THREADS),
        "--force",
        sanitize_path(fasta),
    ]
    rc, out, err = _run(cmd, ann_dir, timeout=7200, env_name="bakta")

    gbk_files = list(ann_dir.glob("*.gbff"))
    if gbk_files:
        db.update_stage(job_id, "annotation", "done", str(gbk_files[0]))
        return gbk_files[0]

    db.update_stage(job_id, "annotation", "skipped",
                    "Bakta did not complete; report continues with AMR findings.")
    return None


# ─── Main Pipeline ────────────────────────────────────────────────────────────

def run_pipeline(job_id: str, fastq_path: Path) -> dict[str, Any]:
    """
    Runs the full pipeline. Returns a dict with all results.
    On failure, sets status to 'failed' in the DB and returns {'error': ...}.
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

        # 1. Read type detection
        db.update_stage(job_id, "detect", "running")
        read_type = detect_read_type(fastq_path)
        db.update_stage(job_id, "detect", "done", read_type)
        db.update_job_status(job_id, "running", read_type=read_type)
        results["read_type"] = read_type
        logger.info("[%s] Read type: %s", job_id, read_type)

        # 2. QC
        results["qc"] = stage_qc(job_id, fastq_path, out_dir, read_type)
        logger.info("[%s] QC done.", job_id)

        # 3. Assembly
        fasta = stage_assembly(job_id, fastq_path, out_dir, read_type)
        results["assembly_fasta"] = str(fasta) if fasta else None
        if not fasta:
            raise RuntimeError("Assembly failed — no output FASTA produced.")
        logger.info("[%s] Assembly: %s", job_id, fasta)

        # 4. MLST
        results["mlst"] = stage_mlst(job_id, fasta, out_dir)
        logger.info("[%s] MLST: ST%s", job_id, results["mlst"].get("st"))

        # 5. AMR
        results["amr"] = stage_amr(job_id, fasta, out_dir)
        logger.info("[%s] AMR: %d genes found.", job_id, results["amr"].get("count", 0))

        # 6. Annotation (optional — continue even if it fails)
        try:
            results["annotation"] = str(
                stage_annotation(job_id, fasta, out_dir, sample_name)
            )
        except Exception as ann_err:
            logger.warning("[%s] Annotation skipped: %s", job_id, ann_err)
            results["annotation"] = None

        db.update_job_status(job_id, "ai_pending")
        return results

    except Exception as exc:
        error_msg = str(exc)
        logger.error("[%s] Pipeline error: %s", job_id, error_msg, exc_info=True)
        db.update_job_status(job_id, "failed", error=error_msg)
        results["error"] = error_msg
        return results
