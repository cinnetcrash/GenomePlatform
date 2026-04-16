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
import shutil
import subprocess
import sys
import threading
import time as _time
from pathlib import Path
from typing import Any

import database as db
from config import (
    ASSEMBLY_QC, ASSEMBLY_THREADS, AUTOCYCLER_BIN, BANDAGE_BIN,
    CHECKM2_DB, CHECKV_DB, CONDA_BASE,
    FLYE_ASM_COVERAGE, HUMAN_DOMINANT_THRESHOLD,
    HOST_TAXIDS, HOST_DEPLETION_THRESHOLD, HOST_DEPLETION_MIN_READS,
    MAX_THREADS, RESULTS_DIR, SCRIPTS_DIR,
    SHOVILL_ASSEMBLER, VIRAL_DOMINANT_THRESHOLD,
)
from security import safe_sample_name, sanitize_path

logger = logging.getLogger("pipeline")

# ─── Cancel Support ───────────────────────────────────────────────────────────

_thread_local = threading.local()       # stores current job_id per thread
_JOB_CANCEL: dict[str, threading.Event] = {}  # job_id → cancel event


class PipelineCancelledError(Exception):
    pass


def register_cancel_event(job_id: str) -> threading.Event:
    ev = threading.Event()
    _JOB_CANCEL[job_id] = ev
    return ev


def get_cancel_event(job_id: str) -> threading.Event | None:
    return _JOB_CANCEL.get(job_id)


def clear_cancel_event(job_id: str) -> None:
    _JOB_CANCEL.pop(job_id, None)


def _check_not_cancelled(job_id: str) -> None:
    ev = _JOB_CANCEL.get(job_id)
    if ev and ev.is_set():
        raise PipelineCancelledError("Cancelled by user.")


# ─── Helpers ──────────────────────────────────────────────────────────────────

def _run(cmd: list[str], cwd: Path, timeout: int = 3600,
         env_name: str = None) -> tuple[int, str, str]:
    """
    Runs a command and returns (returncode, stdout, stderr).
    Polls every 5 s so a cancel event can interrupt a long-running subprocess.
    """
    import os
    job_id    = getattr(_thread_local, "job_id", None)
    cancel_ev = _JOB_CANCEL.get(job_id) if job_id else None

    env = os.environ.copy()
    if env_name:
        conda_bin = CONDA_BASE / "envs" / env_name / "bin"
        env["PATH"] = str(conda_bin) + ":" + env["PATH"]
        env["CONDA_PREFIX"] = str(CONDA_BASE / "envs" / env_name)

    try:
        proc = subprocess.Popen(
            cmd, cwd=str(cwd), env=env,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        )
    except FileNotFoundError as e:
        return -1, "", f"Tool not found: {e}"
    except Exception as e:
        return -1, "", str(e)

    start = _time.monotonic()
    while True:
        try:
            stdout, stderr = proc.communicate(timeout=5)
            return proc.returncode, stdout, stderr
        except subprocess.TimeoutExpired:
            if cancel_ev and cancel_ev.is_set():
                proc.kill()
                proc.communicate()
                raise PipelineCancelledError("Cancelled by user.")
            if _time.monotonic() - start > timeout:
                proc.kill()
                proc.communicate()
                return -1, "", f"Tool timed out after {timeout}s"


# ─── Read Type Detection ──────────────────────────────────────────────────────

def detect_read_type(fastq_path: Path) -> str:
    """
    Determines read type from the length distribution of the first 2000 reads.
    Median ≥ 1000 bp → MinION,  Median < 400 bp → Illumina.
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
                fh.readline(); fh.readline()
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
    return "unknown"


# ─── Assembly Statistics & QC ─────────────────────────────────────────────────

def _assembly_stats(fasta: Path) -> dict[str, Any]:
    """Computes N50, contig count, total length, largest contig, GC% from FASTA."""
    lengths: list[int] = []
    gc_count = total_bases = 0

    with open(fasta) as fh:
        chunks: list[str] = []
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if chunks:
                    s = "".join(chunks)
                    lengths.append(len(s))
                    su = s.upper()
                    gc_count += su.count("G") + su.count("C")
                    total_bases += len(s)
                chunks = []
            else:
                chunks.append(line)
        if chunks:
            s = "".join(chunks)
            lengths.append(len(s))
            su = s.upper()
            gc_count += su.count("G") + su.count("C")
            total_bases += len(s)

    lengths.sort(reverse=True)
    n50 = cumsum = 0
    for ln in lengths:
        cumsum += ln
        if cumsum >= total_bases / 2:
            n50 = ln
            break

    return {
        "total_contigs":     len(lengths),
        "total_length_bp":   total_bases,
        "n50_bp":            n50,
        "largest_contig_bp": lengths[0] if lengths else 0,
        "gc_percent":        round(gc_count / total_bases * 100, 1) if total_bases else 0,
    }


def check_assembly_qc(stats: dict[str, Any]) -> list[dict[str, Any]]:
    """
    Evaluates assembly stats against standard bacterial genome thresholds.
    Returns a list of {metric, value, unit, status, threshold} dicts.
    status: 'pass' | 'warn' | 'fail'
    """
    checks = []

    def _status_minmax(val, cfg):
        p_min = cfg.get("pass_min"); p_max = cfg.get("pass_max")
        w_min = cfg.get("warn_min"); w_max = cfg.get("warn_max")
        if p_min and val < p_min:
            return "fail" if (not w_min or val < w_min) else "warn"
        if p_max and val > p_max:
            return "fail" if (not w_max or val > w_max) else "warn"
        return "pass"

    # N50
    n50 = stats.get("n50_bp", 0)
    cfg = ASSEMBLY_QC["n50_bp"]
    status = "pass" if n50 >= cfg["pass_min"] else ("warn" if n50 >= cfg["warn_min"] else "fail")
    checks.append({"metric": "N50", "value": n50, "unit": "bp", "status": status,
                   "threshold": "≥100 kb pass · ≥20 kb warn"})

    # Contig count
    cnt = stats.get("total_contigs", 0)
    cfg = ASSEMBLY_QC["total_contigs"]
    status = "pass" if cnt <= cfg["pass_max"] else ("warn" if cnt <= cfg["warn_max"] else "fail")
    checks.append({"metric": "Contig count", "value": cnt, "unit": "", "status": status,
                   "threshold": "≤200 pass · ≤500 warn"})

    # Total length
    tlen = stats.get("total_length_bp", 0)
    cfg = ASSEMBLY_QC["total_length_bp"]
    if cfg["pass_min"] <= tlen <= cfg["pass_max"]:
        status = "pass"
    elif cfg["warn_min"] <= tlen <= cfg["warn_max"]:
        status = "warn"
    else:
        status = "fail"
    checks.append({"metric": "Total length", "value": tlen, "unit": "bp", "status": status,
                   "threshold": "2–8 Mb pass · 1–12 Mb warn"})

    # GC%
    gc = stats.get("gc_percent", 0)
    cfg = ASSEMBLY_QC["gc_percent"]
    status = "pass" if cfg["pass_min"] <= gc <= cfg["pass_max"] else "fail"
    checks.append({"metric": "GC content", "value": gc, "unit": "%", "status": status,
                   "threshold": "25–75%"})

    # Largest contig
    lc = stats.get("largest_contig_bp", 0)
    status = "pass" if lc >= ASSEMBLY_QC["largest_contig_bp"]["warn_min"] else "warn"
    checks.append({"metric": "Largest contig", "value": lc, "unit": "bp", "status": status,
                   "threshold": "≥50 kb"})

    return checks


# ─── Kraken2 ──────────────────────────────────────────────────────────────────

def _parse_kraken2_report(report_path: Path) -> dict[str, Any]:
    """Parses a Kraken2 report file, returns classified % + top species."""
    if not report_path.exists():
        return {"classified_percent": 0, "top_taxa": [], "unclassified_percent": 100}

    unclassified_pct = 0.0
    taxa: list[dict] = []

    for line in report_path.read_text().splitlines():
        parts = line.strip().split("\t")
        if len(parts) < 6:
            continue
        pct     = float(parts[0])
        reads   = int(parts[1])
        rank    = parts[3].strip()
        name    = parts[5].strip()

        if name == "unclassified":
            unclassified_pct = pct
            continue
        if rank in ("S", "S1") and pct >= 0.01:
            taxa.append({"name": name, "percent": round(pct, 2),
                         "reads": reads, "rank": rank})

    taxa.sort(key=lambda x: x["percent"], reverse=True)
    return {
        "classified_percent":   round(100 - unclassified_pct, 1),
        "unclassified_percent": round(unclassified_pct, 1),
        "top_taxa":             taxa[:15],
    }


def decide_genomic_context(kraken2: dict | None,
                           sample_type: str) -> dict[str, Any]:
    """
    Decides the genomic context from Kraken2 results + user-supplied sample_type.

    Returns:
        {
          "context": "bacterial" | "viral_dominant" | "human_contaminated" | "metagenomics",
          "skip_bacterial": bool,   # skip MLST / AMR / Abricate / MOB-Suite
          "run_checkv":    bool,
          "reason":        str,     # human-readable explanation
          "viral_pct":     float,
          "human_pct":     float,
        }
    """
    # Default: no Kraken2 data or user explicitly chose bacterial
    if not kraken2 or not kraken2.get("top_taxa"):
        if sample_type == "metagenomics":
            return {"context": "metagenomics", "skip_bacterial": True,
                    "run_checkv": False,
                    "reason": "Sample type set to Metagenomics/Unknown by user.",
                    "viral_pct": 0.0, "human_pct": 0.0}
        return {"context": "bacterial", "skip_bacterial": False,
                "run_checkv": False,
                "reason": "No Kraken2 data — assuming bacterial WGS.",
                "viral_pct": 0.0, "human_pct": 0.0}

    # Sum viral and human percentages from top taxa
    viral_pct = 0.0
    human_pct = 0.0
    top_taxa   = kraken2.get("top_taxa", [])

    # Also check unclassified
    classified_pct = kraken2.get("classified_percent", 100.0)

    for t in top_taxa:
        name = t["name"].lower()
        pct  = t["percent"]
        if "homo sapiens" in name or "human" in name:
            human_pct += pct
        # Broad viral detection: Kraken2 assigns viral reads to virus/phage names
        if any(w in name for w in ["virus", "phage", "viridae", "virales",
                                    "viricota", "viricetes", "bacteriophage",
                                    "rhinovirus", "coronavirus", "influenza",
                                    "adenovirus", "norovirus", "rotavirus"]):
            viral_pct += pct

    # Determine context
    if human_pct >= HUMAN_DOMINANT_THRESHOLD:
        return {
            "context":        "human_contaminated",
            "skip_bacterial": True,
            "run_checkv":     False,
            "reason":         (f"Kraken2 detected {human_pct:.1f}% human reads. "
                               "MLST/AMR/Abricate/MOB-Suite skipped (not applicable to human DNA). "
                               "Assembly and quality checks will still run."),
            "viral_pct": viral_pct, "human_pct": human_pct,
        }

    if viral_pct >= VIRAL_DOMINANT_THRESHOLD or sample_type == "metagenomics":
        reason_parts = []
        if viral_pct >= VIRAL_DOMINANT_THRESHOLD:
            top_virus = next((t["name"] for t in top_taxa
                              if any(w in t["name"].lower()
                                     for w in ["virus","phage","viridae"])), "")
            reason_parts.append(
                f"Kraken2 detected {viral_pct:.1f}% viral reads"
                + (f" (dominant: {top_virus})" if top_virus else "") + "."
            )
        if sample_type == "metagenomics":
            reason_parts.append("Sample type set to Metagenomics/Unknown by user.")
        reason_parts.append(
            "MLST/AMR/Abricate/MOB-Suite skipped (bacterial-specific tools). "
            "Assembly, QUAST, CheckM2, Bandage still running."
        )
        return {
            "context":        "viral_dominant",
            "skip_bacterial": True,
            "run_checkv":     viral_pct >= VIRAL_DOMINANT_THRESHOLD,
            "reason":         " ".join(reason_parts),
            "viral_pct": viral_pct, "human_pct": human_pct,
        }

    if viral_pct > 10 or human_pct > 10:
        return {
            "context":        "metagenomics",
            "skip_bacterial": False,   # still run but with caveat
            "run_checkv":     False,
            "reason":         (f"Mixed sample: {viral_pct:.1f}% viral, {human_pct:.1f}% human. "
                               "Running all steps but results may be less reliable."),
            "viral_pct": viral_pct, "human_pct": human_pct,
        }

    return {
        "context":        "bacterial",
        "skip_bacterial": False,
        "run_checkv":     False,
        "reason":         (f"Kraken2 shows predominantly bacterial reads "
                           f"({viral_pct:.1f}% viral, {human_pct:.1f}% human)."),
        "viral_pct": viral_pct, "human_pct": human_pct,
    }


def stage_kraken2(job_id: str, fastq: Path, out_dir: Path,
                  db_path: Path,
                  fastq_r2: Path | None = None) -> dict[str, Any]:
    """Taxonomic classification with Kraken2 (SE or PE)."""
    db.update_stage(job_id, "kraken2", "running")
    k2_dir = out_dir / "kraken2"
    k2_dir.mkdir(exist_ok=True)

    report_file = k2_dir / "report.txt"
    output_file = k2_dir / "output.txt"

    cmd = [
        "kraken2",
        "--db",      str(db_path),
        "--threads", str(MAX_THREADS),
        "--report",  str(report_file),
        "--output",  str(output_file),
    ]
    if fastq_r2:
        cmd.append("--paired")
        if str(fastq).endswith(".gz"):
            cmd.append("--gzip-compressed")
        cmd += [sanitize_path(fastq), sanitize_path(fastq_r2)]
    else:
        if str(fastq).endswith(".gz"):
            cmd.append("--gzip-compressed")
        cmd.append(sanitize_path(fastq))

    rc, out, err = _run(cmd, k2_dir, timeout=3600)
    if rc != 0:
        logger.warning("[%s] Kraken2 returned %d: %s", job_id, rc, err[:200])

    result = _parse_kraken2_report(report_file)
    db.update_stage(job_id, "kraken2", "done", json.dumps(result))
    logger.info("[%s] Kraken2: %.1f%% classified, top=%s",
                job_id, result["classified_percent"],
                result["top_taxa"][0]["name"] if result["top_taxa"] else "—")
    return result


# ─── Host Depletion ───────────────────────────────────────────────────────────

def _normalize_rid(raw: str) -> str:
    """Strip /1 /2 pair suffixes and any trailing description from a read ID."""
    rid = raw.split()[0]
    if rid.endswith("/1") or rid.endswith("/2"):
        rid = rid[:-2]
    return rid


def _filter_fastq(in_path: Path, out_path: Path, exclude_ids: set) -> int:
    """
    Writes a gzipped FASTQ excluding reads whose ID is in exclude_ids.
    Returns the number of reads kept.
    """
    import gzip as _gz
    open_in = _gz.open if str(in_path).endswith(".gz") else open
    kept = 0
    with open_in(in_path, "rt", errors="ignore") as fin, \
         _gz.open(out_path, "wt") as fout:
        while True:
            header = fin.readline()
            if not header:
                break
            seq  = fin.readline()
            plus = fin.readline()
            qual = fin.readline()
            if not qual:
                break
            if _normalize_rid(header[1:]) not in exclude_ids:
                fout.write(header + seq + plus + qual)
                kept += 1
    return kept


def stage_host_depletion(job_id: str, fastq_r1: Path, out_dir: Path,
                          fastq_r2: Path | None = None,
                          host_taxids: set | None = None,
                          ) -> tuple[Path, Path | None, dict]:
    """
    Removes host reads identified by Kraken2 from the FASTQ file(s).
    Parses kraken2/output.txt (per-read taxid assignments) to build an
    exclusion set, then streams through the FASTQ writing only non-host reads.

    Returns (depleted_r1, depleted_r2_or_None, metrics_dict).
    If no host reads are found, returns the original paths unchanged.
    """
    if host_taxids is None:
        host_taxids = HOST_TAXIDS

    db.update_stage(job_id, "host_depletion", "running",
                    "Scanning Kraken2 classifications for host reads…")

    dep_dir   = out_dir / "host_depletion"
    dep_dir.mkdir(exist_ok=True)
    k2_output = out_dir / "kraken2" / "output.txt"

    # ── Parse Kraken2 per-read output ────────────────────────────────────────
    host_ids: set[str] = set()
    total_reads = 0

    if not k2_output.exists():
        msg = "Kraken2 output.txt not found — host depletion skipped."
        logger.warning("[%s] %s", job_id, msg)
        db.update_stage(job_id, "host_depletion", "skipped", msg)
        return fastq_r1, fastq_r2, {"depleted": False, "reason": msg}

    for line in k2_output.read_text(errors="ignore").splitlines():
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        total_reads += 1
        try:
            taxid = int(parts[2])
        except ValueError:
            continue
        if taxid in host_taxids:
            host_ids.add(_normalize_rid(parts[1]))

    host_count = len(host_ids)
    host_pct   = round(host_count / total_reads * 100, 2) if total_reads else 0.0

    if host_count == 0:
        result = {
            "depleted": False, "total_reads": total_reads,
            "host_reads_removed": 0, "host_percent": 0.0,
            "reads_after": total_reads,
        }
        db.update_stage(job_id, "host_depletion", "done", json.dumps(result))
        logger.info("[%s] Host depletion: 0 host reads found.", job_id)
        return fastq_r1, fastq_r2, result

    logger.info("[%s] Host depletion: removing %d reads (%.1f%%) out of %d",
                job_id, host_count, host_pct, total_reads)

    db.update_stage(job_id, "host_depletion", "running",
                    f"Filtering {host_count:,} host reads from FASTQ…")

    # ── Filter FASTQ files ───────────────────────────────────────────────────
    dep_r1  = dep_dir / "depleted_R1.fastq.gz"
    kept_r1 = _filter_fastq(fastq_r1, dep_r1, host_ids)

    dep_r2: Path | None = None
    if fastq_r2:
        dep_r2  = dep_dir / "depleted_R2.fastq.gz"
        _filter_fastq(fastq_r2, dep_r2, host_ids)

    reads_after     = kept_r1
    reads_after_pct = round(reads_after / total_reads * 100, 1) if total_reads else 0.0

    result = {
        "depleted":            True,
        "total_reads":         total_reads,
        "host_reads_removed":  host_count,
        "host_percent":        host_pct,
        "reads_after":         reads_after,
        "reads_after_percent": reads_after_pct,
        "host_taxids":         sorted(host_taxids),
    }
    db.update_stage(job_id, "host_depletion", "done", json.dumps(result))
    logger.info("[%s] Host depletion done: %d reads kept (%.1f%%)",
                job_id, reads_after, reads_after_pct)
    return dep_r1, dep_r2, result


# ─── Assembly (Illumina + MinION / Autocycler) ────────────────────────────────

def _run_flye(fastq: Path, out_dir: Path, mode: str, timeout: int = 10800,
              genome_size: str = "5m") -> Path | None:
    """Runs Flye with the given mode flag; returns assembly FASTA path or None."""
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        "flye", mode, sanitize_path(fastq),
        "--out-dir",     str(out_dir),
        "--threads",     str(ASSEMBLY_THREADS),
        "--genome-size", genome_size,
        "--no-alt-contigs",          # skip alternative contigs — saves time/memory
    ]
    if FLYE_ASM_COVERAGE > 0:
        # Subsample reads to this coverage for initial assembly.
        # Dramatically faster on high-coverage samples (>100x) with negligible quality loss.
        cmd += ["--asm-coverage", str(FLYE_ASM_COVERAGE)]
    rc, out, err = _run(cmd, out_dir, timeout=timeout, env_name="analiz")
    fasta = out_dir / "assembly.fasta"
    if fasta.exists() and fasta.stat().st_size > 0:
        return fasta
    logger.warning("Flye (%s) produced no assembly: %s", mode, err[-200:])
    return None


def _run_autocycler(assemblies: list[Path], out_dir: Path) -> Path | None:
    """
    Runs the full Autocycler pipeline on a list of assembly FASTAs.
    Returns final consensus FASTA or None on failure.
    """
    if not AUTOCYCLER_BIN.exists():
        logger.warning("Autocycler binary not found at %s", AUTOCYCLER_BIN)
        return None

    ac_dir = out_dir / "autocycler"
    staging = out_dir / "ac_assemblies"
    ac_dir.mkdir(exist_ok=True)
    staging.mkdir(exist_ok=True)

    # Copy assemblies into staging dir
    for i, fasta in enumerate(assemblies):
        shutil.copy2(fasta, staging / f"assembly_{i:02d}.fasta")

    ac = str(AUTOCYCLER_BIN)

    # 1. compress
    rc, _, err = _run([ac, "compress",
                       "--assemblies_dir", str(staging),
                       "--autocycler_dir",  str(ac_dir),
                       "--threads", str(ASSEMBLY_THREADS)],
                      ac_dir, timeout=3600)
    if rc != 0:
        logger.warning("Autocycler compress failed: %s", err[-200:])
        return None

    # 2. cluster
    rc, _, err = _run([ac, "cluster", "--autocycler_dir", str(ac_dir)],
                      ac_dir, timeout=1800)
    if rc != 0:
        logger.warning("Autocycler cluster failed: %s", err[-200:])
        return None

    # 3. trim + resolve each cluster
    cluster_dirs = sorted((ac_dir / "clusters").glob("cluster_*")) if (ac_dir / "clusters").exists() else []
    if not cluster_dirs:
        # Try flat cluster dirs (1/, 2/, ...)
        cluster_dirs = sorted(d for d in ac_dir.iterdir()
                              if d.is_dir() and d.name.isdigit())

    resolved_gfas: list[str] = []
    for cd in cluster_dirs:
        _run([ac, "trim", "--cluster_dir", str(cd),
              "--threads", str(ASSEMBLY_THREADS)], ac_dir, timeout=1800)
        _run([ac, "resolve", "--cluster_dir", str(cd)], ac_dir, timeout=1800)
        # Resolved GFA is the highest-numbered .gfa in the cluster dir
        gfas = sorted(cd.glob("*.gfa"))
        if gfas:
            resolved_gfas.append(str(gfas[-1]))

    if not resolved_gfas:
        logger.warning("Autocycler: no resolved GFAs found.")
        return None

    # 4. combine
    rc, _, err = _run([ac, "combine",
                       "--autocycler_dir", str(ac_dir),
                       "--in_gfas"] + resolved_gfas,
                      ac_dir, timeout=1800)
    if rc != 0:
        logger.warning("Autocycler combine failed: %s", err[-200:])
        return None

    # 5. gfa2fasta — find the combined GFA
    consensus_gfas = list(ac_dir.glob("consensus*.gfa")) + list(ac_dir.glob("*.gfa"))
    if not consensus_gfas:
        logger.warning("Autocycler: no consensus GFA found after combine.")
        return None

    final_fasta = out_dir / "autocycler_assembly.fasta"
    rc, _, err = _run([ac, "gfa2fasta",
                       "--in_gfa",    str(consensus_gfas[0]),
                       "--out_fasta", str(final_fasta)],
                      ac_dir, timeout=600)
    if rc == 0 and final_fasta.exists() and final_fasta.stat().st_size > 0:
        return final_fasta

    logger.warning("Autocycler gfa2fasta failed: %s", err[-200:])
    return None


def _estimate_genome_size(kraken2: dict | None) -> str:
    """
    Rough genome size estimate for Flye based on Kraken2 top hit.
    Flye uses this only for coverage estimation — a rough value is fine.
    Returns a string like "5m", "50k", "3m".
    """
    if not kraken2:
        return "5m"
    top = (kraken2.get("top_taxa") or [{}])[0]
    name = top.get("name", "").lower()
    # Phage / virus → ~50–200 kb
    if any(w in name for w in ("phage", "virus", "viridae", "virales")):
        return "200k"
    # Human / mammal → ~3 Gb (unlikely to assemble well, but set something sane)
    if any(w in name for w in ("homo", "human", "mus ", "rattus")):
        return "3g"
    # Typical bacteria → 2–7 Mb
    return "5m"


def stage_assembly(job_id: str, fastq: Path, out_dir: Path,
                   read_type: str,
                   fastq_r2: Path | None = None,
                   kraken2: dict | None = None) -> tuple[Path | None, list]:
    """
    Assembly stage:
    - MinION  → Flye (--nano-raw) + Flye (--nano-hq) → Autocycler consensus
                Falls back to best single Flye result if Autocycler fails.
    - Illumina → Shovill
    """
    db.update_stage(job_id, "assembly", "running")
    asm_dir = out_dir / "assembly"
    asm_dir.mkdir(exist_ok=True)

    fasta: Path | None = None
    method = "unknown"

    if read_type == "minion":
        gsize = _estimate_genome_size(kraken2)
        logger.info("[%s] Assembly: Flye genome-size estimate: %s", job_id, gsize)

        logger.info("[%s] Assembly: running Flye --nano-raw", job_id)
        flye_raw = _run_flye(fastq, asm_dir / "flye_raw", "--nano-raw", genome_size=gsize)

        logger.info("[%s] Assembly: running Flye --nano-hq", job_id)
        flye_hq  = _run_flye(fastq, asm_dir / "flye_hq",  "--nano-hq",  genome_size=gsize)

        available = [f for f in [flye_raw, flye_hq] if f]
        if len(available) >= 2:
            logger.info("[%s] Assembly: running Autocycler (%d inputs)", job_id, len(available))
            ac_fasta = _run_autocycler(available, asm_dir)
            if ac_fasta:
                fasta  = ac_fasta
                method = "autocycler (flye_raw + flye_hq)"
            else:
                # Fall back: pick the assembly with the best N50
                best = max(available, key=lambda f: _assembly_stats(f).get("n50_bp", 0))
                fasta  = best
                method = f"flye_fallback ({best.parent.name})"
                logger.warning("[%s] Autocycler failed — using best single Flye result", job_id)
        elif available:
            fasta  = available[0]
            method = available[0].parent.name
        else:
            db.update_stage(job_id, "assembly", "failed", "All Flye runs failed.")
            return None, []

    else:  # Illumina
        import psutil
        avail_gb = max(4, int(psutil.virtual_memory().available / 1e9 * 0.90))

        if fastq_r2:
            # ── Paired-end → Shovill ──────────────────────────────────────────
            # --noreadcorr: skip Shovill's lighter/trimmomatic correction step;
            # fastp already handled QC upstream → significant speedup.
            cmd = [
                "shovill",
                "--R1",        sanitize_path(fastq),
                "--R2",        sanitize_path(fastq_r2),
                "--outdir",    str(asm_dir),
                "--cpus",      str(ASSEMBLY_THREADS),
                "--ram",       str(avail_gb),
                "--assembler", SHOVILL_ASSEMBLER,
                "--noreadcorr",
                "--force",
            ]
            rc, out, err = _run(cmd, asm_dir, timeout=10800, env_name="shovill")
            candidates = list(asm_dir.glob("contigs.fa")) + list(asm_dir.glob("assembly.fasta"))
            if candidates:
                fasta  = candidates[0]
                method = f"shovill-{SHOVILL_ASSEMBLER} (PE)"
            else:
                db.update_stage(job_id, "assembly", "failed", err[-500:])
                return None, []

        else:
            # ── Single-end → SPAdes directly ──────────────────────────────────
            # Shovill requires R2; for SE data use SPAdes --s1 mode.
            spades_bin = CONDA_BASE / "envs" / "shovill" / "bin" / "spades.py"
            spades_out = asm_dir / "spades_se"
            spades_out.mkdir(exist_ok=True)
            cmd = [
                str(spades_bin),
                "--s1",         sanitize_path(fastq),
                "-o",           str(spades_out),
                "--threads",    str(ASSEMBLY_THREADS),
                "--memory",     str(avail_gb),
                "--only-assembler",   # fastp already did QC
                "--cov-cutoff",  "auto",
            ]
            rc, out, err = _run(cmd, spades_out, timeout=10800, env_name="shovill")
            # SPAdes outputs contigs.fasta
            candidates = list(spades_out.glob("contigs.fasta")) + \
                         list(spades_out.glob("scaffolds.fasta"))
            if candidates:
                # Copy to standard name expected by downstream stages
                import shutil as _shutil
                fasta = asm_dir / "contigs.fa"
                _shutil.copy2(candidates[0], fasta)
                method = "spades-se"
            else:
                db.update_stage(job_id, "assembly", "failed", err[-500:])
                return None, []

    # Compute and store stats
    try:
        stats = _assembly_stats(fasta)
    except Exception as e:
        logger.warning("Assembly stats failed: %s", e)
        stats = {}

    qc_checks = check_assembly_qc(stats)
    failed_checks = [c for c in qc_checks if c["status"] == "fail"]
    warn_checks   = [c for c in qc_checks if c["status"] == "warn"]

    if failed_checks:
        logger.warning("[%s] Assembly QC FAILED: %s",
                       job_id, [c["metric"] for c in failed_checks])
    elif warn_checks:
        logger.warning("[%s] Assembly QC WARN: %s",
                       job_id, [c["metric"] for c in warn_checks])
    else:
        logger.info("[%s] Assembly QC passed all checks.", job_id)

    detail = {**stats, "fasta_path": str(fasta), "method": method,
              "qc_checks": qc_checks}
    db.update_stage(job_id, "assembly", "done", json.dumps(detail))
    return fasta, qc_checks


# ─── Bandage Graph Visualization ─────────────────────────────────────────────

def _find_assembly_gfa(fasta: Path) -> Path | None:
    """
    Searches for the primary assembly graph GFA near the given FASTA.
    Flye outputs assembly_graph.gfa in the same directory.
    Autocycler leaves a consensus GFA in asm_dir/autocycler/.
    """
    asm_dir = fasta.parent

    # Autocycler: consensus gfa next to the autocycler fasta
    for gfa in asm_dir.glob("*.gfa"):
        return gfa

    # Flye: look in sub-directories
    for subdir in [asm_dir / "flye_raw", asm_dir / "flye_hq"]:
        g = subdir / "assembly_graph.gfa"
        if g.exists():
            return g

    # Shovill: contigs.gfa
    g = asm_dir / "contigs.gfa"
    if g.exists():
        return g

    return None


def stage_bandage(job_id: str, fasta: Path, out_dir: Path) -> dict[str, Any]:
    """
    Generates a PNG image of the assembly graph using Bandage.
    Returns {"image_path": str, "gfa_path": str} or {"error": ...}.
    """
    if not BANDAGE_BIN.exists():
        return {"error": f"Bandage not found at {BANDAGE_BIN}"}

    db.update_stage(job_id, "bandage", "running", "Rendering assembly graph…")

    gfa = _find_assembly_gfa(fasta.parent)
    if not gfa:
        msg = "No assembly graph GFA found."
        db.update_stage(job_id, "bandage", "skipped", msg)
        return {"error": msg}

    img_path = out_dir / "assembly_graph.png"
    import os
    env = os.environ.copy()
    env["QT_QPA_PLATFORM"] = "offscreen"   # headless Qt rendering

    try:
        result = subprocess.run(
            [str(BANDAGE_BIN), "image", str(gfa), str(img_path),
             "--width", "1400", "--height", "1000"],
            capture_output=True, text=True, timeout=300, env=env,
        )
        if img_path.exists() and img_path.stat().st_size > 0:
            db.update_stage(job_id, "bandage", "done",
                            json.dumps({"image_path": str(img_path),
                                        "gfa_path":   str(gfa)}))
            logger.info("[%s] Bandage graph rendered: %s", job_id, img_path)
            return {"image_path": str(img_path), "gfa_path": str(gfa)}

        err = result.stderr[:300]
        logger.warning("[%s] Bandage failed (rc=%d): %s", job_id, result.returncode, err)
        db.update_stage(job_id, "bandage", "failed", err)
        return {"error": err}

    except subprocess.TimeoutExpired:
        db.update_stage(job_id, "bandage", "failed", "Bandage timed out.")
        return {"error": "Bandage timed out."}


# ─── Read Depth / Coverage ────────────────────────────────────────────────────

def stage_coverage(job_id: str, fasta: Path, out_dir: Path,
                   read_type: str,
                   fastq_r1: Path,
                   fastq_r2: Path | None = None) -> dict[str, Any]:
    """
    Maps reads back to the assembly with minimap2, then computes per-base
    coverage with samtools.  Returns mean/median depth and coverage breadth.
    """
    db.update_stage(job_id, "coverage", "running", "Mapping reads to assembly…")
    cov_dir = out_dir / "coverage"
    cov_dir.mkdir(exist_ok=True)

    bam = cov_dir / "reads_vs_assembly.bam"
    preset = "map-ont" if read_type == "minion" else "sr"

    # ── minimap2 → samtools sort ───────────────────────────────────────────
    mm2_cmd = [
        "minimap2", "-ax", preset,
        "-t", str(MAX_THREADS),
        sanitize_path(fasta),
        sanitize_path(fastq_r1),
    ]
    if fastq_r2:
        mm2_cmd.append(sanitize_path(fastq_r2))

    import os
    env = os.environ.copy()
    env["PATH"] = str(CONDA_BASE / "envs" / "analiz" / "bin") + ":" + env["PATH"]

    try:
        mm2 = subprocess.Popen(mm2_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               cwd=str(cov_dir), env=env)
        sort = subprocess.Popen(
            ["samtools", "sort", "-@", str(MAX_THREADS), "-o", str(bam)],
            stdin=mm2.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            cwd=str(cov_dir), env=env,
        )
        mm2.stdout.close()
        _, sort_err = sort.communicate(timeout=3600)
        mm2.wait(timeout=60)
    except Exception as e:
        msg = f"minimap2/samtools sort failed: {e}"
        logger.warning("[%s] Coverage: %s", job_id, msg)
        db.update_stage(job_id, "coverage", "failed", msg)
        return {"error": msg}

    if not bam.exists() or bam.stat().st_size == 0:
        msg = "BAM file not produced"
        db.update_stage(job_id, "coverage", "failed", msg)
        return {"error": msg}

    # ── samtools index ─────────────────────────────────────────────────────
    subprocess.run(["samtools", "index", str(bam)], env=env, cwd=str(cov_dir),
                   capture_output=True)

    # ── samtools coverage → aggregate stats ───────────────────────────────
    rc_cov = subprocess.run(
        ["samtools", "coverage", str(bam)],
        capture_output=True, text=True, env=env, cwd=str(cov_dir),
    )

    result: dict[str, Any] = {}
    if rc_cov.returncode == 0 and rc_cov.stdout:
        depths, covs, lengths = [], [], []
        for line in rc_cov.stdout.splitlines():
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 7:
                continue
            try:
                length    = int(parts[2])
                covbases  = int(parts[3])
                depth     = float(parts[6])
                lengths.append(length)
                depths.append(depth * length)
                covs.append(covbases)
            except (ValueError, IndexError):
                continue

        total_len = sum(lengths)
        if total_len > 0:
            mean_depth   = round(sum(depths) / total_len, 1)
            breadth_1x   = round(sum(covs) / total_len * 100, 1)
            result = {
                "mean_depth":      mean_depth,
                "breadth_1x_pct":  breadth_1x,
                "total_length_bp": total_len,
            }
            logger.info("[%s] Coverage: %.1fx mean depth, %.1f%% at ≥1x",
                        job_id, mean_depth, breadth_1x)
        else:
            result = {"error": "No coverage data computed"}
    else:
        result = {"error": rc_cov.stderr[:300] if rc_cov.stderr else "samtools coverage failed"}

    db.update_stage(job_id, "coverage", "done", json.dumps(result))
    return result


# ─── Serotyping ───────────────────────────────────────────────────────────────

# Maps genus (lowercase) → (tool_name, conda_env)
_SEROTYPE_MAP = {
    "escherichia": ("ectyper",   "analiz"),
    "salmonella":  ("sistr",     "sistr"),
    "klebsiella":  ("kleborate", "kleborate"),
}


def _detect_serotype_tool(kraken2_results: dict | None) -> tuple[str | None, str | None, str]:
    """
    Inspects the Kraken2 top hit to decide which serotyping tool to run.
    Returns (tool_name, conda_env, detected_species) or (None, None, "").
    """
    if not kraken2_results:
        return None, None, ""
    top_taxa = kraken2_results.get("top_taxa") or []
    if not top_taxa:
        return None, None, ""
    top_name = top_taxa[0].get("name", "")
    genus = top_name.split()[0].lower() if top_name else ""
    tool, env = _SEROTYPE_MAP.get(genus, (None, None))
    return tool, env, top_name


def _run_ectyper(job_id: str, fasta: Path, out_dir: Path) -> dict[str, Any]:
    eco_dir = out_dir / "serotyping" / "ectyper"
    eco_dir.mkdir(parents=True, exist_ok=True)
    cmd = ["ectyper", "-i", sanitize_path(fasta), "-o", str(eco_dir),
           "--cores", str(MAX_THREADS)]
    rc, out, err = _run(cmd, eco_dir, timeout=1800, env_name="analiz")

    tsv = eco_dir / "output.tsv"
    if not tsv.exists():
        # ectyper sometimes writes to a subdirectory
        hits = list(eco_dir.rglob("output.tsv"))
        if hits:
            tsv = hits[0]

    result: dict[str, Any] = {"tool": "ECTyper"}
    if tsv.exists():
        lines = tsv.read_text().splitlines()
        if len(lines) >= 2:
            headers = lines[0].split("\t")
            vals    = lines[1].split("\t")
            row = dict(zip(headers, vals))
            result.update({
                "serotype":   row.get("Serotype", row.get("O-type", "—")),
                "o_type":     row.get("O-type", "—"),
                "h_type":     row.get("H-type", "—"),
                "quality":    row.get("QC", row.get("Quality", "—")),
                "evidence":   row.get("Evidence", "—"),
            })
    elif rc != 0:
        result["error"] = err[:300]
    return result


def _run_sistr(job_id: str, fasta: Path, out_dir: Path) -> dict[str, Any]:
    sis_dir = out_dir / "serotyping" / "sistr"
    sis_dir.mkdir(parents=True, exist_ok=True)
    json_out = sis_dir / "sistr_results.json"
    cmd = ["sistr", "-i", sanitize_path(fasta), "assembly",
           "-f", "json", "-o", str(json_out),
           "--qc", "-t", str(MAX_THREADS)]
    rc, out, err = _run(cmd, sis_dir, timeout=1800, env_name="sistr")

    result: dict[str, Any] = {"tool": "SISTR"}
    if json_out.exists():
        try:
            data = json.loads(json_out.read_text())
            if isinstance(data, list):
                data = data[0]
            result.update({
                "serovar":       data.get("serovar", "—"),
                "serogroup":     data.get("serogroup", "—"),
                "h1":            data.get("h1", "—"),
                "h2":            data.get("h2", "—"),
                "o_antigen":     data.get("o_antigen", "—"),
                "cgmlst_st":     data.get("cgmlst_ST", "—"),
                "qc_status":     data.get("qc_status", "—"),
            })
        except Exception as e:
            result["error"] = f"JSON parse error: {e}"
    elif rc != 0:
        result["error"] = err[:300]
    return result


def _run_kleborate(job_id: str, fasta: Path, out_dir: Path) -> dict[str, Any]:
    kleb_dir = out_dir / "serotyping" / "kleborate"
    kleb_dir.mkdir(parents=True, exist_ok=True)
    tsv_out = kleb_dir / "kleborate.txt"
    cmd = ["kleborate", "--all", "-a", sanitize_path(fasta), "-o", str(tsv_out)]
    rc, out, err = _run(cmd, kleb_dir, timeout=1800, env_name="kleborate")

    result: dict[str, Any] = {"tool": "Kleborate"}
    if tsv_out.exists():
        lines = tsv_out.read_text().splitlines()
        if len(lines) >= 2:
            headers = lines[0].split("\t")
            vals    = lines[1].split("\t")
            row = dict(zip(headers, vals))
            result.update({
                "species":        row.get("species", "—"),
                "st":             row.get("ST", "—"),
                "virulence_score":row.get("virulence_score", "—"),
                "resistance_score":row.get("resistance_score", "—"),
                "o_type":         row.get("O_locus", "—"),
                "k_type":         row.get("K_locus", "—"),
            })
    elif rc != 0:
        result["error"] = err[:300]
    return result


def stage_serotyping(job_id: str, fasta: Path, out_dir: Path,
                     kraken2_results: dict | None) -> dict[str, Any]:
    """
    Auto-selects and runs the appropriate serotyping tool based on Kraken2 genus:
      Escherichia → ECTyper   (O:H antigen)
      Salmonella  → SISTR     (serovar prediction)
      Klebsiella  → Kleborate (K/O locus + virulence/resistance scoring)
    Skips gracefully for all other genera.
    """
    tool, env, species = _detect_serotype_tool(kraken2_results)
    if not tool:
        msg = (f"Serotyping skipped — genus not in supported list "
               f"(supported: {', '.join(g.capitalize() for g in _SEROTYPE_MAP)}).")
        db.update_stage(job_id, "serotyping", "skipped", msg)
        return {"skipped": True, "reason": msg}

    db.update_stage(job_id, "serotyping", "running",
                    f"Running {tool} for {species}…")
    logger.info("[%s] Serotyping: tool=%s species=%s", job_id, tool, species)

    try:
        if tool == "ectyper":
            result = _run_ectyper(job_id, fasta, out_dir)
        elif tool == "sistr":
            result = _run_sistr(job_id, fasta, out_dir)
        elif tool == "kleborate":
            result = _run_kleborate(job_id, fasta, out_dir)
        else:
            result = {"skipped": True, "reason": f"Unknown tool: {tool}"}
    except Exception as exc:
        result = {"error": str(exc)}

    result["detected_species"] = species
    db.update_stage(job_id, "serotyping", "done", json.dumps({
        "tool": result.get("tool", tool),
        "detected_species": species,
        "serotype": result.get("serotype") or result.get("serovar") or result.get("st", "—"),
    }))
    logger.info("[%s] Serotyping done: %s", job_id, result)
    return result


# ─── MLST ─────────────────────────────────────────────────────────────────────

def stage_mlst(job_id: str, fasta: Path, out_dir: Path) -> dict[str, Any]:
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


# ─── AMR ──────────────────────────────────────────────────────────────────────

def stage_amr(job_id: str, fasta: Path, out_dir: Path) -> dict[str, Any]:
    db.update_stage(job_id, "amr", "running")
    amr_dir = out_dir / "amr"
    amr_dir.mkdir(exist_ok=True)
    amr_out = amr_dir / "amrfinder.tsv"

    cmd = [
        "amrfinder",
        "--nucleotide", sanitize_path(fasta),
        "--output",     str(amr_out),
        "--threads",    str(MAX_THREADS),
        "--plus",
    ]
    rc, out, err = _run(cmd, amr_dir)

    genes = []
    if amr_out.exists():
        lines   = amr_out.read_text().splitlines()
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


# ─── Abricate ─────────────────────────────────────────────────────────────────

VFDB_SEQUENCES = Path("/home/analysis/miniconda3/envs/analiz/db/vfdb/sequences")

_VFDB_GENERA_CACHE: set[str] | None = None


def _vfdb_genera() -> set[str]:
    """Lazily loads the set of genera present in the VFDB sequences file."""
    global _VFDB_GENERA_CACHE
    if _VFDB_GENERA_CACHE is not None:
        return _VFDB_GENERA_CACHE
    genera: set[str] = set()
    if VFDB_SEQUENCES.exists():
        import re
        # Header format: >vfdb~~~ACC~~~ACC (product) [Genus species str. ...]
        pattern = re.compile(r'\[([A-Z][a-z]+)\s+\w+')
        for line in VFDB_SEQUENCES.read_text(errors="ignore").splitlines():
            if line.startswith(">"):
                m = pattern.search(line)
                if m:
                    genera.add(m.group(1))
    _VFDB_GENERA_CACHE = genera
    logger.info("VFDB genera loaded: %d genera", len(genera))
    return genera


def _genus_in_vfdb(genus: str) -> bool:
    """Returns True if the given genus has entries in the VFDB database."""
    if not genus:
        return False
    return genus.strip().capitalize() in _vfdb_genera()


def _parse_abricate_tsv(tsv_path: Path, db_name: str) -> list[dict]:
    """Parses an Abricate TSV output file, returns list of hit dicts."""
    if not tsv_path.exists():
        return []
    genes = []
    for line in tsv_path.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 15:
            continue
        genes.append({
            "gene":       parts[5],
            "coverage":   parts[9],
            "identity":   parts[10],
            "database":   db_name,
            "accession":  parts[12],
            "product":    parts[13],
            "resistance": parts[14] if len(parts) > 14 else "",
        })
    return genes


def stage_abricate(job_id: str, fasta: Path, out_dir: Path,
                   kraken2_results: dict | None) -> dict:
    """
    Runs Abricate on the assembly:
    - Always: --db card (AMR/resistance genes)
    - Conditionally: --db vfdb (if top Kraken2 genus is in VFDB)
    """
    db.update_stage(job_id, "abricate", "running",
                    f"Running abricate --db card on {fasta.name}")
    abr_dir = out_dir / "abricate"
    abr_dir.mkdir(exist_ok=True)

    all_genes: list[dict] = []

    # Always run CARD
    card_tsv = abr_dir / "card.tsv"
    cmd_card = ["abricate", "--db", "card", "--quiet", sanitize_path(fasta)]
    rc, out_c, err_c = _run(cmd_card, abr_dir, timeout=600, env_name="analiz")
    if rc == 0:
        card_tsv.write_text(out_c)
        card_hits = _parse_abricate_tsv(card_tsv, "CARD")
        all_genes.extend(card_hits)
        logger.info("[%s] Abricate CARD: %d hits", job_id, len(card_hits))
    else:
        logger.warning("[%s] Abricate CARD failed: %s", job_id, err_c[:200])

    # Conditionally run VFDB
    run_vfdb = False
    top_genus = ""
    if kraken2_results and kraken2_results.get("top_taxa"):
        top_name = kraken2_results["top_taxa"][0].get("name", "")
        top_genus = top_name.split()[0] if top_name else ""
        run_vfdb = _genus_in_vfdb(top_genus)
        if run_vfdb:
            logger.info("[%s] %s is in VFDB — running abricate --db vfdb", job_id, top_genus)

    if run_vfdb:
        db.update_stage(job_id, "abricate", "running",
                        f"Running abricate --db vfdb (genus: {top_genus})")
        vfdb_tsv = abr_dir / "vfdb.tsv"
        cmd_vfdb = ["abricate", "--db", "vfdb", "--quiet", sanitize_path(fasta)]
        rc, out_v, err_v = _run(cmd_vfdb, abr_dir, timeout=600, env_name="analiz")
        if rc == 0:
            vfdb_tsv.write_text(out_v)
            vfdb_hits = _parse_abricate_tsv(vfdb_tsv, "VFDB")
            all_genes.extend(vfdb_hits)
            logger.info("[%s] Abricate VFDB: %d hits", job_id, len(vfdb_hits))
        else:
            logger.warning("[%s] Abricate VFDB failed: %s", job_id, err_v[:200])

    # Always run PlasmidFinder
    db.update_stage(job_id, "abricate", "running",
                    "Running abricate --db plasmidfinder")
    pf_tsv = abr_dir / "plasmidfinder.tsv"
    cmd_pf = ["abricate", "--db", "plasmidfinder", "--quiet", sanitize_path(fasta)]
    rc, out_pf, err_pf = _run(cmd_pf, abr_dir, timeout=600, env_name="analiz")
    if rc == 0:
        pf_tsv.write_text(out_pf)
    pf_hits = _parse_abricate_tsv(pf_tsv, "PlasmidFinder") if pf_tsv.exists() else []
    logger.info("[%s] Abricate PlasmidFinder: %d hits", job_id, len(pf_hits))

    result = {
        "genes":               all_genes,
        "count":               len(all_genes),
        "ran_vfdb":            run_vfdb,
        "top_genus":           top_genus,
        "plasmidfinder":       pf_hits,
        "plasmidfinder_count": len(pf_hits),
    }
    db.update_stage(job_id, "abricate", "done",
                    json.dumps({"count": len(all_genes), "ran_vfdb": run_vfdb,
                                "top_genus": top_genus,
                                "plasmidfinder_count": len(pf_hits)}))
    return result


# ─── QUAST ────────────────────────────────────────────────────────────────────

def stage_quast(job_id: str, fasta: Path, out_dir: Path) -> dict[str, Any]:
    """Runs QUAST on the assembly and parses the key metrics."""
    db.update_stage(job_id, "quast", "running", f"Running QUAST on {fasta.name}")
    q_dir = out_dir / "quast"
    q_dir.mkdir(exist_ok=True)

    cmd = [
        "quast.py", sanitize_path(fasta),
        "-o", str(q_dir),
        "--threads", str(MAX_THREADS),
        "--no-html", "--no-plots",
    ]
    rc, out_s, err_s = _run(cmd, q_dir, timeout=1800, env_name="quast5")

    result: dict[str, Any] = {}
    report_tsv = q_dir / "report.tsv"
    if report_tsv.exists():
        for line in report_tsv.read_text().splitlines():
            parts = line.split("\t")
            if len(parts) >= 2:
                key = parts[0].strip()
                val = parts[1].strip()
                result[key] = val

    summary = {
        "contigs":        result.get("# contigs (>= 0 bp)", result.get("# contigs", "—")),
        "total_length":   result.get("Total length (>= 0 bp)", result.get("Total length", "—")),
        "largest_contig": result.get("Largest contig", "—"),
        "n50":            result.get("N50", "—"),
        "n90":            result.get("N90", "—"),
        "l50":            result.get("L50", "—"),
        "gc_pct":         result.get("GC (%)", "—"),
        "ns_per_100k":    result.get("# N's per 100 kbp", "—"),
    }
    db.update_stage(job_id, "quast", "done", json.dumps(summary))
    logger.info("[%s] QUAST: N50=%s, contigs=%s", job_id,
                summary["n50"], summary["contigs"])
    return summary


# ─── CheckM2 ──────────────────────────────────────────────────────────────────

def stage_checkm2(job_id: str, fasta: Path, out_dir: Path) -> dict[str, Any]:
    """Assesses genome completeness and contamination using CheckM2."""
    db.update_stage(job_id, "checkm2", "running", f"Running CheckM2 on {fasta.name}")

    if not CHECKM2_DB.exists():
        msg = f"CheckM2 database not found at {CHECKM2_DB}"
        db.update_stage(job_id, "checkm2", "skipped", msg)
        logger.warning("[%s] %s", job_id, msg)
        return {"error": msg}

    cm_dir   = out_dir / "checkm2"
    fasta_dir = out_dir / "checkm2_input"
    cm_dir.mkdir(exist_ok=True)
    fasta_dir.mkdir(exist_ok=True)

    # CheckM2 expects a directory of FASTA files
    import shutil
    shutil.copy2(fasta, fasta_dir / fasta.name)

    cmd = [
        "checkm2", "predict",
        "--input",            str(fasta_dir),
        "--output-directory", str(cm_dir),
        "--threads",          str(MAX_THREADS),
        "--database_path",    str(CHECKM2_DB),
        "--force",
        "-x", fasta.suffix.lstrip("."),
    ]
    rc, out_s, err_s = _run(cmd, cm_dir, timeout=3600, env_name="checkM")

    result: dict[str, Any] = {}
    report = cm_dir / "quality_report.tsv"
    if report.exists():
        lines = report.read_text().splitlines()
        if len(lines) >= 2:
            headers = lines[0].split("\t")
            vals    = lines[1].split("\t")
            row = dict(zip(headers, vals))
            result = {
                "completeness":   row.get("Completeness", "—"),
                "contamination":  row.get("Contamination", "—"),
                "model_used":     row.get("Completeness_Model_Used", "—"),
            }
    else:
        result = {"error": err_s[:300] if rc != 0 else "No output produced."}

    db.update_stage(job_id, "checkm2", "done", json.dumps(result))
    logger.info("[%s] CheckM2: completeness=%s contamination=%s",
                job_id, result.get("completeness"), result.get("contamination"))
    return result


# ─── CheckV ───────────────────────────────────────────────────────────────────

def stage_checkv(job_id: str, fasta: Path, out_dir: Path) -> dict[str, Any]:
    """Assesses viral genome completeness using CheckV (run when viral dominant)."""
    db.update_stage(job_id, "checkv", "running", f"Running CheckV on {fasta.name}")

    cv_dir = out_dir / "checkv"
    cv_dir.mkdir(exist_ok=True)

    # CheckV will use its internal database if --db not specified
    cmd = ["checkv", "end_to_end",
           sanitize_path(fasta), str(cv_dir),
           "-t", str(MAX_THREADS)]
    if CHECKV_DB.is_dir():
        cmd += ["--db", str(CHECKV_DB)]

    rc, out_s, err_s = _run(cmd, cv_dir, timeout=3600, env_name="checkv")

    result: dict[str, Any] = {}
    summary_tsv = cv_dir / "quality_summary.tsv"
    if summary_tsv.exists():
        lines = summary_tsv.read_text().splitlines()
        if len(lines) >= 2:
            headers = lines[0].split("\t")
            viral_count = len(lines) - 1
            complete = sum(1 for l in lines[1:]
                           if "Complete" in l or "High-quality" in l)
            result = {
                "viral_contigs":      viral_count,
                "complete_or_hq":     complete,
            }
            # Top hits
            tops = []
            for l in lines[1:6]:
                parts = l.split("\t")
                if len(parts) >= len(headers):
                    row = dict(zip(headers, parts))
                    tops.append({
                        "contig":     row.get("contig_id", ""),
                        "checkv_quality": row.get("checkv_quality", ""),
                        "genome_copies":  row.get("genome_copies", ""),
                    })
            result["top_contigs"] = tops
    else:
        result = {"error": err_s[:300] if rc != 0 else "No CheckV output."}

    db.update_stage(job_id, "checkv", "done", json.dumps(result))
    logger.info("[%s] CheckV: %s viral contigs", job_id,
                result.get("viral_contigs", "?"))
    return result


# ─── MOB-Suite (Plasmid typing) ───────────────────────────────────────────────

def stage_mobsuite(job_id: str, fasta: Path, out_dir: Path) -> dict[str, Any]:
    """
    Runs MOB-Recon to classify contigs as chromosome or plasmid.
    Returns summary: plasmid count, replicon types, mobility.
    """
    db.update_stage(job_id, "mobsuite", "running",
                    f"Running mob_recon on {fasta.name}")
    mob_dir = out_dir / "mobsuite"
    mob_dir.mkdir(exist_ok=True)

    cmd = [
        "mob_recon",
        "-i", sanitize_path(fasta),
        "-o", str(mob_dir),
        "-n", str(MAX_THREADS),
        "--force",
    ]
    rc, out_s, err_s = _run(cmd, mob_dir, timeout=3600, env_name="mobsuite")

    result: dict[str, Any] = {"plasmids": [], "plasmid_count": 0, "error": None}

    if rc != 0:
        logger.warning("[%s] mob_recon exited %d: %s", job_id, rc, err_s[:200])
        result["error"] = err_s[:300]
        db.update_stage(job_id, "mobsuite", "done", json.dumps(result))
        return result

    # Parse contig_report.txt
    contig_report = mob_dir / "contig_report.txt"
    mobtyper = mob_dir / "mobtyper_results.txt"

    plasmid_contigs: dict[str, dict] = {}  # plasmid_id -> info

    if contig_report.exists():
        lines = contig_report.read_text().splitlines()
        if lines:
            headers = lines[0].split("\t")
            for line in lines[1:]:
                parts = line.split("\t")
                if len(parts) < len(headers):
                    continue
                row = dict(zip(headers, parts))
                element_type = row.get("molecule_type", "").lower()
                if "plasmid" in element_type:
                    pid = row.get("primary_cluster_id", row.get("plasmid_id", "unknown"))
                    if pid and pid not in plasmid_contigs:
                        plasmid_contigs[pid] = {
                            "id":         pid,
                            "replicons":  row.get("rep_type(s)", "—"),
                            "mobility":   row.get("predicted_mobility", "—"),
                            "mpf":        row.get("mpf_type", "—"),
                            "contigs":    0,
                            "size_bp":    0,
                        }
                    if pid in plasmid_contigs:
                        plasmid_contigs[pid]["contigs"] += 1
                        try:
                            plasmid_contigs[pid]["size_bp"] += int(row.get("contig_length", 0))
                        except (ValueError, TypeError):
                            pass

    # Override with mobtyper data if available
    if mobtyper.exists():
        lines = mobtyper.read_text().splitlines()
        if len(lines) > 1:
            headers = lines[0].split("\t")
            for line in lines[1:]:
                parts = line.split("\t")
                if len(parts) < len(headers):
                    continue
                row = dict(zip(headers, parts))
                pid = row.get("primary_cluster_id", row.get("id", "unknown"))
                if pid in plasmid_contigs:
                    plasmid_contigs[pid]["replicons"] = row.get("rep_type(s)", plasmid_contigs[pid]["replicons"])
                    plasmid_contigs[pid]["mobility"]  = row.get("predicted_mobility", plasmid_contigs[pid]["mobility"])

    plasmid_list = sorted(plasmid_contigs.values(), key=lambda x: -x.get("size_bp", 0))
    result["plasmids"]      = plasmid_list
    result["plasmid_count"] = len(plasmid_list)

    summary = {"plasmid_count": len(plasmid_list),
               "error": result["error"]}
    db.update_stage(job_id, "mobsuite", "done", json.dumps(summary))
    logger.info("[%s] MOB-Suite: %d plasmid(s) detected.", job_id, len(plasmid_list))
    return result


# ─── Annotation ───────────────────────────────────────────────────────────────

def stage_annotation(job_id: str, fasta: Path, out_dir: Path,
                     sample_name: str) -> Path | None:
    db.update_stage(job_id, "annotation", "running")
    ann_dir = out_dir / "annotation"
    ann_dir.mkdir(exist_ok=True)

    cmd = [
        "bakta",
        "--db",     str(CONDA_BASE / "envs" / "bakta" / "db"),
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

def run_pipeline(job_id: str, fastq_path: Path,
                 fastq_r2: Path | None = None,
                 kraken2_db: Path | None = None,
                 sample_type: str = "auto") -> dict[str, Any]:
    """
    Runs the full pipeline. Returns a dict with all results.
    On failure, sets status to 'failed' in the DB and returns {'error': ...}.
    """
    # Thread-local job_id lets _run() find the cancel event without parameter threading
    _thread_local.job_id = job_id
    register_cancel_event(job_id)

    out_dir = RESULTS_DIR / job_id
    out_dir.mkdir(parents=True, exist_ok=True)

    sample_name = fastq_path.stem.replace(".fastq", "").replace(".fq", "")
    sample_name = safe_sample_name(sample_name)

    results: dict[str, Any] = {
        "job_id":      job_id,
        "sample_name": sample_name,
        "paired_end":  fastq_r2 is not None,
    }

    try:
        db.update_job_status(job_id, "running")

        # ── 1. Read type detection ─────────────────────────────────────────
        db.update_stage(job_id, "detect", "running")
        if fastq_r2 is not None:
            # R1 + R2 provided → definitely Illumina paired-end
            read_type = "illumina"
            db.update_stage(job_id, "detect", "done", "illumina (paired-end, R1+R2 supplied)")
        else:
            read_type = detect_read_type(fastq_path)
            db.update_stage(job_id, "detect", "done", read_type)
        db.update_job_status(job_id, "running", read_type=read_type)
        results["read_type"] = read_type
        logger.info("[%s] Read type: %s%s", job_id, read_type,
                    " (PE)" if fastq_r2 else "")

        # ── 2. Kraken2 ────────────────────────────────────────────────────
        if kraken2_db and kraken2_db.is_dir():
            results["kraken2"] = stage_kraken2(
                job_id, fastq_path, out_dir, kraken2_db,
                fastq_r2=fastq_r2,
            )
            logger.info("[%s] Kraken2 done.", job_id)
        else:
            db.update_stage(job_id, "kraken2", "skipped",
                            "No Kraken2 database specified.")
            results["kraken2"] = None

        # ── Genomic context decision ───────────────────────────────────
        gctx = decide_genomic_context(results["kraken2"], sample_type)
        results["genomic_context"] = gctx
        logger.info("[%s] Genomic context: %s (skip_bacterial=%s)",
                    job_id, gctx["context"], gctx["skip_bacterial"])

        # ── 3. QC ─────────────────────────────────────────────────────────
        results["qc"], filtered_r1, filtered_r2 = stage_qc(
            job_id, fastq_path, out_dir, read_type, fastq_r2=fastq_r2
        )
        logger.info("[%s] QC done.", job_id)

        # ── 3b. Host depletion ────────────────────────────────────────────
        human_pct = gctx.get("human_pct", 0.0)
        if results.get("kraken2") and human_pct >= HOST_DEPLETION_THRESHOLD:
            filtered_r1, filtered_r2, depl = stage_host_depletion(
                job_id, filtered_r1, out_dir, fastq_r2=filtered_r2,
            )
            results["host_depletion"] = depl
            # If enough reads remain after depletion, re-enable bacterial analysis
            if depl.get("depleted") and depl.get("reads_after", 0) >= HOST_DEPLETION_MIN_READS:
                gctx["skip_bacterial"] = False
                gctx["context"]        = "bacterial"
                gctx["reason"] = (
                    f"{depl['host_reads_removed']:,} host reads removed "
                    f"({depl['host_percent']}% Homo sapiens). "
                    f"{depl['reads_after']:,} reads remain — continuing bacterial analysis."
                )
                logger.info("[%s] Post-depletion: bacterial analysis re-enabled.", job_id)
        else:
            db.update_stage(job_id, "host_depletion", "skipped",
                            f"No significant host contamination detected "
                            f"({human_pct:.1f}% < {HOST_DEPLETION_THRESHOLD}% threshold).")
            results["host_depletion"] = {"depleted": False, "host_percent": human_pct}

        # ── 4. Assembly ───────────────────────────────────────────────────
        fasta, asm_qc = stage_assembly(
            job_id, filtered_r1, out_dir, read_type,
            fastq_r2=filtered_r2,
            kraken2=results.get("kraken2"),
        )
        results["assembly_fasta"]    = str(fasta) if fasta else None
        results["assembly_qc_checks"] = asm_qc
        if not fasta:
            raise RuntimeError("Assembly failed — no output FASTA produced.")
        logger.info("[%s] Assembly: %s", job_id, fasta)

        skip_reason = gctx["reason"]
        skip_bact   = gctx["skip_bacterial"]

        # 5. Bandage (always)
        try:
            results["bandage"] = stage_bandage(job_id, fasta, out_dir)
        except Exception as band_err:
            logger.warning("[%s] Bandage skipped: %s", job_id, band_err)
            results["bandage"] = {"error": str(band_err)}

        # 5b. Coverage (always — maps filtered reads back to assembly)
        try:
            results["coverage"] = stage_coverage(
                job_id, fasta, out_dir, read_type,
                fastq_r1=filtered_r1, fastq_r2=filtered_r2,
            )
        except Exception as cov_err:
            logger.warning("[%s] Coverage skipped: %s", job_id, cov_err)
            results["coverage"] = {"error": str(cov_err)}

        # 6. QUAST (always)
        try:
            results["quast"] = stage_quast(job_id, fasta, out_dir)
        except Exception as q_err:
            logger.warning("[%s] QUAST skipped: %s", job_id, q_err)
            results["quast"] = {}

        # 7. CheckM2 (always — assesses completeness/contamination)
        try:
            results["checkm2"] = stage_checkm2(job_id, fasta, out_dir)
        except Exception as cm_err:
            logger.warning("[%s] CheckM2 skipped: %s", job_id, cm_err)
            results["checkm2"] = {"error": str(cm_err)}

        # 8. CheckV (only when viral dominant)
        if gctx["run_checkv"]:
            try:
                results["checkv"] = stage_checkv(job_id, fasta, out_dir)
            except Exception as cv_err:
                logger.warning("[%s] CheckV skipped: %s", job_id, cv_err)
                results["checkv"] = {"error": str(cv_err)}
        else:
            db.update_stage(job_id, "checkv", "skipped",
                            "Not run: viral reads not dominant in this sample.")
            results["checkv"] = None

        # 9. MLST (skip if viral/human dominant)
        if skip_bact:
            db.update_stage(job_id, "mlst", "skipped", skip_reason)
            results["mlst"] = {"skipped": True, "reason": skip_reason}
        else:
            results["mlst"] = stage_mlst(job_id, fasta, out_dir)
            logger.info("[%s] MLST: ST%s", job_id, results["mlst"].get("st"))

        # 10. AMR (skip if viral/human dominant)
        if skip_bact:
            db.update_stage(job_id, "amr", "skipped", skip_reason)
            results["amr"] = {"genes": [], "count": 0, "skipped": True, "reason": skip_reason}
        else:
            results["amr"] = stage_amr(job_id, fasta, out_dir)
            logger.info("[%s] AMR: %d genes found.", job_id, results["amr"].get("count", 0))

        # 11. Abricate (skip if viral/human dominant)
        if skip_bact:
            db.update_stage(job_id, "abricate", "skipped", skip_reason)
            results["abricate"] = {"genes": [], "count": 0, "skipped": True, "reason": skip_reason}
        else:
            try:
                results["abricate"] = stage_abricate(
                    job_id, fasta, out_dir, results.get("kraken2")
                )
                logger.info("[%s] Abricate: %d hits.", job_id,
                            results["abricate"].get("count", 0))
            except Exception as abr_err:
                logger.warning("[%s] Abricate error: %s", job_id, abr_err)
                results["abricate"] = {"genes": [], "count": 0, "ran_vfdb": False, "top_genus": ""}

        # 12. MOB-Suite (skip if viral/human dominant)
        if skip_bact:
            db.update_stage(job_id, "mobsuite", "skipped", skip_reason)
            results["mobsuite"] = {"plasmids": [], "plasmid_count": 0,
                                   "skipped": True, "reason": skip_reason}
        else:
            try:
                results["mobsuite"] = stage_mobsuite(job_id, fasta, out_dir)
                logger.info("[%s] MOB-Suite: %d plasmid(s).", job_id,
                            results["mobsuite"].get("plasmid_count", 0))
            except Exception as mob_err:
                logger.warning("[%s] MOB-Suite error: %s", job_id, mob_err)
                results["mobsuite"] = {"plasmids": [], "plasmid_count": 0,
                                       "error": str(mob_err)}

        # 13. Serotyping (auto-detected from Kraken2 genus — skip if not applicable)
        try:
            results["serotyping"] = stage_serotyping(
                job_id, fasta, out_dir, results.get("kraken2")
            )
        except Exception as sero_err:
            logger.warning("[%s] Serotyping skipped: %s", job_id, sero_err)
            results["serotyping"] = {"error": str(sero_err)}

        # 14. Annotation (optional — continue even if it fails)
        try:
            results["annotation"] = str(
                stage_annotation(job_id, fasta, out_dir, sample_name)
            )
        except Exception as ann_err:
            logger.warning("[%s] Annotation skipped: %s", job_id, ann_err)
            results["annotation"] = None

        db.update_job_status(job_id, "ai_pending")
        return results

    except PipelineCancelledError:
        logger.info("[%s] Pipeline cancelled by user.", job_id)
        db.update_job_status(job_id, "cancelled", error="Cancelled by user.")
        results["error"] = "Cancelled by user."
        return results

    except Exception as exc:
        error_msg = str(exc)
        logger.error("[%s] Pipeline error: %s", job_id, error_msg, exc_info=True)
        db.update_job_status(job_id, "failed", error=error_msg)
        results["error"] = error_msg
        return results

    finally:
        clear_cancel_event(job_id)
        _thread_local.job_id = None


# ─── QC (kept here for import by stage_assembly caller) ───────────────────────

def stage_qc(job_id: str, fastq: Path, out_dir: Path,
             read_type: str,
             fastq_r2: Path | None = None) -> tuple[dict[str, Any], Path, Path | None]:
    """
    QC stage:
    - MinION  → NanoPlot (stats) + Filtlong (filter by length/quality)
    - Illumina SE → fastp single-end
    - Illumina PE → fastp paired-end
    Returns (metrics_dict, filtered_r1_path, filtered_r2_path_or_None)
    """
    db.update_stage(job_id, "qc", "running")
    qc_dir = out_dir / "qc"
    qc_dir.mkdir(exist_ok=True)
    result: dict[str, Any] = {}

    if read_type == "minion":
        # ── NanoPlot stats ─────────────────────────────────────────────────
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

        # ── Filtlong — filter reads for downstream assembly ─────────────────
        filt_fastq = qc_dir / "filtered.fastq.gz"
        fl_cmd = [
            "filtlong",
            "--min_length",    "1000",   # discard reads < 1 kb
            "--min_mean_q",    "8",      # discard reads with mean Phred < Q8
            "--target_bases",  "500000000",  # keep top 500 Mb (≈100x for 5 Mb genome)
            sanitize_path(fastq),
        ]
        # Filtlong writes to stdout — capture and gzip-pipe manually
        import os, gzip as _gz
        fl_env = os.environ.copy()
        fl_env["PATH"] = (str(CONDA_BASE / "envs" / "analiz" / "bin") + ":"
                          + fl_env["PATH"])
        try:
            proc = subprocess.Popen(
                fl_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                cwd=str(qc_dir), env=fl_env,
            )
            with _gz.open(filt_fastq, "wb") as fout:
                fout.write(proc.stdout.read())
            proc.wait(timeout=3600)
            if filt_fastq.exists() and filt_fastq.stat().st_size > 0:
                result["filtlong"] = "applied (min_length=1000, min_mean_q=8)"
                logger.info("[%s] Filtlong done: %s bytes", job_id, filt_fastq.stat().st_size)
                filtered_fastq = filt_fastq
            else:
                logger.warning("[%s] Filtlong produced no output — using raw reads", job_id)
                result["filtlong"] = "no output — using raw reads"
                filtered_fastq = fastq
        except Exception as fl_err:
            logger.warning("[%s] Filtlong failed: %s — using raw reads", job_id, fl_err)
            result["filtlong"] = f"skipped ({fl_err})"
            filtered_fastq = fastq

        db.update_stage(job_id, "qc", "done", json.dumps(result))
        return result, filtered_fastq, None

    else:  # Illumina (SE or PE)
        json_out  = qc_dir / "fastp.json"
        filt_r1   = qc_dir / "filtered_R1.fastq.gz"
        filt_r2   = qc_dir / "filtered_R2.fastq.gz"

        if fastq_r2:   # Paired-end
            cmd = [
                "fastp",
                "-i", sanitize_path(fastq),
                "-I", sanitize_path(fastq_r2),
                "-o", str(filt_r1),
                "-O", str(filt_r2),
                "--json", str(json_out),
                "--html", str(qc_dir / "fastp.html"),
                "--thread", str(MAX_THREADS),
                "--detect_adapter_for_pe",
                "-g", "-x",
                "--length_required", "50",
                "-q", "20",
            ]
        else:           # Single-end
            filt_r1 = qc_dir / "filtered.fastq.gz"
            filt_r2 = None
            cmd = [
                "fastp",
                "-i", sanitize_path(fastq),
                "-o", str(filt_r1),
                "--json", str(json_out),
                "--html", str(qc_dir / "fastp.html"),
                "--thread", str(MAX_THREADS),
                "-g", "-x",
                "--length_required", "50",
                "-q", "20",
            ]

        rc, out, err = _run(cmd, qc_dir)
        if json_out.exists():
            data = json.loads(json_out.read_text())
            af = data.get("summary", {}).get("after_filtering", {})
            result = {
                "total_reads":      af.get("total_reads", 0),
                "total_bases":      af.get("total_bases", 0),
                "q30_rate":         round(af.get("q30_rate", 0) * 100, 1),
                "mean_read_length": af.get("read1_mean_length", 0),
                "paired_end":       fastq_r2 is not None,
            }

        # Fallback: if fastp produced no output, use originals
        if not filt_r1.exists():
            filt_r1 = fastq
            filt_r2 = fastq_r2

        db.update_stage(job_id, "qc", "done", json.dumps(result))
        return result, filt_r1, filt_r2 if fastq_r2 else None
