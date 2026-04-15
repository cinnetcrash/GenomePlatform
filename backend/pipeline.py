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
from pathlib import Path
from typing import Any

import database as db
from config import (
    ASSEMBLY_QC, AUTOCYCLER_BIN, CONDA_BASE,
    MAX_THREADS, RESULTS_DIR, SCRIPTS_DIR,
)
from security import safe_sample_name, sanitize_path

logger = logging.getLogger("pipeline")


# ─── Helpers ──────────────────────────────────────────────────────────────────

def _conda_run(env_name: str, cmd: list[str], cwd: Path,
               timeout: int = 7200) -> subprocess.CompletedProcess:
    """Runs a command inside the specified conda environment (shell=False)."""
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


def stage_kraken2(job_id: str, fastq: Path, out_dir: Path,
                  db_path: Path) -> dict[str, Any]:
    """Taxonomic classification with Kraken2."""
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


# ─── Assembly (Illumina + MinION / Autocycler) ────────────────────────────────

def _run_flye(fastq: Path, out_dir: Path, mode: str, timeout: int = 10800) -> Path | None:
    """Runs Flye with the given mode flag; returns assembly FASTA path or None."""
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        "flye", mode, sanitize_path(fastq),
        "--out-dir", str(out_dir),
        "--threads", str(MAX_THREADS),
        "--genome-size", "5m",
    ]
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
                       "--threads", str(MAX_THREADS)],
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
              "--threads", str(MAX_THREADS)], ac_dir, timeout=1800)
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


def stage_assembly(job_id: str, fastq: Path, out_dir: Path,
                   read_type: str) -> tuple[Path | None, list]:
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
        logger.info("[%s] Assembly: running Flye --nano-raw", job_id)
        flye_raw = _run_flye(fastq, asm_dir / "flye_raw", "--nano-raw")

        logger.info("[%s] Assembly: running Flye --nano-hq", job_id)
        flye_hq  = _run_flye(fastq, asm_dir / "flye_hq",  "--nano-hq")

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
        cmd = [
            "shovill",
            "--R1", sanitize_path(fastq),
            "--outdir", str(asm_dir),
            "--cpus", str(MAX_THREADS),
            "--force",
        ]
        rc, out, err = _run(cmd, asm_dir, timeout=10800, env_name="shovill")
        candidates = list(asm_dir.glob("contigs.fa")) + list(asm_dir.glob("assembly.fasta"))
        if candidates:
            fasta  = candidates[0]
            method = "shovill"
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

    result = {
        "genes":      all_genes,
        "count":      len(all_genes),
        "ran_vfdb":   run_vfdb,
        "top_genus":  top_genus,
    }
    db.update_stage(job_id, "abricate", "done",
                    json.dumps({"count": len(all_genes), "ran_vfdb": run_vfdb,
                                "top_genus": top_genus}))
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
                 kraken2_db: Path | None = None) -> dict[str, Any]:
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

        # 2. Kraken2 (optional — runs on raw reads if DB provided)
        if kraken2_db and kraken2_db.is_dir():
            results["kraken2"] = stage_kraken2(job_id, fastq_path, out_dir, kraken2_db)
            logger.info("[%s] Kraken2 done.", job_id)
        else:
            db.update_stage(job_id, "kraken2", "skipped",
                            "No Kraken2 database specified.")
            results["kraken2"] = None

        # 3. QC
        results["qc"] = stage_qc(job_id, fastq_path, out_dir, read_type)
        logger.info("[%s] QC done.", job_id)

        # 4. Assembly
        fasta, asm_qc = stage_assembly(job_id, fastq_path, out_dir, read_type)
        results["assembly_fasta"]    = str(fasta) if fasta else None
        results["assembly_qc_checks"] = asm_qc
        if not fasta:
            raise RuntimeError("Assembly failed — no output FASTA produced.")
        logger.info("[%s] Assembly: %s", job_id, fasta)

        # 5. MLST
        results["mlst"] = stage_mlst(job_id, fasta, out_dir)
        logger.info("[%s] MLST: ST%s", job_id, results["mlst"].get("st"))

        # 6. AMR
        results["amr"] = stage_amr(job_id, fasta, out_dir)
        logger.info("[%s] AMR: %d genes found.", job_id, results["amr"].get("count", 0))

        # 7. Abricate (CARD always; VFDB if Kraken2 top genus is in VFDB)
        try:
            results["abricate"] = stage_abricate(
                job_id, fasta, out_dir, results.get("kraken2")
            )
            logger.info("[%s] Abricate: %d hits.", job_id, results["abricate"].get("count", 0))
        except Exception as abr_err:
            logger.warning("[%s] Abricate skipped: %s", job_id, abr_err)
            results["abricate"] = {"genes": [], "count": 0, "ran_vfdb": False, "top_genus": ""}

        # 8. MOB-Suite (plasmid detection — optional)
        try:
            results["mobsuite"] = stage_mobsuite(job_id, fasta, out_dir)
            logger.info("[%s] MOB-Suite: %d plasmid(s).", job_id,
                        results["mobsuite"].get("plasmid_count", 0))
        except Exception as mob_err:
            logger.warning("[%s] MOB-Suite skipped: %s", job_id, mob_err)
            results["mobsuite"] = {"plasmids": [], "plasmid_count": 0, "error": str(mob_err)}

        # 9. Annotation (optional — continue even if it fails)
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


# ─── QC (kept here for import by stage_assembly caller) ───────────────────────

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
            af = data.get("summary", {}).get("after_filtering", {})
            result = {
                "total_reads":      af.get("total_reads", 0),
                "total_bases":      af.get("total_bases", 0),
                "q30_rate":         round(af.get("q30_rate", 0) * 100, 1),
                "mean_read_length": af.get("read1_mean_length", 0),
            }

    db.update_stage(job_id, "qc", "done", json.dumps(result))
    return result
