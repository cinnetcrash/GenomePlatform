"""
Multi-sample comparison pipeline.
Collects assemblies from completed jobs, runs Mashtree for a
distance-based phylogenetic tree, and aggregates per-sample metadata
(MLST, AMR, serotyping, Kraken2, CheckM2, coverage) into a single report.
"""
import json
import logging
import shutil
import subprocess
from pathlib import Path
from typing import Any

import database as db
from config import MAX_THREADS, CONDA_BASE, RESULTS_DIR
from security import generate_job_id

logger = logging.getLogger("comparison")

COMPARISONS_DIR = RESULTS_DIR.parent / "comparisons"
COMPARISONS_DIR.mkdir(parents=True, exist_ok=True)


# ─── Data extraction helpers ─────────────────────────────────────────────────

def _stage_json(stages: dict, key: str) -> dict:
    """Safely parse a stage's detail JSON."""
    raw = stages.get(key, {}).get("detail", "")
    if not raw:
        return {}
    try:
        return json.loads(raw)
    except Exception:
        return {}


def _collect_sample_data(job: dict) -> dict[str, Any]:
    """Extracts all relevant metadata from a completed job record."""
    stages  = json.loads(job.get("stages", "{}"))
    kraken  = _stage_json(stages, "kraken2")
    mlst    = _stage_json(stages, "mlst")
    amr     = _stage_json(stages, "amr")
    abricate= _stage_json(stages, "abricate")
    checkm2 = _stage_json(stages, "checkm2")
    coverage= _stage_json(stages, "coverage")
    sero    = _stage_json(stages, "serotyping")
    asm     = _stage_json(stages, "assembly")
    quast   = _stage_json(stages, "quast")

    top_organism = "—"
    if kraken.get("top_taxa"):
        top_organism = kraken["top_taxa"][0]["name"]

    amr_genes = [g["gene"] for g in amr.get("genes", [])]
    pf_genes  = [g["gene"] for g in abricate.get("plasmidfinder", [])]

    serotype = "—"
    if sero and not sero.get("skipped") and not sero.get("error"):
        serotype = (sero.get("serotype") or sero.get("serovar")
                    or sero.get("st") or "—")

    fasta_path = asm.get("fasta_path")

    return {
        "job_id":       job["id"],
        "sample_name":  job["filename"].replace(".fastq.gz","").replace(".fq.gz","")
                                        .replace(".fastq","").replace(".fq",""),
        "read_type":    job.get("read_type", "?"),
        "organism":     top_organism,
        "mlst_scheme":  mlst.get("scheme", "—"),
        "mlst_st":      mlst.get("st", "—"),
        "serotype":     serotype,
        "completeness": checkm2.get("completeness", "—"),
        "contamination":checkm2.get("contamination", "—"),
        "mean_depth":   coverage.get("mean_depth", "—"),
        "breadth_1x":   coverage.get("breadth_1x_pct", "—"),
        "n50":          quast.get("n50", asm.get("n50_bp", "—")),
        "total_length": quast.get("total_length", asm.get("total_length_bp", "—")),
        "amr_genes":    amr_genes,
        "amr_count":    len(amr_genes),
        "pf_replicons": pf_genes,
        "fasta_path":   fasta_path,
    }


# ─── Mashtree ────────────────────────────────────────────────────────────────

def _run_mashtree(fasta_paths: list[Path], out_dir: Path,
                  labels: list[str]) -> str | None:
    """
    Runs Mashtree on the given FASTA files.
    Returns the Newick tree string or None on failure.
    """
    if len(fasta_paths) < 2:
        return None

    # Stage FASTAs with readable names (Mashtree uses filename as leaf label)
    stage_dir = out_dir / "fastas"
    stage_dir.mkdir(exist_ok=True)
    staged: list[Path] = []
    for fasta, label in zip(fasta_paths, labels):
        safe_label = label.replace(" ", "_").replace("/", "_")[:40]
        dest = stage_dir / f"{safe_label}.fasta"
        shutil.copy2(fasta, dest)
        staged.append(dest)

    import os
    env = os.environ.copy()
    env["PATH"] = str(CONDA_BASE / "envs" / "mashtree" / "bin") + ":" + env["PATH"]

    nwk_file = out_dir / "tree.nwk"
    cmd = [
        "mashtree",
        "--numcpus", str(MAX_THREADS),
        "--outfile",  str(nwk_file),
    ] + [str(p) for p in staged]

    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True,
            timeout=1800, env=env, cwd=str(out_dir),
        )
        if nwk_file.exists() and nwk_file.stat().st_size > 0:
            return nwk_file.read_text().strip()
        # Some mashtree versions write to stdout
        if result.stdout.strip().startswith("("):
            nwk_file.write_text(result.stdout.strip())
            return result.stdout.strip()
        logger.warning("Mashtree failed: %s", result.stderr[:300])
    except Exception as e:
        logger.warning("Mashtree exception: %s", e)
    return None


# ─── Mash pairwise distances ─────────────────────────────────────────────────

def _mash_distances(fasta_paths: list[Path], labels: list[str],
                    out_dir: Path) -> list[dict]:
    """
    Computes pairwise Mash distances.
    Returns list of {sample_a, sample_b, distance, p_value}.
    """
    import os
    env = os.environ.copy()
    env["PATH"] = str(CONDA_BASE / "envs" / "mashtree" / "bin") + ":" + env["PATH"]

    sketch_file = out_dir / "sketch.msh"
    fasta_strs  = [str(p) for p in fasta_paths]

    try:
        subprocess.run(
            ["mash", "sketch", "-o", str(out_dir / "sketch")] + fasta_strs,
            capture_output=True, env=env, timeout=300, cwd=str(out_dir),
        )
        dist_result = subprocess.run(
            ["mash", "dist", str(sketch_file)] + fasta_strs,
            capture_output=True, text=True, env=env,
            timeout=300, cwd=str(out_dir),
        )
    except Exception as e:
        logger.warning("Mash distance failed: %s", e)
        return []

    # Build label lookup: filename stem → label
    stem_to_label = {p.stem: lbl for p, lbl in zip(fasta_paths, labels)}

    distances = []
    for line in dist_result.stdout.splitlines():
        parts = line.strip().split("\t")
        if len(parts) < 3:
            continue
        try:
            name_a = Path(parts[0]).stem
            name_b = Path(parts[1]).stem
            dist   = float(parts[2])
            pval   = float(parts[3]) if len(parts) > 3 else None
            distances.append({
                "sample_a": stem_to_label.get(name_a, name_a),
                "sample_b": stem_to_label.get(name_b, name_b),
                "distance": round(dist, 6),
                "p_value":  round(pval, 4) if pval is not None else None,
                "identity_pct": round((1 - dist) * 100, 2),
            })
        except (ValueError, IndexError):
            continue
    return distances


# ─── Main comparison runner ───────────────────────────────────────────────────

def run_comparison(comp_id: str, job_ids: list[str]) -> None:
    """
    Runs the full comparison pipeline for the given job IDs.
    Updates the comparisons table with status and report path.
    """
    out_dir = COMPARISONS_DIR / comp_id
    out_dir.mkdir(parents=True, exist_ok=True)

    db.update_comparison(comp_id, "running")
    logger.info("[cmp:%s] Starting comparison for %d jobs", comp_id, len(job_ids))

    try:
        # ── 1. Load job data ───────────────────────────────────────────────
        samples = []
        for jid in job_ids:
            job = db.get_job(jid)
            if not job or job["status"] not in ("completed",):
                logger.warning("[cmp:%s] Job %s not completed, skipping", comp_id, jid)
                continue
            data = _collect_sample_data(job)
            if not data["fasta_path"] or not Path(data["fasta_path"]).exists():
                logger.warning("[cmp:%s] No assembly FASTA for job %s", comp_id, jid)
                data["fasta_path"] = None
            samples.append(data)

        if len(samples) < 2:
            raise RuntimeError(
                "Need at least 2 completed jobs with assemblies for comparison."
            )

        # ── 2. Mashtree + pairwise distances ───────────────────────────────
        fasta_paths  = [Path(s["fasta_path"]) for s in samples if s["fasta_path"]]
        fasta_labels = [s["sample_name"]      for s in samples if s["fasta_path"]]

        newick    = None
        distances = []
        if len(fasta_paths) >= 2:
            logger.info("[cmp:%s] Running Mashtree on %d assemblies", comp_id, len(fasta_paths))
            newick    = _run_mashtree(fasta_paths, out_dir, fasta_labels)
            distances = _mash_distances(fasta_paths, fasta_labels, out_dir)
        else:
            logger.warning("[cmp:%s] Not enough FASTAs for tree — metadata only", comp_id)

        # ── 3. Generate report ─────────────────────────────────────────────
        from comparison_report import generate_comparison_report
        html = generate_comparison_report(comp_id, samples, newick, distances)

        report_path = out_dir / "comparison_report.html"
        report_path.write_text(html, encoding="utf-8")

        db.update_comparison(comp_id, "completed", report_path=str(report_path))
        logger.info("[cmp:%s] Comparison complete: %s", comp_id, report_path)

    except Exception as exc:
        msg = str(exc)
        logger.error("[cmp:%s] Comparison failed: %s", comp_id, msg, exc_info=True)
        db.update_comparison(comp_id, "failed", error=msg)
