"""
Multi-sample comparison pipeline.
Collects assemblies from completed jobs, runs Mashtree for a
distance-based phylogenetic tree, and aggregates per-sample metadata
(MLST, AMR, serotyping, Kraken2, CheckM2, coverage) into a single report.
"""
import json
import logging
import os
import re
import shutil
import subprocess
import urllib.request
import zipfile
from pathlib import Path
from typing import Any

import database as db
from config import MAX_THREADS, CONDA_BASE, RESULTS_DIR

logger = logging.getLogger("comparison")

COMPARISONS_DIR = RESULTS_DIR.parent / "comparisons"
COMPARISONS_DIR.mkdir(parents=True, exist_ok=True)

GTDB_SKETCH = Path("/home/analysis/databases_all/mashdb/GTDBSketch_20231003.msh")
NCBI_DATASETS_API = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc}/download"


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


# ─── NCBI reference genome fetching ──────────────────────────────────────────

def _mash_screen_top5(fasta: Path, out_dir: Path) -> list[dict]:
    """
    Runs `mash screen -w` against the GTDB sketch.
    Returns up to 5 top hits sorted by identity as:
      [{"accession": "GCF_...", "identity": 0.999, "species": "...", "label": "..."}]
    """
    if not GTDB_SKETCH.exists():
        logger.warning("GTDB sketch not found at %s", GTDB_SKETCH)
        return []

    env = os.environ.copy()
    env["PATH"] = str(CONDA_BASE / "envs" / "mashtree" / "bin") + ":" + env["PATH"]

    screen_out = out_dir / "mash_screen.tsv"
    try:
        result = subprocess.run(
            ["mash", "screen", "-w", "-p", str(MAX_THREADS),
             str(GTDB_SKETCH), str(fasta)],
            capture_output=True, text=True, timeout=600,
            env=env, cwd=str(out_dir),
        )
        screen_out.write_text(result.stdout)
    except Exception as e:
        logger.warning("mash screen failed: %s", e)
        return []

    hits = []
    for line in result.stdout.splitlines():
        parts = line.strip().split("\t")
        if len(parts) < 5:
            continue
        try:
            identity = float(parts[0])
            query_id = parts[4].strip()    # GCF_XXXXXXX.X
            comment  = parts[5].strip() if len(parts) > 5 else ""
            # Extract species from GTDB lineage string (last ;s__ field)
            species = ""
            m = re.search(r's__([^;]+)', comment)
            if m:
                species = m.group(1).strip().strip('"')
            accession = query_id.split()[0]   # strip any trailing whitespace
            hits.append({
                "accession": accession,
                "identity":  round(identity, 5),
                "species":   species,
                "label":     f"{accession} ({species})" if species else accession,
            })
        except (ValueError, IndexError):
            continue

    hits.sort(key=lambda h: h["identity"], reverse=True)
    # Deduplicate by accession
    seen = set()
    deduped = []
    for h in hits:
        if h["accession"] not in seen:
            seen.add(h["accession"])
            deduped.append(h)
    return deduped[:5]


def _download_ncbi_fasta(accession: str, dest_dir: Path) -> Path | None:
    """
    Downloads the genomic FASTA for a GCF/GCA accession via the NCBI datasets API.
    Returns the path to the extracted .fna file, or None on failure.
    """
    url = (f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/"
           f"{accession}/download?include_annotation_type=GENOME_FASTA"
           f"&filename={accession}.zip")
    zip_path = dest_dir / f"{accession}.zip"
    fna_path = dest_dir / f"{accession}.fna"

    try:
        logger.info("Downloading %s from NCBI…", accession)
        req = urllib.request.Request(url, headers={"User-Agent": "LycianWay/1.0"})
        with urllib.request.urlopen(req, timeout=120) as resp, \
             open(zip_path, "wb") as fh:
            fh.write(resp.read())

        # Extract the genomic FASTA from the zip
        with zipfile.ZipFile(zip_path, "r") as zf:
            fna_files = [n for n in zf.namelist()
                         if n.endswith("_genomic.fna") or n.endswith(".fna")]
            if not fna_files:
                logger.warning("No .fna in zip for %s", accession)
                return None
            # Write merged content (may be multiple chromosomes/contigs)
            with open(fna_path, "wb") as out:
                for fname in fna_files:
                    out.write(zf.read(fname))

        zip_path.unlink(missing_ok=True)
        logger.info("Downloaded %s → %s (%.1f MB)", accession, fna_path,
                    fna_path.stat().st_size / 1e6)
        return fna_path

    except Exception as e:
        logger.warning("Failed to download %s: %s", accession, e)
        zip_path.unlink(missing_ok=True)
        return None


def fetch_ncbi_references(fasta: Path, out_dir: Path) -> list[dict]:
    """
    Finds the top-5 closest NCBI reference genomes via Mash screen (GTDB),
    downloads their FASTAs, and returns sample-metadata dicts ready for
    inclusion in the comparison.  Downloaded files are tracked in
    out_dir/ncbi_downloads/ for later cleanup.
    """
    dl_dir = out_dir / "ncbi_downloads"
    dl_dir.mkdir(exist_ok=True)

    hits = _mash_screen_top5(fasta, out_dir)
    if not hits:
        return []

    ref_samples = []
    for hit in hits:
        acc   = hit["accession"]
        fna   = _download_ncbi_fasta(acc, dl_dir)
        if not fna:
            continue
        ref_samples.append({
            "job_id":       acc,
            "sample_name":  hit["label"][:45],
            "read_type":    "reference",
            "organism":     hit["species"],
            "mlst_scheme":  "—", "mlst_st":      "—",
            "serotype":     "—",
            "completeness": "—", "contamination": "—",
            "mean_depth":   "—", "breadth_1x":    "—",
            "n50":          "—", "total_length":  "—",
            "amr_genes":    [],  "amr_count":     0,
            "pf_replicons": [],
            "fasta_path":   str(fna),
            "is_reference": True,
            "ncbi_identity": hit["identity"],
        })
        logger.info("NCBI ref added: %s (identity %.4f)", acc, hit["identity"])

    return ref_samples


def cleanup_ncbi_downloads(out_dir: Path) -> None:
    """Deletes downloaded NCBI FASTA files after the report is generated."""
    dl_dir = out_dir / "ncbi_downloads"
    if dl_dir.exists():
        shutil.rmtree(dl_dir, ignore_errors=True)
        logger.info("Cleaned up NCBI downloads in %s", dl_dir)


# ─── Main comparison runner ───────────────────────────────────────────────────

def run_comparison(comp_id: str, job_ids: list[str],
                   include_ncbi: bool = False) -> None:
    """
    Runs the full comparison pipeline for the given job IDs.
    If include_ncbi=True, also fetches the top-5 closest NCBI reference
    genomes for the first sample's assembly and adds them to the tree.
    Downloaded reference FASTAs are deleted after the report is written.
    Updates the comparisons table with status and report path.
    """
    out_dir = COMPARISONS_DIR / comp_id
    out_dir.mkdir(parents=True, exist_ok=True)

    db.update_comparison(comp_id, "running")
    logger.info("[cmp:%s] Starting comparison for %d jobs (ncbi=%s)",
                comp_id, len(job_ids), include_ncbi)

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

        # ── 2. NCBI reference genomes (optional) ───────────────────────────
        ncbi_refs: list[dict] = []
        if include_ncbi and GTDB_SKETCH.exists():
            # Use the first sample's assembly as query
            first_fasta = next(
                (Path(s["fasta_path"]) for s in samples if s["fasta_path"]), None
            )
            if first_fasta:
                logger.info("[cmp:%s] Fetching NCBI references for %s",
                            comp_id, first_fasta.name)
                ncbi_refs = fetch_ncbi_references(first_fasta, out_dir)
                logger.info("[cmp:%s] Got %d NCBI reference(s)", comp_id, len(ncbi_refs))

        all_samples = samples + ncbi_refs

        # ── 3. Mashtree + pairwise distances ───────────────────────────────
        fasta_paths  = [Path(s["fasta_path"]) for s in all_samples if s["fasta_path"]]
        fasta_labels = [s["sample_name"]      for s in all_samples if s["fasta_path"]]

        newick    = None
        distances = []
        if len(fasta_paths) >= 2:
            logger.info("[cmp:%s] Running Mashtree on %d assemblies "
                        "(%d user + %d NCBI refs)",
                        comp_id, len(fasta_paths), len(samples), len(ncbi_refs))
            newick    = _run_mashtree(fasta_paths, out_dir, fasta_labels)
            distances = _mash_distances(fasta_paths, fasta_labels, out_dir)
        else:
            logger.warning("[cmp:%s] Not enough FASTAs for tree — metadata only", comp_id)

        # ── 4. Generate report ─────────────────────────────────────────────
        from comparison_report import generate_comparison_report
        html = generate_comparison_report(comp_id, all_samples, newick, distances)

        report_path = out_dir / "comparison_report.html"
        report_path.write_text(html, encoding="utf-8")

        # ── 5. Delete downloaded NCBI FASTAs ──────────────────────────────
        if ncbi_refs:
            cleanup_ncbi_downloads(out_dir)

        db.update_comparison(comp_id, "completed", report_path=str(report_path))
        logger.info("[cmp:%s] Comparison complete: %s", comp_id, report_path)

    except Exception as exc:
        msg = str(exc)
        logger.error("[cmp:%s] Comparison failed: %s", comp_id, msg, exc_info=True)
        # Still clean up any partial downloads
        cleanup_ncbi_downloads(out_dir)
        db.update_comparison(comp_id, "failed", error=msg)
