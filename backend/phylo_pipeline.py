"""
Global Phylogeny Pipeline
Stages: organism detection → NCBI fetch → Panaroo → IQ-TREE → R annotation → report
"""
import csv
import json
import logging
import os
import subprocess
import urllib.parse
import urllib.request
import zipfile
from collections import Counter
from pathlib import Path
from typing import Any

import database as db
from config import CONDA_BASE, RESULTS_DIR
from errors import (
    PipelineError, ERROR_LABELS,
    PANAROO_FAILED, IQTREE_FAILED, R_ANNOTATION_FAILED,
    NO_ASSEMBLIES, NCBI_FETCH_FAILED,
)
from security import generate_job_id
import scheduler

logger = logging.getLogger("phylo")

PHYLO_DIR = RESULTS_DIR.parent / "phylo"
PHYLO_DIR.mkdir(parents=True, exist_ok=True)

NCBI_DATASETS_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2"
_R_SCRIPT = Path(__file__).parent / "r_scripts" / "annotate_tree.R"


# ─── Stage 1: Organism Detection ──────────────────────────────────────────────

def detect_organism(jobs: list[dict]) -> str | None:
    """Majority vote on top Kraken2 taxon across selected jobs."""
    names = []
    for job in jobs:
        try:
            stages = json.loads(job.get("stages", "{}"))
            k2_detail = stages.get("kraken2", {}).get("detail", "")
            if not k2_detail:
                continue
            k2 = json.loads(k2_detail)
            top = k2.get("top_taxa", [])
            if top:
                names.append(top[0]["name"])
        except (json.JSONDecodeError, KeyError, IndexError):
            continue
    if not names:
        return None
    return Counter(names).most_common(1)[0][0]


# ─── Stage 2: NCBI Reference Genomes (1 per country) ─────────────────────────

def _ncbi_get(url: str) -> dict:
    req = urllib.request.Request(url, headers={"User-Agent": "LycianWay/1.0"})
    with urllib.request.urlopen(req, timeout=30) as resp:
        return json.loads(resp.read())


def _best_assembly_level(assemblies: list[dict]) -> dict:
    order = {"Complete Genome": 0, "Chromosome": 1, "Scaffold": 2, "Contig": 3}
    return min(assemblies, key=lambda a: order.get(
        a.get("assembly_info", {}).get("assembly_level", "Contig"), 4
    ))


def fetch_ncbi_references(organism: str, out_dir: Path,
                          max_countries: int = 50) -> list[dict]:
    """
    Queries NCBI Datasets API for `organism`, groups by country,
    picks 1 best assembly per country, downloads FASTA+GFF.
    """
    dl_dir = out_dir / "ncbi_refs"
    dl_dir.mkdir(exist_ok=True)

    taxon_q = urllib.parse.quote(organism)
    search_url = (
        f"{NCBI_DATASETS_BASE}/genome/taxon/{taxon_q}/dataset_report"
        f"?filters.assembly_source=all&page_size=500&returned_content=COMPLETE"
    )
    try:
        data = _ncbi_get(search_url)
    except Exception as e:
        logger.warning("NCBI search failed for %s: %s", organism, e)
        return []

    reports = data.get("reports", [])
    if not reports:
        logger.warning("No NCBI assemblies found for %s", organism)
        return []

    by_country: dict[str, list[dict]] = {}
    for r in reports:
        geo = (r.get("assembly_info", {})
                .get("biosample", {})
                .get("attributes", []))
        country = next(
            (a["value"] for a in geo if a.get("name") == "geo_loc_name"), None
        )
        if not country:
            continue
        country = country.split(":")[0].strip()
        by_country.setdefault(country, []).append(r)

    logger.info("Found %d countries for %s", len(by_country), organism)

    refs = []
    for country, assemblies in list(by_country.items())[:max_countries]:
        best = _best_assembly_level(assemblies)
        accession = best.get("accession", "")
        if not accession:
            continue
        fasta_path, gff_path = _download_ncbi_assembly(accession, dl_dir)
        if fasta_path is None:
            continue
        refs.append({
            "sample_name":  f"REF_{country.replace(' ', '_')}",
            "country":      country,
            "accession":    accession,
            "fasta_path":   str(fasta_path),
            "gff_path":     str(gff_path) if gff_path else None,
            "is_reference": True,
        })
        logger.info("Downloaded ref: %s (%s)", country, accession)

    return refs


def _download_ncbi_assembly(accession: str, dest_dir: Path) -> tuple[Path | None, Path | None]:
    url = (
        f"{NCBI_DATASETS_BASE}/genome/accession/{accession}/download"
        f"?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF"
        f"&filename={accession}.zip"
    )
    zip_path   = dest_dir / f"{accession}.zip"
    fasta_out  = dest_dir / f"{accession}.fna"
    gff_out    = dest_dir / f"{accession}.gff"

    try:
        req = urllib.request.Request(url, headers={"User-Agent": "LycianWay/1.0"})
        with urllib.request.urlopen(req, timeout=120) as resp, \
             open(zip_path, "wb") as fh:
            fh.write(resp.read())

        with zipfile.ZipFile(zip_path, "r") as zf:
            fnas = [n for n in zf.namelist() if n.endswith("_genomic.fna")]
            gffs = [n for n in zf.namelist() if n.endswith("_genomic.gff")]
            if not fnas:
                return None, None
            with open(fasta_out, "wb") as fh:
                for f in fnas:
                    fh.write(zf.read(f))
            gff_path = None
            if gffs:
                with open(gff_out, "wb") as fh:
                    for g in gffs:
                        fh.write(zf.read(g))
                gff_path = gff_out

        zip_path.unlink(missing_ok=True)
        return fasta_out, gff_path

    except Exception as e:
        logger.warning("Download failed for %s: %s", accession, e)
        zip_path.unlink(missing_ok=True)
        return None, None


# ─── Stage 3: Panaroo Pan-Genome ──────────────────────────────────────────────

def run_panaroo(sample_gffs: list[Path], out_dir: Path) -> Path:
    panaroo_out = out_dir / "panaroo"
    panaroo_out.mkdir(exist_ok=True)

    env = os.environ.copy()
    env["PATH"] = str(CONDA_BASE / "envs" / "panaroo" / "bin") + ":" + env["PATH"]

    cmd = [
        "panaroo",
        "-i", *[str(g) for g in sample_gffs],
        "-o", str(panaroo_out),
        "--clean-mode", "strict",
        "-a", "core",
        "--core_threshold", "0.98",
        "-t", str(scheduler.threads_for_stage()),
    ]
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True,
            timeout=7200, env=env, cwd=str(out_dir),
        )
    except subprocess.TimeoutExpired:
        raise PipelineError(
            code=PANAROO_FAILED, message=ERROR_LABELS[PANAROO_FAILED],
            detail="Panaroo timed out after 2 hours", stage="panaroo",
        )

    if result.returncode != 0:
        raise PipelineError(
            code=PANAROO_FAILED, message=ERROR_LABELS[PANAROO_FAILED],
            detail=result.stderr[:500], stage="panaroo",
        )

    aln = panaroo_out / "core_gene_alignment.aln"
    if not aln.exists():
        raise PipelineError(
            code=PANAROO_FAILED, message=ERROR_LABELS[PANAROO_FAILED],
            detail="core_gene_alignment.aln not produced", stage="panaroo",
        )
    return aln


# ─── Stage 4: IQ-TREE ─────────────────────────────────────────────────────────

def run_iqtree(alignment: Path, out_dir: Path) -> Path:
    iqtree_out = out_dir / "iqtree"
    iqtree_out.mkdir(exist_ok=True)
    prefix = str(iqtree_out / "phylo")

    env = os.environ.copy()
    env["PATH"] = str(CONDA_BASE / "envs" / "iqtree" / "bin") + ":" + env["PATH"]

    cmd = [
        "iqtree2",
        "-s", str(alignment),
        "-m", "GTR+G",
        "-B", "1000",
        "--prefix", prefix,
        "-T", str(scheduler.threads_for_stage()),
        "--redo",
    ]
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True,
            timeout=10800, env=env, cwd=str(iqtree_out),
        )
    except subprocess.TimeoutExpired:
        raise PipelineError(
            code=IQTREE_FAILED, message=ERROR_LABELS[IQTREE_FAILED],
            detail="IQ-TREE timed out after 3 hours", stage="iqtree",
        )

    if result.returncode != 0:
        raise PipelineError(
            code=IQTREE_FAILED, message=ERROR_LABELS[IQTREE_FAILED],
            detail=result.stderr[:500], stage="iqtree",
        )

    treefile = Path(f"{prefix}.treefile")
    if not treefile.exists():
        raise PipelineError(
            code=IQTREE_FAILED, message=ERROR_LABELS[IQTREE_FAILED],
            detail="phylo.treefile not produced", stage="iqtree",
        )
    return treefile


# ─── Stage 5: R Annotation ────────────────────────────────────────────────────

def run_r_annotation(treefile: Path, metadata_csv: Path, out_dir: Path) -> Path | None:
    svg_path = out_dir / "annotated_tree.svg"
    env = os.environ.copy()
    env["PATH"] = str(CONDA_BASE / "envs" / "r_phylo" / "bin") + ":" + env["PATH"]

    cmd = ["Rscript", str(_R_SCRIPT), str(treefile), str(metadata_csv), str(svg_path)]
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True,
            timeout=1800, env=env, cwd=str(out_dir),
        )
        if result.returncode != 0 or not svg_path.exists():
            logger.warning("R annotation failed: %s", result.stderr[:300])
            return None
        return svg_path
    except Exception as e:
        logger.warning("R annotation error: %s", e)
        return None


# ─── Stage 6: Report & Helpers ────────────────────────────────────────────────

def _update_stage(phylo_id: str, stage: str, status: str, detail: str = "") -> None:
    run = db.get_phylo_run(phylo_id)
    stages = json.loads(run["stages"]) if run and run.get("stages") else {}
    stages[stage] = {"status": status, "detail": detail}
    db.update_phylo_run(phylo_id, stages=stages)


def _collect_gffs(jobs: list[dict]) -> list[Path]:
    gffs = []
    for job in jobs:
        jid = job["id"]
        gff = RESULTS_DIR / jid / "annotation" / f"{jid}.gff3"
        if gff.exists():
            gffs.append(gff)
        else:
            logger.warning("No GFF found for job %s", jid)
    return gffs


def _parse_panaroo_stats(panaroo_out: Path) -> dict:
    summary = panaroo_out / "summary_statistics.txt"
    stats: dict[str, Any] = {}
    if not summary.exists():
        return stats
    for line in summary.read_text().splitlines():
        if "Core genes" in line:
            stats["core_genes"] = line.split("\t")[-1].strip()
        elif "Shell genes" in line:
            stats["accessory_genes"] = line.split("\t")[-1].strip()
        elif "Cloud genes" in line:
            stats["unique_genes"] = line.split("\t")[-1].strip()
    return stats


def _build_default_metadata(jobs: list[dict], ncbi_refs: list[dict],
                             out_dir: Path) -> Path:
    csv_path = out_dir / "metadata.csv"
    fields = ["sample_id", "country", "year", "mlst_st",
              "amr_profile", "source", "host"]
    rows = []
    for job in jobs:
        stages = json.loads(job.get("stages", "{}"))
        try:
            mlst = json.loads(stages.get("mlst", {}).get("detail", "{}"))
        except Exception:
            mlst = {}
        try:
            amr = json.loads(stages.get("amr", {}).get("detail", "{}"))
        except Exception:
            amr = {}
        rows.append({
            "sample_id":   job["id"][:8],
            "country":     "",
            "year":        "",
            "mlst_st":     mlst.get("st", ""),
            "amr_profile": ";".join(g["gene"] for g in amr.get("genes", [])),
            "source":      "",
            "host":        "",
        })
    for ref in ncbi_refs:
        rows.append({
            "sample_id":   ref["sample_name"],
            "country":     ref["country"],
            "year":        "",
            "mlst_st":     "",
            "amr_profile": "",
            "source":      "reference",
            "host":        "",
        })
    with open(csv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    return csv_path


def _jobs_to_sample_dicts(jobs: list[dict], meta_csv: Path) -> list[dict]:
    meta: dict[str, dict] = {}
    if meta_csv.exists():
        with open(meta_csv) as fh:
            for row in csv.DictReader(fh):
                meta[row.get("sample_id", "")] = row
    result = []
    for job in jobs:
        jid8 = job["id"][:8]
        m = meta.get(jid8, {})
        result.append({
            "sample_name":  job.get("filename", jid8),
            "country":      m.get("country", ""),
            "year":         m.get("year", ""),
            "mlst_st":      m.get("mlst_st", ""),
            "amr_profile":  m.get("amr_profile", ""),
            "source":       m.get("source", ""),
            "host":         m.get("host", ""),
            "is_reference": False,
        })
    return result


# ─── Orchestrator ─────────────────────────────────────────────────────────────

def run_phylo(phylo_id: str, job_ids: list[str],
              metadata_csv_path: Path | None = None) -> None:
    out_dir = PHYLO_DIR / phylo_id
    out_dir.mkdir(parents=True, exist_ok=True)

    db.update_phylo_run(phylo_id, status="running")
    warnings: list[str] = []
    ncbi_refs: list[dict] = []

    try:
        # ── Stage 1: Organism ───────────────────────────────────────────────
        _update_stage(phylo_id, "organism", "running")
        jobs = [db.get_job(jid) for jid in job_ids if db.get_job(jid)]
        if not jobs:
            raise PipelineError(
                code=NO_ASSEMBLIES, message=ERROR_LABELS[NO_ASSEMBLIES],
                detail="None of the selected jobs exist", stage="organism",
            )
        organism = detect_organism(jobs)
        if not organism:
            warnings.append("Could not determine organism — skipping NCBI references.")
        db.update_phylo_run(phylo_id, organism=organism or "Unknown")
        _update_stage(phylo_id, "organism", "done", organism or "")

        # ── Stage 2: NCBI references ────────────────────────────────────────
        _update_stage(phylo_id, "ncbi", "running")
        if organism:
            try:
                ncbi_refs = fetch_ncbi_references(organism, out_dir)
                if not ncbi_refs:
                    warnings.append(f"No NCBI references found for {organism}.")
            except Exception as e:
                logger.warning("NCBI fetch error: %s", e)
                warnings.append("NCBI reference fetch failed — continuing with user samples only.")
        _update_stage(phylo_id, "ncbi", "done", f"{len(ncbi_refs)} references")

        # ── Collect GFFs ────────────────────────────────────────────────────
        user_gffs = _collect_gffs(jobs)
        ncbi_gffs = [Path(r["gff_path"]) for r in ncbi_refs if r.get("gff_path")]
        all_gffs  = user_gffs + ncbi_gffs

        if len(all_gffs) < 2:
            raise PipelineError(
                code=NO_ASSEMBLIES, message=ERROR_LABELS[NO_ASSEMBLIES],
                detail=f"Only {len(all_gffs)} GFF file(s) available, need at least 2",
                stage="panaroo",
            )

        # ── Stage 3: Panaroo ────────────────────────────────────────────────
        _update_stage(phylo_id, "panaroo", "running")
        aln = run_panaroo(all_gffs, out_dir)
        panaroo_stats = _parse_panaroo_stats(out_dir / "panaroo")
        _update_stage(phylo_id, "panaroo", "done")

        # ── Stage 4: IQ-TREE ────────────────────────────────────────────────
        _update_stage(phylo_id, "iqtree", "running")
        treefile = run_iqtree(aln, out_dir)
        _update_stage(phylo_id, "iqtree", "done")

        # ── Stage 5: R annotation ───────────────────────────────────────────
        _update_stage(phylo_id, "r", "running")
        meta_csv = metadata_csv_path or _build_default_metadata(jobs, ncbi_refs, out_dir)
        svg_path = run_r_annotation(treefile, meta_csv, out_dir)
        if svg_path is None:
            warnings.append("R annotation failed — report includes raw Newick tree.")
        _update_stage(phylo_id, "r", "done" if svg_path else "failed")

        # ── Stage 6: Report ──────────────────────────────────────────────────
        _update_stage(phylo_id, "report", "running")
        from phylo_report import generate_phylo_report
        all_samples = _jobs_to_sample_dicts(jobs, meta_csv) + ncbi_refs
        html = generate_phylo_report(
            phylo_id=phylo_id,
            organism=organism or "Unknown",
            samples=all_samples,
            panaroo_stats=panaroo_stats,
            treefile=treefile,
            svg_path=svg_path,
            ncbi_refs=ncbi_refs,
            warnings=warnings,
        )
        report_path = out_dir / "phylo_report.html"
        report_path.write_text(html, encoding="utf-8")
        _update_stage(phylo_id, "report", "done")
        db.update_phylo_run(phylo_id, status="completed",
                            report_path=str(report_path))

    except PipelineError as exc:
        logger.error("[phylo:%s] %s", phylo_id, exc)
        db.update_phylo_run(phylo_id, status="failed",
                            error=str(exc), error_code=exc.code,
                            error_detail=exc.detail)
    except Exception as exc:
        logger.error("[phylo:%s] Unexpected error: %s", phylo_id, exc, exc_info=True)
        db.update_phylo_run(phylo_id, status="failed", error=str(exc))
