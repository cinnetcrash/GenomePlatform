"""
System dependency checker.
Scans conda environments, system PATH, and databases; reports what is present or missing.
"""
import subprocess
import shutil
import logging
from pathlib import Path

from config import (
    AUTOCYCLER_BIN, BANDAGE_BIN, CHECKM2_DB, CHECKV_DB,
    CONDA_BASE, KRAKEN2_DEFAULT_DB,
)

logger = logging.getLogger("system_check")

# ── Tool definitions ───────────────────────────────────────────────────────────
# (display_name, conda_env_or_None, binary_name, required)
# conda_env=None  → look in system PATH or as an absolute path
TOOLS = [
    # Read QC
    ("fastp",         None,       "fastp",             True),   # /usr/bin/fastp
    ("NanoPlot",      "analiz",   "NanoPlot",          False),  # optional (MinION QC)
    # Assembly
    ("Flye",          "analiz",   "flye",              True),
    ("Shovill",       "shovill",  "shovill",           True),
    ("SPAdes",        "shovill",  "spades.py",         False),
    # Taxonomic
    ("Kraken2",       None,       "kraken2",           False),  # /usr/bin/kraken2
    # Typing & AMR
    ("mlst",          "mlst",     "mlst",              True),
    ("AMRFinderPlus", "analiz",   "amrfinder",         True),
    ("Abricate",      "analiz",   "abricate",          True),
    # Plasmid
    ("MOB-Recon",     "mobsuite", "mob_recon",         True),
    # Quality assessment
    ("QUAST",         "quast5",   "quast.py",          True),
    ("CheckM2",       "checkM",   "checkm2",           True),
    ("CheckV",        "checkv",   "checkv",            False),
    # Annotation
    ("Bakta",         "bakta",    "bakta",             False),
    # Visualisation
    ("Bandage",       None,       str(BANDAGE_BIN),    False),
    ("Autocycler",    None,       str(AUTOCYCLER_BIN), False),
]

# ── Database definitions ───────────────────────────────────────────────────────
DATABASES = [
    ("Kraken2 DB (default)",     KRAKEN2_DEFAULT_DB,    "hash.k2d",            False),
    ("CheckM2 DB",               CHECKM2_DB.parent,     CHECKM2_DB.name,       True),
    ("CheckV DB",                CHECKV_DB,             "genome_db.faa.dmnd",  False),
    ("Bakta DB",                 CONDA_BASE / "envs" / "bakta" / "db",
                                                        "db.db",               False),
    ("VFDB sequences",           Path("/home/analysis/miniconda3/envs/analiz/db/vfdb"),
                                                        "sequences",           False),
    ("Pharokka DB",              Path("/home/analysis/databases_all/database_pharokka"),
                                                        "PHROG_phrogs_profile_db",  False),
]

# ── Standalone binaries ────────────────────────────────────────────────────────
BINARIES = [
    ("Bandage",     BANDAGE_BIN,     True),
    ("Autocycler",  AUTOCYCLER_BIN,  False),
]


def _check_tool(name: str, env: str | None, binary: str) -> dict:
    """Returns {name, status, version, note} for a single tool."""
    if env is None:
        # Try as absolute/relative path first, then fall back to PATH search
        p = Path(binary)
        if not (p.exists() and p.is_file()):
            found = shutil.which(binary)
            p = Path(found) if found else None
        if p and p.exists() and p.is_file():
            return {"name": name, "status": "ok", "env": "system",
                    "note": str(p)}
        return {"name": name, "status": "missing", "env": "system",
                "note": f"Not found: {binary}"}

    conda_bin = CONDA_BASE / "envs" / env / "bin" / binary
    if not conda_bin.exists():
        # Try PATH inside env
        conda_bin = CONDA_BASE / "envs" / env / "bin" / Path(binary).name
    if not conda_bin.exists():
        return {"name": name, "status": "missing", "env": env,
                "note": f"Not in {CONDA_BASE}/envs/{env}/bin/"}

    # Try to get version
    try:
        import os
        e = os.environ.copy()
        e["PATH"] = str(conda_bin.parent) + ":" + e["PATH"]
        r = subprocess.run(
            [str(conda_bin), "--version"],
            capture_output=True, text=True, timeout=10, env=e,
        )
        ver = (r.stdout + r.stderr).strip().split("\n")[0][:80]
    except Exception:
        ver = "installed"

    return {"name": name, "status": "ok", "env": env, "note": ver}


def _check_database(name: str, base: Path, marker: str) -> dict:
    """Returns {name, status, path, note} for a database."""
    target = base / marker
    if target.exists():
        try:
            size_gb = sum(f.stat().st_size for f in base.rglob("*") if f.is_file()) / 1e9
            note = f"{size_gb:.1f} GB"
        except Exception:
            note = "present"
        return {"name": name, "status": "ok",   "path": str(base), "note": note}
    return {"name": name, "status": "missing", "path": str(base),
            "note": f"Expected: {marker}"}


def run_system_check() -> dict:
    """
    Runs all checks and returns a structured report.
    """
    tool_results = []
    for display, env, binary, required in TOOLS:
        r = _check_tool(display, env, binary)
        r["required"] = required
        tool_results.append(r)

    db_results = []
    for display, base, marker, required in DATABASES:
        r = _check_database(display, base, marker)
        r["required"] = required
        db_results.append(r)

    ok_tools    = sum(1 for t in tool_results if t["status"] == "ok")
    ok_dbs      = sum(1 for d in db_results  if d["status"] == "ok")
    miss_req    = [t["name"] for t in tool_results + db_results
                   if t["status"] == "missing" and t.get("required")]

    return {
        "tools":            tool_results,
        "databases":        db_results,
        "summary": {
            "tools_ok":     ok_tools,
            "tools_total":  len(tool_results),
            "dbs_ok":       ok_dbs,
            "dbs_total":    len(db_results),
            "missing_required": miss_req,
            "ready":        len(miss_req) == 0,
        },
    }
