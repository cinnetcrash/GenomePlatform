"""
PCR Primer Design:
1. Extract sequences for AMR/virulence genes
2. Align with MAFFT → find conserved regions
3. Design primers with Primer3
"""
import json
import logging
import subprocess
import tempfile
from pathlib import Path
from typing import Any

from config import CONDA_BASE
from security import sanitize_path

logger = logging.getLogger("primer_designer")

PRIMER3_CORE = "primer3_core"
MAFFT_BIN    = str(CONDA_BASE / "bin" / "mafft")


# ─── Primer3 ─────────────────────────────────────────────────────────────────

def _run_primer3(sequence: str, target_name: str,
                 product_min: int = 100,
                 product_max: int = 800) -> list[dict]:
    """
    Runs Primer3 and returns the designed primer pairs.
    """
    # Primer3 boulder-IO format
    input_block = (
        f"SEQUENCE_ID={target_name}\n"
        f"SEQUENCE_TEMPLATE={sequence}\n"
        f"PRIMER_TASK=generic\n"
        f"PRIMER_PICK_LEFT_PRIMER=1\n"
        f"PRIMER_PICK_RIGHT_PRIMER=1\n"
        f"PRIMER_OPT_SIZE=20\n"
        f"PRIMER_MIN_SIZE=18\n"
        f"PRIMER_MAX_SIZE=25\n"
        f"PRIMER_OPT_TM=60.0\n"
        f"PRIMER_MIN_TM=57.0\n"
        f"PRIMER_MAX_TM=63.0\n"
        f"PRIMER_MIN_GC=40.0\n"
        f"PRIMER_MAX_GC=65.0\n"
        f"PRIMER_MAX_POLY_X=4\n"
        f"PRIMER_SALT_MONOVALENT=50.0\n"
        f"PRIMER_DNA_CONC=250.0\n"
        f"PRIMER_MAX_NS_ACCEPTED=0\n"
        f"PRIMER_PRODUCT_SIZE_RANGE={product_min}-{product_max}\n"
        f"PRIMER_NUM_RETURN=3\n"
        f"PRIMER_EXPLAIN_FLAG=1\n"
        f"=\n"
    )

    try:
        result = subprocess.run(
            [PRIMER3_CORE],
            input=input_block,
            capture_output=True,
            text=True,
            timeout=60,
        )
    except (FileNotFoundError, subprocess.TimeoutExpired) as e:
        logger.error("Failed to run Primer3: %s", e)
        return []

    # Parse output
    pairs = []
    out = result.stdout
    lines = {l.split("=")[0]: l.split("=", 1)[1]
             for l in out.splitlines() if "=" in l}

    num_returned = int(lines.get("PRIMER_PAIR_NUM_RETURNED", "0"))
    for i in range(num_returned):
        left_seq  = lines.get(f"PRIMER_LEFT_{i}_SEQUENCE", "")
        right_seq = lines.get(f"PRIMER_RIGHT_{i}_SEQUENCE", "")
        left_tm   = lines.get(f"PRIMER_LEFT_{i}_TM", "")
        right_tm  = lines.get(f"PRIMER_RIGHT_{i}_TM", "")
        left_gc   = lines.get(f"PRIMER_LEFT_{i}_GC_PERCENT", "")
        right_gc  = lines.get(f"PRIMER_RIGHT_{i}_GC_PERCENT", "")
        product   = lines.get(f"PRIMER_PAIR_{i}_PRODUCT_SIZE", "")
        penalty   = lines.get(f"PRIMER_PAIR_{i}_PENALTY", "")

        if left_seq and right_seq:
            pairs.append({
                "pair":     i + 1,
                "forward":  left_seq,
                "reverse":  right_seq,
                "fwd_tm":   round(float(left_tm), 1) if left_tm else None,
                "rev_tm":   round(float(right_tm), 1) if right_tm else None,
                "fwd_gc":   round(float(left_gc), 1) if left_gc else None,
                "rev_gc":   round(float(right_gc), 1) if right_gc else None,
                "amplicon": int(product) if product else None,
                "penalty":  round(float(penalty), 3) if penalty else None,
            })
    return pairs


# ─── Conserved Region Detection ──────────────────────────────────────────────

def _find_conserved_region(aligned_seqs: list[str],
                           min_length: int = 200) -> str | None:
    """
    Finds the longest conserved block from a MAFFT alignment output.
    Conserved = identical base across all sequences, no gaps.
    """
    if not aligned_seqs:
        return None

    n_seqs = len(aligned_seqs)
    length = len(aligned_seqs[0])

    best_start = best_end = 0
    cur_start = None

    for pos in range(length):
        col = [s[pos] for s in aligned_seqs]
        conserved = (
            col[0] not in "-N" and
            all(b == col[0] for b in col)
        )
        if conserved:
            if cur_start is None:
                cur_start = pos
        else:
            if cur_start is not None:
                if (pos - cur_start) > (best_end - best_start):
                    best_start, best_end = cur_start, pos
                cur_start = None

    if cur_start is not None and (length - cur_start) > (best_end - best_start):
        best_start, best_end = cur_start, length

    block = aligned_seqs[0][best_start:best_end].replace("-", "")
    if len(block) >= min_length:
        return block
    return None


def design_primers_for_gene(gene_name: str,
                             sequences: list[str]) -> dict[str, Any]:
    """
    For the given gene sequences:
    1. Align with MAFFT (if more than one sequence)
    2. Find the conserved region
    3. Design primers with Primer3
    """
    result: dict[str, Any] = {
        "gene": gene_name,
        "n_sequences": len(sequences),
        "conserved_region_bp": None,
        "primers": [],
        "note": "",
    }

    if not sequences:
        result["note"] = "No sequences found."
        return result

    # Single sequence — design primers directly
    if len(sequences) == 1:
        template = sequences[0].upper()
        result["conserved_region_bp"] = len(template)
        result["primers"] = _run_primer3(template, gene_name)
        result["note"] = "Single sequence; MAFFT skipped."
        return result

    # Align with MAFFT
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa",
                                     delete=False) as fh:
        for i, seq in enumerate(sequences):
            fh.write(f">seq{i}\n{seq}\n")
        fasta_path = Path(fh.name)

    try:
        mafft_result = subprocess.run(
            [MAFFT_BIN, "--auto", "--quiet", str(fasta_path)],
            capture_output=True, text=True, timeout=300,
        )
        aligned_raw = mafft_result.stdout.strip()
    except Exception as e:
        logger.error("MAFFT error: %s", e)
        result["note"] = f"MAFFT error: {e}"
        fasta_path.unlink(missing_ok=True)
        return result
    finally:
        fasta_path.unlink(missing_ok=True)

    # Parse aligned sequences
    aligned_seqs = []
    cur_seq = []
    for line in aligned_raw.splitlines():
        if line.startswith(">"):
            if cur_seq:
                aligned_seqs.append("".join(cur_seq).upper())
            cur_seq = []
        else:
            cur_seq.append(line.strip())
    if cur_seq:
        aligned_seqs.append("".join(cur_seq).upper())

    # Find conserved region
    conserved = _find_conserved_region(aligned_seqs)
    if not conserved:
        result["note"] = "No conserved region long enough found (min 200 bp)."
        return result

    result["conserved_region_bp"] = len(conserved)
    result["primers"] = _run_primer3(conserved, gene_name)
    return result


# ─── Main Entry Point ────────────────────────────────────────────────────────

def design_all_primers(amr_results: dict[str, Any],
                        ai_targets: list[dict],
                        assembly_fasta: Path | None) -> list[dict[str, Any]]:
    """
    Designs primers for AI-recommended PCR targets and detected AMR genes.
    """
    primer_results = []

    # Prioritise AI-recommended genes
    target_genes = [t.get("gene", "") for t in ai_targets if t.get("gene")]

    # Add AMR findings (deduplicate)
    for g in amr_results.get("genes", []):
        if g["gene"] not in target_genes:
            target_genes.append(g["gene"])

    if not target_genes:
        logger.info("No target genes found for primer design.")
        return []

    # Extract gene sequences from assembly FASTA
    # Prototype: uses first N contigs as sequence source per gene.
    # Production: would use BLAST to extract the exact gene region.
    gene_sequences: dict[str, list[str]] = {}

    if assembly_fasta and assembly_fasta.exists():
        contigs = _parse_fasta(assembly_fasta)
        for gene in target_genes:
            if contigs:
                gene_sequences[gene] = [contigs[0][1][:3000]]  # first 3000 bp

    for gene in target_genes:
        seqs = gene_sequences.get(gene, [])
        primers = design_primers_for_gene(gene, seqs)
        primer_results.append(primers)
        logger.info("Primer design: %s → %d pair(s)",
                    gene, len(primers.get("primers", [])))

    return primer_results


def _parse_fasta(path: Path) -> list[tuple[str, str]]:
    """Simple FASTA parser."""
    contigs = []
    cur_head, cur_seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if cur_head:
                    contigs.append((cur_head, "".join(cur_seq)))
                cur_head = line[1:]
                cur_seq = []
            else:
                cur_seq.append(line)
    if cur_head:
        contigs.append((cur_head, "".join(cur_seq)))
    return contigs
