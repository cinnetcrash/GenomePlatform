"""
PCR Primer Tasarımı:
1. AMR/virülans genlerinin dizilerini al
2. MAFFT ile hizala → korunmuş bölgeleri bul
3. Primer3 ile primer tasarla
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
    Primer3'ü çalıştırır, primer çiftlerini döndürür.
    """
    # Primer3 boulder-IO formatı
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
        logger.error("Primer3 çalıştırılamadı: %s", e)
        return []

    # Çıktıyı ayrıştır
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
                "pair":      i + 1,
                "forward":   left_seq,
                "reverse":   right_seq,
                "fwd_tm":    round(float(left_tm), 1) if left_tm else None,
                "rev_tm":    round(float(right_tm), 1) if right_tm else None,
                "fwd_gc":    round(float(left_gc), 1) if left_gc else None,
                "rev_gc":    round(float(right_gc), 1) if right_gc else None,
                "amplicon":  int(product) if product else None,
                "penalty":   round(float(penalty), 3) if penalty else None,
            })
    return pairs


# ─── Korunmuş Bölge Tespiti ───────────────────────────────────────────────────

def _find_conserved_region(aligned_seqs: list[str],
                           min_length: int = 200) -> str | None:
    """
    MAFFT hizalama çıktısından en uzun korunmuş bloğu bulur.
    Korunmuş = tüm dizilerde aynı baz, gap yok.
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
    Verilen gen dizileri için:
    1. MAFFT ile hizala (birden fazla dizi varsa)
    2. Korunmuş bölgeyi bul
    3. Primer3 ile primerları tasarla
    """
    result: dict[str, Any] = {
        "gene": gene_name,
        "n_sequences": len(sequences),
        "conserved_region_bp": None,
        "primers": [],
        "note": "",
    }

    if not sequences:
        result["note"] = "Dizi bulunamadı."
        return result

    # Tek dizi varsa doğrudan primer tasarla
    if len(sequences) == 1:
        template = sequences[0].upper()
        result["conserved_region_bp"] = len(template)
        result["primers"] = _run_primer3(template, gene_name)
        result["note"] = "Tek dizi; MAFFT atlandı."
        return result

    # MAFFT ile hizala
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
        logger.error("MAFFT hatası: %s", e)
        result["note"] = f"MAFFT hatası: {e}"
        fasta_path.unlink(missing_ok=True)
        return result
    finally:
        fasta_path.unlink(missing_ok=True)

    # Hizalanmış dizileri ayrıştır
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

    # Korunmuş bölge
    conserved = _find_conserved_region(aligned_seqs)
    if not conserved:
        result["note"] = "Yeterince uzun korunmuş bölge bulunamadı (min 200 bp)."
        return result

    result["conserved_region_bp"] = len(conserved)
    result["primers"] = _run_primer3(conserved, gene_name)
    return result


# ─── Ana Fonksiyon ───────────────────────────────────────────────────────────

def design_all_primers(amr_results: dict[str, Any],
                        ai_targets: list[dict],
                        assembly_fasta: Path | None) -> list[dict[str, Any]]:
    """
    AI'nın önerdiği PCR hedeflerine ve AMR genlerine göre primer tasarlar.
    """
    primer_results = []

    # AI'nın önerdiği genleri önceliklendir
    target_genes = [t.get("gene", "") for t in ai_targets if t.get("gene")]

    # AMR bulunanları da ekle (duplikat önle)
    for g in amr_results.get("genes", []):
        if g["gene"] not in target_genes:
            target_genes.append(g["gene"])

    if not target_genes:
        logger.info("Primer tasarımı için hedef gen bulunamadı.")
        return []

    # Assembly FASTA'sından gen dizilerini çıkar (basit kNN yaklaşımı)
    # Prototip için: assembly'den ilk N kontig'i al, her gen için tek dizi
    gene_sequences: dict[str, list[str]] = {}

    if assembly_fasta and assembly_fasta.exists():
        contigs = _parse_fasta(assembly_fasta)
        # Her hedef gen için tüm contigleri dizi olarak ekle (gerçek uygulamada
        # BLAST ile gen bölgesi çıkarılır; prototipte ilk contig'i kullanıyoruz)
        for gene in target_genes:
            if contigs:
                gene_sequences[gene] = [contigs[0][1][:3000]]  # ilk 3000 bp

    for gene in target_genes:
        seqs = gene_sequences.get(gene, [])
        primers = design_primers_for_gene(gene, seqs)
        primer_results.append(primers)
        logger.info("Primer tasarımı: %s → %d çift",
                    gene, len(primers.get("primers", [])))

    return primer_results


def _parse_fasta(path: Path) -> list[tuple[str, str]]:
    """Basit FASTA ayrıştırıcı."""
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
