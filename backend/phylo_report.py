"""HTML report generator for Global Phylogeny analysis."""
import base64
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def _svg_to_data_uri(svg_path: Path) -> str:
    data = svg_path.read_bytes()
    b64  = base64.b64encode(data).decode()
    return f"data:image/svg+xml;base64,{b64}"


def generate_phylo_report(
    phylo_id: str,
    organism: str,
    samples: list[dict],
    panaroo_stats: dict[str, Any],
    treefile: Path | None,
    svg_path: Path | None,
    ncbi_refs: list[dict],
    warnings: list[str],
) -> str:
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    if svg_path and svg_path.exists():
        tree_block = f'<img src="{_svg_to_data_uri(svg_path)}" style="max-width:100%;height:auto">'
    elif treefile and treefile.exists():
        newick = treefile.read_text().strip()
        tree_block = (
            "<p style='color:#b45309'>⚠ R annotation failed — raw Newick below:</p>"
            f"<pre style='font-size:.75rem;overflow-x:auto;background:#f8fafc;"
            f"padding:1rem;border-radius:8px'>{newick}</pre>"
        )
    else:
        tree_block = "<p style='color:#b45309'>⚠ Phylogenetic tree not available.</p>"

    warn_html = "".join(
        f"<div style='background:#fef9c3;border-left:4px solid #ca8a04;"
        f"padding:.6rem 1rem;margin:.5rem 0;border-radius:4px'>{w}</div>"
        for w in warnings
    )

    pan_html = (
        f"<table style='border-collapse:collapse;width:100%'>"
        f"<tr><td style='padding:.4rem .8rem'>Core genes</td>"
        f"<td style='padding:.4rem .8rem;font-weight:600'>{panaroo_stats.get('core_genes','—')}</td></tr>"
        f"<tr style='background:#f8fafc'><td style='padding:.4rem .8rem'>Accessory genes</td>"
        f"<td style='padding:.4rem .8rem;font-weight:600'>{panaroo_stats.get('accessory_genes','—')}</td></tr>"
        f"<tr><td style='padding:.4rem .8rem'>Unique genes</td>"
        f"<td style='padding:.4rem .8rem;font-weight:600'>{panaroo_stats.get('unique_genes','—')}</td></tr>"
        f"</table>"
    )

    rows = "".join(
        f"<tr style='{'background:#f8fafc' if i%2 else ''}'>"
        f"<td style='padding:.35rem .7rem'>{s.get('sample_name','')}</td>"
        f"<td style='padding:.35rem .7rem'>{s.get('country','')}</td>"
        f"<td style='padding:.35rem .7rem'>{s.get('mlst_st','')}</td>"
        f"<td style='padding:.35rem .7rem'>{s.get('amr_profile','')}</td>"
        f"<td style='padding:.35rem .7rem'>{s.get('source','')}</td>"
        f"<td style='padding:.35rem .7rem'>{s.get('host','')}</td>"
        f"<td style='padding:.35rem .7rem'>{s.get('year','')}</td>"
        f"</tr>"
        for i, s in enumerate(samples)
    )

    ref_rows = "".join(
        f"<tr style='{'background:#f8fafc' if i%2 else ''}'>"
        f"<td style='padding:.35rem .7rem'>{r.get('country','')}</td>"
        f"<td style='padding:.35rem .7rem;font-family:monospace'>{r.get('accession','')}</td>"
        f"</tr>"
        for i, r in enumerate(ncbi_refs)
    )

    return f"""<!DOCTYPE html>
<html lang="en">
<head><meta charset="UTF-8">
<title>Global Phylogeny — {organism}</title>
<style>
body{{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;
  margin:0;padding:2rem;background:#f8fafc;color:#1e293b}}
h1{{font-size:1.5rem;font-weight:700;margin-bottom:.25rem}}
h2{{font-size:1.1rem;font-weight:600;margin:2rem 0 .75rem;color:#334155}}
.card{{background:#fff;border-radius:12px;padding:1.5rem;margin-bottom:1.5rem;
  box-shadow:0 1px 3px rgba(0,0,0,.08)}}
table{{border-collapse:collapse;width:100%;font-size:.875rem}}
th{{background:#f1f5f9;padding:.5rem .8rem;text-align:left;font-weight:600;color:#475569}}
td{{border-top:1px solid #e2e8f0}}
</style>
</head>
<body>
<div style="max-width:1200px;margin:0 auto">
<h1>🌍 Global Phylogeny Report</h1>
<p style="color:#64748b">Organism: <strong>{organism}</strong> &nbsp;|&nbsp; Generated: {now}</p>
<p><a href="/" style="color:#6366f1">← New Analysis</a> &nbsp;|&nbsp;
   <a href="/phylo" style="color:#6366f1">Global Phylogeny</a></p>

{warn_html}

<div class="card">
<h2>Phylogenetic Tree</h2>
{tree_block}
</div>

<div class="card">
<h2>Pan-Genome Statistics (Panaroo)</h2>
{pan_html}
</div>

<div class="card">
<h2>Samples ({len(samples)})</h2>
<table>
<tr><th>Sample</th><th>Country</th><th>MLST ST</th><th>AMR Profile</th>
    <th>Source</th><th>Host</th><th>Year</th></tr>
{rows}
</table>
</div>

<div class="card">
<h2>NCBI Reference Genomes ({len(ncbi_refs)} countries)</h2>
<table><tr><th>Country</th><th>Accession</th></tr>
{ref_rows}
</table>
</div>

</div>
</body></html>"""
