"""
HTML report generator.
Combines pipeline + AI results into a single HTML report.
"""
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


import base64

RISK_COLORS = {
    "LOW":      ("#22c55e", "#f0fdf4"),
    "MEDIUM":   ("#f59e0b", "#fffbeb"),
    "HIGH":     ("#ef4444", "#fef2f2"),
    "CRITICAL": ("#7f1d1d", "#fef2f2"),
    "UNKNOWN":  ("#6b7280", "#f9fafb"),
}


def _kraken2_section(k2: dict) -> str:
    """Renders Kraken2 results as a horizontal CSS bar chart."""
    if not k2:
        return "<p style='color:#6b7280'>Kraken2 was not run for this analysis.</p>"

    classified = k2.get("classified_percent", 0)
    unclassified = k2.get("unclassified_percent", 100)
    taxa = k2.get("top_taxa", [])

    summary = (
        f"<p style='margin-bottom:1rem'>"
        f"<strong>{classified}%</strong> of reads classified &nbsp;|&nbsp; "
        f"{unclassified}% unclassified</p>"
    )

    if not taxa:
        return summary + "<p style='color:#6b7280'>No species-level hits above 0.01%.</p>"

    bars = []
    max_pct = taxa[0]["percent"] if taxa else 1
    for t in taxa:
        bar_w = round(t["percent"] / max_pct * 100, 1)
        color  = "#6366f1" if t["percent"] >= 5 else "#a5b4fc"
        bars.append(f"""
      <div style="display:flex;align-items:center;gap:.6rem;margin-bottom:.55rem">
        <div style="width:200px;font-size:.8rem;color:#334155;white-space:nowrap;
                    overflow:hidden;text-overflow:ellipsis" title="{t['name']}">
          <em>{t['name']}</em>
        </div>
        <div style="flex:1;background:#e2e8f0;border-radius:9999px;height:14px;overflow:hidden">
          <div style="width:{bar_w}%;background:{color};height:100%;border-radius:9999px"></div>
        </div>
        <div style="width:48px;text-align:right;font-size:.8rem;font-weight:600;
                    color:#1e40af;font-family:monospace">{t['percent']}%</div>
        <div style="width:70px;text-align:right;font-size:.75rem;color:#94a3b8;
                    font-family:monospace">{t['reads']:,} reads</div>
      </div>""")

    return summary + "\n".join(bars)


def _host_depletion_section(depl: dict) -> str:
    """Renders host depletion metrics — shown only when host reads were detected."""
    if not depl or not depl.get("depleted"):
        return ""

    total   = depl.get("total_reads", 0)
    removed = depl.get("host_reads_removed", 0)
    host_pct = depl.get("host_percent", 0.0)
    after   = depl.get("reads_after", 0)
    after_pct = depl.get("reads_after_percent", 0.0)

    def _fmt(n):
        try:
            return f"{int(n):,}"
        except (ValueError, TypeError):
            return str(n)

    # Colour the host bar red, the remaining bar green
    return f"""
    <div class="section">
      <div class="section-title">&#x1F9F9; Host Read Depletion</div>
      <div class="section-body">
        <p style="margin-bottom:1.2rem;font-size:.9rem;color:#374151">
          <strong>{host_pct}%</strong> of reads were classified as
          <em>Homo sapiens</em> (taxid&nbsp;9606) by Kraken2 and removed before assembly.
        </p>

        <!-- Before -->
        <div style="margin-bottom:.8rem">
          <div style="display:flex;justify-content:space-between;font-size:.8rem;
                      color:#64748b;margin-bottom:.25rem">
            <span>Before depletion</span>
            <span><strong>{_fmt(total)}</strong> reads · 100%</span>
          </div>
          <div style="background:#e2e8f0;border-radius:9999px;height:14px;overflow:hidden">
            <div style="width:100%;background:#94a3b8;height:100%;border-radius:9999px"></div>
          </div>
        </div>

        <!-- After -->
        <div style="margin-bottom:1.2rem">
          <div style="display:flex;justify-content:space-between;font-size:.8rem;
                      color:#64748b;margin-bottom:.25rem">
            <span>After depletion</span>
            <span><strong>{_fmt(after)}</strong> reads · {after_pct}%</span>
          </div>
          <div style="background:#e2e8f0;border-radius:9999px;height:14px;overflow:hidden">
            <div style="width:{after_pct}%;background:#22c55e;height:100%;border-radius:9999px"></div>
          </div>
        </div>

        <!-- Summary chips -->
        <div style="display:flex;gap:.75rem;flex-wrap:wrap">
          <span style="background:#fef2f2;border:1px solid #fecaca;color:#991b1b;
                       border-radius:9999px;padding:.25rem .85rem;font-size:.8rem;font-weight:600">
            &#x2796; {_fmt(removed)} host reads removed ({host_pct}%)
          </span>
          <span style="background:#f0fdf4;border:1px solid #bbf7d0;color:#166534;
                       border-radius:9999px;padding:.25rem .85rem;font-size:.8rem;font-weight:600">
            &#x2714; {_fmt(after)} reads used for assembly ({after_pct}%)
          </span>
        </div>
      </div>
    </div>"""


def _assembly_qc_section(qc_checks: list[dict]) -> str:
    """Renders assembly QC pass/warn/fail indicators."""
    if not qc_checks:
        return "<p style='color:#6b7280'>Assembly QC data not available.</p>"

    STATUS_STYLE = {
        "pass": ("&#x2705;", "#166534", "#f0fdf4", "#bbf7d0"),
        "warn": ("&#x26A0;&#xFE0F;", "#92400e", "#fffbeb", "#fde68a"),
        "fail": ("&#x274C;",  "#991b1b", "#fef2f2", "#fecaca"),
    }

    cards = []
    for c in qc_checks:
        icon, text_c, bg, border = STATUS_STYLE.get(c["status"], STATUS_STYLE["warn"])
        val = c["value"]
        if isinstance(val, (int, float)):
            if val >= 1_000_000:
                display = f"{val/1_000_000:.2f} Mb"
            elif val >= 1_000:
                display = f"{val/1_000:.1f} kb"
            else:
                display = str(val)
        else:
            display = str(val)
        if c["unit"] and c["unit"] not in display:
            display += f" {c['unit']}"

        cards.append(f"""
      <div style="background:{bg};border:1px solid {border};border-radius:10px;
                  padding:.9rem 1rem;display:flex;align-items:center;gap:.7rem">
        <span style="font-size:1.1rem">{icon}</span>
        <div>
          <div style="font-size:.72rem;color:#64748b;text-transform:uppercase;
                      letter-spacing:.04em">{c['metric']}</div>
          <div style="font-size:1.1rem;font-weight:700;color:{text_c}">{display}</div>
          <div style="font-size:.7rem;color:#94a3b8;margin-top:.1rem">{c['threshold']}</div>
        </div>
      </div>""")

    return f'<div style="display:grid;grid-template-columns:repeat(auto-fill,minmax(185px,1fr));gap:.75rem">' \
           + "\n".join(cards) + "</div>"


def _error_banner(error: str) -> str:
    """Shows a red error banner with the failure reason at the top of the report."""
    if not error:
        return ""
    # Truncate very long tracebacks for display
    display = error if len(error) <= 600 else error[:600] + "…"
    return (
        "<div class='section'><div class='section-body' style='"
        "background:#fef2f2;border-left:4px solid #dc2626;padding:1rem 1.2rem;"
        "border-radius:0 8px 8px 0'>"
        "<strong style='color:#dc2626'>&#x274C; Analysis Failed — Partial Results</strong><br>"
        "<span style='font-size:.85rem;color:#7f1d1d'>The pipeline encountered an error. "
        "Partial results collected before the failure are shown below.</span>"
        f"<pre style='margin-top:.8rem;background:#fff5f5;padding:.8rem;border-radius:6px;"
        f"font-size:.78rem;color:#991b1b;overflow-x:auto;white-space:pre-wrap'>{display}</pre>"
        "</div></div>"
    )


def _context_banner(gctx: dict) -> str:
    """Shows a coloured banner when bacterial steps were skipped."""
    if not gctx or not gctx.get("skip_bacterial"):
        return ""
    ctx = gctx.get("context", "")
    color = "#f59e0b" if ctx == "metagenomics" else "#ef4444"
    bg    = "#fffbeb" if ctx == "metagenomics" else "#fef2f2"
    icon  = "⚠️" if ctx == "metagenomics" else "🛑"
    return (
        f"<div class='section'><div class='section-body' style='"
        f"background:{bg};border-left:4px solid {color};padding:1rem 1.2rem;"
        f"border-radius:0 8px 8px 0'>"
        f"<strong style='color:{color}'>{icon} Genomic Context: {ctx.replace('_',' ').title()}</strong><br>"
        f"<span style='font-size:.88rem'>{gctx.get('reason','')}</span>"
        f"</div></div>"
    )


def _quast_section(quast: dict) -> str:
    """Renders QUAST assembly quality metrics."""
    if not quast:
        return "<p style='color:#6b7280'>QUAST was not run.</p>"
    label_map = {
        "contigs":        ("# Contigs",       ""),
        "total_length":   ("Total Length",     "bp"),
        "largest_contig": ("Largest Contig",   "bp"),
        "n50":            ("N50",              "bp"),
        "n90":            ("N90",              "bp"),
        "l50":            ("L50",              ""),
        "gc_pct":         ("GC Content",       "%"),
        "ns_per_100k":    ("Ns per 100 kbp",   ""),
    }
    cards = []
    for key, (label, unit) in label_map.items():
        val = quast.get(key, "—")
        if val == "—" or val is None:
            continue
        # Format numbers
        try:
            n = float(str(val).replace(",", ""))
            if unit == "bp":
                if n >= 1e6:  val = f"{n/1e6:.2f} Mb"
                elif n >= 1e3: val = f"{n/1e3:.1f} kb"
                else:          val = f"{int(n)} bp"
            elif unit == "%":
                val = f"{n:.1f}%"
            else:
                val = f"{int(n):,}" if n == int(n) else f"{n:.2f}"
        except (ValueError, TypeError):
            pass
        cards.append(
            f"<div style='background:#f8fafc;border:1px solid #e2e8f0;border-radius:10px;"
            f"padding:.9rem 1rem'>"
            f"<div style='font-size:.72rem;color:#64748b;text-transform:uppercase;"
            f"letter-spacing:.04em'>{label}</div>"
            f"<div style='font-size:1.2rem;font-weight:700;color:#1e40af;margin-top:.1rem'>{val}</div>"
            f"</div>"
        )
    if not cards:
        return "<p style='color:#6b7280'>No QUAST data available.</p>"
    return (f"<div style='display:grid;grid-template-columns:repeat(auto-fill,minmax(160px,1fr));gap:.75rem'>"
            + "".join(cards) + "</div>")


def _checkm2_section(cm: dict) -> str:
    """Renders CheckM2 genome completeness results."""
    if not cm:
        return "<p style='color:#6b7280'>CheckM2 was not run.</p>"
    if cm.get("error"):
        return f"<p style='color:#6b7280'>CheckM2: {cm['error']}</p>"
    comp = cm.get("completeness", "—")
    cont = cm.get("contamination", "—")
    model = cm.get("model_used", "—")
    try:
        comp_f = float(comp)
        comp_color = "#22c55e" if comp_f >= 90 else ("#f59e0b" if comp_f >= 70 else "#ef4444")
    except (ValueError, TypeError):
        comp_color = "#6b7280"
    try:
        cont_f = float(cont)
        cont_color = "#22c55e" if cont_f <= 5 else ("#f59e0b" if cont_f <= 10 else "#ef4444")
    except (ValueError, TypeError):
        cont_color = "#6b7280"
    return (
        f"<div style='display:grid;grid-template-columns:repeat(3,1fr);gap:.75rem'>"
        f"<div style='background:#f0fdf4;border:1px solid #bbf7d0;border-radius:10px;padding:1rem'>"
        f"<div style='font-size:.72rem;color:#64748b;text-transform:uppercase'>Completeness</div>"
        f"<div style='font-size:1.5rem;font-weight:700;color:{comp_color}'>{comp}%</div></div>"
        f"<div style='background:#fef2f2;border:1px solid #fecaca;border-radius:10px;padding:1rem'>"
        f"<div style='font-size:.72rem;color:#64748b;text-transform:uppercase'>Contamination</div>"
        f"<div style='font-size:1.5rem;font-weight:700;color:{cont_color}'>{cont}%</div></div>"
        f"<div style='background:#f8fafc;border:1px solid #e2e8f0;border-radius:10px;padding:1rem'>"
        f"<div style='font-size:.72rem;color:#64748b;text-transform:uppercase'>Model</div>"
        f"<div style='font-size:.9rem;font-weight:600;color:#1e40af;margin-top:.2rem'>{model}</div></div>"
        f"</div>"
    )


def _bandage_section(bandage: dict) -> str:
    """Embeds the Bandage assembly graph image as base64 HTML, or returns a placeholder."""
    if not bandage:
        return ""
    img_path = bandage.get("image_path")
    if not img_path or not Path(img_path).exists():
        err = bandage.get("error", "Graph image not available.")
        return (
            "<div class='section'>"
            "<div class='section-title'>&#x1F5FA; Assembly Graph (Bandage)</div>"
            "<div class='section-body'>"
            f"<p style='color:#6b7280'>{err}</p>"
            "</div></div>"
        )
    img_b64 = base64.b64encode(Path(img_path).read_bytes()).decode()
    return (
        "<div class='section'>"
        "<div class='section-title'>&#x1F5FA; Assembly Graph (Bandage)</div>"
        "<div class='section-body' style='text-align:center'>"
        f"<img src='data:image/png;base64,{img_b64}' "
        "style='max-width:100%;border-radius:8px;border:1px solid #e2e8f0' "
        "alt='Assembly graph'>"
        "<p style='font-size:.75rem;color:#94a3b8;margin-top:.5rem'>"
        "Each node is a contig; edges represent overlaps. "
        "Colour intensity reflects contig depth of coverage.</p>"
        "</div></div>"
    )


def _mobsuite_section(mob: dict) -> str:
    """Renders MOB-Suite plasmid detection results."""
    if not mob:
        return "<p style='color:#6b7280'>MOB-Suite was not run.</p>"
    if mob.get("error") and not mob.get("plasmids"):
        return f"<p style='color:#ef4444'>MOB-Suite error: {mob['error']}</p>"

    plasmids = mob.get("plasmids", [])
    if not plasmids:
        return "<p style='color:#22c55e'>&#x2705; No plasmids detected — likely a purely chromosomal isolate.</p>"

    rows = [
        "<table><thead><tr>"
        "<th>Plasmid ID</th><th>Replicon(s)</th><th>Mobility</th><th>MPF Type</th>"
        "<th>Contigs</th><th>Size</th>"
        "</tr></thead><tbody>"
    ]
    for p in plasmids:
        size = p.get("size_bp", 0)
        if size >= 1_000_000:
            size_str = f"{size/1_000_000:.2f} Mb"
        elif size >= 1_000:
            size_str = f"{size/1_000:.1f} kb"
        else:
            size_str = str(size) + " bp" if size else "—"
        rows.append(
            f"<tr>"
            f"<td><code>{p.get('id','—')}</code></td>"
            f"<td>{p.get('replicons','—')}</td>"
            f"<td>{p.get('mobility','—')}</td>"
            f"<td>{p.get('mpf','—')}</td>"
            f"<td>{p.get('contigs','—')}</td>"
            f"<td>{size_str}</td>"
            f"</tr>"
        )
    rows.append("</tbody></table>")
    return "\n".join(rows)


def _abricate_section(abricate: dict) -> str:
    """Renders Abricate (CARD + VFDB) results as a table."""
    if not abricate:
        return "<p style='color:#6b7280'>Abricate was not run.</p>"

    genes = abricate.get("genes", [])
    ran_vfdb  = abricate.get("ran_vfdb", False)
    top_genus = abricate.get("top_genus", "")

    vfdb_note = (
        f"<p style='margin-bottom:.8rem;font-size:.85rem;color:#6b7280'>"
        f"VFDB also searched (detected genus: <em>{top_genus}</em>).</p>"
        if ran_vfdb else
        f"<p style='margin-bottom:.8rem;font-size:.85rem;color:#6b7280'>"
        f"CARD database searched. VFDB skipped (genus not in VFDB or Kraken2 not run).</p>"
    )

    if not genes:
        return vfdb_note + "<p style='color:#22c55e'>&#x2705; No hits found in searched databases.</p>"

    # Group by database
    by_db: dict[str, list] = {}
    for g in genes:
        by_db.setdefault(g["database"], []).append(g)

    sections = [vfdb_note]
    for db_name, hits in by_db.items():
        db_color = "#1e40af" if db_name == "CARD" else "#7c3aed"
        sections.append(
            f"<h4 style='margin:.6rem 0 .4rem;color:{db_color}'>"
            f"{db_name} &mdash; {len(hits)} hit{'s' if len(hits)!=1 else ''}</h4>"
        )
        rows = ["<table><thead><tr>"
                "<th>Gene</th><th>Product</th><th>% Coverage</th>"
                "<th>% Identity</th><th>Resistance</th>"
                "</tr></thead><tbody>"]
        for g in hits:
            rows.append(
                f"<tr>"
                f"<td><strong>{g.get('gene','')}</strong></td>"
                f"<td style='font-size:.82rem'>{g.get('product','')}</td>"
                f"<td>{g.get('coverage','')}</td>"
                f"<td>{g.get('identity','')}</td>"
                f"<td style='font-size:.8rem;color:#7c3aed'>{g.get('resistance','')}</td>"
                f"</tr>"
            )
        rows.append("</tbody></table>")
        sections.append("\n".join(rows))

    return "\n".join(sections)


def _coverage_section(cov: dict) -> str:
    """Renders read-depth / coverage breadth metrics."""
    if not cov:
        return "<p style='color:#6b7280'>Coverage analysis was not run.</p>"
    if cov.get("error"):
        return f"<p style='color:#6b7280'>Coverage: {cov['error']}</p>"

    depth  = cov.get("mean_depth", "—")
    bread  = cov.get("breadth_1x_pct", "—")
    tlen   = cov.get("total_length_bp", 0)

    try:
        depth_f = float(depth)
        depth_color = "#22c55e" if depth_f >= 20 else ("#f59e0b" if depth_f >= 10 else "#ef4444")
        depth_label = f"{depth_f:.1f}x"
    except (ValueError, TypeError):
        depth_color, depth_label = "#6b7280", str(depth)

    try:
        bread_f = float(bread)
        bread_color = "#22c55e" if bread_f >= 90 else ("#f59e0b" if bread_f >= 70 else "#ef4444")
        bread_label = f"{bread_f:.1f}%"
    except (ValueError, TypeError):
        bread_color, bread_label = "#6b7280", str(bread)

    try:
        tlen_label = f"{int(tlen)/1e6:.2f} Mb"
    except (ValueError, TypeError):
        tlen_label = str(tlen)

    return (
        f"<div style='display:grid;grid-template-columns:repeat(3,1fr);gap:.75rem'>"
        f"<div style='background:#f0fdf4;border:1px solid #bbf7d0;border-radius:10px;padding:1rem'>"
        f"<div style='font-size:.72rem;color:#64748b;text-transform:uppercase'>Mean Depth</div>"
        f"<div style='font-size:1.5rem;font-weight:700;color:{depth_color}'>{depth_label}</div>"
        f"<div style='font-size:.72rem;color:#94a3b8;margin-top:.2rem'>≥20x recommended</div></div>"
        f"<div style='background:#eff6ff;border:1px solid #bfdbfe;border-radius:10px;padding:1rem'>"
        f"<div style='font-size:.72rem;color:#64748b;text-transform:uppercase'>Breadth ≥1x</div>"
        f"<div style='font-size:1.5rem;font-weight:700;color:{bread_color}'>{bread_label}</div>"
        f"<div style='font-size:.72rem;color:#94a3b8;margin-top:.2rem'>% assembly covered</div></div>"
        f"<div style='background:#f8fafc;border:1px solid #e2e8f0;border-radius:10px;padding:1rem'>"
        f"<div style='font-size:.72rem;color:#64748b;text-transform:uppercase'>Assembly Length</div>"
        f"<div style='font-size:1.5rem;font-weight:700;color:#1e40af'>{tlen_label}</div>"
        f"<div style='font-size:.72rem;color:#94a3b8;margin-top:.2rem'>mapped target</div></div>"
        f"</div>"
    )


def _serotyping_section(sero: dict) -> str:
    """Renders serotyping results (ECTyper / SISTR / Kleborate)."""
    if not sero:
        return ""
    if sero.get("skipped"):
        return f"<p style='color:#6b7280'>{sero.get('reason','Serotyping skipped.')}</p>"
    if sero.get("error"):
        return f"<p style='color:#ef4444'>Serotyping error: {sero['error']}</p>"

    tool    = sero.get("tool", "")
    species = sero.get("detected_species", "")

    header = (f"<p style='margin-bottom:1rem;font-size:.88rem;color:#374151'>"
              f"Tool: <strong>{tool}</strong> &nbsp;|&nbsp; Detected: <em>{species}</em></p>")

    cards = []

    if tool == "ECTyper":
        items = [
            ("Serotype",  sero.get("serotype", "—")),
            ("O Antigen", sero.get("o_type",   "—")),
            ("H Antigen", sero.get("h_type",   "—")),
            ("Quality",   sero.get("quality",  "—")),
        ]
    elif tool == "SISTR":
        items = [
            ("Serovar",     sero.get("serovar",   "—")),
            ("Serogroup",   sero.get("serogroup", "—")),
            ("H1 Antigen",  sero.get("h1",        "—")),
            ("H2 Antigen",  sero.get("h2",        "—")),
            ("O Antigen",   sero.get("o_antigen", "—")),
            ("cgMLST ST",   sero.get("cgmlst_st", "—")),
            ("QC Status",   sero.get("qc_status", "—")),
        ]
    elif tool == "Kleborate":
        items = [
            ("Species",          sero.get("species",         "—")),
            ("ST",               sero.get("st",              "—")),
            ("K Locus",          sero.get("k_type",          "—")),
            ("O Locus",          sero.get("o_type",          "—")),
            ("Virulence Score",  sero.get("virulence_score", "—")),
            ("Resistance Score", sero.get("resistance_score","—")),
        ]
    else:
        items = [(k, v) for k, v in sero.items()
                 if k not in ("tool", "detected_species") and not k.startswith("_")]

    for label, val in items:
        if val and val != "—":
            cards.append(
                f"<div style='background:#f8fafc;border:1px solid #e2e8f0;border-radius:10px;"
                f"padding:.9rem 1rem'>"
                f"<div style='font-size:.72rem;color:#64748b;text-transform:uppercase;"
                f"letter-spacing:.04em'>{label}</div>"
                f"<div style='font-size:1.1rem;font-weight:700;color:#1e40af;margin-top:.15rem'>{val}</div>"
                f"</div>"
            )

    grid = (f"<div style='display:grid;grid-template-columns:repeat(auto-fill,minmax(160px,1fr));"
            f"gap:.75rem'>" + "".join(cards) + "</div>") if cards else ""

    return header + grid


def _plasmidfinder_section(pf_hits: list[dict]) -> str:
    """Renders PlasmidFinder (Abricate) results."""
    if not pf_hits:
        return "<p style='color:#22c55e'>&#x2705; No plasmid replicons detected (PlasmidFinder).</p>"
    rows = [
        "<table><thead><tr>"
        "<th>Replicon</th><th>Product</th><th>% Coverage</th><th>% Identity</th><th>Accession</th>"
        "</tr></thead><tbody>"
    ]
    for g in pf_hits:
        rows.append(
            f"<tr>"
            f"<td><strong>{g.get('gene','')}</strong></td>"
            f"<td style='font-size:.82rem'>{g.get('product','')}</td>"
            f"<td>{g.get('coverage','')}</td>"
            f"<td>{g.get('identity','')}</td>"
            f"<td><code style='font-size:.78rem'>{g.get('accession','')}</code></td>"
            f"</tr>"
        )
    rows.append("</tbody></table>")
    return "\n".join(rows)


def _amr_table(genes: list[dict]) -> str:
    if not genes:
        return "<p style='color:#22c55e'>&#x2705; No AMR genes detected.</p>"

    rows = ["<table><thead><tr>"
            "<th>Gene</th><th>Class</th><th>Subclass</th>"
            "<th>% Identity</th><th>Method</th>"
            "</tr></thead><tbody>"]
    for g in genes:
        rows.append(
            f"<tr>"
            f"<td><strong>{g.get('gene','')}</strong></td>"
            f"<td>{g.get('class','')}</td>"
            f"<td>{g.get('subclass','')}</td>"
            f"<td>{g.get('identity','')}</td>"
            f"<td>{g.get('method','')}</td>"
            f"</tr>"
        )
    rows.append("</tbody></table>")
    return "\n".join(rows)


def generate_html_report(job_id: str,
                          pipeline_results: dict[str, Any],
                          ai_results: dict[str, Any],
                          job_error: str | None = None) -> str:
    """Returns the complete HTML report as a string."""

    sample      = pipeline_results.get("sample_name", "Unknown")
    read_type   = pipeline_results.get("read_type", "?")
    qc          = pipeline_results.get("qc", {})
    mlst        = pipeline_results.get("mlst", {})
    amr         = pipeline_results.get("amr", {})
    abricate    = pipeline_results.get("abricate") or {}
    mobsuite    = pipeline_results.get("mobsuite") or {}
    bandage     = pipeline_results.get("bandage") or {}
    quast       = pipeline_results.get("quast") or {}
    checkm2     = pipeline_results.get("checkm2") or {}
    kraken2     = pipeline_results.get("kraken2") or {}
    gctx        = pipeline_results.get("genomic_context") or {}
    host_depl   = pipeline_results.get("host_depletion") or {}
    coverage    = pipeline_results.get("coverage") or {}
    serotyping  = pipeline_results.get("serotyping") or {}
    pf_hits     = abricate.get("plasmidfinder", [])

    risk        = ai_results.get("risk_level", "UNKNOWN").upper()
    risk_color, risk_bg = RISK_COLORS.get(risk, RISK_COLORS["UNKNOWN"])

    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>LycianWay &mdash; {sample}</title>
  <style>
    * {{ box-sizing: border-box; margin: 0; padding: 0; }}
    body {{
      font-family: 'Segoe UI', system-ui, sans-serif;
      background: #f8fafc; color: #1e293b; line-height: 1.6;
    }}
    .header {{
      background: linear-gradient(135deg, #1e40af, #7c3aed);
      color: white; padding: 2rem 3rem;
    }}
    .header h1 {{ font-size: 1.8rem; font-weight: 700; }}
    .header p  {{ opacity: .8; margin-top: .3rem; }}
    .risk-badge {{
      display: inline-block; padding: .4rem 1.2rem;
      border-radius: 9999px; font-weight: 700;
      background: {risk_color}; color: white;
      font-size: 1rem; margin-top: .8rem;
    }}
    .container {{ max-width: 1100px; margin: 2rem auto; padding: 0 1.5rem; }}
    .section {{
      background: white; border-radius: 12px;
      box-shadow: 0 1px 4px rgba(0,0,0,.08);
      margin-bottom: 1.5rem; overflow: hidden;
    }}
    .section-title {{
      background: #f1f5f9; padding: .8rem 1.5rem;
      font-weight: 700; font-size: 1rem;
      border-bottom: 1px solid #e2e8f0;
      display: flex; align-items: center; gap: .5rem;
    }}
    .section-body {{ padding: 1.5rem; }}
    .grid-2 {{ display: grid; grid-template-columns: 1fr 1fr; gap: 1rem; }}
    .metric {{
      background: #f8fafc; border-radius: 8px;
      padding: 1rem 1.2rem; border: 1px solid #e2e8f0;
    }}
    .metric-label {{ font-size: .8rem; color: #64748b; text-transform: uppercase; }}
    .metric-value {{ font-size: 1.5rem; font-weight: 700; color: #1e40af; }}
    table {{
      width: 100%; border-collapse: collapse; margin-top: .5rem;
      font-size: .9rem;
    }}
    th {{ background: #f1f5f9; padding: .6rem .8rem; text-align: left;
          border-bottom: 2px solid #e2e8f0; }}
    td {{ padding: .55rem .8rem; border-bottom: 1px solid #f1f5f9; }}
    tr:hover td {{ background: #f8fafc; }}
    code {{
      background: #f1f5f9; padding: .15rem .4rem;
      border-radius: 4px; font-family: monospace; font-size: .85rem;
    }}
    .gene-block {{
      border: 1px solid #e2e8f0; border-radius: 8px;
      padding: 1rem; margin-bottom: 1rem;
    }}
    .gene-block h3 {{ margin-bottom: .5rem; color: #1e40af; }}
    .ai-text {{
      background: {risk_bg}; border-left: 4px solid {risk_color};
      padding: 1rem 1.2rem; border-radius: 0 8px 8px 0;
      margin-bottom: .8rem; line-height: 1.7;
    }}
    .footer {{
      text-align: center; padding: 2rem; color: #94a3b8;
      font-size: .85rem;
    }}
    @media print {{
      .header {{ -webkit-print-color-adjust: exact; print-color-adjust: exact; }}
    }}
  </style>
</head>
<body>

<div class="header">
  <div style="display:flex;align-items:flex-start;justify-content:space-between;gap:1rem;flex-wrap:wrap;">
    <div>
      <h1>&#x1F9EC; LycianWay &mdash; Genomic Analysis Report</h1>
      <p>Sample: <strong>{sample}</strong> &nbsp;|&nbsp; Job ID: {job_id} &nbsp;|&nbsp; {now}</p>
      <p>Read type: {read_type.upper()}</p>
      <div class="risk-badge">Risk: {risk}</div>
    </div>
    <a href="/" style="
      display:inline-block; margin-top:.4rem; padding:.55rem 1.2rem;
      background:rgba(255,255,255,0.15); color:white; border-radius:8px;
      font-size:.85rem; font-weight:600; text-decoration:none;
      border:1px solid rgba(255,255,255,0.3); white-space:nowrap;
      backdrop-filter:blur(4px);
    ">&#x2190; New Analysis</a>
  </div>
</div>

<div class="container">

  {_error_banner(job_error)}

  {_context_banner(gctx)}

  <!-- Kraken2 -->
  <div class="section">
    <div class="section-title">&#x1F9AB; Taxonomic Classification (Kraken2)</div>
    <div class="section-body">
      {_kraken2_section(kraken2)}
    </div>
  </div>

  <!-- Host Depletion (only rendered when host reads were removed) -->
  {_host_depletion_section(host_depl)}

  <!-- QC -->
  <div class="section">
    <div class="section-title">&#x1F4CA; Quality Control (QC)</div>
    <div class="section-body">
      <div class="grid-2">
        {''.join(
            f'<div class="metric"><div class="metric-label">{k.replace("_"," ").title()}</div>'
            f'<div class="metric-value">{v}</div></div>'
            for k, v in qc.items()
        ) or '<p style="color:#6b7280">No QC data available.</p>'}
      </div>
    </div>
  </div>

  <!-- Assembly QC -->
  <div class="section">
    <div class="section-title">&#x1F527; Assembly Quality Control</div>
    <div class="section-body">
      {_assembly_qc_section(pipeline_results.get("assembly_qc_checks", []))}
    </div>
  </div>

  <!-- Assembly Graph -->
  {_bandage_section(bandage)}

  <!-- Coverage -->
  <div class="section">
    <div class="section-title">&#x1F4CA; Read Depth &amp; Coverage</div>
    <div class="section-body">{_coverage_section(coverage)}</div>
  </div>

  <!-- QUAST -->
  <div class="section">
    <div class="section-title">&#x1F4D0; Assembly Quality (QUAST)</div>
    <div class="section-body">{_quast_section(quast)}</div>
  </div>

  <!-- CheckM2 -->
  <div class="section">
    <div class="section-title">&#x2705; Genome Completeness (CheckM2)</div>
    <div class="section-body">{_checkm2_section(checkm2)}</div>
  </div>

  <!-- MLST -->
  <div class="section">
    <div class="section-title">&#x1F52C; MLST Typing</div>
    <div class="section-body">
      <div class="grid-2">
        <div class="metric">
          <div class="metric-label">Scheme</div>
          <div class="metric-value">{mlst.get('scheme','&mdash;')}</div>
        </div>
        <div class="metric">
          <div class="metric-label">Sequence Type (ST)</div>
          <div class="metric-value">{mlst.get('st','&mdash;')}</div>
        </div>
      </div>
      {f'<p style="margin-top:.8rem">Alleles: <code>{", ".join(mlst.get("alleles",[]))}</code></p>' if mlst.get("alleles") else ''}
    </div>
  </div>

  <!-- AMR -->
  <div class="section">
    <div class="section-title">&#x1F48A; Antimicrobial Resistance Genes ({amr.get('count',0)} genes)</div>
    <div class="section-body">
      {_amr_table(amr.get('genes', []))}
    </div>
  </div>

  <!-- Abricate -->
  <div class="section">
    <div class="section-title">&#x1F9AB; Virulence &amp; Resistance Screen (Abricate)</div>
    <div class="section-body">
      {_abricate_section(abricate)}
    </div>
  </div>

  <!-- MOB-Suite -->
  <div class="section">
    <div class="section-title">&#x1F48A; Mobile Genetic Elements (MOB-Suite)</div>
    <div class="section-body">
      {_mobsuite_section(mobsuite)}
    </div>
  </div>

  <!-- PlasmidFinder -->
  <div class="section">
    <div class="section-title">&#x1F9EC; Plasmid Replicons (PlasmidFinder / Abricate)</div>
    <div class="section-body">
      {_plasmidfinder_section(pf_hits)}
    </div>
  </div>

  <!-- Serotyping -->
  <div class="section">
    <div class="section-title">&#x1F9EC; Serotyping</div>
    <div class="section-body">
      {_serotyping_section(serotyping)}
    </div>
  </div>

  <!-- AI Interpretation -->
  <div class="section">
    <div class="section-title">&#x1F916; AI Clinical Interpretation</div>
    <div class="section-body">
      <div class="ai-text">
        <strong>Species Prediction:</strong> {ai_results.get('species_prediction','&mdash;')}
      </div>
      <div class="ai-text">
        <strong>Clinical Significance:</strong><br>{ai_results.get('clinical_significance','&mdash;')}
      </div>
      <div class="ai-text">
        <strong>Resistance Profile:</strong><br>{ai_results.get('resistance_profile','&mdash;')}
      </div>
      <div class="ai-text">
        <strong>Treatment Implications:</strong><br>{ai_results.get('treatment_implications','&mdash;')}
      </div>
      <div class="ai-text">
        <strong>Epidemiology:</strong><br>{ai_results.get('epidemiology','&mdash;')}
      </div>
    </div>
  </div>

  <!-- Summary -->
  <div class="section">
    <div class="section-title">&#x1F4DD; Overall Summary</div>
    <div class="section-body">
      <div class="ai-text">{ai_results.get('summary','&mdash;')}</div>
      <p style="margin-top:1rem;font-size:.8rem;color:#94a3b8">
        &#x26A0;&#xFE0F; This report was generated by an automated analysis system.
        Expert review is required for clinical decision-making.
        All data will be automatically deleted after 24 hours.
      </p>
    </div>
  </div>

</div>

<div class="footer">
  LycianWay &bull; {now} &bull; Job: {job_id}
</div>

</body>
</html>"""

    return html


def save_report(job_id: str, html_content: str, out_dir: Path) -> Path:
    """Saves the report to disk and returns the path."""
    report_path = out_dir / "report.html"
    report_path.write_text(html_content, encoding="utf-8")
    return report_path
