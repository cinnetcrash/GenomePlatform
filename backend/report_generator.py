"""
HTML report generator.
Combines pipeline + AI results into a single HTML report.
"""
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


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
                          ai_results: dict[str, Any]) -> str:
    """Returns the complete HTML report as a string."""

    sample      = pipeline_results.get("sample_name", "Unknown")
    read_type   = pipeline_results.get("read_type", "?")
    qc          = pipeline_results.get("qc", {})
    mlst        = pipeline_results.get("mlst", {})
    amr         = pipeline_results.get("amr", {})
    abricate    = pipeline_results.get("abricate") or {}
    mobsuite    = pipeline_results.get("mobsuite") or {}
    kraken2     = pipeline_results.get("kraken2") or {}

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

  <!-- Kraken2 -->
  <div class="section">
    <div class="section-title">&#x1F9AB; Taxonomic Classification (Kraken2)</div>
    <div class="section-body">
      {_kraken2_section(kraken2)}
    </div>
  </div>

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
