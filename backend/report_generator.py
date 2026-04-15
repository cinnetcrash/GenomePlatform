"""
HTML rapor üretici.
Pipeline + AI + Primer sonuçlarından tek bir HTML rapor oluşturur.
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


def _primer_rows(primer_list: list[dict]) -> str:
    if not primer_list:
        return "<p style='color:#6b7280'>Primer tasarımlayacak veri bulunamadı.</p>"

    rows = []
    for gene_result in primer_list:
        gene   = gene_result.get("gene", "?")
        note   = gene_result.get("note", "")
        n_seqs = gene_result.get("n_sequences", 0)
        cons   = gene_result.get("conserved_region_bp")
        pairs  = gene_result.get("primers", [])

        rows.append(f"""
        <div class="gene-block">
          <h3>🧬 {gene}</h3>
          <p><small>Dizi sayısı: {n_seqs} |
             Korunmuş bölge: {cons or '—'} bp
             {f'| <em>{note}</em>' if note else ''}</small></p>
        """)

        if pairs:
            rows.append("""
          <table>
            <thead>
              <tr>
                <th>#</th><th>Yön</th><th>Dizi (5'→3')</th>
                <th>Tm (°C)</th><th>GC%</th><th>Amplikon (bp)</th>
              </tr>
            </thead><tbody>""")
            for p in pairs:
                rows.append(f"""
              <tr>
                <td rowspan="2">{p['pair']}</td>
                <td>Forward</td>
                <td><code>{p['forward']}</code></td>
                <td>{p['fwd_tm'] or '—'}</td>
                <td>{p['fwd_gc'] or '—'}</td>
                <td rowspan="2">{p['amplikon'] or '—'}</td>
              </tr>
              <tr>
                <td>Reverse</td>
                <td><code>{p['reverse']}</code></td>
                <td>{p['rev_tm'] or '—'}</td>
                <td>{p['rev_gc'] or '—'}</td>
              </tr>""")
            rows.append("</tbody></table>")
        else:
            rows.append("<p><em>Bu gen için primer tasarlanamadı.</em></p>")

        rows.append("</div>")
    return "\n".join(rows)


def _amr_table(genes: list[dict]) -> str:
    if not genes:
        return "<p style='color:#22c55e'>✅ AMR geni tespit edilmedi.</p>"

    rows = ["<table><thead><tr>"
            "<th>Gen</th><th>Sınıf</th><th>Alt Sınıf</th>"
            "<th>%Kimlik</th><th>Yöntem</th>"
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
                          primer_results: list[dict]) -> str:
    """Tam HTML raporu string olarak döndürür."""

    sample      = pipeline_results.get("sample_name", "Bilinmeyen")
    read_type   = pipeline_results.get("read_type", "?")
    qc          = pipeline_results.get("qc", {})
    mlst        = pipeline_results.get("mlst", {})
    amr         = pipeline_results.get("amr", {})

    risk        = ai_results.get("risk_level", "UNKNOWN").upper()
    risk_color, risk_bg = RISK_COLORS.get(risk, RISK_COLORS["UNKNOWN"])

    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    # AI PCR hedefleri
    pcr_targets = ai_results.get("pcr_targets", [])
    pcr_rows = ""
    for t in pcr_targets:
        pcr_rows += (
            f"<tr><td><strong>{t.get('gene','')}</strong></td>"
            f"<td>{t.get('rationale','')}</td>"
            f"<td>{t.get('clinical_use','')}</td></tr>"
        )

    html = f"""<!DOCTYPE html>
<html lang="tr">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>GenomePlatform — {sample}</title>
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
  <h1>🧬 Genomik Analiz Raporu</h1>
  <p>Örnek: <strong>{sample}</strong> &nbsp;|&nbsp; İş ID: {job_id} &nbsp;|&nbsp; {now}</p>
  <p>Okuma tipi: {read_type.upper()}</p>
  <div class="risk-badge">Risk: {risk}</div>
</div>

<div class="container">

  <!-- QC -->
  <div class="section">
    <div class="section-title">📊 Kalite Kontrol (QC)</div>
    <div class="section-body">
      <div class="grid-2">
        {''.join(
            f'<div class="metric"><div class="metric-label">{k.replace("_"," ").title()}</div>'
            f'<div class="metric-value">{v}</div></div>'
            for k, v in qc.items()
        ) or '<p style="color:#6b7280">QC verisi bulunamadı.</p>'}
      </div>
    </div>
  </div>

  <!-- MLST -->
  <div class="section">
    <div class="section-title">🔬 MLST Tiplemesi</div>
    <div class="section-body">
      <div class="grid-2">
        <div class="metric">
          <div class="metric-label">Şema</div>
          <div class="metric-value">{mlst.get('scheme','—')}</div>
        </div>
        <div class="metric">
          <div class="metric-label">Sekans Tipi (ST)</div>
          <div class="metric-value">{mlst.get('st','—')}</div>
        </div>
      </div>
      {f'<p style="margin-top:.8rem">Alleller: <code>{", ".join(mlst.get("alleles",[]))}</code></p>' if mlst.get("alleles") else ''}
    </div>
  </div>

  <!-- AMR -->
  <div class="section">
    <div class="section-title">💊 Antimikrobiyal Direnç Genleri ({amr.get('count',0)} gen)</div>
    <div class="section-body">
      {_amr_table(amr.get('genes', []))}
    </div>
  </div>

  <!-- AI Yorum -->
  <div class="section">
    <div class="section-title">🤖 AI Klinik Yorumu</div>
    <div class="section-body">
      <div class="ai-text">
        <strong>Tür Tahmini:</strong> {ai_results.get('species_prediction','—')}
      </div>
      <div class="ai-text">
        <strong>Klinik Önem:</strong><br>{ai_results.get('clinical_significance','—')}
      </div>
      <div class="ai-text">
        <strong>Direnç Profili:</strong><br>{ai_results.get('resistance_profile','—')}
      </div>
      <div class="ai-text">
        <strong>Tedavi Önerileri:</strong><br>{ai_results.get('treatment_implications','—')}
      </div>
      <div class="ai-text">
        <strong>Epidemiyoloji:</strong><br>{ai_results.get('epidemiology','—')}
      </div>
      {f'''
      <h4 style="margin-top:1rem">Önerilen PCR Hedefleri</h4>
      <table>
        <thead><tr><th>Gen</th><th>Neden?</th><th>Klinik Kullanım</th></tr></thead>
        <tbody>{pcr_rows}</tbody>
      </table>''' if pcr_rows else ''}
    </div>
  </div>

  <!-- Primerler -->
  <div class="section">
    <div class="section-title">🧪 PCR Primer Tasarımı</div>
    <div class="section-body">
      {_primer_rows(primer_results)}
    </div>
  </div>

  <!-- Özet -->
  <div class="section">
    <div class="section-title">📝 Genel Özet</div>
    <div class="section-body">
      <div class="ai-text">{ai_results.get('summary','—')}</div>
      <p style="margin-top:1rem;font-size:.8rem;color:#94a3b8">
        ⚠️ Bu rapor otomatik analiz sistemi tarafından üretilmiştir.
        Klinik kararlar için uzman değerlendirmesi gereklidir.
        Veriler 24 saat sonra otomatik silinecektir.
      </p>
    </div>
  </div>

</div>

<div class="footer">
  GenomePlatform • {now} • Job: {job_id}
</div>

</body>
</html>"""

    return html


def save_report(job_id: str, html_content: str, out_dir: Path) -> Path:
    """Raporu diske kaydeder, yolunu döndürür."""
    report_path = out_dir / "report.html"
    report_path.write_text(html_content, encoding="utf-8")
    return report_path
