"""
HTML report generator for multi-sample comparisons.
"""
from datetime import datetime, timezone
from typing import Any


def _newick_to_ascii(newick: str) -> str:
    """
    Minimal Newick → indented text tree.
    Falls back to raw Newick if parsing fails.
    """
    try:
        # Strip branch lengths for display: (a:0.1,b:0.2) → (a,b)
        import re
        clean = re.sub(r':[0-9.e+-]+', '', newick).strip().rstrip(";")

        def _parse(s: str) -> list:
            """Parse nested parentheses into a list tree."""
            s = s.strip()
            if not s.startswith("("):
                return [s]
            inner = s[1:-1]
            depth, start, parts = 0, 0, []
            for i, ch in enumerate(inner):
                if ch == "(": depth += 1
                elif ch == ")": depth -= 1
                elif ch == "," and depth == 0:
                    parts.append(inner[start:i].strip())
                    start = i + 1
            parts.append(inner[start:].strip())
            return [_parse(p) for p in parts]

        def _render(node, prefix="", is_last=True) -> list[str]:
            connector = "└─ " if is_last else "├─ "
            if isinstance(node, str):
                return [prefix + connector + node]
            lines = [prefix + connector + "+"]
            child_prefix = prefix + ("   " if is_last else "│  ")
            for i, child in enumerate(node):
                lines += _render(child, child_prefix, i == len(node) - 1)
            return lines

        tree = _parse(clean)
        lines = []
        for i, child in enumerate(tree):
            lines += _render(child, "", i == len(tree) - 1)
        return "\n".join(lines)
    except Exception:
        return newick


def _distance_heatmap(distances: list[dict], samples: list[dict]) -> str:
    """Renders a CSS pairwise distance matrix."""
    if not distances:
        return "<p style='color:#6b7280'>Distance data not available.</p>"

    labels = [s["sample_name"] for s in samples]
    n = len(labels)

    # Build lookup
    dist_map: dict[tuple, float] = {}
    for d in distances:
        dist_map[(d["sample_a"], d["sample_b"])] = d["distance"]
        dist_map[(d["sample_b"], d["sample_a"])] = d["distance"]

    def _color(dist: float | None) -> str:
        if dist is None:
            return "#f1f5f9"
        if dist == 0:
            return "#e0f2fe"
        if dist < 0.001:
            return "#bbf7d0"   # very close
        if dist < 0.01:
            return "#fde68a"   # moderate
        if dist < 0.05:
            return "#fca5a5"   # distant
        return "#f87171"       # very distant

    col_w = max(80, min(150, 700 // (n + 1)))
    cell_style = (f"width:{col_w}px;height:36px;text-align:center;font-size:.75rem;"
                  f"font-family:monospace;border:1px solid #e2e8f0;")

    rows = [
        f"<div style='overflow-x:auto'><table style='border-collapse:collapse'>"
        f"<thead><tr><th style='{cell_style}background:#f1f5f9'></th>"
    ]
    for lbl in labels:
        short = lbl[:12] + "…" if len(lbl) > 12 else lbl
        rows.append(
            f"<th style='{cell_style}background:#f1f5f9;font-size:.7rem;"
            f"white-space:nowrap' title='{lbl}'>{short}</th>"
        )
    rows.append("</tr></thead><tbody>")

    for row_lbl in labels:
        short_row = row_lbl[:12] + "…" if len(row_lbl) > 12 else row_lbl
        rows.append(
            f"<tr><td style='{cell_style}background:#f1f5f9;font-size:.7rem;"
            f"white-space:nowrap;font-weight:600' title='{row_lbl}'>{short_row}</td>"
        )
        for col_lbl in labels:
            if row_lbl == col_lbl:
                rows.append(f"<td style='{cell_style}background:#e0f2fe;color:#1e40af;font-weight:700'>—</td>")
            else:
                d = dist_map.get((row_lbl, col_lbl))
                bg = _color(d)
                val = f"{d:.4f}" if d is not None else "?"
                rows.append(f"<td style='{cell_style}background:{bg}' title='Mash distance'>{val}</td>")
        rows.append("</tr>")

    rows.append("</tbody></table></div>")
    rows.append(
        "<div style='display:flex;gap:1rem;flex-wrap:wrap;margin-top:.75rem;font-size:.75rem'>"
        "<span><span style='display:inline-block;width:12px;height:12px;background:#bbf7d0;border-radius:2px;margin-right:4px'></span>&lt; 0.001 (very close)</span>"
        "<span><span style='display:inline-block;width:12px;height:12px;background:#fde68a;border-radius:2px;margin-right:4px'></span>0.001–0.01</span>"
        "<span><span style='display:inline-block;width:12px;height:12px;background:#fca5a5;border-radius:2px;margin-right:4px'></span>0.01–0.05</span>"
        "<span><span style='display:inline-block;width:12px;height:12px;background:#f87171;border-radius:2px;margin-right:4px'></span>&gt; 0.05 (distant)</span>"
        "</div>"
    )
    return "\n".join(rows)


def _amr_matrix(samples: list[dict]) -> str:
    """Renders AMR gene presence/absence as a colour matrix."""
    all_genes: list[str] = []
    for s in samples:
        for g in s.get("amr_genes", []):
            if g not in all_genes:
                all_genes.append(g)

    if not all_genes:
        return "<p style='color:#22c55e'>✅ No AMR genes detected in any sample.</p>"

    col_w = max(60, min(120, 700 // (len(all_genes) + 1)))
    cell  = (f"width:{col_w}px;height:34px;text-align:center;font-size:.72rem;"
             f"border:1px solid #e2e8f0;")

    rows = [
        f"<div style='overflow-x:auto'><table style='border-collapse:collapse'>"
        f"<thead><tr><th style='{cell}background:#f1f5f9;text-align:left;padding:0 .5rem'>Sample</th>"
    ]
    for gene in all_genes:
        short = gene[:10] + "…" if len(gene) > 10 else gene
        rows.append(
            f"<th style='{cell}background:#f1f5f9;font-size:.68rem;writing-mode:vertical-rl;"
            f"transform:rotate(180deg);height:80px;white-space:nowrap' title='{gene}'>{short}</th>"
        )
    rows.append("</tr></thead><tbody>")

    for s in samples:
        sample_genes = set(s.get("amr_genes", []))
        short_name = s["sample_name"][:20] + "…" if len(s["sample_name"]) > 20 else s["sample_name"]
        sname = s["sample_name"]
        rows.append(
            f"<tr><td style='{cell}background:#f8fafc;text-align:left;padding:0 .5rem;"
            f"font-weight:600;white-space:nowrap' title='{sname}'>{short_name}</td>"
        )
        for gene in all_genes:
            if gene in sample_genes:
                rows.append(f"<td style='{cell}background:#fee2e2;color:#b91c1c;font-weight:700'>✓</td>")
            else:
                rows.append(f"<td style='{cell}background:#f0fdf4;color:#86efac'>·</td>")
        rows.append("</tr>")

    rows.append("</tbody></table></div>")
    return "\n".join(rows)


def _summary_table(samples: list[dict]) -> str:
    """Renders per-sample metadata summary table."""
    def _badge(val: str, good_fn, warn_fn=None) -> str:
        try:
            v = float(str(val).replace(",", ""))
        except (ValueError, TypeError):
            return f"<span style='color:#6b7280'>{val}</span>"
        if good_fn(v):
            color = "#166534"; bg = "#f0fdf4"; border = "#bbf7d0"
        elif warn_fn and warn_fn(v):
            color = "#92400e"; bg = "#fffbeb"; border = "#fde68a"
        else:
            color = "#991b1b"; bg = "#fef2f2"; border = "#fecaca"
        return (f"<span style='background:{bg};border:1px solid {border};color:{color};"
                f"border-radius:4px;padding:.1rem .4rem;font-size:.82rem;font-weight:600'>{val}</span>")

    rows = [
        "<div style='overflow-x:auto'><table style='width:100%;border-collapse:collapse;font-size:.88rem'>",
        "<thead><tr style='background:#f1f5f9'>",
        "<th style='padding:.6rem .8rem;text-align:left;border-bottom:2px solid #e2e8f0'>Sample</th>",
        "<th style='padding:.6rem .8rem;text-align:left'>Organism (Kraken2)</th>",
        "<th style='padding:.6rem .8rem;text-align:left'>MLST</th>",
        "<th style='padding:.6rem .8rem;text-align:left'>Serotype</th>",
        "<th style='padding:.6rem .8rem;text-align:center'>Completeness</th>",
        "<th style='padding:.6rem .8rem;text-align:center'>Contam.</th>",
        "<th style='padding:.6rem .8rem;text-align:center'>Depth</th>",
        "<th style='padding:.6rem .8rem;text-align:center'>N50</th>",
        "<th style='padding:.6rem .8rem;text-align:center'>AMR genes</th>",
        "<th style='padding:.6rem .8rem;text-align:center'>Plasmids</th>",
        "</tr></thead><tbody>",
    ]
    for i, s in enumerate(samples):
        is_ref = s.get("is_reference", False)
        bg = "#f0f9ff" if is_ref else ("#ffffff" if i % 2 == 0 else "#f8fafc")
        td = f"padding:.55rem .8rem;border-bottom:1px solid #f1f5f9;background:{bg}"

        st_str = f"ST{s['mlst_st']}" if s['mlst_st'] not in ("—", "", None) else "—"
        scheme  = s['mlst_scheme'] if s['mlst_scheme'] not in ("—", "", None) else ""
        mlst_str = f"{scheme} {st_str}".strip() if scheme else st_str

        try:
            n50_val = int(str(s['n50']).replace(",",""))
            if n50_val >= 1_000_000: n50_disp = f"{n50_val/1_000_000:.2f} Mb"
            elif n50_val >= 1_000:   n50_disp = f"{n50_val/1_000:.0f} kb"
            else:                    n50_disp = str(n50_val)
        except (ValueError, TypeError):
            n50_disp = str(s['n50'])

        amr_count = s.get("amr_count", 0)
        amr_color = "#991b1b" if amr_count > 5 else ("#92400e" if amr_count > 0 else "#166534")
        amr_bg    = "#fef2f2" if amr_count > 5 else ("#fffbeb" if amr_count > 0 else "#f0fdf4")
        pf_count  = len(s.get("pf_replicons", []))

        ref_badge = (" <span style='background:#dbeafe;color:#1d4ed8;border-radius:4px;"
                     "padding:.05rem .35rem;font-size:.68rem;font-weight:600'>REF</span>"
                     if is_ref else "")
        rows.append(
            f"<tr>"
            f"<td style='{td};font-weight:600'>{s['sample_name']}{ref_badge}</td>"
            f"<td style='{td};font-style:italic;color:#374151'>{s['organism']}</td>"
            f"<td style='{td}'><code style='font-size:.8rem'>{mlst_str}</code></td>"
            f"<td style='{td}'>{s['serotype']}</td>"
            f"<td style='{td};text-align:center'>"
            f"{_badge(s['completeness'], lambda v: v>=90, lambda v: v>=70)}%</td>"
            f"<td style='{td};text-align:center'>"
            f"{_badge(s['contamination'], lambda v: v<=5, lambda v: v<=10)}%</td>"
            f"<td style='{td};text-align:center'>"
            f"{_badge(s['mean_depth'], lambda v: v>=20, lambda v: v>=10)}x</td>"
            f"<td style='{td};text-align:center;font-family:monospace'>{n50_disp}</td>"
            f"<td style='{td};text-align:center'>"
            f"<span style='background:{amr_bg};color:{amr_color};border-radius:9999px;"
            f"padding:.15rem .6rem;font-weight:700;font-size:.82rem'>{amr_count}</span></td>"
            f"<td style='{td};text-align:center'>{pf_count if pf_count else '—'}</td>"
            f"</tr>"
        )
    rows.append("</tbody></table></div>")
    return "\n".join(rows)


def generate_comparison_report(comp_id: str,
                                samples: list[dict[str, Any]],
                                newick: str | None,
                                distances: list[dict]) -> str:
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    n   = len(samples)

    tree_html = ""
    if newick:
        ascii_tree = _newick_to_ascii(newick)
        tree_html = (
            "<div class='section'>"
            "<div class='section-title'>🌳 Phylogenetic Tree (Mashtree / Mash distances)</div>"
            "<div class='section-body'>"
            f"<pre style='background:#f8fafc;padding:1.2rem;border-radius:8px;"
            f"font-size:.82rem;overflow-x:auto;line-height:1.6'>{ascii_tree}</pre>"
            "<p style='font-size:.75rem;color:#94a3b8;margin-top:.5rem'>"
            "Neighbour-joining tree based on Mash (MinHash) distances. "
            "Branch lengths are proportional to Mash distance (not scaled above).</p>"
            "</div></div>"
        )

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>LycianWay — Comparison {comp_id[:8]}</title>
  <style>
    * {{ box-sizing:border-box; margin:0; padding:0; }}
    body {{ font-family:'Segoe UI',system-ui,sans-serif; background:#f8fafc; color:#1e293b; line-height:1.6; }}
    .header {{ background:linear-gradient(135deg,#0f172a,#1e40af); color:white; padding:2rem 3rem; }}
    .header h1 {{ font-size:1.6rem; font-weight:700; }}
    .header p  {{ opacity:.8; margin-top:.3rem; font-size:.9rem; }}
    .container {{ max-width:1200px; margin:2rem auto; padding:0 1.5rem; }}
    .section {{ background:white; border-radius:12px; box-shadow:0 1px 4px rgba(0,0,0,.08);
                margin-bottom:1.5rem; overflow:hidden; }}
    .section-title {{ background:#f1f5f9; padding:.8rem 1.5rem; font-weight:700; font-size:1rem;
                      border-bottom:1px solid #e2e8f0; }}
    .section-body {{ padding:1.5rem; }}
    table {{ width:100%; border-collapse:collapse; }}
    th {{ background:#f1f5f9; padding:.6rem .8rem; text-align:left; border-bottom:2px solid #e2e8f0; }}
    td {{ padding:.55rem .8rem; border-bottom:1px solid #f1f5f9; }}
    code {{ background:#f1f5f9; padding:.15rem .4rem; border-radius:4px; font-family:monospace; font-size:.85rem; }}
    .footer {{ text-align:center; padding:2rem; color:#94a3b8; font-size:.85rem; }}
    a {{ color:#1e40af; text-decoration:none; }}
  </style>
</head>
<body>

<div class="header">
  <div style="display:flex;align-items:flex-start;justify-content:space-between;flex-wrap:wrap;gap:1rem">
    <div>
      <h1>🧬 LycianWay — Multi-Sample Comparison</h1>
      <p>Comparison ID: {comp_id[:8]} &nbsp;|&nbsp; {n} samples &nbsp;|&nbsp; {now}</p>
    </div>
    <a href="/" style="display:inline-block;padding:.55rem 1.2rem;background:rgba(255,255,255,.15);
       color:white;border-radius:8px;font-size:.85rem;font-weight:600;border:1px solid rgba(255,255,255,.3)">
      ← New Analysis</a>
  </div>
</div>

<div class="container">

  <!-- Summary table -->
  <div class="section">
    <div class="section-title">📋 Sample Summary</div>
    <div class="section-body">
      {_summary_table(samples)}
    </div>
  </div>

  <!-- Phylogenetic tree -->
  {tree_html}

  <!-- Distance matrix -->
  <div class="section">
    <div class="section-title">📏 Pairwise Mash Distance Matrix</div>
    <div class="section-body">
      {_distance_heatmap(distances, samples)}
      <p style="font-size:.75rem;color:#94a3b8;margin-top:.75rem">
        Mash distance ≈ 1 − Average Nucleotide Identity (ANI).
        Values &lt; 0.001 indicate highly similar isolates (likely same outbreak cluster).
      </p>
    </div>
  </div>

  <!-- AMR matrix -->
  <div class="section">
    <div class="section-title">💊 AMR Gene Presence / Absence</div>
    <div class="section-body">
      {_amr_matrix(samples)}
    </div>
  </div>

  <!-- Per-sample links -->
  <div class="section">
    <div class="section-title">🔗 Individual Reports</div>
    <div class="section-body">
      <div style="display:flex;flex-wrap:wrap;gap:.75rem">
        {''.join(
            f'<a href="/report/{s["job_id"]}" target="_blank" '
            f'style="background:#eff6ff;border:1px solid #bfdbfe;color:#1e40af;'
            f'border-radius:8px;padding:.5rem 1rem;font-size:.85rem;font-weight:600">'
            f'📄 {s["sample_name"]}</a>'
            for s in samples
        )}
      </div>
    </div>
  </div>

</div>

<div class="footer">LycianWay &bull; {now} &bull; Comparison: {comp_id[:8]}</div>
</body>
</html>"""
