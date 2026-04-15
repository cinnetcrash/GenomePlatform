# 🧬 GenomePlatform

**Automated bacterial genomic analysis · AI clinical interpretation · PCR primer design**

[![CI](https://github.com/cinnetcrash/GenomePlatform/actions/workflows/ci.yml/badge.svg)](https://github.com/cinnetcrash/GenomePlatform/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.11](https://img.shields.io/badge/Python-3.11-blue)](https://python.org)

> Upload a single FASTQ file. Get a full clinical genomic report in ~15–45 minutes.

🌐 **Landing page:** https://cinnetcrash.github.io/GenomePlatform

---

## What it does

1. **Auto-detects** MinION (long reads) vs Illumina (short reads)
2. **QC** — NanoPlot or FastP
3. **Assembly** — Flye (MinION) or Shovill/SPAdes (Illumina)
4. **MLST typing** — 100+ schemes via CGE MLST
5. **AMR profiling** — AMRFinderPlus
6. **Annotation** — Bakta
7. **AI interpretation** — Claude synthesises findings into a clinical narrative
8. **PCR primer design** — MAFFT conserved regions → Primer3
9. **HTML report** — downloadable, self-contained

All uploaded data is **automatically deleted after 24 hours**.

---

## Quick start

### Option A — Docker (recommended)

```bash
git clone https://github.com/cinnetcrash/GenomePlatform
cd GenomePlatform
cp .env.example .env          # add your ANTHROPIC_API_KEY
docker compose build
docker compose up -d
open http://localhost:8000
```

### Option B — Direct (requires conda environments)

```bash
git clone https://github.com/cinnetcrash/GenomePlatform
cd GenomePlatform
export ANTHROPIC_API_KEY='sk-ant-...'
bash start.sh
```

---

## System requirements

| | Minimum | Recommended |
|---|---|---|
| CPU | 8 cores | 16+ cores |
| RAM | 16 GB | 32 GB |
| Disk | 50 GB | 500 GB SSD |

**Required tools** (auto-installed via Docker):
`fastp` · `NanoPlot` · `Flye` · `Shovill` · `mlst` · `AMRFinderPlus` · `Bakta` · `MAFFT` · `Primer3`

---

## Project structure

```
GenomePlatform/
├── backend/
│   ├── main.py             # FastAPI app, all routes
│   ├── config.py           # Central configuration
│   ├── security.py         # File validation, path safety
│   ├── database.py         # SQLite job tracking
│   ├── cleanup.py          # 24h auto-deletion
│   ├── pipeline.py         # Read detection + analysis stages
│   ├── ai_interpreter.py   # Claude API integration
│   ├── primer_designer.py  # MAFFT + Primer3
│   └── report_generator.py # HTML report
├── frontend/
│   └── templates/index.html
├── docker/
│   ├── Dockerfile
│   └── environment.yml
├── docs/                   # GitHub Pages landing page
├── docker-compose.yml
├── start.sh
└── requirements.txt
```

---

## Security

- File magic byte validation (no fake FASTQ.gz)
- Path traversal protection
- Rate limiting (5 uploads/min per IP)
- Max 3 concurrent jobs per IP
- `subprocess` always uses list args (`shell=False`)
- IP addresses stored as SHA-256 hashes (GDPR)
- All data auto-deleted after 24 hours

---

## License

MIT © Gültekin Ünal — [gultekinnunal@gmail.com](mailto:gultekinnunal@gmail.com)

> ⚠️ Research use only. Not validated for clinical diagnostics.
