"""
Claude API ile klinik yorum üretimi.
Pipeline sonuçlarını alır, yapılandırılmış bir yorum döndürür.
"""
import json
import logging
from typing import Any

import anthropic

from config import ANTHROPIC_API_KEY, CLAUDE_MODEL

logger = logging.getLogger("ai_interpreter")


def _build_prompt(results: dict[str, Any]) -> str:
    """Pipeline sonuçlarından klinik yorum için sistem promptu oluşturur."""

    mlst    = results.get("mlst", {})
    amr     = results.get("amr", {})
    qc      = results.get("qc", {})
    rtype   = results.get("read_type", "unknown")
    sample  = results.get("sample_name", "Bilinmeyen")

    amr_genes = amr.get("genes", [])
    amr_summary = "\n".join(
        f"  - {g['gene']} | Sınıf: {g['class']} | Kimlik: {g['identity']}%"
        for g in amr_genes
    ) or "  AMR geni tespit edilmedi."

    prompt = f"""Sen bir klinik mikrobiyoloji uzmanısın.
Aşağıdaki genomik analiz sonuçlarını incele ve kapsamlı bir klinik yorum yaz.

## Örnek Bilgisi
- Örnek adı: {sample}
- Okuma tipi: {rtype}

## QC Metrikleri
{json.dumps(qc, indent=2, ensure_ascii=False)}

## MLST Tiplemesi
- Şema: {mlst.get('scheme', 'Tespit edilemedi')}
- ST tipi: {mlst.get('st', 'Tespit edilemedi')}
- Alleller: {', '.join(mlst.get('alleles', []))}

## AMR Genleri ({amr.get('count', 0)} gen)
{amr_summary}

## İstenen Çıktı (JSON formatında yanıt ver):
{{
  "species_prediction": "Tahmin edilen tür veya yakın akraba",
  "clinical_significance": "Bu izolatın klinik önemi (2-3 cümle)",
  "resistance_profile": "Direnç profilinin özeti",
  "treatment_implications": "Tedavi seçenekleri ve öneriler",
  "epidemiology": "Epidemiyolojik önem (lineage, klonal kompleks vb.)",
  "pcr_targets": [
    {{
      "gene": "Hedef gen adı",
      "rationale": "Neden bu gen PCR için önerildi",
      "clinical_use": "Bu PCR'ın klinik/epidemiyolojik kullanımı"
    }}
  ],
  "risk_level": "LOW | MEDIUM | HIGH | CRITICAL",
  "summary": "Genel özet (İngilizce, yayına uygun dil)"
}}

Sadece JSON döndür, başka metin ekleme."""

    return prompt


def interpret(results: dict[str, Any]) -> dict[str, Any]:
    """
    Claude API'ye analiz sonuçlarını gönderir, yapılandırılmış yorum alır.
    API anahtarı yoksa fallback mesaj döndürür.
    """
    if not ANTHROPIC_API_KEY:
        logger.warning("ANTHROPIC_API_KEY ayarlanmamış — AI yorumu atlanıyor.")
        return {
            "species_prediction": "API anahtarı gerekli",
            "clinical_significance": "AI yorumu için ANTHROPIC_API_KEY ortam değişkenini ayarlayın.",
            "resistance_profile": "",
            "treatment_implications": "",
            "epidemiology": "",
            "pcr_targets": [],
            "risk_level": "UNKNOWN",
            "summary": "AI yorum servisi devre dışı.",
        }

    client = anthropic.Anthropic(api_key=ANTHROPIC_API_KEY)
    prompt = _build_prompt(results)

    try:
        message = client.messages.create(
            model=CLAUDE_MODEL,
            max_tokens=2048,
            messages=[{"role": "user", "content": prompt}],
        )
        raw = message.content[0].text.strip()

        # JSON bloğunu temizle (```json ... ``` varsa)
        if raw.startswith("```"):
            raw = raw.split("```")[1]
            if raw.startswith("json"):
                raw = raw[4:]
        raw = raw.strip()

        interpretation = json.loads(raw)
        logger.info("AI yorumu başarıyla alındı.")
        return interpretation

    except json.JSONDecodeError as e:
        logger.error("AI yanıtı JSON ayrıştırılamadı: %s", e)
        return {"summary": raw, "risk_level": "UNKNOWN", "pcr_targets": []}

    except anthropic.APIError as e:
        logger.error("Claude API hatası: %s", e)
        return {
            "summary": f"API hatası: {e}",
            "risk_level": "UNKNOWN",
            "pcr_targets": [],
        }
