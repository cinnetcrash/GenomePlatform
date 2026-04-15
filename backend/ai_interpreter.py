"""
Clinical interpretation using the Claude API.
Receives pipeline results and returns a structured interpretation.
"""
import json
import logging
from typing import Any

import anthropic

from config import ANTHROPIC_API_KEY, CLAUDE_MODEL

logger = logging.getLogger("ai_interpreter")


def _build_prompt(results: dict[str, Any]) -> str:
    """Builds the clinical interpretation prompt from pipeline results."""

    mlst    = results.get("mlst", {})
    amr     = results.get("amr", {})
    qc      = results.get("qc", {})
    rtype   = results.get("read_type", "unknown")
    sample  = results.get("sample_name", "Unknown")

    amr_genes = amr.get("genes", [])
    amr_summary = "\n".join(
        f"  - {g['gene']} | Class: {g['class']} | Identity: {g['identity']}%"
        for g in amr_genes
    ) or "  No AMR genes detected."

    prompt = f"""You are a clinical microbiology expert.
Review the following genomic analysis results and provide a comprehensive clinical interpretation.

## Sample Information
- Sample name: {sample}
- Read type: {rtype}

## QC Metrics
{json.dumps(qc, indent=2, ensure_ascii=False)}

## MLST Typing
- Scheme: {mlst.get('scheme', 'Not determined')}
- Sequence type: {mlst.get('st', 'Not determined')}
- Alleles: {', '.join(mlst.get('alleles', []))}

## AMR Genes ({amr.get('count', 0)} genes)
{amr_summary}

## Required Output (respond in JSON format):
{{
  "species_prediction": "Predicted species or closest relative",
  "clinical_significance": "Clinical significance of this isolate (2-3 sentences)",
  "resistance_profile": "Summary of resistance profile",
  "treatment_implications": "Treatment options and recommendations",
  "epidemiology": "Epidemiological significance (lineage, clonal complex, etc.)",
  "pcr_targets": [
    {{
      "gene": "Target gene name",
      "rationale": "Why this gene is recommended for PCR",
      "clinical_use": "Clinical/epidemiological application of this PCR"
    }}
  ],
  "risk_level": "LOW | MEDIUM | HIGH | CRITICAL",
  "summary": "Overall summary (English, publication-ready language)"
}}

Return only JSON, no additional text."""

    return prompt


def interpret(results: dict[str, Any]) -> dict[str, Any]:
    """
    Sends analysis results to the Claude API and returns a structured interpretation.
    Returns a fallback message if the API key is not configured.
    """
    if not ANTHROPIC_API_KEY:
        logger.warning("ANTHROPIC_API_KEY not set — skipping AI interpretation.")
        return {
            "species_prediction": "API key required",
            "clinical_significance": "Set the ANTHROPIC_API_KEY environment variable to enable AI interpretation.",
            "resistance_profile": "",
            "treatment_implications": "",
            "epidemiology": "",
            "pcr_targets": [],
            "risk_level": "UNKNOWN",
            "summary": "AI interpretation service is disabled.",
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

        # Strip markdown code fences if present (```json ... ```)
        if raw.startswith("```"):
            raw = raw.split("```")[1]
            if raw.startswith("json"):
                raw = raw[4:]
        raw = raw.strip()

        interpretation = json.loads(raw)
        logger.info("AI interpretation received successfully.")
        return interpretation

    except json.JSONDecodeError as e:
        logger.error("Failed to parse AI response as JSON: %s", e)
        return {"summary": raw, "risk_level": "UNKNOWN", "pcr_targets": []}

    except anthropic.APIError as e:
        logger.error("Claude API error: %s", e)
        return {
            "summary": f"API error: {e}",
            "risk_level": "UNKNOWN",
            "pcr_targets": [],
        }
