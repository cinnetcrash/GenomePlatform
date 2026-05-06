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

        if not message.content:
            logger.error("Claude API returned empty content (stop_reason=%s)", message.stop_reason)
            return {
                "species_prediction": "N/A",
                "clinical_significance": f"AI returned no content (stop_reason: {message.stop_reason}). The model may have refused or been rate-limited.",
                "resistance_profile": "",
                "treatment_implications": "",
                "epidemiology": "",
                "risk_level": "UNKNOWN",
                "summary": "AI interpretation returned empty content.",
            }

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
        return {"summary": raw, "risk_level": "UNKNOWN"}

    except anthropic.APIError as e:
        logger.error("Claude API error: %s", e)
        return {"summary": f"API error: {e}", "risk_level": "UNKNOWN"}


_EXPLAIN_ERROR_SYSTEM = (
    "You are a bioinformatics platform assistant. "
    "A genomic analysis pipeline has failed. "
    "Given an error code, stage name, and raw error detail, "
    "explain in 2-3 plain sentences what went wrong and what the user can try. "
    "Be concise and practical. Do not repeat the error code."
)


def explain_error(error_code: str, stage: str, detail: str) -> str | None:
    if not ANTHROPIC_API_KEY:
        return None
    client = anthropic.Anthropic(api_key=ANTHROPIC_API_KEY)
    user_msg = (
        f"Error code: {error_code}\n"
        f"Pipeline stage: {stage}\n"
        f"Error detail: {detail[:500]}"
    )
    try:
        msg = client.messages.create(
            model=CLAUDE_MODEL, max_tokens=200,
            system=[{"type": "text", "text": _EXPLAIN_ERROR_SYSTEM,
                     "cache_control": {"type": "ephemeral"}}],
            messages=[{"role": "user", "content": user_msg}],
        )
        return msg.content[0].text.strip() if msg.content else None
    except Exception as exc:
        logger.warning("explain_error failed: %s", exc)
        return None


_ASSEMBLY_QC_SYSTEM = (
    "You are a clinical microbiology expert. "
    "Given bacterial genome assembly statistics and the most likely organism, "
    "write 1-2 sentences of clinical context: note any quality concerns "
    "(fragmented assembly, unusual size, GC anomaly) and what they might mean. "
    "Be concise. If all metrics look normal, say so briefly."
)


def interpret_assembly_qc(stats: dict, top_organism: str | None) -> str | None:
    if not ANTHROPIC_API_KEY:
        return None
    client = anthropic.Anthropic(api_key=ANTHROPIC_API_KEY)
    organism_line = f"Most likely organism: {top_organism}" if top_organism else "Organism: not determined"
    user_msg = (
        f"{organism_line}\n"
        f"N50: {stats.get('n50_bp', 0):,} bp\n"
        f"Contigs: {stats.get('total_contigs', 0)}\n"
        f"Total length: {stats.get('total_length_bp', 0):,} bp\n"
        f"GC content: {stats.get('gc_percent', 0)}%"
    )
    try:
        msg = client.messages.create(
            model=CLAUDE_MODEL, max_tokens=150,
            system=[{"type": "text", "text": _ASSEMBLY_QC_SYSTEM,
                     "cache_control": {"type": "ephemeral"}}],
            messages=[{"role": "user", "content": user_msg}],
        )
        return msg.content[0].text.strip() if msg.content else None
    except Exception as exc:
        logger.warning("interpret_assembly_qc failed: %s", exc)
        return None


_COMPARISON_SYSTEM = (
    "You are a public health microbiologist specialising in outbreak surveillance. "
    "Given a summary of a multi-sample genomic comparison (MLST types, AMR genes, "
    "and the closest sample pair by Mash distance), write one paragraph of "
    "epidemiological interpretation: note clonal clusters, shared AMR profiles, "
    "and any outlier samples. Be concise and actionable."
)


def summarise_comparison(
    sample_names: list[str],
    mlst_types: list[str],
    amr_matrix: dict[str, list[str]],
    closest_pair: tuple[str, str, float] | None,
) -> str | None:
    if not ANTHROPIC_API_KEY:
        return None
    client = anthropic.Anthropic(api_key=ANTHROPIC_API_KEY)
    amr_lines = "\n".join(
        f"  {name}: {', '.join(genes) or 'none'}"
        for name, genes in amr_matrix.items()
    )
    pair_line = (
        f"Closest pair: {closest_pair[0]} & {closest_pair[1]} "
        f"(Mash distance {closest_pair[2]:.5f})"
        if closest_pair else "Distance data not available"
    )
    user_msg = (
        f"Samples ({len(sample_names)}): {', '.join(sample_names)}\n"
        f"MLST types: {', '.join(mlst_types)}\n"
        f"AMR genes per sample:\n{amr_lines}\n"
        f"{pair_line}"
    )
    try:
        msg = client.messages.create(
            model=CLAUDE_MODEL, max_tokens=300,
            system=[{"type": "text", "text": _COMPARISON_SYSTEM,
                     "cache_control": {"type": "ephemeral"}}],
            messages=[{"role": "user", "content": user_msg}],
        )
        return msg.content[0].text.strip() if msg.content else None
    except Exception as exc:
        logger.warning("summarise_comparison failed: %s", exc)
        return None
