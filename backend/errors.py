from __future__ import annotations
from dataclasses import dataclass

# ── Code constants ────────────────────────────────────────────────────────────
INVALID_FILE_FORMAT = "GP-101"
FILE_TOO_LARGE      = "GP-102"
NANOPLOT_FAILED     = "GP-201"
FASTP_FAILED        = "GP-202"
INSUFFICIENT_READS  = "GP-203"
FLYE_FAILED         = "GP-301"
SHOVILL_FAILED      = "GP-302"
ASSEMBLY_QC_FAILED  = "GP-303"
MLST_FAILED         = "GP-401"
AMRFINDER_FAILED    = "GP-402"
BAKTA_FAILED        = "GP-403"
DISK_INSUFFICIENT   = "GP-501"
TOOL_NOT_FOUND      = "GP-502"
STAGE_TIMEOUT       = "GP-503"
PANAROO_FAILED      = "GP-601"
IQTREE_FAILED       = "GP-602"
R_ANNOTATION_FAILED = "GP-603"
NO_ASSEMBLIES       = "GP-604"
NCBI_FETCH_FAILED   = "GP-605"

ERROR_LABELS: dict[str, str] = {
    INVALID_FILE_FORMAT: "Invalid file format",
    FILE_TOO_LARGE:      "File too large",
    NANOPLOT_FAILED:     "NanoPlot QC failed",
    FASTP_FAILED:        "FastP QC failed",
    INSUFFICIENT_READS:  "Insufficient reads after host depletion",
    FLYE_FAILED:         "Flye assembly failed",
    SHOVILL_FAILED:      "Shovill/Megahit assembly failed",
    ASSEMBLY_QC_FAILED:  "Assembly QC thresholds not met",
    MLST_FAILED:         "MLST typing failed",
    AMRFINDER_FAILED:    "AMRFinder profiling failed",
    BAKTA_FAILED:        "Bakta annotation failed",
    DISK_INSUFFICIENT:   "Insufficient disk space",
    TOOL_NOT_FOUND:      "Required tool not found",
    STAGE_TIMEOUT:       "Pipeline stage timed out",
    PANAROO_FAILED:      "Panaroo pan-genome analysis failed",
    IQTREE_FAILED:       "IQ-TREE phylogeny failed",
    R_ANNOTATION_FAILED: "R tree annotation failed",
    NO_ASSEMBLIES:       "No assemblies found for selected jobs",
    NCBI_FETCH_FAILED:   "NCBI reference genome fetch failed",
}


@dataclass
class PipelineError(Exception):
    code: str
    message: str
    detail: str
    stage: str
    recoverable: bool = False

    def __post_init__(self) -> None:
        super().__init__(self.code, self.message, self.detail, self.stage)

    def __str__(self) -> str:
        return f"[{self.code}] {self.message}: {self.detail[:200]}"
