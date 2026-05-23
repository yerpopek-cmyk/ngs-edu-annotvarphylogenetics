#!/usr/bin/env bash
# =============================================================================
# pipelines/variants/run_all.sh — Variant calling pipeline orchestrator
#
# WHAT THIS PIPELINE DOES
# ═══════════════════════
# End-to-end variant discovery for a single sample, in five steps:
#
#   01_align.sh   — Align FASTQ reads to the reference genome (BWA + samtools)
#   02_call.sh    — Call SNVs and small indels (FreeBayes)
#   03_filter.sh  — Normalize and soft-filter variants (bcftools)
#   04_annotate.sh — Predict functional consequences (Ensembl VEP)
#   05_report.py  — Prioritize and generate a Markdown report
#
# BACKGROUND: WHAT IS VARIANT CALLING?
# ─────────────────────────────────────
# We start with raw sequencing reads (FASTQ), align them to a reference genome,
# and then statistically determine which positions in the genome differ from the
# reference in our sample. Each position is represented by:
#
#   QUAL = −10 × log₁₀(P(variant is wrong))
#
# So QUAL = 30 means we expect 1 false call per 1000 variants.
#
# USAGE
# ═════
#   conda activate variants
#   bash run_all.sh [--from STEP] [--to STEP] [--config FILE] [--run-dir DIR]
#
# OPTIONS
#   --from N   Start from step N (e.g., --from 03 to re-run from filtering)
#   --to   N   Stop after step N
#   --config   Alternative config.sh path
#   --run-dir  Existing or custom run directory (needed when resuming a prior run)
#   -h, --help Show this help
#
# EXAMPLE (re-run just annotation + report):
#   bash run_all.sh --from 04 --run-dir outputs/run_YYYYMMDD_HHMMSS
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

source "${REPO_ROOT}/shared/logging.sh"

# --- Defaults ---
CONFIG_FILE="${SCRIPT_DIR}/config.sh"
FROM_STEP="01"
TO_STEP="05"
RUN_DIR_OVERRIDE=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --from)    FROM_STEP="$2";         shift 2 ;;
        --to)      TO_STEP="$2";           shift 2 ;;
        --config)  CONFIG_FILE="$2";       shift 2 ;;
        --run-dir) RUN_DIR_OVERRIDE="$2";  shift 2 ;;
        -h|--help)
            sed -n '/^# USAGE/,/^# =============================================================================/p' "$0" \
                | sed '/^# =============================================================================/d; s/^# \?//'
            exit 0 ;;
        *) die "Unknown argument: $1  (use --help)" ;;
    esac
done

[[ -f "$CONFIG_FILE" ]] || die "Config file not found: $CONFIG_FILE"
source "$CONFIG_FILE"

# Export for child scripts
export SCRIPT_DIR CONFIG_FILE SAMPLE_ID REF_FASTA READS_R1 READS_R2
export BAM_INPUT THREADS MEM_SORT MIN_AB MIN_QUAL MIN_DP MAX_DP MAX_AB
export VEP_CACHE_DIR VEP_ASSEMBLY REPORT_TOP_N

# Create a timestamped run directory, or reuse one when resuming.
if [[ -n "$RUN_DIR_OVERRIDE" ]]; then
    export RUN_DIR="$RUN_DIR_OVERRIDE"
else
    TIMESTAMP=$(date +%Y%m%d_%H%M%S)
    export RUN_DIR="${OUTROOT}/run_${TIMESTAMP}"
fi
mkdir -p "${RUN_DIR}"/{1_bams,2_vcf_raw,3_vcf_filtered,4_reports}

# Tee all output to a log file
LOG_FILE="${RUN_DIR}/pipeline.log"
exec > >(tee -i "$LOG_FILE") 2>&1

banner "NGS Variant Calling Pipeline"
log_info "Sample ID    : $SAMPLE_ID"
log_info "Reference    : $REF_FASTA"
[[ -n "$BAM_INPUT" ]] && log_info "BAM input    : $BAM_INPUT (skipping alignment)"
log_info "Threads      : $THREADS"
log_info "Run dir      : $RUN_DIR"
log_info "Log file     : $LOG_FILE"
log_info "Steps        : ${FROM_STEP} → ${TO_STEP}"

PIPELINE_START=$(date +%s)

# Helper: run a step script if its number is within [FROM_STEP, TO_STEP]
run_step() {
    local num="$1"     # e.g. "02"
    local script="$2"  # e.g. "02_call.sh"
    if [[ "$num" -ge "${FROM_STEP#0}" && "$num" -le "${TO_STEP#0}" ]]; then
        bash "${SCRIPT_DIR}/${script}"
    else
        log_warn "Skipping step ${num} (outside range ${FROM_STEP}–${TO_STEP})"
    fi
}

run_report() {
    local num=5
    local annotated_vcf="${RUN_DIR}/3_vcf_filtered/${SAMPLE_ID}.annotated.vcf.gz"
    if [[ "$num" -ge "${FROM_STEP#0}" && "$num" -le "${TO_STEP#0}" ]]; then
        require_cmd python3
        require_file "$annotated_vcf"
        python3 "${SCRIPT_DIR}/05_report.py" \
            "$annotated_vcf" \
            "${RUN_DIR}/4_reports" \
            --top "$REPORT_TOP_N" \
            --sample-id "$SAMPLE_ID"
    else
        log_warn "Skipping step 05 (outside range ${FROM_STEP}–${TO_STEP})"
    fi
}

run_step 1 "01_align.sh"
run_step 2 "02_call.sh"
run_step 3 "03_filter.sh"
run_step 4 "04_annotate.sh"
run_report

banner "Pipeline Complete"
log_success "All steps finished in $(elapsed $PIPELINE_START)"
echo ""
echo "  Output directory: $RUN_DIR"
echo "  ├── 1_bams/                    — aligned, deduplicated BAM"
echo "  ├── 2_vcf_raw/                 — raw FreeBayes calls"
echo "  ├── 3_vcf_filtered/            — filtered + VEP-annotated VCF"
echo "  └── 4_reports/"
echo "      ├── ${SAMPLE_ID}.prioritized.tsv"
echo "      └── ${SAMPLE_ID}.report.md"
