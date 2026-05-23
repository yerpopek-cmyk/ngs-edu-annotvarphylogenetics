#!/usr/bin/env bash
# =============================================================================
# 03_filter.sh — Variant normalization and soft-filtering
#
# STAGE OVERVIEW
# ──────────────
# Raw variant calls contain a mix of true variants and technical artifacts.
# This step applies two layers of processing:
#
#   1. Normalization (bcftools norm)
#      Converts variants to a canonical representation so that the same
#      biological variant is never described two different ways. Two key steps:
#
#      a) Left-alignment of indels:
#         AATTT / A-TTT  →  AAT-T / A--- (move gap as far LEFT as possible)
#         This matters because ATTTG:ATTG and ATTTG:ATG are the same deletion
#         but look different without normalization.
#
#      b) Split multi-allelic sites (-m -any):
#         chr1 100 . A T,G  →  two records: A→T and A→G
#         Downstream tools (VEP, annotation databases) generally expect
#         one ALT per record.
#
#   2. Soft-filtering (bcftools filter --soft-filter)
#      Rather than REMOVING variants that fail quality criteria, we MARK them
#      with a label in the FILTER column. The PASS value means a variant
#      passed all filters. This approach is preferred in clinical settings
#      because it preserves the full dataset for audit and re-analysis.
#
#      Filter tags applied:
#        LowQual  — QUAL < MIN_QUAL OR depth outside [MIN_DP, MAX_DP]
#        LowAB    — allele balance outside [MIN_AB, MAX_AB] for heterozygotes
#
# WHY SOFT-FILTER INSTEAD OF HARD-FILTER?
# ─────────────────────────────────────────
# Hard filtering deletes variants permanently. If you later discover your
# depth threshold was too strict, you must re-run the caller. Soft filtering
# adds a label; downstream tools can apply `bcftools view -f PASS` when they
# only want high-confidence variants, while the full file remains available.
#
# OUTPUT
#   ${RUN_DIR}/3_vcf_filtered/${SAMPLE_ID}.filtered.vcf.gz  — normalized, soft-filtered
#   ${RUN_DIR}/3_vcf_filtered/${SAMPLE_ID}.stats.txt        — bcftools stats output
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
source "${REPO_ROOT}/shared/logging.sh"

banner "Step 03 — Normalization and Soft-Filtering"

VCF_IN="${RUN_DIR}/2_vcf_raw/${SAMPLE_ID}.raw.vcf.gz"
VCF_OUT="${RUN_DIR}/3_vcf_filtered/${SAMPLE_ID}.filtered.vcf.gz"
STATS_OUT="${RUN_DIR}/3_vcf_filtered/${SAMPLE_ID}.stats.txt"

require_file "$VCF_IN"
require_file "$REF_FASTA"
require_cmd  bcftools
require_cmd  tabix

STEP_START=$(date +%s)

log_info "Input VCF    : $VCF_IN"
log_info "QUAL cutoff  : $MIN_QUAL"
log_info "Depth range  : $MIN_DP – $MAX_DP"
log_info "Allele bal.  : $MIN_AB – $MAX_AB"

# ── Single-pass normalization + soft-filtering ────────────────────────────────
# We chain two bcftools filter commands in a pipe to keep it a single pass.
# Mode + (--mode +) means "add a tag" rather than replacing PASS if the
# variant already has another filter tag.
#
# The allele balance expression uses FreeBayes FORMAT fields:
#   AO = ALT observation count per sample
#   RO = REF observation count per sample
# The ratio AO / (RO + AO) is the allele balance (AB).
# For a true heterozygous diploid variant we expect AB ≈ 0.5.
# Very low (< MIN_AB) or very high (> MAX_AB) values suggest:
#   - Low AB: not enough reads support the ALT; possible sequencing artifact
#   - High AB: suspiciously close to 1.0; may be a strand-bias artifact

bcftools norm \
        --fasta-ref "$REF_FASTA" \
        -m -any \
        "$VCF_IN" \
    | bcftools filter \
        --soft-filter LowQual \
        --mode + \
        --exclude "QUAL < ${MIN_QUAL} || INFO/DP < ${MIN_DP} || INFO/DP > ${MAX_DP}" \
    | bcftools filter \
        --soft-filter LowAB \
        --mode + \
        --exclude "GT[*]=\"het\" && (FORMAT/AO[*:0] / (FORMAT/RO[*:0] + FORMAT/AO[*:0]) < ${MIN_AB} || \
                                     FORMAT/AO[*:0] / (FORMAT/RO[*:0] + FORMAT/AO[*:0]) > ${MAX_AB})" \
        -Oz \
        -o "$VCF_OUT"

tabix -p vcf "$VCF_OUT"

# ── Statistics ────────────────────────────────────────────────────────────────
step 1 "Generating VCF statistics"

# bcftools stats reports:
#   SN lines: summary counts (SNPs, MNPs, indels, multiallelic, transitions, transversions)
#   TSTV: Ti/Tv ratio — for a well-calibrated human WGS call set this should be ≈ 2.0–2.2
#         A lower value (e.g., 1.5) often indicates excess false positives
#   IDD: indel size distribution
#   ST: per-chromosome counts

bcftools stats "$VCF_OUT" | tee "$STATS_OUT" | grep '^SN'

TOTAL=$(count_vcf "$VCF_OUT")
PASS=$(bcftools view -f PASS -H "$VCF_OUT" | wc -l | tr -d ' ')
SNP=$(bcftools view -f PASS -v snps   -H "$VCF_OUT" | wc -l | tr -d ' ')
INDEL=$(bcftools view -f PASS -v indels -H "$VCF_OUT" | wc -l | tr -d ' ')

echo ""
log_info "Total variants : $TOTAL"
log_info "PASS variants  : $PASS  (SNPs: $SNP, Indels: $INDEL)"
log_info "Filtered out   : $(( TOTAL - PASS ))"

# Quick Ti/Tv check — grab it from the stats file
TSTV=$(grep '^TSTV' "$STATS_OUT" | awk '{printf "%.2f", $5}' || echo "N/A")
log_info "Ti/Tv ratio    : $TSTV  (expected ≈ 2.0–2.2 for WGS)"

log_success "Filtering done in $(elapsed $STEP_START)"
log_info    "Output VCF: $VCF_OUT"
