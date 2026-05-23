#!/usr/bin/env bash
# =============================================================================
# 02_call.sh — Variant calling with FreeBayes
#
# WHAT IS A VARIANT CALLER?
# ─────────────────────────
# A variant caller reads the pileup of aligned sequencing reads at every
# position of the genome and asks: "Is the difference from the reference
# here due to a real genetic variant, or just a sequencing error?"
#
# HOW FREEBAYES WORKS
# ───────────────────
# FreeBayes is a haplotype-based caller: instead of looking at each position
# independently, it constructs small local haplotypes from the reads and
# evaluates the likelihood of observing the data under each possible genotype.
#
# It uses a Bayesian framework:
#   P(genotype | reads) ∝ P(reads | genotype) × P(genotype)
#
# The output is a Phred-scaled quality score:
#   QUAL = −10 × log₁₀( P(variant is wrong) )
#   QUAL 30  ⟹  1 error per 1000 calls
#   QUAL 60  ⟹  1 error per 1 000 000 calls
#
# VCF FORMAT FIELDS (what FreeBayes writes):
#   QUAL      — overall variant quality
#   GT        — genotype: 0/0 (ref), 0/1 (het), 1/1 (hom-alt)
#   DP        — total read depth at this position
#   AO / RO   — alternative / reference observation counts
#   GQ        — genotype quality (Phred)
#
# KEY PARAMETERS (all configurable in config.sh)
#   --min-alternate-count   — minimum number of ALT reads required
#   --min-alternate-fraction — minimum ALT / (REF + ALT); i.e., allele balance
#   --min-base-quality      — ignore bases with Phred quality below this value
#   --ploidy 2              — assumed ploidy of the sample
#
# OUTPUT
#   ${RUN_DIR}/2_vcf_raw/${SAMPLE_ID}.raw.vcf.gz
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
source "${REPO_ROOT}/shared/logging.sh"

banner "Step 02 — Variant Calling with FreeBayes"

# Use the pre-made BAM if provided, otherwise use the output of step 01
BAM_IN="${BAM_INPUT:-${RUN_DIR}/1_bams/${SAMPLE_ID}.bam}"
VCF_OUT="${RUN_DIR}/2_vcf_raw/${SAMPLE_ID}.raw.vcf.gz"

require_file "$BAM_IN"
require_file "$REF_FASTA"
require_cmd  freebayes
require_cmd  bcftools
require_cmd  bgzip
require_cmd  tabix

STEP_START=$(date +%s)

log_info "BAM input         : $BAM_IN"
log_info "Reference         : $REF_FASTA"
log_info "Min allele balance: $MIN_AB"

# ── Call variants ─────────────────────────────────────────────────────────────
step 1 "Running FreeBayes"

# freebayes writes an uncompressed VCF to stdout.
# We pipe it through:
#   bcftools sort  — sort by position (FreeBayes may emit records out of order)
#   bgzip          — block-gzip compression (indexed by tabix)
#
# --min-base-quality 20: only count bases with Phred ≥ 20 (1% error rate).
#   This filters low-quality base calls that would inflate false-positive rates.
# --min-alternate-fraction $MIN_AB: the minimum fraction of reads supporting
#   the ALT allele. Filtering at this stage (caller level) is faster than
#   filtering post-hoc; we apply soft filters at the next step anyway.
# --ploidy 2: tell FreeBayes to assume diploid genotypes.

freebayes \
    -f "$REF_FASTA" \
    -b "$BAM_IN" \
    --ploidy 2 \
    --min-base-quality     20 \
    --min-alternate-count  2 \
    --min-alternate-fraction "$MIN_AB" \
    | bcftools sort --max-mem 1G \
    | bgzip -@ "$THREADS" > "$VCF_OUT"

# tabix builds a .tbi index for fast random access by genomic region.
# Without it, tools like bcftools annotate and VEP cannot work with this file.
tabix -p vcf "$VCF_OUT"

TOTAL_VARIANTS=$(count_vcf "$VCF_OUT")
log_success "Variant calling done in $(elapsed $STEP_START)"
log_info    "Total raw variants: $TOTAL_VARIANTS"
log_info    "Output VCF: $VCF_OUT"
