#!/usr/bin/env bash
# =============================================================================
# 01_align.sh — Align paired-end reads to the reference genome
#
# TOOLS USED
# ──────────
#   bwa mem     — Burrows-Wheeler Aligner, MEM algorithm.
#                 MEM (Maximal Exact Matches) finds the longest exact match
#                 between a read and the reference, then extends it.
#                 Better than classic BWA-backtrack for reads > 70 bp.
#
#   samtools fixmate
#               — Fills in mate-pair information (RNEXT, PNEXT, TLEN fields)
#                 and marks "orphaned" read pairs. Required before markdup.
#
#   samtools sort
#               — Sorts reads by coordinate. Most downstream tools (BAM indexing,
#                 variant callers, IGV) expect coordinate-sorted BAM.
#
#   samtools markdup
#               — Marks reads that are likely PCR or optical duplicates.
#                 Duplicates arise because the same original DNA molecule was
#                 amplified and sequenced multiple times. They inflate apparent
#                 coverage and create false allele counts. markdup tags them
#                 in the FLAG field (bit 1024); most callers then ignore them.
#
# OUTPUT
#   ${RUN_DIR}/1_bams/${SAMPLE_ID}.bam  — coordinate-sorted, markdup BAM
#   ${RUN_DIR}/1_bams/${SAMPLE_ID}.bam.bai — BAM index (auto-created by --write-index)
#   ${RUN_DIR}/1_bams/${SAMPLE_ID}.flagstat — alignment statistics
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
source "${REPO_ROOT}/shared/logging.sh"

banner "Step 01 — Alignment with BWA"

# If the user already has a BAM, skip this step entirely
if [[ -n "${BAM_INPUT:-}" ]]; then
    log_info "BAM_INPUT is set → skipping alignment"
    log_info "Using: $BAM_INPUT"
    exit 0
fi

require_file "$READS_R1"
require_file "$READS_R2"
require_file "$REF_FASTA"

OUT_BAM="${RUN_DIR}/1_bams/${SAMPLE_ID}.bam"

STEP_START=$(date +%s)

# ── Build BWA and FASTA indices if they are missing ──────────────────────────
# These are only built once; on subsequent runs the existing index is reused.

if [[ ! -f "${REF_FASTA}.fai" ]]; then
    step 1 "Building FASTA index (samtools faidx)"
    # The .fai index allows random access into the FASTA file.
    # It is also required by some downstream tools (bcftools, GATK).
    samtools faidx "$REF_FASTA"
    log_success "FASTA index created"
fi

if [[ ! -f "${REF_FASTA}.bwt" ]]; then
    step 2 "Building BWA index (bwa index)"
    # BWA creates five index files (.amb, .ann, .bwt, .pac, .sa).
    # For large genomes (> 2 Gb) add -a bwtsw to use a more memory-efficient algorithm.
    bwa index "$REF_FASTA"
    log_success "BWA index created"
fi

# ── Main alignment pipeline ───────────────────────────────────────────────────
step 3 "Aligning reads: bwa mem | samtools fixmate | sort | markdup"

# Read Group (@RG) tag — required by many variant callers (GATK demands it).
# Fields:
#   ID   — read group identifier (must be unique per sequencing run)
#   SM   — sample name (used in VCF FORMAT column header)
#   PL   — platform (ILLUMINA, PACBIO, ONT, …)
#   LB   — library preparation identifier
RG="@RG\tID:${SAMPLE_ID}\tSM:${SAMPLE_ID}\tPL:ILLUMINA\tLB:lib1"

log_info "Read group: $RG"
log_info "Threads: $THREADS"

# This is a single streaming pipeline to avoid writing large intermediate files:
#
#   bwa mem  →  fixmate  →  sort  →  markdup  →  BAM
#
# -M flag (bwa): mark split-read secondary alignments as supplementary
#                (compatibility with Picard-based tools)
# -u flag (fixmate/sort): use uncompressed BAM in the pipe for speed
# -m flag (fixmate): add mate score (used by markdup to choose which
#                    duplicate in a cluster to keep)
# --write-index (markdup): write the .bai index without a separate samtools index step

bwa mem \
        -t  "$THREADS" \
        -R  "$RG" \
        "$REF_FASTA" "$READS_R1" "$READS_R2" \
    | samtools fixmate \
        -m \
        -u \
        - \
        - \
    | samtools sort \
        -@ "$THREADS" \
        -m "$MEM_SORT" \
        -u \
        - \
    | samtools markdup \
        -@ "$THREADS" \
        --write-index \
        - \
        "$OUT_BAM"

# ── Alignment statistics ──────────────────────────────────────────────────────
step 4 "Generating alignment statistics (samtools flagstat)"

# flagstat counts reads in each alignment category.
# Key metrics to check:
#   - Mapped reads: ideally ≥ 95% for a human sample
#   - Properly paired: reads where both mates mapped to the same chromosome
#     in the expected orientation and distance
#   - Duplicates: typically 5–30% for standard WGS; high values suggest over-amplification

samtools flagstat "$OUT_BAM" | tee "${RUN_DIR}/1_bams/${SAMPLE_ID}.flagstat"

require_file "$OUT_BAM"
log_success "Alignment done in $(elapsed $STEP_START)"
log_info    "Output BAM: $OUT_BAM"
