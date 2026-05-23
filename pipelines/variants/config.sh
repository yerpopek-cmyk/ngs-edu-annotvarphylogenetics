#!/usr/bin/env bash
# =============================================================================
# pipelines/variants/config.sh — User-editable configuration
#
# ╔══════════════════════════════════════════════════════════════════════╗
# ║  EDIT THIS FILE — do NOT edit any step script directly.            ║
# ║  Every path and parameter the pipeline needs lives here.           ║
# ╚══════════════════════════════════════════════════════════════════════╝
#
# Copy this file from config.sh.example:
#   cp pipelines/variants/config.sh.example pipelines/variants/config.sh
# =============================================================================

# --- Input: raw reads (FASTQ) or pre-aligned BAM ---
# If BAM_INPUT is set, steps 01_align.sh is skipped entirely.
# BAM_INPUT="/path/to/existing/sample.bam"   # uncomment to use a pre-made BAM
BAM_INPUT=""

READS_R1="/path/to/sample_R1.fastq.gz"
READS_R2="/path/to/sample_R2.fastq.gz"
SAMPLE_ID="SAMPLE"                           # used in output file names

# --- Reference genome ---
# Must be accompanied by a BWA index (.bwt, .amb, .ann, .pac, .sa)
# and a FASTA index (.fai).  Both are built automatically on first run.
REF_FASTA="/path/to/reference.fasta"

# --- Output root ---
# A timestamped run directory is created inside OUTROOT on each run.
OUTROOT="${SCRIPT_DIR}/outputs"

# --- Performance ---
THREADS=$(( $(nproc 2>/dev/null || echo 4) - 2 ))
THREADS=$(( THREADS < 1 ? 1 : THREADS ))
MEM_SORT="2G"        # memory per thread for samtools sort (-m flag)

# --- FreeBayes variant calling thresholds ---
# MIN_AB: minimum allele balance (ALT reads / total reads) to call a variant.
# Lower values increase sensitivity but also increase false positives.
# Typical values: 0.2 for diploid, 0.05–0.1 for somatic/low-frequency variants.
MIN_AB=0.2

# --- Filtering thresholds ---
# Variants NOT meeting these criteria are soft-filtered (marked, not removed).
# Soft filtering preserves all variants in the VCF while flagging suspect ones
# — important for clinical pipelines that need an audit trail.
MIN_QUAL=30       # Phred-scaled variant quality; QUAL=30 → P(error) = 0.001
MIN_DP=10         # Minimum total read depth at the variant site
MAX_DP=1000       # Maximum depth; extremely high depth often signals repeats
MAX_AB=0.8        # Maximum allele balance (filters out apparent homozygous ALT calls
                  # that sit suspiciously close to 0.5, suggesting strand bias)

# --- VEP annotation ---
# VEP_CACHE_DIR: path to a local Ensembl VEP cache (~15 GB for GRCh38).
# Download once with:
#   vep_install -a cf -s homo_sapiens -y GRCh38 \
#               --CACHEDIR data/db/vep_cache --NO_HTSLIB
VEP_CACHE_DIR="${SCRIPT_DIR}/data/db/vep_cache"
VEP_ASSEMBLY="GRCh38"

# --- Report settings ---
REPORT_TOP_N=20   # Number of top-priority variants shown in the Markdown report
