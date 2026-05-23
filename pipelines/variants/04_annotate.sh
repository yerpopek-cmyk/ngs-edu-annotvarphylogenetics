#!/usr/bin/env bash
# =============================================================================
# 04_annotate.sh — Functional consequence annotation with Ensembl VEP
#
# WHAT DOES VEP DO?
# ─────────────────
# The Variant Effect Predictor (McLaren et al. 2016, Genome Biology) takes
# each variant and asks: what is the molecular consequence of this change?
#
# It annotates every variant with:
#   Consequence    — SO term: missense_variant, stop_gained, splice_donor_variant …
#   IMPACT         — HIGH / MODERATE / LOW / MODIFIER
#   SYMBOL         — gene name (e.g., BRCA2)
#   HGVSc          — coding sequence notation (e.g., c.2T>A)
#   HGVSp          — protein notation (e.g., p.Met1Lys)
#   gnomADe/g_AF   — allele frequency in gnomAD exomes / genomes
#   CLIN_SIG       — ClinVar clinical significance
#
# IMPACT HIERARCHY (most to least severe):
#   HIGH     — frameshift, stop_gained, splice_donor/acceptor, start_lost
#   MODERATE — missense, in-frame indel, protein_altering
#   LOW      — synonymous, splice_region
#   MODIFIER — UTR, intron, intergenic
#
# VEP CACHE
# ─────────
# VEP can run in "online" mode (API calls) or "offline" mode (local cache).
# Offline mode is MUCH faster and reproducible. The cache for GRCh38 is ~15 GB
# and only needs to be downloaded once.
#
# If the cache directory does not exist, this script will download it
# automatically using `vep_install`.
#
# OUTPUT
#   ${RUN_DIR}/3_vcf_filtered/${SAMPLE_ID}.annotated.vcf.gz
#   Annotations are added to the INFO/CSQ field (one block per transcript).
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
source "${REPO_ROOT}/shared/logging.sh"

banner "Step 04 — Functional Annotation with VEP"

VCF_IN="${RUN_DIR}/3_vcf_filtered/${SAMPLE_ID}.filtered.vcf.gz"
VCF_OUT="${RUN_DIR}/3_vcf_filtered/${SAMPLE_ID}.annotated.vcf.gz"

require_file "$VCF_IN"
require_cmd  vep
require_cmd  tabix
require_cmd  bcftools

STEP_START=$(date +%s)

# ── VEP cache setup ───────────────────────────────────────────────────────────
if [[ ! -d "$VEP_CACHE_DIR" ]]; then
    log_warn "VEP cache not found at: $VEP_CACHE_DIR"
    log_info "Downloading cache (homo_sapiens, ${VEP_ASSEMBLY}) — this may take 20–60 min..."
    mkdir -p "$VEP_CACHE_DIR"
    vep_install \
        -a cf \
        -s homo_sapiens \
        -y "$VEP_ASSEMBLY" \
        --CACHEDIR "$VEP_CACHE_DIR" \
        --NO_HTSLIB
    log_success "VEP cache installed"
fi

log_info "Cache dir  : $VEP_CACHE_DIR"
log_info "Assembly   : $VEP_ASSEMBLY"
log_info "Threads    : $THREADS"

# ── Run VEP ──────────────────────────────────────────────────────────────────
step 1 "Running Ensembl VEP"

# Key flags explained:
#   --offline / --cache   : use the local cache; no internet needed
#   --canonical           : mark the canonical transcript for each gene
#                           (usually the longest protein-coding transcript)
#   --hgvs                : compute HGVS notation for coding and protein changes
#   --check_existing      : look up existing dbSNP / ClinVar entries
#   --af_gnomade          : add gnomAD exome allele frequencies
#   --af_gnomadg          : add gnomAD genome allele frequencies
#   --no_stats            : skip the HTML stats page (faster, less disk)
#   --force_overwrite     : overwrite output if it already exists
#   --fields              : restrict the CSQ tag to only the fields we need
#                           (keeps the VCF manageable)
#   --synonyms            : chromosome name synonym file so VEP can match
#                           chr1 (UCSC) to 1 (Ensembl) and vice-versa

vep \
    --input_file     "$VCF_IN" \
    --output_file    "$VCF_OUT" \
    --format         vcf \
    --vcf \
    --compress_output bgzip \
    --offline \
    --cache \
    --dir_cache      "$VEP_CACHE_DIR" \
    --assembly       "$VEP_ASSEMBLY" \
    --fork           "$THREADS" \
    --canonical \
    --hgvs \
    --check_existing \
    --af_gnomade \
    --af_gnomadg \
    --no_stats \
    --quiet \
    --force_overwrite \
    --fields "Consequence,IMPACT,SYMBOL,Gene,Feature,HGVSc,HGVSp,gnomADe_AF,gnomADg_AF,CLIN_SIG"

tabix -p vcf "$VCF_OUT"

# ── Quick consequence summary ─────────────────────────────────────────────────
step 2 "Consequence summary"

# Extract the Consequence field from CSQ and count occurrences
bcftools query -f '%INFO/CSQ\n' "$VCF_OUT" \
    | tr ',' '\n' \
    | cut -d'|' -f1 \
    | sort | uniq -c | sort -rn \
    | head -15 \
    | awk '{printf "  %-10s %s\n", $1, $2}'

log_success "Annotation done in $(elapsed $STEP_START)"
log_info    "Output VCF: $VCF_OUT"
