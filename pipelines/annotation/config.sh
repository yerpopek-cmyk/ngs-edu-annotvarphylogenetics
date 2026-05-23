#!/usr/bin/env bash
# =============================================================================
# pipelines/annotation/config.sh — User-editable configuration
#
# ╔══════════════════════════════════════════════════════════════════════╗
# ║  EDIT THIS FILE — do NOT edit run_annotation.sh directly.          ║
# ║  Every path and parameter the pipeline needs lives here.           ║
# ╚══════════════════════════════════════════════════════════════════════╝
# =============================================================================

# --- Input ---
# ASSEMBLY: path to your genome assembly in FASTA format.
# For the test run this is downloaded automatically (B. subtilis 168, ~4 MB).
# For your own genome: copy it to data/input/assembly.fasta and set OFFLINE=true.
ASSEMBLY="${SCRIPT_DIR}/data/input/assembly.fasta"

# --- Output ---
# A timestamped subdirectory is created automatically inside OUTDIR so that
# successive runs never overwrite each other.
OUTDIR="${SCRIPT_DIR}/outputs"

# --- Reference databases ---
# These are downloaded once on the first run (OFFLINE=false).
# On subsequent runs set OFFLINE=true to skip the download.
DB_DIR="${SCRIPT_DIR}/data/db"
PROTEIN_DB_FASTA="${DB_DIR}/swissprot_subset.fasta"   # DIAMOND query target
DIAMOND_DB="${DB_DIR}/swissprot_subset"               # Path without .dmnd extension
HMM_DB="${DB_DIR}/pfam_subset.hmm"                   # Pfam domain profiles

# --- Prokka settings ---
# PREFIX: prefix for all Prokka output files (e.g., "BSUB" → BSUB.faa, BSUB.gff …)
# KINGDOM: Bacteria | Archaea | Viruses | Mitochondria
# GENUS / SPECIES: optional but improve annotation quality when set
PROKKA_PREFIX="GENOME"
PROKKA_KINGDOM="Bacteria"
PROKKA_GENUS=""
PROKKA_SPECIES=""

# --- Performance ---
# CPUS: threads used by Prokka, DIAMOND, and HMMER.
# Defaults to (total CPUs − 2), leaving headroom for the OS.
TOTAL_CPUS=$(nproc 2>/dev/null || echo 4)
CPUS=$(( TOTAL_CPUS > 2 ? TOTAL_CPUS - 2 : TOTAL_CPUS ))

# --- Download behavior ---
# Set to true if you already have the data and want to run fully offline.
OFFLINE=false

# --- DIAMOND thresholds ---
# E-value: maximum allowed expected number of random hits.
# 1e-5 is the standard threshold for confident homology.
DIAMOND_EVALUE="1e-5"
DIAMOND_MAX_HITS=3           # How many database hits to keep per query protein

# --- HMMER thresholds ---
# Sequence-level E-value for the hmmscan search.
HMMER_SEQ_EVALUE="1e-5"
# Domain-level E-value (usually less stringent than sequence-level).
HMMER_DOM_EVALUE="1e-3"

# --- Subset for testing ---
# Set MAX_PROTEINS > 0 to run hmmscan on only the first N proteins.
# Useful for a quick sanity check without waiting for a full run.
MAX_PROTEINS=0
