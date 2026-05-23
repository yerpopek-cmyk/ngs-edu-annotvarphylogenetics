#!/usr/bin/env bash
# =============================================================================
# shared/logging.sh — Logging and utility functions used by all pipelines
#
# HOW TO USE IN YOUR SCRIPT:
#   SHARED_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
#   source "${SHARED_DIR}/../../shared/logging.sh"
#
# WHAT THIS FILE PROVIDES:
#   banner   — Print a large section header (start of a new pipeline stage)
#   step     — Print a smaller step header within a stage
#   log_info / log_warn / log_error / log_success — Timestamped messages
#   die      — Print an error message and exit with code 1
#   require_cmd / require_file / require_dir — Guard clauses
#   elapsed  — Human-readable elapsed time
# =============================================================================

# --- ANSI Color Codes ---
# These make console output much easier to scan visually.
# \033[...m sets a "Select Graphic Rendition" mode. 0 = reset.
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
BOLD='\033[1m'
RESET='\033[0m'

# Return the current timestamp in a log-friendly format
timestamp() { date '+%Y-%m-%d %H:%M:%S'; }

# =============================================================================
# banner — Print a prominent section separator.
#
# Use at the START of each major logical block (e.g., before alignment,
# before variant calling). Makes pipeline logs much easier to grep and read.
#
# Usage: banner "Alignment with BWA-MEM2"
# =============================================================================
banner() {
    local title="${1:-}"
    local width=68
    # Build a line of '═' characters of the right length
    local bar
    printf -v bar '%*s' "$width" ''
    bar="${bar// /═}"
    echo -e "\n${BLUE}${BOLD}╔${bar}╗${RESET}"
    printf "${BLUE}${BOLD}║  %-${width}s║${RESET}\n" "  ${title}"
    printf "${BLUE}${BOLD}║  %-${width}s║${RESET}\n" "  $(timestamp)"
    echo -e "${BLUE}${BOLD}╚${bar}╝${RESET}\n"
}

# =============================================================================
# step — Print a numbered step header within a section.
#
# Usage: step 3 "Normalize and soft-filter variants"
# =============================================================================
step() {
    local num="$1"
    local title="$2"
    echo -e "\n${CYAN}${BOLD}▶ Step ${num}: ${title}${RESET}"
    echo -e "${CYAN}  $(timestamp)${RESET}"
    echo -e "${CYAN}  ─────────────────────────────────────────────────${RESET}"
}

# =============================================================================
# Timestamped log functions
# =============================================================================
log_info()    { echo -e "$(timestamp)  ${GREEN}INFO${RESET}    $*"; }
log_warn()    { echo -e "$(timestamp)  ${YELLOW}WARN${RESET}    $*" >&2; }
log_error()   { echo -e "$(timestamp)  ${RED}ERROR${RESET}   $*" >&2; }
log_success() { echo -e "$(timestamp)  ${GREEN}✓ OK${RESET}    $*"; }

# =============================================================================
# die — Print an error message and terminate immediately.
#
# Called in guard clauses and anywhere a fatal condition is detected.
# Combined with `set -euo pipefail` this ensures no silent failures.
# =============================================================================
die() { log_error "$*"; exit 1; }

# =============================================================================
# require_cmd — Assert that a command is present on PATH.
#
# If it is missing, we print a helpful message pointing to the environment
# file rather than a cryptic "command not found" later in the script.
#
# Usage: require_cmd samtools
# =============================================================================
require_cmd() {
    local cmd="$1"
    command -v "$cmd" >/dev/null 2>&1 \
        || die "'${cmd}' not found in PATH. Activate the correct conda environment."
}

# =============================================================================
# require_file — Assert that a file exists and is non-empty.
# require_dir  — Assert that a directory exists.
# =============================================================================
require_file() { [[ -f "$1" && -s "$1" ]] || die "Required file missing or empty: $1"; }
require_dir()  { [[ -d "$1" ]]            || die "Required directory missing: $1"; }

# =============================================================================
# elapsed — Pretty-print elapsed wall-clock time.
#
# Usage:
#   START=$(date +%s)
#   ... do work ...
#   log_info "Finished in $(elapsed $START)"
# =============================================================================
elapsed() {
    local start="$1"
    local secs=$(( $(date +%s) - start ))
    printf '%dh %02dm %02ds' $(( secs/3600 )) $(( (secs%3600)/60 )) $(( secs%60 ))
}

# =============================================================================
# count_fasta — Count sequences in a FASTA file (by counting '>' lines).
# count_vcf   — Count non-header variant records in a VCF/BCF.
# =============================================================================
count_fasta() { grep -c '^>' "${1:?}" 2>/dev/null || echo 0; }
count_vcf()   { bcftools view -H "${1:?}" 2>/dev/null | wc -l | tr -d ' ' || echo 0; }
