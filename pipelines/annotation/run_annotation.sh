#!/usr/bin/env bash
# =============================================================================
# pipelines/annotation/run_annotation.sh
#
# WHAT THIS PIPELINE DOES
# ═══════════════════════
# Bacterial genome annotation in three stages:
#
#   1. Prokka — Structural annotation
#      Finds ORFs (Open Reading Frames), rRNA, tRNA, and assigns initial
#      functional labels based on curated databases. Output: .gff, .faa, .ffn
#
#   2. DIAMOND blastp — Homology-based functional annotation
#      Aligns predicted proteins against Swiss-Prot using DIAMOND (a fast
#      reimplementation of BLASTX/BLASTP). Proteins with hits receive
#      a functional description from the database.
#      Key concept: E-value = K × m × n × e^(-λ × S)
#      where m=query length, n=DB size, S=alignment score.
#      We use E-value ≤ 1e-5 as the significance threshold.
#
#   3. hmmscan — Domain-based annotation (Pfam)
#      Searches for conserved protein domains using Hidden Markov Models.
#      HMMs capture the statistical signature of a protein family across
#      many aligned sequences — much more sensitive than BLAST alone.
#
# USAGE
# ═════
#   conda activate annotation
#   bash run_annotation.sh [OPTIONS]
#
# OPTIONS
#   --offline            Skip downloading; use cached data in data/db/ and data/input/
#   --max-proteins N     Run hmmscan on only the first N proteins (testing)
#   --cpus N             Override thread count from config.sh
#   --config FILE        Load an alternative config file
#   -h, --help           Show this help
#
# QUICK TEST (downloads B. subtilis genome automatically):
#   bash run_annotation.sh
# =============================================================================

set -euo pipefail

# --- Resolve script location so config.sh can reference $SCRIPT_DIR reliably ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# --- Load shared logging utilities ---
source "${REPO_ROOT}/shared/logging.sh"

# --- Load default configuration, then allow CLI overrides ---
CONFIG_FILE="${SCRIPT_DIR}/config.sh"

# ── CLI Argument Parsing ────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --offline)            OFFLINE_OVERRIDE=true;       shift ;;
        --max-proteins)       MAX_PROTEINS_OVERRIDE="$2";  shift 2 ;;
        --cpus)               CPUS_OVERRIDE="$2";          shift 2 ;;
        --config)             CONFIG_FILE="$2";            shift 2 ;;
        -h|--help)
            sed -n '/^# USAGE/,/^# =============================================================================/p' "$0" \
                | sed '/^# =============================================================================/d; s/^# \?//'
            exit 0 ;;
        *) die "Unknown argument: $1  (use --help for usage)" ;;
    esac
done

# Source the config file (must happen after --config is parsed)
[[ -f "$CONFIG_FILE" ]] || die "Config file not found: $CONFIG_FILE"
source "$CONFIG_FILE"

# Apply CLI overrides (take precedence over config.sh values)
[[ -n "${OFFLINE_OVERRIDE:-}"      ]] && OFFLINE=true
[[ -n "${MAX_PROTEINS_OVERRIDE:-}" ]] && MAX_PROTEINS="$MAX_PROTEINS_OVERRIDE"
[[ -n "${CPUS_OVERRIDE:-}"         ]] && CPUS="$CPUS_OVERRIDE"

# ── Output directory ─────────────────────────────────────────────────────────
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
RUN_DIR="${OUTDIR}/run_${TIMESTAMP}"
mkdir -p \
    "${SCRIPT_DIR}/data/input" \
    "${SCRIPT_DIR}/data/db" \
    "${RUN_DIR}"/{prokka,diamond,hmmer}

# ── Redirect all output to a log file ───────────────────────────────────────
LOG_FILE="${RUN_DIR}/pipeline.log"
exec > >(tee -i "$LOG_FILE") 2>&1

# ── Print pipeline header ────────────────────────────────────────────────────
banner "Genome Annotation Pipeline"
log_info "Config file   : $CONFIG_FILE"
log_info "Output dir    : $RUN_DIR"
log_info "Threads       : $CPUS"
log_info "Offline mode  : $OFFLINE"
log_info "Log file      : $LOG_FILE"

PIPELINE_START=$(date +%s)

# ============================================================================
# STAGE 0 — Dependency check
# ============================================================================
step 0 "Checking required tools"

# We check for every tool upfront so the user gets a clear, complete list of
# what is missing rather than discovering each tool absence one at a time.
MISSING_TOOLS=()
for tool in prokka diamond hmmscan hmmpress seqkit wget bgzip tabix; do
    if ! command -v "$tool" &>/dev/null; then
        MISSING_TOOLS+=("$tool")
    else
        log_success "$tool  $(${tool} --version 2>&1 | head -1 || true)"
    fi
done

if [[ ${#MISSING_TOOLS[@]} -gt 0 ]]; then
    log_error "The following tools are missing: ${MISSING_TOOLS[*]}"
    log_error "Run:  conda env create -f environments/annotation.yml"
    log_error "Then: conda activate annotation"
    exit 1
fi

# ============================================================================
# STAGE 1 — Download reference data (skipped in offline mode)
# ============================================================================
step 1 "Downloading reference data"

if [[ "$OFFLINE" == false ]]; then
    # ── B. subtilis 168 reference genome ────────────────────────────────────
    # This is a well-annotated model organism, ideal for testing annotation
    # pipelines (genome size ~4.2 Mb, ~4100 genes).
    if [[ ! -f "$ASSEMBLY" ]]; then
        log_info "Downloading B. subtilis 168 genome (~4 MB)..."
        wget -qc --show-progress \
            -O "${ASSEMBLY}.gz" \
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/045/GCF_000009045.1_ASM904v1/GCF_000009045.1_ASM904v1_genomic.fna.gz"
        gunzip -f "${ASSEMBLY}.gz"
        log_success "Genome saved: $ASSEMBLY"
    else
        log_info "Genome: using cached file $ASSEMBLY"
    fi

    # ── Swiss-Prot protein database (first 100 k records) ───────────────────
    # Swiss-Prot is the manually reviewed section of UniProt — very high
    # confidence annotations. The full file is ~600 MB compressed; we take
    # the first 100 k entries as a manageable, representative subset.
    if [[ ! -f "$PROTEIN_DB_FASTA" ]]; then
        log_info "Downloading Swiss-Prot subset (~40 MB)..."
        wget -qc --show-progress \
            -O "${DB_DIR}/uniprot_sprot.fasta.gz" \
            "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
        # awk trick: count '>' header lines, stop after 100 k sequences
        zcat "${DB_DIR}/uniprot_sprot.fasta.gz" \
            | awk '/^>/{n++} n>100000{exit} {print}' \
            > "$PROTEIN_DB_FASTA"
        rm -f "${DB_DIR}/uniprot_sprot.fasta.gz"
        log_success "Swiss-Prot subset: $(count_fasta "$PROTEIN_DB_FASTA") sequences"
    else
        log_info "Swiss-Prot: using cached file"
    fi

    # ── Pfam domain profiles (10-domain subset for testing) ─────────────────
    # Pfam stores protein domain families as HMM profiles. Each profile is a
    # statistical model of the conservation pattern across many aligned sequences.
    # We extract 10 commonly occurring domain families for the test run:
    #   PF00005 — ABC transporter ATPase domain
    #   PF00009 — Elongation factor Tu, GTP-binding domain
    #   PF00012 — Hsp70 chaperone domain
    #   PF00013 — RRM (RNA Recognition Motif)
    #   PF00023 — Ankyrin repeat
    #   PF00027 — Cyclic nucleotide binding domain
    #   PF00028 — Cadherin domain
    #   PF00043 — Glutathione S-transferase C-terminal domain
    #   PF00044 — Glyceraldehyde 3-phosphate dehydrogenase NAD binding
    #   PF00072 — Response regulator receiver domain
    if [[ ! -f "${HMM_DB}.h3i" ]]; then
        log_info "Downloading Pfam-A HMM database (~500 MB compressed)..."
        wget -qc --show-progress \
            -O "${DB_DIR}/Pfam-A.hmm.gz" \
            "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
        gunzip -kf "${DB_DIR}/Pfam-A.hmm.gz"

        log_info "Extracting 10-domain subset from Pfam-A..."
        DOMAINS="PF00005 PF00009 PF00012 PF00013 PF00023 PF00027 PF00028 PF00043 PF00044 PF00072"
        python3 - "$DOMAINS" "${DB_DIR}/Pfam-A.hmm" "$HMM_DB" <<'PY'
import sys

domains   = set(sys.argv[1].split())
input_hmm = sys.argv[2]
output    = sys.argv[3]

# Each HMM block starts with "NAME" and ends with "//".
# We stream through and only write blocks whose ACC field is in our set.
capture, buffer = False, []
with open(input_hmm, encoding="utf-8") as fin, \
     open(output,    "w", encoding="utf-8") as fout:
    for line in fin:
        buffer.append(line)
        if line.startswith("ACC "):
            acc = line.split()[1].split(".")[0]   # strip version number (e.g. .3)
            capture = acc in domains
        elif line.startswith("//"):
            if capture:
                fout.writelines(buffer)
            capture = False
            buffer  = []
PY
        # hmmpress builds binary index files (.h3f, .h3i, .h3m, .h3p) that
        # allow hmmscan to do random access into the profile database.
        hmmpress -f "$HMM_DB"
        rm -f "${DB_DIR}/Pfam-A.hmm" "${DB_DIR}/Pfam-A.hmm.gz"
        log_success "Pfam subset pressed and ready: $HMM_DB"
    else
        log_info "Pfam: using cached profiles"
    fi
else
    log_info "Offline mode — skipping all downloads"
    require_file "$ASSEMBLY"
    require_file "$PROTEIN_DB_FASTA"
    require_file "${HMM_DB}.h3i"
fi

# ============================================================================
# STAGE 2 — Build DIAMOND database
# ============================================================================
step 2 "Building DIAMOND protein database"

# DIAMOND uses a binary database format (.dmnd) for fast alignment.
# diamond makedb converts a FASTA protein file into this format once;
# subsequent searches reuse the binary without rebuilding.
#
# WHY DIAMOND instead of BLAST?
#   DIAMOND is 500–20,000× faster than BLASTP at comparable sensitivity.
#   It achieves this by: (1) double-indexing both query and database,
#   (2) a spaced seed strategy, and (3) SIMD vectorization.

if [[ ! -f "${DIAMOND_DB}.dmnd" ]]; then
    log_info "Building DIAMOND database (this may take a few minutes)..."
    diamond makedb \
        --in  "$PROTEIN_DB_FASTA" \
        --db  "$DIAMOND_DB" \
        --threads "$CPUS" \
        --quiet
    log_success "DIAMOND database ready: ${DIAMOND_DB}.dmnd"
else
    log_info "DIAMOND database: using cached ${DIAMOND_DB}.dmnd"
fi

# ============================================================================
# STAGE 3 — Prokka structural annotation
# ============================================================================
step 3 "Structural annotation with Prokka"

# Prokka (Seemann 2014) is a prokaryote annotation pipeline that integrates:
#   1. Prodigal  — Gene prediction (finds ORFs using hexamer statistics)
#   2. RNAmmer   — rRNA prediction
#   3. Aragorn    — tRNA prediction
#   4. BLAST/DIAMOND — Functional annotation against curated databases
#   5. Infernal   — ncRNA prediction via Rfam profiles
#
# --force    : overwrite the output directory if it already exists
# --quiet    : suppress verbose tool output (logged separately)
# --kingdom  : affects which databases are searched (Bacteria uses Uniprot-sprot
#              + AMRfinder + Resfams)

PROKKA_ARGS=(
    --outdir "${RUN_DIR}/prokka"
    --prefix "$PROKKA_PREFIX"
    --kingdom "$PROKKA_KINGDOM"
    --cpus    "$CPUS"
    --force
    --quiet
)
[[ -n "$PROKKA_GENUS"   ]] && PROKKA_ARGS+=(--genus   "$PROKKA_GENUS")
[[ -n "$PROKKA_SPECIES" ]] && PROKKA_ARGS+=(--species "$PROKKA_SPECIES")

prokka "${PROKKA_ARGS[@]}" "$ASSEMBLY"

FAA="${RUN_DIR}/prokka/${PROKKA_PREFIX}.faa"
require_file "$FAA"
log_success "Prokka done: $(count_fasta "$FAA") proteins predicted"

# ============================================================================
# STAGE 4 — DIAMOND blastp: homology-based functional annotation
# ============================================================================
step 4 "Homology annotation with DIAMOND blastp"

# outfmt 6 fields (tab-separated):
#   qseqid  — query sequence ID (our predicted protein)
#   sseqid  — subject (database) sequence ID
#   pident  — percentage of identical positions in the alignment
#   length  — alignment length
#   evalue  — E-value: expected number of random hits with this score or better
#   bitscore — normalized alignment score (higher = better)
#   stitle  — full subject title (functional description from Swiss-Prot)
#
# --max-target-seqs 3  : report up to 3 hits per query protein; more hits
#                        allow downstream tools to pick the best description.
# --very-sensitive     : increases sensitivity at the cost of speed;
#                        recommended when DB coverage is lower (our small subset).

diamond blastp \
    --db             "$DIAMOND_DB" \
    --query          "$FAA" \
    --out            "${RUN_DIR}/diamond/blastp.tsv" \
    --outfmt 6 qseqid sseqid pident length evalue bitscore stitle \
    --evalue         "$DIAMOND_EVALUE" \
    --max-target-seqs "$DIAMOND_MAX_HITS" \
    --threads        "$CPUS" \
    --block-size     1.0 \
    --very-sensitive \
    --quiet

DIAMOND_HITS=$(cut -f1 "${RUN_DIR}/diamond/blastp.tsv" | sort -u | wc -l | tr -d ' ')
log_success "DIAMOND: ${DIAMOND_HITS}/$(count_fasta "$FAA") proteins have homology hits"

# ============================================================================
# STAGE 5 — hmmscan: domain-based annotation (Pfam)
# ============================================================================
step 5 "Domain annotation with HMMER hmmscan"

# hmmscan searches each PROTEIN against all profiles in the HMM database.
# (contrast: hmmsearch searches each PROFILE against a sequence database)
#
# --domtblout  : per-domain tabular output — one line per domain hit;
#                easier to parse than the human-readable main output.
# --noali      : omit the alignment display (saves space; we only need the table)
# -E           : sequence-level E-value threshold
# --domE       : domain-level E-value threshold (more permissive, catches
#                multi-domain proteins where some domains score just below the
#                sequence threshold alone)
# -o           : human-readable output (saved to log for debugging)

QUERY_FAA="$FAA"
if [[ "$MAX_PROTEINS" -gt 0 ]]; then
    log_info "Limiting hmmscan to first $MAX_PROTEINS proteins (testing mode)"
    seqkit head -n "$MAX_PROTEINS" "$FAA" > "${RUN_DIR}/hmmer/query.faa"
    QUERY_FAA="${RUN_DIR}/hmmer/query.faa"
fi

hmmscan \
    --cpu        "$CPUS" \
    --domtblout  "${RUN_DIR}/hmmer/domtblout.txt" \
    --noali \
    -E           "$HMMER_SEQ_EVALUE" \
    --domE       "$HMMER_DOM_EVALUE" \
    -o           "${RUN_DIR}/hmmer/hmmscan.log" \
    "$HMM_DB" \
    "$QUERY_FAA"

DOMAIN_HITS=$(grep -cv '^#' "${RUN_DIR}/hmmer/domtblout.txt" || echo 0)
log_success "HMMER: ${DOMAIN_HITS} domain hits found"

# ============================================================================
# STAGE 6 — Summary report
# ============================================================================
step 6 "Pipeline summary"

TOTAL_PROTEINS=$(count_fasta "$FAA")

echo ""
echo -e "${BOLD}┌─────────────────────────────────────────────────────┐${RESET}"
echo -e "${BOLD}│              Annotation Summary                     │${RESET}"
echo -e "${BOLD}├─────────────────────────────────────────────────────┤${RESET}"

if [[ -f "${RUN_DIR}/prokka/${PROKKA_PREFIX}.txt" ]]; then
    # The .txt file contains one statistic per line: "field: value"
    grep -E "CDS|tRNA|rRNA|tmRNA|CRISPR" \
        "${RUN_DIR}/prokka/${PROKKA_PREFIX}.txt" \
        | awk '{printf "│  %-50s│\n", $0}'
fi

printf "│  %-50s│\n" ""
printf "│  %-50s│\n" "DIAMOND  hits  : ${DIAMOND_HITS} / ${TOTAL_PROTEINS}"
printf "│  %-50s│\n" "HMMER   domains: ${DOMAIN_HITS}"
printf "│  %-50s│\n" ""
printf "│  %-50s│\n" "Log file: ${LOG_FILE}"
printf "│  %-50s│\n" "Output dir: ${RUN_DIR}"
echo -e "${BOLD}└─────────────────────────────────────────────────────┘${RESET}"

echo ""
log_success "Annotation pipeline finished in $(elapsed $PIPELINE_START)"
echo ""
echo "  Key output files:"
echo "  ├── ${RUN_DIR}/prokka/${PROKKA_PREFIX}.gff  — gene coordinates (GFF3)"
echo "  ├── ${RUN_DIR}/prokka/${PROKKA_PREFIX}.faa  — protein sequences (FASTA)"
echo "  ├── ${RUN_DIR}/diamond/blastp.tsv           — homology hits (TSV)"
echo "  └── ${RUN_DIR}/hmmer/domtblout.txt          — domain annotations (TSV)"
