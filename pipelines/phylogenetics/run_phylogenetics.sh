#!/usr/bin/env bash
# =============================================================================
# pipelines/phylogenetics/run_phylogenetics.sh
#
# WHAT THIS PIPELINE DOES
# ═══════════════════════
# Builds a maximum-likelihood phylogenetic tree from unaligned sequences.
# Seven steps:
#
#   1. MAFFT     — Multiple sequence alignment (MSA)
#   2. ClipKit   — Trim uninformative alignment columns
#   3. FastTree  — Rapid approximate ML tree (for exploration)
#   4. IQ-TREE   — Full ML tree with model selection and bootstrap
#   5. Rerooting — Midpoint rooting of the best tree (gotree)
#   6. PhyloBayes prep — Convert to PHYLIP format (for Bayesian analysis)
#   7. BEAST2 prep — Generate XML alignment draft (for dated phylogenies)
#
# BACKGROUND: WHY DO WE NEED THESE STEPS?
# ────────────────────────────────────────
#
#   ALIGNMENT (MAFFT):
#     Raw sequences have different lengths — insertions and deletions shift
#     positions relative to each other. Alignment inserts gap characters (-)
#     to make homologous positions line up in the same column.
#     MAFFT uses Fast Fourier Transforms (FFT) to quickly find regions of
#     similarity and then refines the alignment iteratively.
#
#   TRIMMING (ClipKit):
#     Not every aligned column contains phylogenetic signal. Highly gapped
#     columns and invariant columns (same base in all sequences) add noise
#     without contributing information about evolutionary relationships.
#     kpic = keep parsimony-informative and constant sites.
#
#   TREE CONSTRUCTION:
#     Maximum Likelihood (ML) finds the tree topology and branch lengths
#     that maximize P(alignment | tree, model).
#
#     The substitution model describes how bases change over time.
#     For DNA, the most general model is GTR+Γ+I:
#       GTR: 6 independent substitution rates (A↔C, A↔G, A↔T, C↔G, C↔T, G↔T)
#       Γ : gamma distribution of rate variation across sites
#       I : proportion of invariant sites
#
#   BOOTSTRAP:
#     Resample alignment columns with replacement, rebuild the tree N times.
#     Bootstrap support at a node = % of resampled trees containing that clade.
#     ≥70% is commonly used as the threshold for a reliable clade.
#
# USAGE
# ═════
#   conda activate phylogenetics_env
#   bash run_phylogenetics.sh -i sequences.fasta -t DNA
#
# OPTIONS
#   -i, --input   Input FASTA file (required)
#   -t, --type    DNA (default) or AA
#   -o, --outdir  Custom output directory
#   --config FILE Alternative config.sh
#   --skip-STEP   Skip a step (mafft|clipkit|fasttree|iqtree|phylobayes|beast2)
#   -h, --help    Show this help
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
source "${REPO_ROOT}/shared/logging.sh"

# --- Load config then allow CLI overrides ---
CONFIG_FILE="${SCRIPT_DIR}/config.sh"
INPUT_FASTA_OVERRIDE=""
SEQ_TYPE_OVERRIDE=""
OUTDIR_OVERRIDE=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input)      INPUT_FASTA_OVERRIDE="$2";  shift 2 ;;
        -t|--type)       SEQ_TYPE_OVERRIDE="${2^^}"; shift 2 ;;
        -o|--outdir)     OUTDIR_OVERRIDE="$2";       shift 2 ;;
        --config)        CONFIG_FILE="$2";            shift 2 ;;
        --skip-mafft)    SKIP_MAFFT=1;               shift ;;
        --skip-clipkit)  SKIP_CLIPKIT=1;             shift ;;
        --skip-fasttree) SKIP_FASTTREE=1;            shift ;;
        --skip-iqtree)   SKIP_IQTREE=1;              shift ;;
        --skip-phylobayes) SKIP_PHYLOBAYES=1;        shift ;;
        --skip-beast2)   SKIP_BEAST2=1;              shift ;;
        -h|--help)
            sed -n '/^# USAGE/,/^# =============================================================================/p' "$0" \
                | sed '/^# =============================================================================/d; s/^# \?//'
            exit 0 ;;
        *) die "Unknown argument: $1  (use --help)" ;;
    esac
done

[[ -f "$CONFIG_FILE" ]] || die "Config file not found: $CONFIG_FILE"
source "$CONFIG_FILE"

# CLI overrides
[[ -n "$INPUT_FASTA_OVERRIDE" ]] && INPUT_FASTA="$INPUT_FASTA_OVERRIDE"
[[ -n "$SEQ_TYPE_OVERRIDE"    ]] && SEQ_TYPE="$SEQ_TYPE_OVERRIDE"
[[ -n "$OUTDIR_OVERRIDE"      ]] && OUTROOT="$OUTDIR_OVERRIDE"

[[ -f "$INPUT_FASTA" ]] || die "Input FASTA not found: ${INPUT_FASTA:-<not set>}. Use -i <file>."
[[ "$SEQ_TYPE" == "DNA" || "$SEQ_TYPE" == "AA" ]] || die "SEQ_TYPE must be DNA or AA."

# Derive a clean basename for output files (strip .fasta/.fa/.fas)
BASENAME="$(basename "$INPUT_FASTA")"
BASENAME="${BASENAME%.fasta}"; BASENAME="${BASENAME%.fa}"; BASENAME="${BASENAME%.fas}"

RUN_DIR="${OUTROOT}/${BASENAME}_results"
mkdir -p "$RUN_DIR"

LOG_FILE="${RUN_DIR}/pipeline.log"
exec > >(tee -i "$LOG_FILE") 2>&1

# ── File paths used across steps ─────────────────────────────────────────────
ALIGNED="${RUN_DIR}/${BASENAME}_aligned.fasta"
TRIMMED="${RUN_DIR}/${BASENAME}_trimmed.fasta"
TREE_FASTTREE="${RUN_DIR}/${BASENAME}_fasttree.tree"
TREE_ML_PRE="${RUN_DIR}/${BASENAME}_iqtree"
TREE_ML="${TREE_ML_PRE}.treefile"
TREE_ROOTED="${RUN_DIR}/${BASENAME}_rooted.tree"
PHYLIP_OUT="${RUN_DIR}/${BASENAME}_trimmed.phy"
PHYLIP_MAP="${RUN_DIR}/${BASENAME}_phylip_name_map.tsv"
BEAST_XML="${RUN_DIR}/${BASENAME}_beast_alignment.xml"

alignment_columns() {
    awk '
        BEGIN { len = 0; in_first = 0; printed = 0 }
        /^>/ {
            if (in_first) {
                print len
                printed = 1
                exit
            }
            in_first = 1
            next
        }
        in_first { len += length($0) }
        END {
            if (in_first && !printed) {
                print len
            }
        }
    ' "$1"
}

banner "Phylogenetic Analysis Pipeline"
log_info "Input    : $INPUT_FASTA  ($(count_fasta "$INPUT_FASTA") sequences)"
log_info "Type     : $SEQ_TYPE"
log_info "Threads  : $THREADS"
log_info "Run dir  : $RUN_DIR"

PIPELINE_START=$(date +%s)

# ── Validate input ────────────────────────────────────────────────────────────
step 0 "Validating input FASTA"

N_SEQS=$(count_fasta "$INPUT_FASTA")
[[ "$N_SEQS" -ge 3 ]] \
    || die "Need ≥ 3 sequences for tree inference; found $N_SEQS in $INPUT_FASTA"

# Check for empty sequence headers
if grep -qP '^>(\s*)$' "$INPUT_FASTA" 2>/dev/null; then
    die "FASTA contains empty sequence headers — please fix before running."
fi

log_success "Input is valid: $N_SEQS sequences"

# ============================================================================
# STEP 1 — MAFFT: multiple sequence alignment
# ============================================================================
if [[ "${SKIP_MAFFT:-0}" -eq 0 ]]; then
    step 1 "Multiple sequence alignment with MAFFT"
    require_cmd mafft

    # --auto: MAFFT automatically chooses between L-INS-i (accurate, slow,
    #         for < 200 sequences) and FFT-NS-2 (fast, for large datasets).
    # --thread: use multiple CPU cores for the alignment computation.
    #
    # The output is a FASTA file with '-' gap characters inserted so that
    # all sequences have the same length and homologous positions are aligned.

    mafft \
        ${MAFFT_EXTRA_ARGS} \
        --thread "$THREADS" \
        "$INPUT_FASTA" \
        > "$ALIGNED"

    log_success "Alignment: $(count_fasta "$ALIGNED") sequences aligned"
    log_info    "Alignment length: $(alignment_columns "$ALIGNED") columns"
else
    log_warn "Skipping MAFFT (--skip-mafft)"
    [[ -f "$ALIGNED" ]] || die "Aligned file expected at $ALIGNED but not found."
fi

# ============================================================================
# STEP 2 — ClipKit: trim uninformative alignment columns
# ============================================================================
if [[ "${SKIP_CLIPKIT:-0}" -eq 0 ]]; then
    step 2 "Trimming with ClipKit (method: ${CLIPKIT_METHOD})"
    require_cmd clipkit

    # Parsimony-informative sites: sites where at least 2 different bases
    # each appear in at least 2 sequences. Only these sites can resolve
    # tree topology; invariant and gappy sites just add noise.
    #
    # Constant sites (same base in all taxa) are kept by kpic because they
    # are needed for some model parameters (proportion of invariant sites).

    ALIGNED_LEN=$(alignment_columns "$ALIGNED")
    clipkit "$ALIGNED" -m "$CLIPKIT_METHOD" -o "$TRIMMED"
    TRIMMED_LEN=$(alignment_columns "$TRIMMED")

    log_success "Trimmed: ${ALIGNED_LEN} → ${TRIMMED_LEN} columns (removed $(( ALIGNED_LEN - TRIMMED_LEN )))"
else
    log_warn "Skipping ClipKit (--skip-clipkit)"
    cp "$ALIGNED" "$TRIMMED"
fi

# ============================================================================
# STEP 3 — FastTree: rapid approximate ML tree
# ============================================================================
if [[ "${SKIP_FASTTREE:-0}" -eq 0 ]]; then
    step 3 "Quick ML tree with FastTree (exploratory)"
    require_cmd fasttree

    # FastTree uses a hybrid approach:
    #   - Neighbor-Joining for the initial topology
    #   - Nearest Neighbor Interchange (NNI) and Subtree-Pruning-Regrafting (SPR)
    #     to optimize the ML score
    # It is approximately 100–1000× faster than IQ-TREE but less accurate.
    # Use this tree for a quick look; rely on IQ-TREE for publication.
    #
    # -nt    : nucleotide data
    # -gtr   : use the GTR substitution model (most general for DNA)
    # -gamma : use a gamma distribution for rate variation across sites
    # -noml  : skip the final ML optimization pass (faster)
    # -boot  : number of bootstrap resamples (local bootstrap, not ML)

    if [[ "$SEQ_TYPE" == "DNA" ]]; then
        fasttree -nt -gtr -gamma -noml -boot "$FASTTREE_BOOTS" \
            < "$TRIMMED" > "$TREE_FASTTREE" 2>"${RUN_DIR}/fasttree.log"
    else
        # For amino acids, FastTree uses the LG model by default
        fasttree -gamma -noml -boot "$FASTTREE_BOOTS" \
            < "$TRIMMED" > "$TREE_FASTTREE" 2>"${RUN_DIR}/fasttree.log"
    fi

    log_success "FastTree done: $TREE_FASTTREE"
else
    log_warn "Skipping FastTree (--skip-fasttree)"
fi

# ============================================================================
# STEP 4 — IQ-TREE: maximum likelihood tree with model selection
# ============================================================================
if [[ "${SKIP_IQTREE:-0}" -eq 0 ]]; then
    step 4 "Maximum likelihood tree with IQ-TREE"
    require_cmd iqtree

    # ModelFinder Plus (MFP) evaluates dozens of substitution models and
    # selects the best one by BIC. This takes extra time but produces a more
    # accurate tree than assuming a fixed model.
    #
    # -B 1000  : ultrafast bootstrap with 1000 resamples.
    #            Ultrafast bootstrap ≠ standard bootstrap — it is much faster
    #            but values > 95 are generally comparable to standard > 70.
    # -T AUTO  : let IQ-TREE choose the number of threads automatically
    #            (it tests performance on your specific alignment)
    # --seqtype : explicit declaration avoids auto-detection errors
    # -pre     : prefix for all output files
    # -quiet   : suppress verbose per-iteration output

    iqtree \
        -s        "$TRIMMED" \
        -m        "$IQTREE_MODEL" \
        -B        "$IQTREE_BOOTSTRAP" \
        --seqtype "$SEQ_TYPE" \
        -T        AUTO \
        -pre      "$TREE_ML_PRE" \
        -quiet

    # Extract the best model chosen by ModelFinder for logging
    BEST_MODEL=$(grep "Best-fit model:" "${TREE_ML_PRE}.log" 2>/dev/null | awk '{print $NF}' || echo "N/A")
    log_success "IQ-TREE done"
    log_info    "Best substitution model: $BEST_MODEL"
    log_info    "ML tree: $TREE_ML"
else
    log_warn "Skipping IQ-TREE (--skip-iqtree)"
fi

# ============================================================================
# STEP 5 — Midpoint rerooting
# ============================================================================
step 5 "Rerooting tree (midpoint)"

# Phylogenetic trees are initially unrooted — the position of the root
# (common ancestor) is not determined by the sequence data alone.
#
# Two rerooting strategies:
#   a) Outgroup rooting: include a known distant relative and place the root
#      on the branch connecting it to the ingroup.
#   b) Midpoint rooting: place the root at the midpoint of the longest
#      path between any two leaves. Assumes a molecular clock.
#
# We use midpoint rooting as a default since no outgroup is specified.

# Choose the best available tree
TARGET_TREE=""
[[ -f "$TREE_ML"       ]] && TARGET_TREE="$TREE_ML"
[[ -z "$TARGET_TREE" && -f "$TREE_FASTTREE" ]] && TARGET_TREE="$TREE_FASTTREE"

if [[ -n "$TARGET_TREE" ]]; then
    if command -v gotree &>/dev/null; then
        gotree reroot midpoint -i "$TARGET_TREE" -o "$TREE_ROOTED"
        log_success "Midpoint-rooted tree: $TREE_ROOTED"
    else
        cp "$TARGET_TREE" "$TREE_ROOTED"
        log_warn "'gotree' not found — tree copied without rerooting."
        log_warn "Install it via: conda install -c bioconda gotree"
    fi
else
    log_warn "No tree file found for rerooting."
fi

# ============================================================================
# STEP 6 — Convert to PHYLIP format (for PhyloBayes)
# ============================================================================
if [[ "${SKIP_PHYLOBAYES:-0}" -eq 0 ]]; then
    step 6 "Converting to PHYLIP format for PhyloBayes"
    require_cmd python3

    # PHYLIP format requires taxon names ≤ 10 characters (strict mode).
    # We replace original names with short IDs (T000000001, T000000002, …)
    # and write a mapping table so you can recover the original names later.
    #
    # PhyloBayes implements Bayesian phylogenetics using MCMC sampling
    # of the posterior distribution P(tree, model | data).

    python3 - "$TRIMMED" "$PHYLIP_OUT" "$PHYLIP_MAP" <<'PY'
import sys
from Bio import SeqIO

trimmed, phylip_out, name_map = sys.argv[1], sys.argv[2], sys.argv[3]
records = list(SeqIO.parse(trimmed, "fasta"))
if not records:
    raise SystemExit("[ERROR] No sequences in trimmed alignment.")

with open(name_map, "w") as mf:
    mf.write("phylip_id\toriginal_id\n")
    for idx, rec in enumerate(records, start=1):
        new_id = f"T{idx:09d}"
        mf.write(f"{new_id}\t{rec.id}\n")
        rec.id = rec.name = new_id
        rec.description = ""

SeqIO.write(records, phylip_out, "phylip-sequential")
print(f"Written {len(records)} sequences to {phylip_out}")
PY

    log_success "PHYLIP file: $PHYLIP_OUT"
    log_success "Name map:    $PHYLIP_MAP"
else
    log_warn "Skipping PhyloBayes prep (--skip-phylobayes)"
fi

# ============================================================================
# STEP 7 — Generate BEAST2 XML alignment draft
# ============================================================================
if [[ "${SKIP_BEAST2:-0}" -eq 0 ]]; then
    step 7 "Generating BEAST2 alignment XML draft"
    require_cmd python3

    # BEAST2 (Bayesian Evolutionary Analysis Sampling Trees) performs
    # Bayesian phylogenetic inference with molecular clock models.
    # It estimates both the tree topology AND divergence times.
    # Input is an XML configuration file that embeds the alignment.
    # This script generates the alignment block; the user still needs to
    # configure the clock model, priors, and MCMC settings in BEAUti.

    python3 - "$TRIMMED" "$BEAST_XML" <<'PY'
import sys
from html import escape
from Bio import SeqIO

trimmed, beast_xml = sys.argv[1], sys.argv[2]
records = list(SeqIO.parse(trimmed, "fasta"))
if not records:
    raise SystemExit("[ERROR] No sequences in trimmed alignment.")

with open(beast_xml, "w") as out:
    out.write('<?xml version="1.0" encoding="UTF-8"?>\n')
    out.write('<beast version="2.7">\n')
    out.write('  <!-- Alignment block — paste into BEAUti or a complete BEAST2 XML -->\n')
    out.write('  <data id="alignment" spec="Alignment" name="alignment" dataType="nucleotide">\n')
    for rec in records:
        out.write(
            f'    <sequence id="seq_{escape(rec.id)}" '
            f'taxon="{escape(rec.id)}" '
            f'value="{escape(str(rec.seq))}"/>\n'
        )
    out.write("  </data>\n")
    out.write("</beast>\n")
print(f"Written {len(records)} sequences to {beast_xml}")
PY

    log_success "BEAST2 XML: $BEAST_XML"
else
    log_warn "Skipping BEAST2 prep (--skip-beast2)"
fi

# ============================================================================
# Summary
# ============================================================================
banner "Pipeline Complete"

echo ""
echo "  Input sequences : $N_SEQS"
echo "  Run directory   : $RUN_DIR"
echo ""
echo "  Key output files:"
echo "  ├── $(basename "$ALIGNED")       — MAFFT alignment"
echo "  ├── $(basename "$TRIMMED")       — Trimmed alignment (ClipKit)"
[[ -f "$TREE_FASTTREE" ]] && echo "  ├── $(basename "$TREE_FASTTREE") — FastTree quick tree"
[[ -f "$TREE_ML"       ]] && echo "  ├── $(basename "$TREE_ML")        — IQ-TREE ML tree"
[[ -f "$TREE_ROOTED"   ]] && echo "  ├── $(basename "$TREE_ROOTED")    — Midpoint-rooted tree"
[[ -f "$PHYLIP_OUT"    ]] && echo "  ├── $(basename "$PHYLIP_OUT")     — PHYLIP alignment"
[[ -f "$BEAST_XML"     ]] && echo "  └── $(basename "$BEAST_XML")      — BEAST2 XML draft"
echo ""

log_success "Total runtime: $(elapsed $PIPELINE_START)"
