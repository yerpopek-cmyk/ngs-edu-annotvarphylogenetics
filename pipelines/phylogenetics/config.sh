#!/usr/bin/env bash
# =============================================================================
# pipelines/phylogenetics/config.sh — User-editable configuration
#
# ╔══════════════════════════════════════════════════════════════════════╗
# ║  EDIT THIS FILE — do NOT edit run_phylogenetics.sh directly.       ║
# ╚══════════════════════════════════════════════════════════════════════╝
# =============================================================================

# --- Input ---
# Unaligned sequences in FASTA format.
# Must contain ≥ 3 sequences for tree inference.
INPUT_FASTA=""   # e.g., "data/my_sequences.fasta"

# --- Sequence type ---
# DNA : nucleotide sequences (16S rRNA, whole-genome SNPs, …)
# AA  : amino acid sequences (protein families)
SEQ_TYPE="DNA"

# --- Output root directory ---
OUTROOT="${SCRIPT_DIR}/outputs"

# --- Performance ---
THREADS=$(nproc 2>/dev/null || echo 4)

# --- MAFFT alignment ---
# auto: MAFFT selects the strategy based on sequence number and length.
# For proteins with distinct domains, consider --localpair --maxiterate 1000.
# For large datasets (> 1000 sequences), use --retree 2 for speed.
MAFFT_EXTRA_ARGS="--auto"

# --- ClipKit trimming ---
# Method kpic (keep parsimony-informative and constant sites) is recommended
# for phylogenetics — it removes sites that carry no phylogenetic signal.
# Other methods: kpi (keep parsimony-informative only), gappy (remove gappy columns)
CLIPKIT_METHOD="kpic"

# --- FastTree quick tree ---
# Used as a fast exploratory tree before the full ML analysis.
# Bootstrap resamples for FastTree support values:
FASTTREE_BOOTS=100

# --- IQ-TREE maximum likelihood tree ---
# MFP = ModelFinder Plus: IQ-TREE tests many substitution models and
# selects the best one by BIC (Bayesian Information Criterion).
# You can hard-code a model (e.g., "GTR+G+I" for DNA, "LG+G" for proteins)
# to skip model selection and save time.
IQTREE_MODEL="MFP"
IQTREE_BOOTSTRAP=1000

# --- Skip flags ---
# Set to 1 to skip a particular stage (useful if re-running from a checkpoint)
SKIP_MAFFT=0
SKIP_CLIPKIT=0
SKIP_FASTTREE=0
SKIP_IQTREE=0
SKIP_PHYLOBAYES=0
SKIP_BEAST2=0
