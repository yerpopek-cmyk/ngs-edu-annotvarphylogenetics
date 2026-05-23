# =============================================================================
# Makefile — NGS Education Hub
#
# Usage: make <target>
#   make help       — show this menu
#   make setup      — create all conda environments
#   make test-data  — generate synthetic test datasets
#   make test-annot — run annotation on test data
#   make test-phyl  — run phylogenetics on test data
#   make clean      — remove all pipeline outputs
# =============================================================================

.PHONY: help setup test-data test-annot test-phyl clean lint check-envs

# ── Color codes for pretty output ────────────────────────────────────────────
BLUE  := \033[0;34m
CYAN  := \033[0;36m
GREEN := \033[0;32m
BOLD  := \033[1m
RESET := \033[0m

# ── Help ─────────────────────────────────────────────────────────────────────
help:
	@echo ""
	@echo "$(BOLD)$(BLUE)╔══════════════════════════════════════════════════════╗$(RESET)"
	@echo "$(BOLD)$(BLUE)║           NGS Education Hub — Makefile              ║$(RESET)"
	@echo "$(BOLD)$(BLUE)╚══════════════════════════════════════════════════════╝$(RESET)"
	@echo ""
	@echo "$(BOLD)SETUP$(RESET)"
	@echo "  $(CYAN)make setup$(RESET)          Create all three conda environments"
	@echo "  $(CYAN)make check-envs$(RESET)     Verify all environments are installed"
	@echo ""
	@echo "$(BOLD)TEST DATA$(RESET)"
	@echo "  $(CYAN)make test-data$(RESET)      Generate synthetic datasets (no internet needed)"
	@echo ""
	@echo "$(BOLD)PIPELINE TESTS$(RESET)"
	@echo "  $(CYAN)make test-annot$(RESET)     Annotation pipeline: B. subtilis download + annotate"
	@echo "  $(CYAN)make test-phyl$(RESET)      Phylogenetics: align + trim + tree (test 16S data)"
	@echo ""
	@echo "$(BOLD)MAINTENANCE$(RESET)"
	@echo "  $(CYAN)make clean$(RESET)          Remove all generated outputs/"
	@echo "  $(CYAN)make lint$(RESET)           Check bash scripts for common errors (shellcheck)"
	@echo ""
	@echo "$(BOLD)LEARNING PATH$(RESET)"
	@echo "  1. make setup          — install environments"
	@echo "  2. make test-data      — generate test sequences"
	@echo "  3. make test-annot     — see annotation in action"
	@echo "  4. make test-phyl      — see phylogenetics in action"
	@echo "  5. read docs/          — understand the theory"
	@echo ""

# ── Environment setup ─────────────────────────────────────────────────────────
setup:
	@echo "$(BOLD)Creating conda environments (this may take 5–15 minutes)...$(RESET)"
	@echo ""
	@echo "$(CYAN)[1/3] annotation$(RESET)"
	conda env create -f environments/annotation.yml || conda env update -f environments/annotation.yml
	@echo ""
	@echo "$(CYAN)[2/3] variants$(RESET)"
	conda env create -f environments/variants.yml || conda env update -f environments/variants.yml
	@echo ""
	@echo "$(CYAN)[3/3] phylogenetics_env$(RESET)"
	conda env create -f environments/phylogenetics_env.yml || conda env update -f environments/phylogenetics_env.yml
	@echo ""
	@echo "$(GREEN)✓ All environments installed$(RESET)"
	@echo "Activate with: conda activate annotation | variants | phylogenetics_env"

check-envs:
	@echo "Checking conda environments..."
	@conda env list | grep -E "annotation|variants|phylogenetics_env" \
		|| echo "Some environments are missing — run: make setup"

# ── Test data generation ──────────────────────────────────────────────────────
test-data:
	@echo "$(BOLD)Generating synthetic test datasets...$(RESET)"
	mkdir -p test_data
	python3 scripts/generate_test_data.py --outdir test_data/ --seed 42
	@echo ""
	@echo "$(GREEN)✓ Test data ready in test_data/$(RESET)"
	@echo ""
	@echo "  test_data/phylogenetics/test_16S.fasta   — 20 synthetic 16S sequences"
	@echo "  test_data/annotation/test_genome.fasta   — synthetic bacterial genome (10 kb)"
	@echo "  test_data/variants/test_ref.fasta        — synthetic reference (50 kb)"
	@echo "  test_data/variants/test_R1/R2.fastq.gz   — synthetic reads (5000 × 150 bp)"

# ── Annotation test ───────────────────────────────────────────────────────────
test-annot:
	@echo "$(BOLD)Running annotation pipeline test (downloads ~44 MB on first run)...$(RESET)"
	@echo "Activate environment first: conda activate annotation"
	@echo ""
	conda run -n annotation bash pipelines/annotation/run_annotation.sh

# ── Phylogenetics test ────────────────────────────────────────────────────────
test-phyl: test-data
	@echo "$(BOLD)Running phylogenetics pipeline test...$(RESET)"
	conda run -n phylogenetics_env \
		bash pipelines/phylogenetics/run_phylogenetics.sh \
		-i test_data/phylogenetics/test_16S.fasta \
		-t DNA
	@echo "$(GREEN)✓ Tree files written to pipelines/phylogenetics/outputs/$(RESET)"

# ── Lint bash scripts ─────────────────────────────────────────────────────────
lint:
	@echo "Running shellcheck on all bash scripts..."
	@if command -v shellcheck >/dev/null 2>&1; then \
		find pipelines/ shared/ -name "*.sh" -exec shellcheck {} +; \
		echo "$(GREEN)✓ All scripts passed shellcheck$(RESET)"; \
	else \
		echo "shellcheck not found — install with: conda install -c conda-forge shellcheck"; \
	fi

# ── Clean ─────────────────────────────────────────────────────────────────────
clean:
	@echo "Removing all pipeline outputs..."
	rm -rf pipelines/annotation/outputs/
	rm -rf pipelines/variants/outputs/
	rm -rf pipelines/phylogenetics/outputs/
	@echo "$(GREEN)✓ Clean done (environments and test_data/ preserved)$(RESET)"
