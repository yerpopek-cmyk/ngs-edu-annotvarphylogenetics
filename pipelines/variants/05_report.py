#!/usr/bin/env python3
"""
05_report.py — Variant prioritization and Markdown report generation
====================================================================

WHAT THIS SCRIPT DOES
─────────────────────
Reads the VEP-annotated VCF produced by step 04, parses the CSQ INFO field,
computes a composite priority score for each PASS variant, and writes:

  1. A TSV table with all PASS variants and their scores
  2. A Markdown report showing the top-N variants with explanations

SCORING LOGIC
─────────────
The priority score is a simple additive model:

  SCORE = IMPACT_score + gnomAD_AF_score + ClinVar_score

  IMPACT_score:
    HIGH     = 4  (frameshift, stop_gained, splice±2 …)
    MODERATE = 2  (missense, in-frame indel …)
    LOW      = 1  (synonymous, splice_region …)
    MODIFIER = 0  (UTR, intron, intergenic …)

  gnomAD_AF_score (population allele frequency):
    AF < 0.001  → +3  (very rare; more likely pathogenic)
    AF < 0.01   → +1  (rare)
    AF > 0.05   → −1  (common; unlikely to be pathogenic)
    not in gnomAD → +2  (ultra-rare or de novo)

  ClinVar_score:
    Pathogenic           → +5
    Likely_pathogenic    → +3
    Uncertain_significance → 0
    Likely_benign        → −1
    Benign               → −2

NOTE: This is an educational heuristic model.  Clinical interpretation
requires additional evidence (ACMG/AMP criteria, functional studies, etc.).

USAGE
─────
  python3 05_report.py <annotated.vcf.gz> <output_dir> [--top N] [--sample-id ID]
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path
from datetime import datetime


# =============================================================================
# Score tables
# =============================================================================

IMPACT_SCORE: dict[str, int] = {
    "HIGH": 4, "MODERATE": 2, "LOW": 1, "MODIFIER": 0,
}

# ClinVar values arrive lowercase from VEP; normalise spaces → underscores
CLNSIG_SCORE: dict[str, int] = {
    "pathogenic":            +5,
    "likely_pathogenic":     +3,
    "uncertain_significance": 0,
    "vus":                    0,
    "likely_benign":         -1,
    "benign":                -2,
}


# =============================================================================
# CLI
# =============================================================================

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("vcf",       help="Annotated VCF (VEP, bgzip-compressed)")
    p.add_argument("outdir",    help="Output directory")
    p.add_argument("--top",     type=int, default=20,
                   help="Number of top variants in the report (default: 20)")
    p.add_argument("--sample-id", default="",
                   help="Override sample ID (default: derived from filename)")
    return p.parse_args()


# =============================================================================
# VCF parsing helpers
# =============================================================================

def csq_fields_from_header(header_lines: list[str]) -> list[str]:
    """
    Extract the CSQ field names from the VEP header line.

    VEP writes a line like:
      ##INFO=<ID=CSQ,…,Description="…Format: Consequence|IMPACT|SYMBOL|…">
    We parse the Format: section to get the field names in order.
    """
    for line in header_lines:
        if "ID=CSQ" in line and "Format:" in line:
            fragment = line.split("Format:")[1].rstrip('">').strip()
            return [f.strip() for f in fragment.split("|")]
    return []


def info_value(info: str, key: str) -> str:
    """Return the value of KEY= from a VCF INFO string, or empty string."""
    prefix = f"{key}="
    for field in info.split(";"):
        if field.startswith(prefix):
            return field[len(prefix):]
    return ""


def compute_score(impact: str, clnsig_raw: str, af_raw: str) -> tuple[int, float | None]:
    """
    Calculate the priority score and parsed allele frequency.

    Parameters
    ----------
    impact    : VEP IMPACT value (HIGH/MODERATE/LOW/MODIFIER)
    clnsig_raw: raw ClinVar CLIN_SIG value (may contain '&'-joined terms)
    af_raw    : raw gnomAD AF string (may be '.', empty, or '&'-joined)

    Returns
    -------
    (score, af) where af is float or None if not available
    """
    score = IMPACT_SCORE.get(impact, 0)

    # ClinVar — a variant can have multiple significance terms joined by '&'
    if clnsig_raw:
        for sig in clnsig_raw.lower().replace(" ", "_").split("&"):
            score += CLNSIG_SCORE.get(sig.strip(), 0)

    # gnomAD allele frequency
    af: float | None = None
    if af_raw and af_raw not in ("", "."):
        try:
            # Take the first value if multiple are '&'-joined
            af = float(af_raw.split("&")[0])
        except ValueError:
            pass

    if af is None:
        score += 2   # not in gnomAD → ultra-rare bonus
    elif af < 0.001:
        score += 3
    elif af < 0.01:
        score += 1
    elif af > 0.05:
        score -= 1

    return score, af


# =============================================================================
# Main parsing logic
# =============================================================================

def parse_vcf(vcf_path: str) -> list[dict]:
    """
    Parse VEP-annotated VCF and return a list of variant dicts, sorted by
    descending priority score.  Only PASS (or unfiltered) variants are kept.
    """
    open_fn = gzip.open if vcf_path.endswith(".gz") else open
    header_lines: list[str] = []
    csq_fields:   list[str] = []
    variants:     list[dict] = []

    with open_fn(vcf_path, "rt") as fh:
        for raw_line in fh:
            line = raw_line.rstrip()

            # Collect header; extract CSQ field definitions on the fly
            if line.startswith("##"):
                header_lines.append(line)
                if not csq_fields and "ID=CSQ" in line:
                    csq_fields = csq_fields_from_header([line])
                continue

            if line.startswith("#"):
                continue   # column-names row; skip

            cols = line.split("\t")
            if len(cols) < 8:
                continue

            chrom, pos, _id, ref, alt, qual, filt, info = cols[:8]

            # Skip soft-filtered variants (LowQual, LowAB, …)
            if filt not in (".", "PASS", ""):
                continue

            csq_raw = info_value(info, "CSQ")
            if not csq_raw:
                continue

            # VEP lists transcripts in decreasing severity order;
            # take the first (most severe) transcript annotation.
            first_csq = csq_raw.split(",")[0].split("|")
            csq: dict[str, str] = (
                dict(zip(csq_fields, first_csq)) if csq_fields else {}
            )

            impact = csq.get("IMPACT", "MODIFIER")
            gene   = csq.get("SYMBOL", ".")
            conseq = csq.get("Consequence", ".").replace("_variant", "")
            hgvsp  = csq.get("HGVSp", ".")

            # Prefer gnomAD genome AF; fall back to exome AF
            af_raw = (
                csq.get("gnomADg_AF")
                or csq.get("gnomADe_AF")
                or ""
            )
            clnsig = csq.get("CLIN_SIG", "")

            score, af = compute_score(impact, clnsig, af_raw)

            variants.append({
                "CHROM":       chrom,
                "POS":         pos,
                "REF":         ref,
                "ALT":         alt,
                "GENE":        gene,
                "CONSEQUENCE": conseq,
                "IMPACT":      impact,
                "HGVSp":       hgvsp,
                "AF":          f"{af:.5f}" if af is not None else "novel",
                "CLIN_SIG":    clnsig or ".",
                "SCORE":       score,
            })

    variants.sort(key=lambda v: v["SCORE"], reverse=True)
    return variants


# =============================================================================
# Output writers
# =============================================================================

def write_tsv(variants: list[dict], path: Path) -> None:
    """Write all variants as a tab-separated table."""
    if not variants:
        path.write_text("# No PASS variants found\n")
        return
    keys = list(variants[0].keys())
    lines = ["\t".join(keys)]
    for v in variants:
        lines.append("\t".join(str(v[k]) for k in keys))
    path.write_text("\n".join(lines) + "\n")


def write_report(variants: list[dict], path: Path, top_n: int, sample_id: str) -> None:
    """Write a Markdown report with a summary table and top-N variant table."""
    now  = datetime.now().strftime("%Y-%m-%d %H:%M")
    top  = variants[:top_n]

    n_high   = sum(1 for v in variants if v["IMPACT"] == "HIGH")
    n_mod    = sum(1 for v in variants if v["IMPACT"] == "MODERATE")
    n_path   = sum(1 for v in variants if "pathogenic" in v["CLIN_SIG"].lower())
    n_novel  = sum(1 for v in variants if v["AF"] == "novel")

    # Markdown table header
    header = ("| # | Gene | Position | Change | Consequence "
               "| IMPACT | gnomAD AF | ClinVar | Score |")
    sep    = ("|---|------|----------|--------|-------------|"
               "--------|-----------|---------|-------|")
    rows = []
    for i, v in enumerate(top, start=1):
        rows.append(
            f"| {i} | **{v['GENE']}** | {v['CHROM']}:{v['POS']} "
            f"| {v['REF']}→{v['ALT']} | `{v['CONSEQUENCE']}` "
            f"| **{v['IMPACT']}** | {v['AF']} | {v['CLIN_SIG']} | {v['SCORE']} |"
        )

    path.write_text(f"""\
# Variant Prioritization Report

**Sample:** {sample_id or "unknown"}  
**Generated:** {now}  

---

## Summary

| Metric | Count |
|--------|-------|
| Total PASS variants | {len(variants)} |
| HIGH impact | {n_high} |
| MODERATE impact | {n_mod} |
| Pathogenic / Likely pathogenic (ClinVar) | {n_path} |
| Novel (not in gnomAD) | {n_novel} |

---

## Top {top_n} Priority Variants

Sorted by composite score (IMPACT + rarity + ClinVar significance).

{header}
{sep}
{chr(10).join(rows)}

---

## Scoring Model

> **This is an educational heuristic model — not for clinical use.**

```
SCORE = IMPACT_score + gnomAD_AF_score + ClinVar_score

IMPACT:       HIGH=4  MODERATE=2  LOW=1  MODIFIER=0
gnomAD AF:    novel → +2 | <0.001 → +3 | <0.01 → +1 | >0.05 → −1
ClinVar:      Pathogenic → +5 | Likely pathogenic → +3
              Likely benign → −1 | Benign → −2
```

### ACMG/AMP Classification Framework

Real clinical variant interpretation follows the ACMG/AMP 2015 guidelines,
which categorise evidence as:

| Code | Category | Example criteria |
|------|----------|-----------------|
| PVS1 | Very strong pathogenic | Null variant in gene where LOF is disease mechanism |
| PS1–PS4 | Strong pathogenic | Same AA change as known pathogenic; de novo; functional evidence |
| PM1–PM6 | Moderate pathogenic | Located in hotspot; absent from controls (PM2) |
| PP1–PP5 | Supporting pathogenic | Co-segregation; computational prediction |
| BA1 | Stand-alone benign | AF > 5% in gnomAD |
| BS1–BS4 | Strong benign | AF too high for disorder; no disease in carriers |
| BP1–BP7 | Supporting benign | Synonymous in non-splice gene; no functional effect |

Variants are classified as: **Pathogenic · Likely Pathogenic · VUS ·  
Likely Benign · Benign** based on the combination of evidence codes.

---

## Further Reading

- [ACMG/AMP 2015 standards](https://www.nature.com/articles/gim201530)
- [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/)
- [gnomAD browser](https://gnomad.broadinstitute.org/)
- [VEP documentation](https://www.ensembl.org/info/docs/tools/vep/)
""")


# =============================================================================
# Entry point
# =============================================================================

def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    sample_id = args.sample_id or Path(args.vcf).name.split(".")[0]

    print(f"[05_report] Parsing: {args.vcf}")
    variants = parse_vcf(args.vcf)
    print(f"[05_report] PASS variants found: {len(variants)}")

    if not variants:
        print("[05_report] No PASS variants to report — check filters in step 03.",
              file=sys.stderr)
        return

    tsv_path = outdir / f"{sample_id}.prioritized.tsv"
    md_path  = outdir / f"{sample_id}.report.md"

    write_tsv(variants, tsv_path)
    write_report(variants, md_path, args.top, sample_id)

    print(f"[05_report] TSV    → {tsv_path}")
    print(f"[05_report] Report → {md_path}")


if __name__ == "__main__":
    main()
