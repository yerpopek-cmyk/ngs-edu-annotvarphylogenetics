# 0. Introduction to Next-Generation Sequencing

> **Start here if you are new to bioinformatics.**  
> This document explains what NGS is, how it produces data, and why we need
> computational analysis — before any pipeline runs.

---

## 0.1 What Is DNA Sequencing?

DNA sequencing determines the precise order of nucleotide bases (A, C, G, T) in a DNA
molecule. The first practical method (Sanger sequencing, 1977) could read ~1000 bases per
reaction. The human genome (~3.2 billion bases) took 13 years and $2.7 billion to sequence
using this approach.

**Next-Generation Sequencing (NGS)** — also called high-throughput sequencing or
massively parallel sequencing — changed everything by sequencing millions of fragments
simultaneously.

---

## 0.2 How Illumina Sequencing Works (The Most Common NGS Platform)

### Step 1 — Library preparation

```
Genomic DNA
    │
    ▼ Fragmentation (sonication or enzymatic)
Short fragments (~150–500 bp)
    │
    ▼ Adapter ligation
Fragments with adapters on both ends
    │
    ▼ PCR amplification
Library (many copies of each fragment)
```

### Step 2 — Cluster generation (bridge amplification)

Each fragment attaches to a solid surface (flow cell) and is amplified in place,
creating a "cluster" of ~1000 identical copies. This amplification is necessary
because the signal from a single molecule would be too weak to detect.

### Step 3 — Sequencing by synthesis

Fluorescently labeled, reversibly terminated nucleotides are incorporated one at a time.
After each cycle, the color of each cluster is imaged:

```
Cycle 1:  add labeled A,C,G,T → one incorporates → image → remove label → repeat
Cycle 2:  same
...
Cycle 150: done → each cluster has yielded 150 base-calls
```

### Step 4 — Base calling and quality scores

The fluorescence intensity at each cycle is converted to a base call and a **Phred quality
score** Q, where:

```
Q = −10 × log₁₀(P_error)
```

A base with Q30 has a 0.1% chance of being called incorrectly.

---

## 0.3 The Output: FASTQ Files

Every NGS experiment produces **FASTQ** files — one or two per sample (R1 and R2 for
paired-end sequencing). Each read is four lines:

```
@read_name              ← Read identifier (starts with @)
ACGTACGTACGTACGT...     ← DNA sequence (150 bases for standard Illumina)
+                       ← Separator (always +)
IIIIFFFFFHHHHHJJ...     ← Phred+33 encoded quality scores (one char per base)
```

**Phred+33 encoding:** the quality score Q is stored as the ASCII character `chr(Q + 33)`.
So Q=40 → `chr(73)` → `I`.

```python
# Decode quality string
def decode_quality(qual_string: str) -> list[int]:
    return [ord(c) - 33 for c in qual_string]

decode_quality("IIIII")  # → [40, 40, 40, 40, 40] — excellent quality
decode_quality("#####")  # → [2, 2, 2, 2, 2]  — terrible quality
```

---

## 0.4 Paired-End Sequencing

Modern Illumina instruments sequence both ends of each DNA fragment:

```
Fragment:   ←───────────────── 350 bp ────────────────→
              │                                      │
              R1 (150 bp →)              (← 150 bp) R2
```

This produces two FASTQ files: `sample_R1.fastq.gz` and `sample_R2.fastq.gz`.
Every read in R1 has a corresponding mate in R2 (same line number, same read name
with `/1` or `/2` suffix).

**Why paired-end?** Knowing that R1 and R2 came from the same fragment gives you:
- Better alignment accuracy (two anchors instead of one)
- Fragment size information (useful for detecting structural variants)
- Higher sensitivity for detecting indels

---

## 0.5 The Reference Genome

A **reference genome** is a representative DNA sequence assembly used as a coordinate system.
Reads from any individual sample are aligned to the reference to identify where they come from
and what differs from it.

```
Reference:  ACGTAGCTAGCTAGCTA
Read 1:     ACGTAGCTAGCTAGCTA   ← perfect match
Read 2:     ACGTAGCTAG T TAGCTA ← SNP (C → T at position 11)
Read 3:     ACGTA--TAGCTAGCTA  ← deletion of 2 bases (GC)
```

**Important:** The reference is not "normal" — it is a mosaic of many individuals and
contains all the variants of those individuals. Absence of a variant in your sample
relative to the reference does not mean it is wild-type.

Human reference genome assemblies:
- **GRCh37 / hg19** — released 2009, still widely used in clinical settings
- **GRCh38 / hg38** — released 2013, current gold standard
- **T2T-CHM13 / hg38** — 2022, first truly complete assembly including centromeres

---

## 0.6 Read Alignment (Mapping)

To compare a sample to the reference, we must find where each read originated. This is
**read alignment** or **mapping**.

### The challenge

- A human genome has 3.2 billion positions
- A typical WGS experiment produces 600 million reads
- Each read must be matched to its source position in seconds

### How BWA-MEM works

BWA (Li & Durbin 2009) uses the **Burrows-Wheeler Transform (BWT)** — a reversible string
transformation that makes suffix searching extremely fast. Key idea:

1. Transform the reference so that identical suffixes group together
2. Build a suffix array index
3. Use FM-index (a compressed representation) for rapid exact matching of read seeds
4. Extend seeds using Smith-Waterman alignment

**MEM (Maximal Exact Match) algorithm:** find the longest exact match between the read and
the reference, then extend outward. Far better than older "backtrack" algorithm for reads ≥ 70 bp.

### BAM format

The alignment output is stored in **BAM** (Binary Alignment Map), the compressed binary
version of SAM (Sequence Alignment Map). Each line describes one aligned read:

```
QNAME     FLAG  RNAME  POS    MAPQ  CIGAR   RNEXT  PNEXT  TLEN  SEQ            QUAL
read001   99    chr1   10001  60    150M    =      10354  503   ACGT...ACGT    IIII...
```

**CIGAR string** describes the alignment: `150M` = 150 matches; `10M2I138M` = 10 matches,
2-base insertion, 138 matches.

**FLAG field** encodes properties as bitwise flags:
- Bit 1 (1): read is paired
- Bit 2 (2): pair is properly mapped
- Bit 64 (64): this is R1
- Bit 128 (128): this is R2
- Bit 1024 (1024): read is a PCR or optical duplicate

---

## 0.7 PCR Duplicates and markdup

Library preparation amplifies each DNA fragment with PCR. Sometimes the same molecule is
amplified more than once, and both copies get sequenced. These **PCR duplicates** look like
real reads but are just copies of the same original molecule.

```
Original fragment: 5'────────────────3' (position 1000–1150)
After PCR:
  Copy 1: 5'────────────────3' (same sequence, same start/end)
  Copy 2: 5'────────────────3' (identical)
  Copy 3: 5'────────────────3' (identical)
```

If you count all three copies as independent evidence for a variant, you overestimate
confidence. `samtools markdup` identifies pairs with identical start positions (likely
duplicates) and tags them in the FLAG field. Variant callers then ignore tagged reads.

---

## 0.8 Types of Genetic Variation

| Type | Abbreviation | Description | Example |
|------|-------------|-------------|---------|
| Single nucleotide variant | SNV / SNP | One base changed | `A → G` at position 12345 |
| Insertion | ins | Extra bases inserted | `ACG → ACTTG` |
| Deletion | del | Bases removed | `ACGT → AT` |
| Indel | indel | Insertion or deletion < 50 bp | Either of above |
| Multi-nucleotide variant | MNV / MNP | Multiple adjacent SNVs | `AC → GT` |
| Structural variant | SV | Large rearrangement > 50 bp | 10 kb deletion, inversion |
| Copy number variant | CNV | Change in number of copies | Gene duplicated 3× |

---

## 0.9 Germline vs Somatic Variants

| | Germline | Somatic |
|--|---------|---------|
| Origin | Inherited from parents | Acquired during life |
| Present in | Every cell in the body | Only in a subset of cells (e.g., tumour) |
| Allele frequency | ~50% (het) or ~100% (hom) | 1–50% (depends on tumour purity) |
| Application | Mendelian disease diagnosis | Cancer genomics |
| Typical caller | FreeBayes, GATK HaplotypeCaller | Mutect2, Strelka2 |

---

## 0.10 The Full NGS Analysis Journey

```
┌─────────────────────────────────────────────────────────────────────────┐
│                    FROM BLOOD TUBE TO CLINICAL REPORT                   │
└─────────────────────────────────────────────────────────────────────────┘

  🧪 WET LAB                          💻 BIOINFORMATICS
  ──────────────                      ─────────────────────────────────────
  DNA extraction                       ↓
  Library preparation                  ↓
  Sequencing                           ↓ FASTQ files
                                       ↓
                              Quality control (FastQC, MultiQC)
                                       ↓
                              Read alignment (BWA, STAR)
                                       ↓ BAM files
                                       ↓
                              Mark PCR duplicates (samtools markdup)
                                       ↓
                              Variant calling (FreeBayes, GATK, Mutect2)
                                       ↓ raw VCF
                                       ↓
                              Normalization + Filtering (bcftools)
                                       ↓ filtered VCF
                                       ↓
                              Functional annotation (VEP, snpEff)
                                       ↓ annotated VCF
                                       ↓
                              Clinical interpretation (ACMG criteria)
                                       ↓
                              Report
```

---

## 0.11 Key Data Volume to Expect

| Data type | Typical size per sample |
|-----------|------------------------|
| Raw FASTQ (WGS 30×) | 30–50 GB (compressed) |
| BAM (aligned) | 30–50 GB |
| VCF (all variants) | 50–500 MB |
| VCF (filtered, PASS only) | 2–20 MB |
| Annotation report | < 1 MB |

This is why large-scale genomics requires computing clusters and object storage —
a cohort of 1000 patients easily exceeds 50 TB.

---

## 🔗 Where to Go Next

1. **Run the annotation pipeline** → [`docs/01_genome_annotation.md`](01_genome_annotation.md)
2. **Learn variant calling details** → [`docs/02_variant_calling.md`](02_variant_calling.md)
3. **Understand phylogenetics** → [`docs/03_phylogenetics.md`](03_phylogenetics.md)
4. **Look up any term** → [`docs/glossary.md`](glossary.md)
