# 2. Variant Calling (SNV / Indel)

## The Core Question

A variant caller reads the pile of sequencing reads that align to each genomic position and
statistically asks: **"Is the observed difference from the reference a real genetic variant,
or a sequencing error?"**

This is not trivial. A typical short-read experiment has a per-base error rate of ~0.1–1%.
If a variant exists at 50% allele frequency (heterozygous diploid), roughly half the reads
at that position carry the alternate allele. If a variant exists at 5% frequency (somatic
or low-frequency), only 1 in 20 reads looks different — yet the error rate might be 1 in 100.

---

## 2.1 The VCF Format

All variant callers write output in **VCF** (Variant Call Format). Understanding VCF is
essential before you can filter, annotate, or interpret calls.

### Mandatory columns

```
#CHROM   POS   ID   REF   ALT   QUAL   FILTER   INFO   FORMAT   SAMPLE_NAME
chr1     12345  .    A     G     50.2   PASS     DP=100  GT:DP:AD:GQ:PL  0/1:100:45,55:99:0,50,500
```

| Column | Description | Example |
|--------|-------------|---------|
| CHROM | Chromosome / contig | `chr1` |
| POS | 1-based position | `12345` |
| ID | rsID or `.` if unknown | `.` |
| REF | Reference allele | `A` |
| ALT | Alternate allele(s), comma-separated | `G` |
| QUAL | Phred-scaled variant quality | `50.2` |
| FILTER | `PASS` or filter tag(s) | `PASS` |
| INFO | Semicolon-separated annotations | `DP=100;AF=0.55` |
| FORMAT | Colon-separated field descriptors for samples | `GT:DP:AD` |
| Sample columns | Per-sample values matching FORMAT order | `0/1:100:45,55` |

### Key FORMAT fields

| Field | Full name | Meaning |
|-------|-----------|---------|
| `GT` | Genotype | `0/1` = heterozygous; `1/1` = homozygous alt; `/` = unphased; `|` = phased |
| `DP` | Depth | Total reads at this position in this sample |
| `AD` | Allele Depth | `45,55` = 45 REF reads, 55 ALT reads |
| `GQ` | Genotype Quality | Phred confidence in the genotype call |
| `PL` | Phred-scaled Likelihoods | Log-scaled genotype probabilities (PL=0 for best genotype) |
| `AO` | Alt Observation count | FreeBayes-specific: ALT reads |
| `RO` | Ref Observation count | FreeBayes-specific: REF reads |

---

## 2.2 Phred-Scale Quality Scores

Both sequencing quality and variant quality use the **Phred scale**, named after the
Phred software (Ewing & Green 1998). Phred converts a probability of error to a convenient
positive integer:

```
Q = −10 × log₁₀(P_error)

Inverted:
P_error = 10^(−Q / 10)
```

| Q score | P(error) | Accuracy | In practice |
|---------|----------|----------|-------------|
| Q10 | 0.10 | 90.0% | Terrible — discard |
| Q20 | 0.01 | 99.0% | Minimum acceptable |
| Q30 | 0.001 | 99.9% | Standard threshold for variant calling |
| Q40 | 0.0001 | 99.99% | Good quality |
| Q60 | 0.000001 | 99.9999% | Exceptional |

**QUAL in VCF** = Phred probability that the variant call is wrong (not a real variant).
`QUAL = 30` means a 0.1% chance the called variant is a false positive.

---

## 2.3 Genotype Likelihoods (PL and GQ)

For a diploid sample at a biallelic site, there are three possible genotypes:
`0/0` (ref/ref), `0/1` (ref/alt), `1/1` (alt/alt).

**PL (Phred Likelihood)** scores each genotype:

```
PL(genotype) = −10 × log₁₀[ P(reads | genotype) ]
```

PL values are normalized so the most likely genotype has PL = 0.

**Example** (heterozygous call):

| Genotype | Raw likelihood | PL |
|----------|---------------|-----|
| `0/0` | 0.00001 | **50** |
| `0/1` | 1.0 | **0** ← best |
| `1/1` | 0.01 | **20** |

VCF stores: `PL=50,0,20` in the order: ref/ref, ref/alt, alt/alt.

**GQ (Genotype Quality):**

```
GQ = min(99, PL_second_best − PL_best) = min(99, 20 − 0) = 20
```

GQ tells you how much better the best genotype is than the second-best. GQ = 20 means
the probability of the wrong genotype is 100× higher than the called one.

---

## 2.4 How FreeBayes Works

FreeBayes (Garrison & Marth 2012) is a **haplotype-based** variant caller, which means
instead of looking at single positions independently, it constructs local haplotypes from
the reads and evaluates genotype likelihoods over those haplotypes.

### Key steps

1. **Window construction:** Collect all reads overlapping a region into a local window.

2. **Candidate haplotype discovery:** Identify all unique sequences seen in the reads.
   This naturally handles MNPs (multi-nucleotide polymorphisms) and nearby indels.

3. **Likelihood computation:** For each possible genotype (combination of haplotypes),
   compute the probability of observing the read data:
   ```
   P(reads | genotype) = ∏_read P(read | genotype)
   ```
   where each read likelihood is computed from base qualities.

4. **Genotype calling:** Select the genotype that maximizes the posterior probability.
   Apply a Bayesian prior for allele frequencies in the population.

5. **Output:** Write VCF with QUAL = −10 log₁₀(1 − P_MAP_genotype).

### Key parameters

```bash
freebayes \
    -f reference.fasta \
    -b aligned.bam \
    --ploidy 2 \                        # diploid organism
    --min-base-quality 20 \             # ignore bases with Phred < 20
    --min-alternate-count 2 \           # at least 2 reads must support the ALT
    --min-alternate-fraction 0.2 \      # ALT / (REF+ALT) must be ≥ 0.2
    --pooled-continuous \               # for pooled samples or non-diploid data
    > variants.vcf
```

---

## 2.5 Variant Normalization

The same biological deletion can be written multiple ways in VCF:

```
chr1  100  .  ATT  A   .  PASS  .     # deletion of TT, starting at position 100
chr1  101  .  TT   .   .  PASS  .     # same deletion, starting at position 101
```

**Left-alignment** is the convention: shift indels as far left (toward lower genomic
coordinates) as possible while keeping the sequence equivalent. This ensures that
different callers produce the same representation of the same variant.

**Multi-allelic splitting:** Two ALT alleles at one site should be split into two records
for most downstream tools:

```
chr1  100  .  A  T,G  .  PASS  .     # multi-allelic: split into two records:
chr1  100  .  A  T    .  PASS  .
chr1  100  .  A  G    .  PASS  .
```

```bash
bcftools norm \
    --fasta-ref reference.fasta \   # needed for left-alignment
    -m -any \                        # split multi-allelic sites
    input.vcf.gz \
    -Oz -o normalized.vcf.gz
```

---

## 2.6 Soft vs Hard Filtering

### Hard filtering (removes records)
```bash
bcftools filter -i 'QUAL >= 30 && INFO/DP >= 10' -Oz -o filtered.vcf.gz input.vcf.gz
```
Variants below threshold are deleted. Simple but irreversible — if the threshold was wrong,
you must re-run the caller.

### Soft filtering (marks records)
```bash
bcftools filter \
    --soft-filter LowQual \    # adds "LowQual" to FILTER column (not "PASS")
    --mode + \                  # keep existing FILTER tags, add to them
    --exclude 'QUAL < 30 || INFO/DP < 10' \
    -Oz -o softfiltered.vcf.gz input.vcf.gz
```

After soft-filtering, the VCF still contains all variants. Tools that want only
high-confidence calls use `-f PASS` to select them:

```bash
bcftools view -f PASS softfiltered.vcf.gz   # only PASS variants
bcftools view -f ''   softfiltered.vcf.gz   # all variants including filtered
```

**Why soft-filter?** Clinical pipelines require an audit trail. Soft filtering allows
you to change thresholds later without re-running the caller, and preserves low-frequency
variants that may be interesting for other analyses.

---

## 2.7 Quality Metrics and Red Flags

### Allele Balance (AB)

```
AB = ALT_reads / (REF_reads + ALT_reads)
```

For a heterozygous diploid variant, we expect AB ≈ 0.5. Extreme values suggest artifacts:

| AB range | Interpretation |
|----------|---------------|
| 0.45–0.55 | Expected heterozygote |
| 0.10–0.45 | Heterozygote with strand bias, or low-frequency variant |
| < 0.10 | Very low frequency — may be real or artifact, needs scrutiny |
| > 0.90 | Should be homozygous ALT; if called het, check depth |

### Strand Bias (FS — Fisher Strand)

Sequencing artifacts often affect reads from only one strand. The Fisher Strand score tests
whether REF and ALT reads are equally distributed across forward and reverse strands:

```
FS = −10 × log₁₀(p_Fisher)
```

where `p_Fisher` is the p-value of Fisher's exact test on the 2×2 contingency table:

```
              Forward reads   Reverse reads
REF allele:       30              35
ALT allele:       28               2     ← suspicious: ALT only on forward strand!
```

| FS value | Interpretation |
|----------|---------------|
| < 10 | No strand bias ✅ |
| 10–30 | Mild bias, investigate |
| 30–60 | Moderate bias ⚠️ |
| > 60 | Strong bias — likely artifact 🚨 |

### Five Red Flags for False Positive Variants 🚩

1. **Low mapping quality (MQ < 30):** Reads that align to multiple places — often in repeats.
2. **Strong strand bias (FS > 60):** ALT reads only on one strand.
3. **Only 1–2 ALT reads:** Even at low depth, this is not enough statistical evidence.
4. **Variant at read edge:** Alignment errors are more common at the ends of reads.
5. **Abnormal allele balance:** AB very far from 0.5 for a supposed heterozygote.

---

## 2.8 Quality Control Metrics

### Ti/Tv Ratio (Transition/Transversion)

DNA substitutions are not equally likely. **Transitions** (A↔G, C↔T — purine↔purine or
pyrimidine↔pyrimidine) occur ~2–3× more often than **transversions** (A↔C, A↔T, G↔C, G↔T).

```
Ti/Tv = count(transitions) / count(transversions)
```

| Sample type | Expected Ti/Tv |
|-------------|---------------|
| Human WGS (SNPs) | 2.0–2.2 |
| Human WES (exome) | 2.8–3.0 |
| Too low (< 1.5) | Excess false positives (artifacts have random Ti/Tv ≈ 0.5) |
| Too high (> 3.5) | Possible systematic error, check read trimming |

```bash
bcftools stats sample.vcf.gz | grep "^TSTV"
# TSTV   id  ts  tv  ts/tv  ts (1st ALT)  tv (1st ALT)  ts/tv (1st ALT)
# TSTV   0   1234  567  2.18   1234         567           2.18
```

### Coverage and Depth

```
Mean depth = (N_reads × read_length) / genome_length
```

| Depth | Use case |
|-------|---------|
| 5–15× | Germline calling in research (minimum) |
| 30× | Standard WGS germline |
| 60–100× | Somatic variant calling (tumour) |
| > 200× | Low-frequency variant detection |

---

## 2.9 Functional Annotation with VEP

After calling and filtering variants, we want to know: **what is the biological consequence
of each variant?** Ensembl VEP (Variant Effect Predictor) answers this.

### Consequence hierarchy (IMPACT levels)

| IMPACT | Consequence terms | Effect |
|--------|------------------|--------|
| **HIGH** | stop_gained, frameshift_variant, splice_donor_variant, splice_acceptor_variant, start_lost | Likely disrupts protein function |
| **MODERATE** | missense_variant, inframe_insertion, inframe_deletion, protein_altering_variant | May affect protein function |
| **LOW** | synonymous_variant, splice_region_variant | Unlikely to affect function |
| **MODIFIER** | intron_variant, intergenic_variant, upstream_gene_variant, UTR variant | Regulatory or unknown effect |

### VCF CSQ field

VEP adds annotations in a semicolon-separated `CSQ` INFO field. Each transcript that
overlaps the variant gets one pipe-separated block:

```
CSQ=missense_variant|MODERATE|BRCA2|ENSG00000139618|ENST00000380152|c.2T>A|p.Met1Lys|0.00001|.|Pathogenic
```

Parse with:
```bash
# Count variants by consequence
bcftools query -f '%INFO/CSQ\n' sample.annotated.vcf.gz \
    | tr ',' '\n' | cut -d'|' -f1 | sort | uniq -c | sort -rn
```

---

## 2.10 Full Pipeline Overview

```bash
# 1. Index reference (once)
bwa index reference.fasta
samtools faidx reference.fasta

# 2. Align reads
bwa mem -t 8 -R "@RG\tID:sample\tSM:sample\tPL:ILLUMINA" \
    reference.fasta R1.fastq.gz R2.fastq.gz \
  | samtools fixmate -m - - \
  | samtools sort -@4 -m 2G \
  | samtools markdup --write-index - sample.bam

# 3. Alignment QC
samtools flagstat sample.bam

# 4. Call variants
freebayes -f reference.fasta -b sample.bam \
    --min-base-quality 20 --min-alternate-fraction 0.2 \
  | bcftools sort | bgzip > raw.vcf.gz
tabix -p vcf raw.vcf.gz

# 5. Normalize and filter
bcftools norm -f reference.fasta -m -any raw.vcf.gz \
  | bcftools filter --soft-filter LowQual -e 'QUAL < 30 || INFO/DP < 10' \
  -Oz -o filtered.vcf.gz
tabix -p vcf filtered.vcf.gz

# 6. QC: Ti/Tv ratio
bcftools stats -f PASS filtered.vcf.gz | grep "^TSTV"

# 7. Annotate with VEP
vep --input_file filtered.vcf.gz --output_file annotated.vcf.gz \
    --format vcf --vcf --compress_output bgzip \
    --offline --cache --dir_cache vep_cache/ --assembly GRCh38 \
    --canonical --hgvs --af_gnomade --af_gnomadg --check_existing

# 8. View high-impact variants
bcftools view -f PASS annotated.vcf.gz \
  | bcftools query -f '[%SAMPLE] %CHROM:%POS %REF>%ALT %INFO/CSQ\n'
```

---

## Further Reading

- Garrison E. & Marth G. (2012) "Haplotype-based variant detection from short-read sequencing" — arXiv
- McLaren W. et al. (2016) "The Ensembl Variant Effect Predictor" — Genome Biology
- [bcftools manual](https://samtools.github.io/bcftools/bcftools.html)
- [VCF specification](https://samtools.github.io/hts-specs/VCFv4.3.pdf)
- [gnomAD browser](https://gnomad.broadinstitute.org/)

---

## 🔗 Related Files

- [`pipelines/variants/`](../pipelines/variants/) — implementation
- [`docs/04_variant_interpretation.md`](04_variant_interpretation.md) — clinical significance
- [`docs/formulas_reference.md`](formulas_reference.md) — QUAL, PL, GQ, FS formulas
