# 4. Variant Interpretation

## From Calls to Meaning

Variant calling tells you **where** a genome differs from the reference.
Interpretation tells you **what it means** — Is this difference harmless? Disease-causing?
Clinically actionable? Only a minority of variants in any individual are clinically significant.

A typical WGS experiment on a human sample produces ~4–5 million variants relative to
GRCh38. Of these:
- ~99.5% are common population variants (allele frequency > 1%)
- ~4,000–10,000 are rare (AF < 0.1%)
- ~100–500 are protein-altering
- ~10–50 are in disease-associated genes
- **1–5** may be genuinely clinically significant

The task of interpretation is filtering and evaluating this funnel.

---

## 4.1 The Interpretation Ladder

Work through these levels in order. Stop when you have enough evidence.

```
Level 1 ─── Molecular consequence
             What does the variant do to the transcript?
             (frameshift, missense, synonymous, splice …)
             Tools: VEP, SnpEff, bcftools csq

Level 2 ─── Transcript context
             Which isoform? Is this the canonical transcript?
             Is the affected exon in all isoforms?
             Sources: Ensembl, RefSeq, MANE transcripts

Level 3 ─── Population frequency
             How common is this in healthy people?
             Very common = probably benign
             Sources: gnomAD, TOPMed, 1000 Genomes

Level 4 ─── Known disease links
             Has this exact variant been seen before?
             Is it in a disease gene?
             Sources: ClinVar, OMIM, HGMD

Level 5 ─── Functional predictions
             Do computational tools predict damage?
             (SIFT, PolyPhen-2, CADD, REVEL, SpliceAI)
             Note: weak evidence alone; must converge with other levels

Level 6 ─── Actionability
             Is there clinical guidance? Therapeutic options?
             Sources: CIViC, OncoKB, ClinGen, PharmGKB
```

---

## 4.2 Consequence Vocabulary (SO Terms)

VEP uses **Sequence Ontology (SO)** terms to describe molecular consequences.
These are grouped by expected functional impact:

### HIGH impact — likely disrupts protein function

| Consequence | Example | Mechanism |
|-------------|---------|-----------|
| `stop_gained` | CAA → TAA (Gln → Stop) | Truncates protein |
| `frameshift_variant` | Insertion/deletion not divisible by 3 | Shifts reading frame → likely PTC |
| `splice_donor_variant` | Variant at GT+2 of intron | Destroys 5' splice site |
| `splice_acceptor_variant` | Variant at AG−2 of intron | Destroys 3' splice site |
| `start_lost` | ATG → TTG | Prevents translation initiation |
| `stop_lost` | TAA → CAA | Read-through; extended protein |

### MODERATE impact — may affect protein function

| Consequence | Example |
|-------------|---------|
| `missense_variant` | Single amino acid change (non-conservative) |
| `inframe_insertion` | Extra codon(s) inserted |
| `inframe_deletion` | Codon(s) deleted (frame maintained) |
| `protein_altering_variant` | Non-standard coding change |

### LOW impact — probably harmless

| Consequence | Example |
|-------------|---------|
| `synonymous_variant` | Codon changes but amino acid stays the same |
| `splice_region_variant` | Within 1–3 bp of splice site (not canonical) |
| `stop_retained_variant` | Variant in stop codon but still stop codon |

### MODIFIER — regulatory, UTR, intergenic

| Consequence | Region |
|-------------|--------|
| `5_prime_UTR_variant` | 5' untranslated region |
| `3_prime_UTR_variant` | 3' untranslated region |
| `intron_variant` | Intronic |
| `upstream_gene_variant` | Within 5 kb upstream of gene |
| `intergenic_variant` | Between genes |

---

## 4.3 Population Frequency Databases

### Why frequency matters

A variant present at 5% in the general population is almost certainly not causing a rare
Mendelian disease. If the disease affects 1 in 10,000 people, a disease-causing variant
cannot be this common in controls.

### gnomAD — The Gold Standard

[gnomAD](https://gnomad.broadinstitute.org/) (Genome Aggregation Database) is the largest
public collection of population variant frequencies, aggregated from:
- > 125,000 exomes
- > 15,000 genomes
- Diverse global populations

**Population groups in gnomAD v3:**
African/African American · Latino/Admixed American · Ashkenazi Jewish ·
East Asian · Finnish · Non-Finnish European · South Asian · Other

**Key thresholds:**

| gnomAD AF | Interpretation |
|-----------|---------------|
| > 5% (0.05) | Common — likely benign (BA1 criterion in ACMG) |
| 1–5% | Fairly common — strong evidence against pathogenicity |
| 0.1–1% | Uncommon — could be associated with common disease |
| < 0.1% | Rare — worth investigating |
| Absent (AF = 0) | Ultra-rare or novel — highest priority |

**Important caveat:** gnomAD includes some individuals with disease. A variant absent from
gnomAD is not automatically pathogenic — it might just be rare and benign.

---

## 4.4 Clinical Databases

### ClinVar

[ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) is NCBI's archive of variants submitted
by clinical laboratories with their interpretations. Each submission includes:
- Classification (pathogenic, likely pathogenic, VUS, likely benign, benign)
- The evidence used
- The submitting laboratory
- Review status (number of stars):
  - ★★★★ — Practice guideline
  - ★★★ — Expert panel reviewed
  - ★★ — Multiple submitters with no conflict
  - ★ — Single submitter
  - No star — No assertion criteria

**Searching ClinVar programmatically:**
```bash
# Download all variants for a gene
esearch -db clinvar -query "BRCA2[gene]" | efetch -format tabular

# Or use the API
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term=BRCA2[gene]&retmax=100"
```

### OMIM

[OMIM](https://omim.org/) (Online Mendelian Inheritance in Man) catalogs:
- Gene–disease relationships
- Inheritance patterns (autosomal dominant/recessive, X-linked, etc.)
- Molecular mechanisms
- Clinical descriptions

Every OMIM entry has a numeric ID (e.g., BRCA2 = OMIM:600185).

---

## 4.5 The ACMG/AMP Classification Framework

The American College of Medical Genetics and Genomics (ACMG) and Association for Molecular
Pathology (AMP) published a landmark 2015 guideline for classifying variants into five tiers:

| Class | Name | P(pathogenic) | Meaning |
|-------|------|--------------|---------|
| 5 | **Pathogenic** | > 99% | Causes disease |
| 4 | **Likely pathogenic** | 90–99% | Probably causes disease |
| 3 | **VUS** | 10–89% | Uncertain significance — insufficient evidence |
| 2 | **Likely benign** | 1–10% | Probably harmless |
| 1 | **Benign** | < 1% | Does not cause this disease |

### Evidence criteria

The framework assigns evidence codes that combine to support or refute pathogenicity:

#### Pathogenic evidence codes

| Code | Strength | Criterion |
|------|---------|-----------|
| **PVS1** | Very strong | Null variant (nonsense, frameshift, canonical ±1,2 splice) in a gene where LoF is the known mechanism |
| **PS1** | Strong | Same amino acid change as a previously established pathogenic variant |
| **PS2** | Strong | De novo variant (confirmed parentage) in a patient with the disease |
| **PS3** | Strong | Well-established functional studies demonstrate damage |
| **PS4** | Strong | Variant significantly more frequent in cases than controls |
| **PM1** | Moderate | In a mutational hotspot or critical functional domain |
| **PM2** | Moderate | Absent from controls (gnomAD) — see caveats |
| **PM3** | Moderate | In trans with a known pathogenic variant (recessive) |
| **PM4** | Moderate | In-frame indel affecting protein length |
| **PM5** | Moderate | Novel missense at the same codon as a known pathogenic missense |
| **PM6** | Moderate | Assumed de novo (parentage not confirmed) |
| **PP1** | Supporting | Segregates with disease in multiple affected family members |
| **PP2** | Supporting | Missense in a gene with low missense tolerance and known missense mechanism |
| **PP3** | Supporting | Multiple computational predictions of pathogenicity |
| **PP4** | Supporting | Patient phenotype highly specific for a single-gene disorder |
| **PP5** | Supporting | Reputable source reported pathogenic (without full criteria shared) |

#### Benign evidence codes

| Code | Strength | Criterion |
|------|---------|-----------|
| **BA1** | Stand-alone | AF > 5% in gnomAD → benign (with exceptions) |
| **BS1** | Strong | AF higher than expected for the disorder |
| **BS2** | Strong | Observed in healthy adults for a fully penetrant childhood-onset disease |
| **BS3** | Strong | Functional studies demonstrate no damage |
| **BS4** | Strong | Lack of segregation in affected family members |
| **BP1** | Supporting | Missense in a gene where only truncating variants cause disease |
| **BP2** | Supporting | Observed in trans with a known pathogenic variant (AD disease) |
| **BP3** | Supporting | In-frame indel in a repetitive region without known function |
| **BP4** | Supporting | Multiple computational predictions of benignity |
| **BP5** | Supporting | Variant found in a case where an alternative molecular explanation exists |
| **BP6** | Supporting | Reputable source reported benign |
| **BP7** | Supporting | Synonymous with no predicted splice effect |

### Combining evidence

The combination of codes determines the final classification:
- **Pathogenic:** (1 Very Strong + 1 Strong) or (2 Strong) or (1 Strong + 3 Moderate) etc.
- **Likely Pathogenic:** (1 Very Strong + 1 Moderate) or (1 Strong + 2 Moderate) etc.

Use the [ClinGen Pathogenicity Calculator](https://www.clinicalgenome.org) for the exact rules.

---

## 4.6 Computational Prediction Scores

These tools predict whether a variant is damaging based on sequence features.
They provide **supporting evidence only** — never classify a variant based on
computational predictions alone.

| Tool | Type | Score range | Pathogenic threshold |
|------|------|-------------|---------------------|
| **SIFT** | Evolutionary conservation | 0–1 | < 0.05 |
| **PolyPhen-2** | Structural + evolutionary | 0–1 | > 0.908 (probably damaging) |
| **CADD** | Integrated (many features) | 0–99 (Phred) | > 20 |
| **REVEL** | Ensemble (for missense) | 0–1 | > 0.75 |
| **SpliceAI** | Deep learning (splice effect) | 0–1 per splice type | > 0.5 |
| **AlphaMissense** | Protein structure (DeepMind) | 0–1 | > 0.564 |

**Why PP3 alone is weak evidence:**
These tools are trained on known pathogenic and benign variants. They work well for
well-characterized protein families but can be circular (they may already know about
your variant) and often disagree with each other.

---

## 4.7 Special Scenarios

### Mosaicism

A mosaic variant is present in only a fraction of cells because the mutation occurred
**after fertilization** (postzygotic). Features:
- Allele balance far below 0.5 (e.g., 5–30% AF)
- Often missed by standard filtering thresholds
- Can cause disease even at low allele frequency

```bash
# Lower filtering thresholds for mosaic detection
# In config.sh, set:
MIN_AB=0.05      # detect down to 5% allele frequency
MIN_DP=100       # need much higher depth for reliable low-frequency calls
```

### Compound Heterozygosity

Two different variants in the **same gene**, on **opposite chromosomes** (in trans).
Together they may cause autosomal recessive disease even though neither variant alone
is homozygous.

```
Chromosome 1: ─────── variant A ─────────────
Chromosome 2: ────────────────── variant B ───
                (from mother)    (from father)
```

Detection requires either:
- **Trio sequencing** (parents + child) to assign variants to chromosomes
- **Long-read sequencing** (reads long enough to span both variants)
- Statistical phasing (less reliable)

### Structural Variants in Interpretation

SVs are harder to interpret than SNVs because:
- Breakpoints may be imprecise (CIPOS/CIEND uncertainty intervals)
- The same SV can affect multiple genes
- Population databases for SVs are less complete than for SNVs

Key questions for SV interpretation:
1. Does it delete, disrupt, or duplicate a coding region?
2. Does it separate a gene's promoter from its coding sequence (regulatory)?
3. Is the SV documented in gnomAD-SV or DGV (common structural variation)?

---

## 4.8 Variant Interpretation Workflow

```bash
# 1. Filter to rare, protein-altering variants
bcftools view -f PASS annotated.vcf.gz \
  | bcftools filter -i 'INFO/gnomADg_AF < 0.01 || INFO/gnomADg_AF = "."' \
  | bcftools filter -i 'INFO/IMPACT = "HIGH" || INFO/IMPACT = "MODERATE"' \
  > candidate_variants.vcf

# 2. Count by consequence
bcftools query -f '%INFO/CSQ\n' candidate_variants.vcf \
  | cut -d'|' -f1 | sort | uniq -c | sort -rn

# 3. Extract candidates with ClinVar pathogenic classifications
bcftools filter -i 'INFO/CLIN_SIG ~ "pathogenic"' candidate_variants.vcf

# 4. Manual review of top candidates in IGV
igv.sh -g hg38 sample.bam candidate_variants.vcf
```

---

## 4.9 Resources for Clinical Interpretation

| Resource | URL | Use |
|----------|-----|-----|
| ClinVar | clinvar.ncbi.nlm.nih.gov | Variant classifications |
| gnomAD | gnomad.broadinstitute.org | Population frequencies |
| OMIM | omim.org | Gene–disease relationships |
| ClinGen | clinicalgenome.org | Gene–disease validity, ACMG calculator |
| LOVD | lovd.nl | Locus-specific databases |
| CIViC | civicdb.org | Cancer variant evidence |
| OncoKB | oncokb.org | Oncology biomarkers |
| VarSome | varsome.com | Integrated variant interpretation |
| Franklin | franklin.genoox.com | ACMG-based classification |

---

## Further Reading

- Richards S. et al. (2015) "Standards and guidelines for the interpretation of sequence variants" — Genetics in Medicine (the ACMG/AMP guideline)
- Karczewski K.J. et al. (2020) "The mutational constraint spectrum quantified from variation in 141,456 humans" — Nature (gnomAD v2 paper)
- McLaren W. et al. (2016) "The Ensembl Variant Effect Predictor" — Genome Biology

---

## 🔗 Related Files

- [`pipelines/variants/05_report.py`](../pipelines/variants/05_report.py) — scoring implementation
- [`docs/02_variant_calling.md`](02_variant_calling.md) — how variants are called
- [`docs/formulas_reference.md`](formulas_reference.md) — quality score formulas
