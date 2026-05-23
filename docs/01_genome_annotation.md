# 1. Genome Annotation

## What and Why

Genome annotation is the process of identifying the locations and functions of genes and other
functional elements in a DNA sequence. A raw genome assembly is just a string of A/C/G/T — annotation
adds the biological layer: "this region is a gene, it encodes a ribosomal protein, it belongs to
protein family X."

There are two complementary types:

| Type | Question answered | Output |
|------|-------------------|--------|
| **Structural annotation** | Where are the genes? | GFF3 coordinates |
| **Functional annotation** | What do those genes do? | Database IDs, descriptions |

> Every automated annotation is a **hypothesis**, not a fact. Confidence increases as more
> independent lines of evidence converge on the same conclusion.

---

## 1.1 The Evidence Ladder

Think of annotation as a ladder — each rung adds a more confident interpretation:

```
Level 1 ─── ORF
             Just a stretch of DNA between a start codon (ATG) and a stop codon (TAA/TAG/TGA).
             Every genome has thousands; most are not real genes.

Level 2 ─── Coding potential score
             Statistical models (Prodigal) score how much a sequence "looks like" a real gene
             using hexamer frequencies (codon usage biases built from the genome itself).

Level 3 ─── Homology hit
             BLAST / DIAMOND finds a similar sequence in a curated database (Swiss-Prot, NCBI).
             Strong hits (low E-value, high identity, long alignment) are informative.

Level 4 ─── Domain hit
             HMMER finds a conserved protein domain (Pfam, InterPro) within the sequence.
             Domains are "evolutionary atoms" — they survive even when the overall sequence
             diverges significantly.

Level 5 ─── Curated rules / taxonomy
             Genus- or species-specific databases, known gene names, manually reviewed entries.

Level 6 ─── Human sanity check
             A biologist reviews flagged or unusual annotations manually.
```

---

## 1.2 Structural Annotation: Prokka and Prodigal

### Prokka

[Prokka](https://github.com/tseemann/prokka) (Seemann 2014) is a prokaryote annotation pipeline
that runs in minutes on a laptop. It wraps several tools internally:

| Internal tool | Task |
|--------------|------|
| **Prodigal** | ORF prediction (CDS) |
| **RNAmmer / Barrnap** | rRNA prediction |
| **Aragorn / tRNAscan-SE** | tRNA / tmRNA prediction |
| **Infernal** | ncRNA via Rfam profiles |
| **BLAST / DIAMOND** | Functional assignment |

**Key command and flags:**

```bash
prokka \
    --outdir   ./annotation   \  # output directory
    --prefix   GENOME          \  # prefix for all output files
    --kingdom  Bacteria        \  # affects which reference databases are searched
    --cpus     8               \  # parallel threads
    --evalue   1e-9            \  # E-value cutoff for BLAST searches
    --genus    Bacillus        \  # (optional) improves genus-specific annotations
    --rfam                     \  # enable ncRNA search via Rfam
    assembly.fasta
```

**Output files:**

| File | Format | Contents |
|------|--------|----------|
| `.gff` | GFF3 | Gene coordinates with attributes |
| `.faa` | FASTA | Protein sequences of all predicted CDS |
| `.ffn` | FASTA | Nucleotide sequences of all predicted features |
| `.gbk` | GenBank | Full annotation in GenBank format |
| `.txt` | Text | Summary statistics (CDS count, tRNA count, …) |
| `.tsv` | TSV | Feature table (easy to import into Excel / pandas) |

**Quick post-run analysis:**

```bash
# How many CDS were predicted?
grep -c ">" GENOME.faa

# What fraction are "hypothetical proteins" (no known function)?
grep -c "hypothetical protein" GENOME.faa

# What is the ratio? (a well-annotated genome is typically < 30% hypothetical)
python3 -c "
hyp=$(grep -c 'hypothetical' GENOME.faa)
total=$(grep -c '>' GENOME.faa)
print(f'{hyp}/{total} = {hyp/total*100:.1f}% hypothetical')
"
```

---

### How Prodigal Finds Genes

Prodigal (Hyatt et al. 2010) uses dynamic programming to find the set of ORFs that maximizes
the total score across the genome. The score for each candidate ORF has two parts:

```
Total gene score = Start score + Coding score
```

**Start score:** evaluates whether this ATG is a genuine translation initiation site.
Considers the upstream Shine-Dalgarno (ribosome binding) sequence:

```
   ─────── AGGAG ────────── ATG ─────────────────────
    Shine-Dalgarno motif     Start     Coding region
    (5–10 bp upstream)       codon
```

**Coding score:** log-likelihood ratio of the hexamer frequencies observed in the ORF
compared to a background (non-coding) model:

```
Coding score = log[ P(hexamers | coding) / P(hexamers | non-coding) ]
```

This is why Prodigal needs ≥ 100 kb of sequence to train itself — it builds the hexamer model
from the genome being annotated.

For **metagenomes** (mixed species, no single codon usage), use:

```bash
prodigal -p meta -i metagenome.fasta -a proteins.faa -f gff -o genes.gff
```

---

## 1.3 Homology-Based Annotation: BLAST and E-value

### The E-value Formula

When you run a BLAST or DIAMOND search, the key statistic is the **E-value** (Expect value):

```
E-value = K × m × n × e^(−λ × S)
```

| Symbol | Meaning | Note |
|--------|---------|------|
| `K` | Karlin-Altschul constant | Depends on scoring matrix |
| `m` | Query length (aa or bp) | Your sequence |
| `n` | Total database length | Grows over time → E-values change! |
| `λ` | Decay parameter | ~0.267 for BLOSUM62 |
| `S` | Raw alignment score | Sum of substitution matrix values |
| `e` | Euler's number (2.718) | Natural exponential base |

**Intuition:** E-value = expected number of database hits with score ≥ S purely by chance.
E-value = 1e-5 means we expect one false hit per 100,000 searches.

**Practical thresholds:**

| E-value | Interpretation |
|---------|---------------|
| < 1e-50 | Almost certainly homologous |
| 1e-20 to 1e-50 | Probably homologous, same family |
| 1e-5 to 1e-20 | Likely related, possibly same superfamily |
| 0.001 to 1e-5 | Weak signal — treat with caution |
| > 0.001 | Probably random match |

**Common mistake:** Comparing E-values across searches with different database sizes is misleading.
Use **bit-score** for database-size-independent comparisons:

```
S' (bit-score) = (λ × S − ln K) / ln 2
```

---

### Why DIAMOND Instead of BLASTP?

DIAMOND (Buchfink et al. 2015) reimplements BLASTP with two algorithmic improvements:

1. **Double indexing:** indexes both the query and the database with spaced seeds
2. **SIMD vectorization:** uses CPU vector instructions for the alignment calculation

Result: **500 – 20,000× faster** than BLASTP at similar sensitivity. For annotating a typical
bacterial genome (~4000 proteins) against Swiss-Prot (~500,000 proteins), DIAMOND takes
seconds instead of hours.

```bash
# Build the database once (only needed once per database file)
diamond makedb --in swissprot.fasta --db swissprot --threads 8

# Search (outfmt 6 = tab-separated BLAST format)
diamond blastp \
    --db         swissprot \
    --query      GENOME.faa \
    --out        hits.tsv \
    --outfmt 6 qseqid sseqid pident length evalue bitscore stitle \
    --evalue     1e-5 \
    --max-target-seqs 3 \
    --very-sensitive \
    --threads    8
```

**Key parameters:**

| Flag | Meaning |
|------|---------|
| `--max-target-seqs 3` | Keep only top 3 hits per query (saves disk) |
| `--very-sensitive` | Highest sensitivity mode (slower but catches distant homologs) |
| `--block-size` | RAM vs speed trade-off (default 2.0 GB per thread) |

---

## 1.4 Domain-Based Annotation: HMMER and Pfam

### Why Domains Matter

Evolution recombines protein domains like Lego bricks. The same domain (e.g., a kinase domain,
a DNA-binding domain) appears in thousands of different proteins across distant lineages.

**Domains are more evolutionarily conserved than full-protein sequences.**
HMMER finds domains even when the overall protein identity is too low for BLAST to detect.

### Profile Hidden Markov Models

A Pfam HMM profile is a statistical model built from a multiple alignment of dozens or hundreds
of known family members. It captures:

- Which positions are conserved (high information content)
- Which positions tolerate variation (low information content)
- Typical insertion and deletion patterns

```
                 Position in the domain
                 1  2  3  4  5  6  7  8
                 ─────────────────────
High conservation → A  G  K  T  x  x  x  x   ← ATP-binding P-loop
Any amino acid  → x  x  x  x  A  x  x  x   ← variable loop
```

### Running hmmscan

```bash
# Build binary index (only once per HMM database)
hmmpress Pfam-A.hmm

# Search protein sequences against Pfam profiles
hmmscan \
    --cpu        8 \
    --domtblout  hits.domtbl \    # per-domain tabular output
    --noali \                      # omit alignment text (faster)
    -E           1e-5 \            # sequence-level E-value cutoff
    --domE       1e-3 \            # domain-level E-value cutoff (more permissive)
    Pfam-A.hmm \
    proteins.faa
```

**domtblout format (key columns):**

```
#                                          --- full sequence --- ----------- best 1 domain -----
# target name    accession  query name    E-value   score   E-value   score   c-Evalue
PF00005.30       PF00005    PROTEIN_0001  1.2e-20   72.1    3.4e-08   32.5    0.0011
```

| Column | Description |
|--------|-------------|
| target name | HMM profile name (Pfam family) |
| accession | Pfam accession (PF00005) |
| query name | Protein being searched |
| E-value (full) | E-value considering the whole sequence |
| c-Evalue | Conditional E-value for this domain instance |

---

## 1.5 Functional Databases: KEGG and GO

### KEGG Orthology (KO)

KEGG groups genes from different organisms into **orthologous groups** (KO numbers) linked to
metabolic pathways. KO assignment lets you ask: "Does this genome encode a complete TCA cycle?"

```bash
# Use BlastKOALA webserver (web upload):
# https://www.kegg.jp/blastkoala/

# Or command-line via KAAS:
# Reconstruct pathways at https://www.kegg.jp/kegg/mapper/reconstruct.html
```

**KAAS modified identity formula** (corrects for alignment length bias):

```
modified_identity = pident × min(1, (alignment_length × 2) / (query_length + target_length))
```

Why the correction? A 100-aa protein with 95% identity over 50 aa is NOT the same as a 500-aa
protein with 95% identity over 500 aa. The formula penalizes partial alignments.

### Gene Ontology (GO)

GO provides a controlled vocabulary for gene function across all species, organized into three
independent hierarchies:

| Ontology | Question answered | Example term |
|----------|------------------|--------------|
| **Molecular Function** | What biochemical activity? | GO:0003723 RNA binding |
| **Biological Process** | What pathway/process? | GO:0006281 DNA repair |
| **Cellular Component** | Where in the cell? | GO:0005694 chromosome |

GO terms form a **directed acyclic graph** (DAG): each term can have multiple parent terms.
A gene annotated with "DNA repair" is automatically also annotated with "cellular response to
DNA damage stimulus" (more general), and so on up the hierarchy.

---

## 1.6 Assessing Assembly Quality: N50

Before annotating, check whether the assembly is good enough.

**N50** = the length L such that contigs of length ≥ L collectively cover 50% of the genome.

**Algorithm:**

```python
def n50(lengths: list[int]) -> int:
    lengths.sort(reverse=True)           # 1. Sort descending
    total = sum(lengths)                 # 2. Sum all lengths
    cumsum = 0
    for L in lengths:
        cumsum += L
        if cumsum >= total * 0.5:        # 3. Find the 50% threshold
            return L                     # 4. Return that contig's length
```

**Worked example:**

```
Contig lengths: 5000, 3000, 2000, 1500, 1000   (total = 12500 bp)
50% threshold: 6250 bp

Cumulative sum:
  5000      → 5000  (not yet ≥ 6250)
  5000+3000 → 8000  (≥ 6250 ✓)

N50 = 3000 bp
```

**Interpreting N50 in context:**

| Genome type | Good N50 |
|-------------|---------|
| Bacteria (short-read Illumina) | ≥ 100 kb |
| Bacteria (long-read ONT/PacBio) | ≥ 2 Mb (often chromosome-complete) |
| Human (draft assembly) | ≥ 10 Mb |

---

## 1.7 Quick Annotation Workflow

```bash
# 1. Structural annotation
prokka --outdir annotation --prefix GENOME --kingdom Bacteria --cpus 8 assembly.fasta

# 2. Check annotation statistics
cat annotation/GENOME.txt

# 3. How many proteins?
grep -c ">" annotation/GENOME.faa

# 4. Build DIAMOND database and search Swiss-Prot
diamond makedb --in swissprot.fasta --db swissprot --threads 8
diamond blastp --db swissprot --query annotation/GENOME.faa --out hits.tsv \
    --outfmt 6 qseqid sseqid pident length evalue bitscore stitle \
    --evalue 1e-5 --max-target-seqs 3 --threads 8

# 5. Search Pfam domains
hmmpress Pfam-A.hmm
hmmscan --domtblout domains.txt --noali -E 1e-5 Pfam-A.hmm annotation/GENOME.faa

# 6. Summarize domain hits
grep -v "^#" domains.txt | awk '{print $1}' | sort | uniq -c | sort -rn | head -20
```

---

## Further Reading

- Seemann T. (2014) "Prokka: rapid prokaryotic genome annotation" — Bioinformatics
- Buchfink B. et al. (2015) "Fast and sensitive protein alignment using DIAMOND" — Nature Methods
- Hyatt D. et al. (2010) "Prodigal: prokaryotic gene recognition and translation initiation site identification" — BMC Bioinformatics
- [Pfam user guide](https://pfam.xfam.org/help)
- [KEGG pathway reconstruction tutorial](https://www.kegg.jp/kegg/mapper/)

---

## 🔗 Related Files

- [`pipelines/annotation/run_annotation.sh`](../pipelines/annotation/run_annotation.sh) — implementation
- [`pipelines/annotation/config.sh`](../pipelines/annotation/config.sh) — configuration
- [`docs/formulas_reference.md`](formulas_reference.md) — E-value formula, N50 examples
- [`docs/glossary.md`](glossary.md) — ORF, CDS, GFF3, HMM
