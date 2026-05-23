# 📐 Formulas Reference

Every formula used across this repository, with worked examples and parameter explanations.

---

## 1. Sequence Alignment and Database Search

### 1.1 E-value (BLAST / DIAMOND)

**Definition:** Expected number of alignments with a score ≥ S occurring by chance in a
database of the given size.

```
E-value = K × m × n × e^(−λ × S)
```

| Symbol | Meaning | Typical value |
|--------|---------|--------------|
| K | Karlin-Altschul constant (scale) | ~0.041 for BLOSUM62 |
| m | Query sequence length (aa or bp) | 100–10,000 |
| n | Total database length (aa or bp) | 10⁸–10¹¹ |
| λ | Decay rate parameter | ~0.267 for BLOSUM62 |
| S | Raw alignment score | depends on matrix |
| e | Euler's number | 2.71828… |

**Significance thresholds:**

| E-value | Interpretation |
|---------|---------------|
| < 1e-50 | Extremely significant — same protein family |
| 1e-20 to 1e-50 | Strong homology |
| 1e-5 to 1e-20 | Moderate homology |
| 0.001 to 1e-5 | Weak signal — treat carefully |
| > 0.001 | Likely random match |

**Important:** E-value depends on database size `n`. A search against a small database
yields higher (worse) E-values than the same search against a large one for the same
alignment. Use **bit-score** for database-size-independent comparisons.

---

### 1.2 Bit-Score

**Definition:** Normalized alignment score, independent of database size.

```
S' = (λ × S − ln K) / ln 2
```

Relationship to E-value:

```
E-value = m × n × 2^(−S')
```

Higher bit-score = better alignment. Comparable across databases.

---

### 1.3 KEGG Modified Identity (KAAS / BlastKOALA)

Penalizes partial alignments that would inflate identity by only aligning the best
matching portion of two proteins:

```
modified_identity = pident × min(1, (L_aln × 2) / (L_query + L_target))
```

| Symbol | Meaning |
|--------|---------|
| pident | Percentage identity of the alignment |
| L_aln | Alignment length |
| L_query | Query protein length |
| L_target | Database protein length |

**Example:**
- pident = 95%, L_aln = 50, L_query = 500, L_target = 490
- `min(1, (50×2)/(500+490)) = min(1, 100/990) = 0.101`
- modified_identity = 95% × 0.101 = 9.6% — correctly penalized for partial alignment

---

## 2. Genome Assembly Quality

### 2.1 N50

**Definition:** The length L such that contigs of length ≥ L collectively cover ≥ 50% of
the total assembly length.

**Algorithm:**

```
1. Sort contig lengths in descending order: L₁ ≥ L₂ ≥ … ≥ Lₙ
2. Compute total length: T = Σᵢ Lᵢ
3. Compute cumulative sum; stop when cumsum ≥ 0.5 × T
4. N50 = the length Lₖ at which the threshold was crossed
```

**Worked example:**

| Contig | Length | Cumulative sum |
|--------|--------|---------------|
| 1 | 5000 | 5000 |
| 2 | 3000 | 8000 ← 8000 ≥ 6250 (50% of 12500) |
| 3 | 2000 | 10000 |
| 4 | 1500 | 11500 |
| 5 | 1000 | 12500 |

Total = 12500 bp. 50% threshold = 6250 bp. **N50 = 3000 bp.**

---

## 3. Variant Calling

### 3.1 Phred Quality Score

**Definition:** Logarithmic transformation of the probability of an error.

```
Q = −10 × log₁₀(P_error)

Inverse: P_error = 10^(−Q/10)
```

| Q | P_error | Accuracy |
|---|---------|---------|
| 10 | 10⁻¹ = 0.1 | 90.0% |
| 20 | 10⁻² = 0.01 | 99.0% |
| 30 | 10⁻³ = 0.001 | 99.9% ✅ |
| 40 | 10⁻⁴ = 0.0001 | 99.99% |
| 60 | 10⁻⁶ = 0.000001 | 99.9999% |

**QUAL in VCF:** `QUAL = 30` means 0.1% probability the variant is a false positive.

---

### 3.2 Genotype Likelihood (PL)

**Definition:** Phred-scaled probability of observing the read data given a specific genotype.

```
PL(genotype) = −10 × log₁₀[ P(reads | genotype) ]
```

PL values are normalized so the **best (most probable) genotype has PL = 0**.

**Full example for a heterozygous call (REF=A, ALT=G):**

```
Reads: AAAAAAGGGGG  (5 A, 5 G; depth = 10)

P(reads | 0/0) = 0.5^10 × (error_rate)^5 ≈ 3.1e-11 → PL = 100
P(reads | 0/1) = 0.5^10 ≈ 9.8e-4          → PL = 30  ← best
P(reads | 1/1) = 0.5^10 × (error_rate)^5 ≈ 3.1e-11 → PL = 100

Normalized: PL = [70, 0, 70]   (subtract PL_best = 30 from all)
VCF writes: PL=70,0,70
```

VCF encodes PL in genotype order: `ref/ref, ref/alt, alt/alt`.

---

### 3.3 Genotype Quality (GQ)

**Definition:** Confidence in the assigned genotype — the difference between the best and
second-best PL values.

```
GQ = min(99, PL_second_best − PL_best)
```

**Example:**
- PL values: `[70, 0, 70]` — best = `0/1` (PL=0), second best = `0/0` or `1/1` (PL=70)
- `GQ = min(99, 70 − 0) = 70`
- Interpretation: the probability of the wrong genotype is 10^(−70/10) = 10⁻⁷ = 1 in 10,000,000

---

### 3.4 Allele Balance (AB)

**Definition:** Fraction of reads supporting the alternate allele.

```
AB = ALT_count / (REF_count + ALT_count) = AO / (RO + AO)
```

| Genotype | Expected AB |
|----------|------------|
| Homozygous REF (0/0) | ≈ 0.0 |
| Heterozygous (0/1) | 0.4 – 0.6 |
| Homozygous ALT (1/1) | ≈ 1.0 |

AB far from 0.5 in a called heterozygote → suspect strand bias or misalignment.

---

### 3.5 Strand Bias — Fisher Strand (FS)

Tests whether REF and ALT alleles are equally distributed between forward and reverse strands.

```
FS = −10 × log₁₀(p_Fisher)
```

where `p_Fisher` = p-value from Fisher's exact test on:

```
             Forward   Reverse
REF allele:    a          b
ALT allele:    c          d
```

**Example:**
```
             Forward   Reverse
REF allele:    30         35
ALT allele:    28          2    ← ALT mostly on forward strand
```
Fisher p = 0.0001 → FS = 40 → moderate strand bias ⚠️

| FS | p-value | Interpretation |
|----|---------|---------------|
| < 10 | > 0.1 | No bias ✅ |
| 10 – 30 | 0.001 – 0.1 | Mild |
| 30 – 60 | 1e⁻⁶ – 0.001 | Moderate ⚠️ |
| > 60 | < 1e⁻⁶ | Strong bias 🚨 |

---

### 3.6 Ti/Tv Ratio

**Transitions (Ti):** A↔G (purine↔purine) and C↔T (pyrimidine↔pyrimidine)
**Transversions (Tv):** A↔C, A↔T, G↔C, G↔T (purine↔pyrimidine)

```
Ti/Tv = N(transitions) / N(transversions)
```

| Sample | Expected Ti/Tv |
|--------|---------------|
| Human WGS | 2.0 – 2.2 |
| Human WES (exome) | 2.8 – 3.0 |
| Random mutations | ~0.5 (Ti/Tv for purely random base changes) |
| Too low (< 1.5) | Excess false positives |
| Too high (> 3.5) | Systematic bias or poor trimming |

---

### 3.7 Coverage (Depth)

**Definition:** Average number of reads covering each genome position.

```
Depth = (N_reads × L_read) / L_genome
```

**Example:**
- N_reads = 10,000,000 reads
- L_read = 150 bp
- L_genome = 3,000,000,000 bp (human genome)
- `Depth = (10⁷ × 150) / (3 × 10⁹) = 0.5×`

Only 0.5× — need more reads. For 30× coverage:
- `N_reads_needed = (30 × 3×10⁹) / 150 = 600,000,000 reads`

---

## 4. Phylogenetics

### 4.1 Neighbor-Joining: Q-matrix

**Purpose:** Correct for rate variation when choosing which taxa to join.
The pair with the lowest Q value is joined next.

```
Q(i,j) = (n − 2) × d(i,j) − Σₖ d(i,k) − Σₖ d(j,k)
```

| Symbol | Meaning |
|--------|---------|
| n | Number of current taxa |
| d(i,j) | Distance between taxa i and j |
| Σₖ d(i,k) | Sum of all distances from taxon i |

**Worked example (n=4, taxa A,B,C,D):**

Distance matrix and row sums:
```
d(A,B)=5  d(A,C)=7  d(A,D)=10   ΣA=22
          d(B,C)=8  d(B,D)=11   ΣB=24
                    d(C,D)=6    ΣC=21
                                ΣD=27
```

Q values:
```
Q(A,B) = (4−2)×5 − 22 − 24 = 10 − 46 = −36   ← minimum → join A,B
Q(A,C) = (4−2)×7 − 22 − 21 = 14 − 43 = −29
Q(C,D) = (4−2)×6 − 21 − 27 = 12 − 48 = −36   (tie; choose either)
```

---

### 4.2 NJ: Branch Lengths When Joining

When taxa f and g are joined into new node u:

```
d(f,u) = ½·d(f,g) + [1/(2n−2)] × (Σₖ d(f,k) − Σₖ d(g,k))
d(g,u) = d(f,g) − d(f,u)
```

**Example (join A,B → u; n=4):**
```
d(A,u) = ½×5 + [1/(2×4−2)] × (22−24)
        = 2.5 + (1/6) × (−2)
        = 2.5 − 0.333
        = 2.167

d(B,u) = 5 − 2.167 = 2.833
```

---

### 4.3 NJ: Update Distance Matrix

After joining f and g into u, compute distance from u to every remaining taxon x:

```
d(u,x) = ½ × [d(f,x) + d(g,x) − d(f,g)]
```

**Example (distance from new node u to C):**
```
d(u,C) = ½ × [d(A,C) + d(B,C) − d(A,B)]
        = ½ × [7 + 8 − 5]
        = ½ × 10
        = 5.0
```

---

### 4.4 Maximum Likelihood

**Likelihood:** probability of observing the alignment given the tree and model.

```
L(T, θ | data) = P(data | T, θ) = ∏_sites P(column_i | T, θ)
```

In practice, use log-likelihood (sum of log-probabilities):

```
ln L = Σ_sites ln P(column_i | T, θ)
```

**Model selection by BIC:**
```
BIC = k × ln(n) − 2 × ln(L)
```

Lower BIC = better model. k = free parameters, n = alignment length (columns).

---

### 4.5 Bayesian Inference: Posterior Probability

```
P(Tree | Data) = P(Data | Tree) × P(Tree) / P(Data)
```

| Term | Name | Meaning |
|------|------|---------|
| P(Tree \| Data) | Posterior | Probability of this tree given the observed data |
| P(Data \| Tree) | Likelihood | Probability of the data given this tree |
| P(Tree) | Prior | Prior belief about the tree (usually uniform) |
| P(Data) | Marginal | Normalizing constant (usually intractable → use MCMC) |

---

### 4.6 Robinson-Foulds Distance

**Definition:** The number of clades (bipartitions) present in one tree but not the other.

```
RF(T₁, T₂) = |C(T₁) \ C(T₂)| + |C(T₂) \ C(T₁)|
```

**Normalized RF distance:**
```
RF_norm = RF / (2n − 6)     (for bifurcating trees with n leaves)
```

| RF_norm | Meaning |
|---------|---------|
| 0 | Trees are identical |
| 1 | Trees are completely different (no shared clades) |

---

## 5. Quick-Reference Thresholds

| Metric | Threshold | Meaning |
|--------|-----------|---------|
| E-value | ≤ 1e-5 | Significant homology |
| QUAL | ≥ 30 | Good variant quality |
| GQ | ≥ 20 | Reliable genotype |
| DP | ≥ 10 | Minimum read depth |
| MAPQ | ≥ 30 | Read mapped uniquely |
| AB (het) | 0.40 – 0.60 | Expected allele balance |
| FS | < 60 | Acceptable strand bias |
| Ti/Tv | 2.0 – 2.2 | Expected for human WGS |
| Bootstrap | ≥ 70% | Reliable phylogenetic clade |
| N50 bacteria | ≥ 100 kb | Good short-read assembly |

---

## 🔗 Related Theory Documents

- [`docs/01_genome_annotation.md`](01_genome_annotation.md) — E-value, DIAMOND, N50
- [`docs/02_variant_calling.md`](02_variant_calling.md) — QUAL, PL, GQ, FS, Ti/Tv
- [`docs/03_phylogenetics.md`](03_phylogenetics.md) — NJ, ML, bootstrap, RF distance
- [`docs/04_variant_interpretation.md`](04_variant_interpretation.md) — ACMG scoring
