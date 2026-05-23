# 3. Phylogenetics

## What Is Phylogenetics?

Phylogenetics is the study of evolutionary relationships among organisms or genes, inferred
from molecular sequence data (DNA, RNA, or protein). The result is a **phylogenetic tree**:
a branching diagram where tips (leaves) are the sequences you gave it, and internal nodes
represent their inferred common ancestors.

Phylogenetics answers questions like:
- Are two bacterial strains in the same outbreak?
- Did this gene family arise by duplication or speciation?
- When did two species last share a common ancestor?
- Which SARS-CoV-2 lineage is circulating in my city?

---

## 3.1 Key Concepts

### Homology, Analogy, and Paralogy

These are frequently confused terms:

| Term | Definition | Example |
|------|------------|---------|
| **Homology** | Similarity due to shared ancestry | Human and mouse Hox genes |
| **Analogy** | Similarity due to convergent evolution, not ancestry | Wings of birds and insects |
| **Orthologues** | Homologues separated by a **speciation** event | Human BRCA1 vs mouse Brca1 |
| **Paralogues** | Homologues separated by a **duplication** event | Human α-globin vs β-globin |
| **Homoplasy** | Same character state evolved independently | White fur in polar bears and arctic foxes |

> **Why it matters for analysis:** When building a gene tree, you want **orthologues**
> (same gene in different species) to reflect species history. **Paralogues** reflect
> gene duplication history, not speciation — mixing them produces incorrect trees.

### Tree components

```
        0.05   ┌──── Species A    ← leaf (tip, taxon)
       ┌───────┤
       │       └──── Species B
───────┤              branch length = evolutionary distance
       │       ┌──── Species C
       └───────┤
               └──── Species D

       ↑                ↑
     root           internal node
   (common ancestor)  (divergence event)
```

| Component | Description |
|-----------|-------------|
| **Leaf (tip)** | The input sequences — one per taxon |
| **Internal node** | Hypothetical common ancestor of the descending clade |
| **Branch** | Lineage connecting two nodes |
| **Branch length** | Evolutionary distance (substitutions per site) |
| **Root** | Most ancient common ancestor of all leaves |
| **Clade** | A node + all its descendants (monophyletic group) |
| **Topology** | The branching pattern (ignoring branch lengths) |

---

## 3.2 Multiple Sequence Alignment (MSA)

Before building a tree, sequences must be **aligned** — columns must contain homologous
positions (bases that trace back to the same ancestral base).

### The alignment problem

Given two sequences:
```
Seq 1: ACGTACGT
Seq 2: ACGACGT
```

Are these aligned as:
```
Option A:                 Option B:
Seq 1: ACGTACGT           Seq 1: ACG-TACGT
Seq 2: ACG-ACGT           Seq 2: ACGA-CGT
       ↑ gap                      ↑ gap in different position
```

Option A implies a 1-base deletion in Seq 2; Option B implies a different deletion.
The correct alignment requires understanding the evolutionary history — which we don't
have directly. Alignment is itself an inference.

### Pairwise alignment algorithms

**Needleman-Wunsch (1970):** Global alignment — finds the best alignment of the
entire sequences. Uses dynamic programming to fill a scoring matrix:

```
Score(i,j) = max(
    Score(i-1, j-1) + substitution_score(a[i], b[j]),  # match/mismatch
    Score(i-1, j)   + gap_penalty,                       # gap in sequence b
    Score(i,   j-1) + gap_penalty                        # gap in sequence a
)
```

**Smith-Waterman (1981):** Local alignment — finds the best alignment of any
subsequences. Same algorithm but resets negative scores to zero, allowing the
alignment to "start fresh" in regions of similarity.

### Multiple sequence alignment with MAFFT

Direct extension of pairwise DP to N sequences would require O(L^N) time — impossible
for more than a few sequences. MAFFT (Katoh et al. 2002) uses:

1. **Fast Fourier Transform** to find regions of similarity (seeds) quickly
2. **Progressive alignment** — align the most similar pair first, then add more
   sequences incrementally (guided by an initial distance tree)
3. **Iterative refinement** — repeatedly realign subsets to improve the overall score

```bash
# Auto mode: MAFFT chooses strategy based on data size
mafft --auto --thread 8 input.fasta > aligned.fasta

# High accuracy for < 200 proteins:
mafft --localpair --maxiterate 1000 --thread 8 proteins.fasta > aligned.fasta

# Fast mode for > 1000 sequences:
mafft --retree 2 --thread 8 large.fasta > aligned.fasta
```

**What to check in the alignment:**
- Are known conserved positions (e.g., catalytic residues) in the same column?
- Are there excessively gappy columns (> 50% gaps)?
- Do sequences with very different lengths align sensibly?

---

## 3.3 Alignment Trimming with ClipKit

Not all alignment columns contain phylogenetic signal:

| Column type | Phylogenetic value | Action |
|-------------|-------------------|--------|
| **Parsimony-informative** | Useful: supports at least one tree topology over another | Keep |
| **Constant** | Neutral: all sequences identical | Keep (needed for model parameters) |
| **Gappy** | Noise: mostly gaps, unreliable homology | Remove |
| **Singleton** | Uninformative: only one sequence differs | Remove |

```bash
# kpic = keep parsimony-informative and constant sites (recommended)
clipkit aligned.fasta -m kpic -o trimmed.fasta

# Check what was removed
echo "Before: $(grep -v '>' aligned.fasta | head -1 | wc -c) columns"
echo "After:  $(grep -v '>' trimmed.fasta | head -1 | wc -c) columns"
```

---

## 3.4 Tree Building Methods

### Method 1 — Distance-based (Neighbor-Joining)

**Concept:** Convert the alignment into a pairwise distance matrix, then build a tree
that best explains those distances.

**Neighbor-Joining algorithm (Saitou & Nei 1987):**

The key insight: simply joining the two closest taxa is wrong because it doesn't account
for the fact that some lineages evolve faster than others. NJ corrects this with the
**Q-matrix:**

```
Q(i,j) = (n−2) × d(i,j) − Σₖ d(i,k) − Σₖ d(j,k)
```

The pair with the **lowest Q** is joined — this corrects for "long branch attraction"
(fast-evolving sequences falsely appearing close together).

**Branch length after joining taxa f and g into node u:**

```
d(f,u) = ½·d(f,g) + [1/(2n−2)] × [Σₖ d(f,k) − Σₖ d(g,k)]
d(g,u) = d(f,g) − d(f,u)
```

**Update distance matrix for new node u:**

```
d(u,x) = ½ × [d(f,x) + d(g,x) − d(f,g)]   for all remaining taxa x
```

**Full worked example (4 taxa, A B C D):**

Distance matrix:
```
      A    B    C    D
A     —    5    7   10
B     5    —    8   11
C     7    8    —    6
D    10   11    6    —

Row sums: A=22, B=24, C=21, D=27
```

Q-matrix (n=4):
```
Q(A,B) = (4−2)×5 − 22 − 24 = 10 − 46 = −36   ← minimum → join A and B
Q(A,C) = (4−2)×7 − 22 − 21 = 14 − 43 = −29
Q(A,D) = (4−2)×10 − 22 − 27 = 20 − 49 = −29
...
```

A and B are joined into node u:
```
d(A,u) = ½×5 + [1/6]×(22−24) = 2.5 − 0.333 = 2.167
d(B,u) = 5 − 2.167 = 2.833

d(u,C) = ½×(7+8−5) = 5.0
d(u,D) = ½×(10+11−5) = 8.0
```

Continue with (u, C, D) until only 3 nodes remain.

**Pros:** Very fast (O(n³)), produces reasonable trees.
**Cons:** Sensitive to rate variation across lineages; no statistical model.

---

### Method 2 — Maximum Likelihood (ML)

**Concept:** Find the tree topology and branch lengths that make the observed alignment
data most probable under an explicit evolutionary model.

```
ML tree = argmax_T  P(alignment | T, model)
```

**Likelihood calculation:** For each column of the alignment, compute the probability
of observing that pattern by summing over all possible ancestral states at internal nodes.
The total likelihood is the product over all columns (or sum of log-likelihoods).

**Substitution models for DNA:**

| Model | Parameters | Description |
|-------|-----------|-------------|
| **JC69** | 1 | All substitutions equally likely |
| **K80** | 2 | Transitions ≠ transversions |
| **HKY85** | 4 | Unequal base frequencies + Ti ≠ Tv |
| **GTR** | 9 | All 6 substitution rates free + base frequencies |
| **GTR+Γ** | 10 | GTR + gamma-distributed rate variation across sites |
| **GTR+Γ+I** | 11 | GTR+Γ + proportion of invariant sites |

**Rate variation (Γ model):** In real data, some alignment columns evolve fast (exposed
surface residues, non-essential regions) and others evolve slowly (active site residues,
structural cores). The gamma distribution models this:

```
Rate variation across sites ~ Γ(α, α)
α < 1: highly variable rates (most sites slow, some very fast)
α > 1: relatively uniform rates
α → ∞: all sites evolve at the same rate
```

**ModelFinder (IQ-TREE -m MFP):** Tests dozens of models and selects the best by
Bayesian Information Criterion (BIC):

```
BIC = k × ln(n) − 2 × ln(L)
```
where k = number of free parameters, n = alignment length, L = likelihood.
Lower BIC = better model (rewards fit but penalizes complexity).

---

### Method 3 — Bayesian Inference

**Concept:** Instead of finding one best tree, sample the posterior distribution of trees
using MCMC (Markov Chain Monte Carlo). The result is a probability distribution over
tree topologies and branch lengths.

```
P(tree | data) = P(data | tree) × P(tree) / P(data)
```

**MCMC sampling:** Start with a random tree, propose small changes (nearest-neighbor
interchange, branch length modification), accept changes with probability:

```
acceptance_ratio = min(1, P(data | new_tree) × P(new_tree) / P(data | old_tree) × P(old_tree))
```

After many iterations (typically millions), the sampled trees form the posterior distribution.
**Posterior probability of a clade** = fraction of sampled trees containing that clade.

**Programs:** MrBayes, BEAST2 (adds molecular clock for dating), PhyloBayes (complex models).

---

## 3.5 Bootstrap: Measuring Tree Reliability

### The Problem

A phylogenetic tree will always produce an answer — but how reliable is each branch?

### Bootstrap Algorithm

1. Generate B **bootstrap replicates** of the alignment by randomly sampling columns
   with replacement (same length, different column composition)
2. Build a tree for each replicate
3. **Bootstrap support** = percentage of B trees in which a given clade appears

```
Original:  ACGT ACGT TTTT GGGG CCCC
Bootstrap: ACGT CCCC ACGT TTTT ACGT   ← same columns, different sample
```

### Interpretation

| Support | Meaning | Recommendation |
|---------|---------|---------------|
| ≥ 95% | Very strongly supported | Report with high confidence |
| 70–94% | Well supported | Generally reliable |
| 50–69% | Weakly supported | Treat with caution |
| < 50% | Unreliable | Branch may collapse |

**Standard in publications:** ≥ 70% bootstrap is considered reliable.

### IQ-TREE Ultrafast Bootstrap (UFBoot2)

UFBoot2 (Hoang et al. 2018) achieves similar results to standard bootstrap but 10–40×
faster, by reusing computations across replicates. Note: UFBoot values are not directly
comparable to standard bootstrap values — UFBoot ≥ 95 ≈ standard ≥ 70.

```bash
iqtree -s alignment.fasta -m MFP -B 1000 -T AUTO --seqtype DNA
```

---

## 3.6 Tree Rooting

An unrooted tree shows relationships but not direction (which end is ancestral).
Rooting places the common ancestor:

### Outgroup rooting

Include a sequence **known to be outside** the group of interest. The root falls on
the branch connecting the outgroup to all ingroup sequences.

**Choosing a good outgroup:**
- Similar enough to the ingroup to align well
- Diverged enough to be clearly external
- Exactly one outgroup taxon is ideal (multiple outgroups can introduce artifacts)

### Midpoint rooting

Place the root at the midpoint of the longest path between any two leaves.
Assumption: all lineages evolve at the same rate (molecular clock).

```bash
gotree reroot midpoint -i tree.nwk -o rooted_tree.nwk
```

---

## 3.7 Common Artifacts

### Long Branch Attraction (LBA)

**Problem:** Sequences that evolved very fast (long branches) accumulate so many
substitutions that they saturate — all four bases become equally likely. Two
unrelated fast-evolving sequences then appear similar by chance, and distance-based
or parsimony methods incorrectly group them together.

```
True tree:        LBA artifact:
    A                 A
   / \               / \
  B   C             B   D    ← D wrongly attracted to B
  |   |             |   |
  D   E             C   E
(D evolved fast)  (D and B both fast)
```

**Solutions:**
- Use models that account for rate variation (+Γ)
- Add more taxa to "break up" long branches
- Remove fastest-evolving sites
- Use Bayesian methods (more robust than parsimony)

### Horizontal Gene Transfer (HGT)

In bacteria, genes can be transferred between distantly related organisms without
reproduction. A gene tree for an HGT-affected gene will conflict with the species tree.

**Signs of HGT:**
- Unexpected GC content or codon usage (different from the rest of the genome)
- Gene tree topology strongly contradicts a well-established species tree
- Gene is flanked by mobile elements (transposons, integrons)

---

## 3.8 Comparing Trees: Robinson-Foulds Distance

To quantify how similar two trees are, use the **Robinson-Foulds (RF) distance** (1981):

```
RF(T₁, T₂) = |C(T₁) \ C(T₂)| + |C(T₂) \ C(T₁)|
```

where C(T) = set of all non-trivial clades (bipartitions) in tree T.
`\` = set difference (clades present in one tree but not the other).

**Normalized RF distance** (0 = identical, 1 = completely different):

```
RF_norm = RF / (2n − 6)     for bifurcating trees with n leaves
```

| RF_norm | Interpretation |
|---------|---------------|
| 0 | Trees are identical |
| 0–0.1 | Very similar |
| 0.1–0.3 | Minor topological differences |
| > 0.5 | Major disagreements |

```bash
# Compare two trees with gotree
gotree compare trees -i tree1.nwk -c tree2.nwk
```

---

## 3.9 File Formats

| Format | Extension | Use |
|--------|-----------|-----|
| **FASTA** | `.fasta`, `.fa` | Input sequences |
| **Aligned FASTA** | `.fasta` | MSA output (same as FASTA, just aligned) |
| **PHYLIP** | `.phy` | MSA; required by PhyloBayes, RAxML |
| **NEXUS** | `.nex` | Alignment + tree + metadata; used by MrBayes, BEAST |
| **Newick** | `.nwk`, `.tre`, `.treefile` | Tree topology + branch lengths |

**Newick format example:**

```
((TaxonA:0.12, TaxonB:0.09):0.05, (TaxonC:0.21, TaxonD:0.18):0.08);
```

Read as: A and B share an ancestor (branch lengths 0.12 and 0.09);
C and D share an ancestor; those two clades share a root.

---

## 3.10 Quick Command Reference

```bash
# 1. Align with MAFFT
mafft --auto --thread 8 sequences.fasta > aligned.fasta

# 2. Trim with ClipKit
clipkit aligned.fasta -m kpic -o trimmed.fasta

# 3. Quick tree with FastTree
fasttree -nt -gtr -gamma < trimmed.fasta > quick_tree.nwk

# 4. Full ML tree with IQ-TREE (model selection + bootstrap)
iqtree -s trimmed.fasta -m MFP -B 1000 -T AUTO --seqtype DNA -pre ml_tree

# 5. Check bootstrap values in the tree file
grep -o ":[0-9.]*)" ml_tree.treefile | head

# 6. Midpoint root with gotree
gotree reroot midpoint -i ml_tree.treefile -o rooted.nwk

# 7. Visualize (iTOL web, FigTree desktop, or ggtree in R)
# Upload rooted.nwk to https://itol.embl.de/
```

---

## Further Reading

- Saitou N. & Nei M. (1987) "The neighbor-joining method" — Molecular Biology and Evolution
- Stamatakis A. (2014) "RAxML version 8: a tool for phylogenetic analysis" — Bioinformatics
- Nguyen L-T et al. (2015) "IQ-TREE: a fast and effective stochastic algorithm for estimating ML phylogenies" — Molecular Biology and Evolution
- Hoang D.T. et al. (2018) "UFBoot2: Improving Ultrafast Bootstrap Approximation" — Molecular Biology and Evolution
- Katoh K. & Standley D.M. (2013) "MAFFT Multiple Sequence Alignment Software" — Molecular Biology and Evolution

---

## 🔗 Related Files

- [`pipelines/phylogenetics/run_phylogenetics.sh`](../pipelines/phylogenetics/run_phylogenetics.sh)
- [`pipelines/phylogenetics/config.sh`](../pipelines/phylogenetics/config.sh)
- [`docs/formulas_reference.md`](formulas_reference.md) — NJ formulas, RF distance
