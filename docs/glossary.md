# Glossary of Bioinformatics Terms

> Plain-English definitions for terms used across this repository.
> Alphabetically sorted. Click any term in other docs to find it here.

---

## A

**Adapter** — Short synthetic DNA sequences ligated to both ends of every library fragment
during NGS library preparation. Adapters allow fragments to attach to the flow cell and
provide primer binding sites for sequencing. Must be removed (trimmed) before alignment
because they are not part of the biological sequence.

**Allele** — One of two or more alternative versions of a DNA sequence at a given position.
For example, at position 12345, 60% of chromosomes might carry `A` (the reference allele)
and 40% might carry `G` (an alternate allele).

**Allele Balance (AB)** — The fraction of reads supporting the alternate allele:
`AB = ALT_reads / (REF_reads + ALT_reads)`. For a heterozygous diploid variant, AB ≈ 0.5.

**Alignment** — (1) In sequence analysis: the process of inserting gap characters into
sequences so that homologous positions line up in the same column (multiple sequence
alignment, MSA). (2) In read mapping: matching each sequencing read to its position in the
reference genome.

**ALT** — Alternate allele — the non-reference nucleotide(s) observed at a variant site.
Written in the ALT column of a VCF file.

**Assembly** — The computational reconstruction of a genome sequence from many short
sequencing reads, using overlap or de Bruijn graph methods. The output is a set of
contigs or scaffolds.

---

## B

**BAM** — Binary Alignment Map. The compressed binary version of SAM format. Stores the
alignment of every sequencing read to the reference genome. Requires a `.bai` index file
for random access.

**Base calling** — The process of converting raw fluorescence signal from a sequencer into
nucleotide sequence (A/C/G/T) and quality scores. Performed by the sequencer software.

**bcftools** — A suite of command-line tools for manipulating VCF/BCF files: filtering,
merging, normalizing, annotating, and computing statistics.

**BIC** — Bayesian Information Criterion. A model selection criterion that penalizes model
complexity: `BIC = k × ln(n) − 2 × ln(L)`. Lower BIC is better. Used by ModelFinder
(IQ-TREE) to choose the best substitution model.

**Bit-score** — A database-size-normalized alignment score, making alignments comparable
across searches against databases of different sizes. Higher is better.
`S' = (λ × S − ln K) / ln 2`.

**Bootstrap** — A resampling method to estimate confidence in phylogenetic clades. Alignment
columns are sampled with replacement to create pseudo-replicate datasets; a new tree is built
for each. Bootstrap support = % of trees that contain a given clade.

**BWA** — Burrows-Wheeler Aligner. A widely used tool for aligning short reads to a reference
genome. Uses the Burrows-Wheeler Transform to enable fast exact string matching.

---

## C

**CADD** — Combined Annotation Dependent Depletion. A computational score that integrates
dozens of annotations to predict the deleteriousness of a variant. Higher CADD score = more
likely damaging.

**CDS** — Coding Sequence. The portion of an mRNA (and by extension, gene) that encodes the
protein, from start codon to stop codon. Does not include UTRs or introns.

**CIGAR string** — Compact Idiosyncratic Gapped Alignment Report. Encodes how a read aligns
to the reference in SAM/BAM format. `150M` = 150 matches; `10M2I138M` = 10 matches, 2-base
insertion, 138 matches.

**ClinVar** — NCBI database of human genetic variants and their clinical interpretations.
Contains pathogenicity classifications submitted by labs and researchers, along with
supporting evidence.

**CNV** — Copy Number Variant. A structural variant in which a region of the genome is
duplicated or deleted, changing the number of copies. Can affect gene expression.

**Consensus sequence** — A sequence representing the "most common" base at each position
across a set of aligned sequences or reads.

**Contig** — A contiguous DNA sequence assembled from overlapping reads. Unlike scaffolds,
contigs contain no gaps.

**Coverage** — See **Depth**.

---

## D

**De Bruijn graph** — A directed graph where nodes represent k-mers (substrings of length k)
and edges represent k-mer overlaps. Used in genome assemblers (SPAdes, Velvet) to avoid
explicit pairwise read comparisons.

**Depth** (also coverage) — The average number of sequencing reads covering each position
in the genome. `Depth = (N_reads × read_length) / genome_length`. Higher depth = more
statistical confidence in variant calls.

**DIAMOND** — A fast sequence aligner that reimplements BLASTP/BLASTX with 500–20,000×
speedup using double indexing and SIMD vectorization.

**Duplication (PCR)** — Reads that arise from PCR amplification of the same DNA molecule
during library preparation. Duplicates have identical start and end positions and inflate
apparent depth without adding independent evidence.

---

## E

**E-value** — The expected number of alignments with a score ≥ S that would occur by chance
in a database of a given size. Lower E-value = more significant hit. Formula:
`E = K × m × n × e^(−λ × S)`.

**Exome** — The protein-coding portion of the genome (all exons), comprising ~1–2% of the
human genome but ~85% of disease-causing mutations. Whole-exome sequencing (WES) targets
this region with capture probes.

**Exon** — A segment of a gene that is present in the mature mRNA and translated into protein
(or retained in functional non-coding RNA). Contrast with **intron**.

---

## F

**FASTA** — A text format for nucleotide or protein sequences. Each record starts with
`>identifier description` on one line, followed by the sequence on subsequent lines.

**FASTQ** — A text format that combines a FASTA sequence with per-base Phred quality scores.
The standard output of most sequencers.

**Filter (VCF)** — The FILTER column in VCF records the quality status of each variant.
`PASS` = passed all filters; anything else = a tag explaining why it was flagged.

**FLAG (BAM)** — A bitwise integer encoding properties of a read alignment: paired, properly
paired, unmapped, mate unmapped, reverse strand, R1/R2, duplicate, etc.

**FreeBayes** — A haplotype-based variant caller that builds local haplotypes from reads and
uses Bayesian inference to assign genotype likelihoods.

**FS (Fisher Strand)** — A Phred-scaled p-value from Fisher's exact test for strand bias.
`FS = −10 × log₁₀(p)`. High FS (> 60) suggests the ALT allele appears predominantly on one
strand, often indicating a technical artifact.

---

## G

**GATK** — Genome Analysis Toolkit. A widely used software suite from the Broad Institute
for variant calling and genomic analysis. Includes HaplotypeCaller, Mutect2, and many other
tools.

**GFF / GFF3** — General Feature Format. A tab-separated text format for storing genomic
feature annotations (gene positions, exon boundaries, etc.) relative to a reference sequence.

**gnomAD** — Genome Aggregation Database. A large public database of genetic variants
observed in >125,000 exomes and >15,000 genomes from diverse populations. Used to
determine population allele frequencies.

**GQ (Genotype Quality)** — Phred-scaled confidence that the assigned genotype is correct.
`GQ = min(99, PL_second_best − PL_best)`. GQ = 20 means 1% probability of wrong genotype.

**GT (Genotype)** — VCF FORMAT field encoding the called genotype. `0/0` = homozygous
reference; `0/1` = heterozygous; `1/1` = homozygous alternate. `/` = unphased; `|` = phased.

**GTR model** — General Time Reversible model for nucleotide substitution. The most general
reversible model, with 6 free substitution rate parameters and 4 base frequency parameters.
Often extended: GTR+Γ (rate variation) or GTR+Γ+I (invariant sites).

---

## H

**Haplotype** — A specific combination of alleles along a chromosome that tend to be
inherited together. In variant calling, local haplotypes are the short sequences observed
in a window of reads.

**Hard filter** — A variant filter that permanently removes records below a threshold.
Contrast with **soft filter**, which marks but preserves them.

**HMM (Hidden Markov Model)** — A statistical model with hidden states and observed outputs.
In sequence analysis, used to model protein families (Pfam HMM profiles) where hidden states
correspond to conserved vs variable positions.

**HMMER** — A software suite for protein sequence analysis using profile HMMs. Includes
`hmmbuild` (build a profile from an alignment), `hmmscan` (search sequences against profiles),
and `hmmsearch` (search a profile against sequences).

**HGVS notation** — Human Genome Variation Society nomenclature for describing variants.
`c.2T>A` = coding sequence change at position 2 (T → A); `p.Met1Lys` = protein change
(methionine → lysine at position 1).

---

## I

**IGV** — Integrative Genomics Viewer. A desktop tool for visualizing BAM, VCF, and
annotation files. Essential for manual review of variant calls.

**IMPACT** — VEP impact category: HIGH / MODERATE / LOW / MODIFIER. Describes the expected
functional severity of a variant's consequence on the protein.

**Indel** — Insertion or deletion of bases. Indels < 50 bp are called small indels (handled
by standard variant callers). Indels > 50 bp are classified as structural variants.

**INFO field (VCF)** — The eighth column of a VCF record, containing semicolon-separated
key=value annotations for the variant site (e.g., `DP=100;AF=0.55;CSQ=...`).

**IQ-TREE** — A fast maximum-likelihood tree builder with automated model selection
(ModelFinder) and ultrafast bootstrap (UFBoot2).

---

## L

**Left-alignment** — The convention of shifting indels as far toward lower genomic coordinates
as possible while keeping the sequence equivalent. Ensures consistent representation of the
same variant from different callers.

**Library** — In NGS: the collection of DNA fragments prepared for sequencing, with adapters
attached. Each library corresponds to one biological sample.

---

## M

**MAFFT** — Multiple Alignment using Fast Fourier Transform. A widely used tool for multiple
sequence alignment of DNA or protein sequences.

**MAPQ** — Mapping Quality. A Phred-scaled score indicating the confidence that a read is
correctly aligned to its reported position. MAPQ = 0 means the read maps equally well to
multiple locations. Filter on MAPQ ≥ 30 for high-confidence alignments.

**markdup** — The process of identifying and tagging PCR duplicate reads in a BAM file
(`samtools markdup`). Tagged reads are typically excluded from variant calling.

**MSA** — Multiple Sequence Alignment. An alignment of three or more sequences so that
homologous positions are in the same column.

---

## N

**N50** — A statistic for genome assembly quality. The N50 length L is defined such that
contigs of length ≥ L cover 50% of the total assembly length. Higher N50 = better
contiguity.

**Neighbor-Joining (NJ)** — A distance-based phylogenetic tree-building algorithm. Fast
(O(n³)) and suitable for exploratory analysis. Corrects for rate variation using the
Q-matrix transformation.

**NGS** — Next-Generation Sequencing. All massively parallel sequencing technologies
(Illumina, Ion Torrent, PacBio, Oxford Nanopore) as opposed to first-generation Sanger
sequencing.

---

## O

**ORF** — Open Reading Frame. A stretch of DNA between a start codon (ATG) and a stop codon
(TAA, TAG, or TGA) in the same reading frame. ORFs are the primary candidates for protein-
coding genes, though not all ORFs encode real proteins.

**Orthologues** — Genes in different species that descended from the same gene in the common
ancestor (separated by speciation). Used in phylogenetic analysis to trace species history.

---

## P

**Paired-end sequencing** — An NGS library design where both ends of each DNA fragment are
sequenced, producing two reads (R1 and R2) per fragment. Provides better alignment accuracy
and enables structural variant detection.

**Paralogues** — Genes within the same genome that arose by gene duplication. They share a
common ancestral gene but diverged after duplication. Paralogues should not be used to infer
species relationships.

**Pfam** — A database of protein family alignments and HMM profiles, maintained by EMBL-EBI.
Each Pfam entry covers one conserved domain or protein family.

**Phred score** — See **Q score**.

**PL (Phred-scaled Likelihood)** — VCF FORMAT field. `PL = −10 × log₁₀(L)` for each
possible genotype. Normalized so the most likely genotype has PL = 0.

**Prokka** — A command-line tool for rapid prokaryote genome annotation. Integrates Prodigal
(gene prediction), Barrnap (rRNA), Aragorn (tRNA), and BLAST/DIAMOND (function).

**Prodigal** — PROkaryotic DYnamic programming Gene-finding ALogithm. Predicts protein-coding
genes in bacterial and archaeal genomes using statistical models trained on the genome itself.

---

## Q

**Q score** (Phred score, quality score) — A logarithmic measure of the probability of an
incorrect base call: `Q = −10 × log₁₀(P_error)`. Stored in FASTQ files as ASCII characters
(Phred+33 encoding).

**QUAL (VCF)** — Phred-scaled probability that the variant call is incorrect.
`QUAL = −10 × log₁₀(P_variant_is_wrong)`. QUAL = 30 means 0.1% chance of a false call.

---

## R

**Read** — A single short DNA sequence produced by an NGS instrument, typically 50–300 bp
for Illumina. Each read also comes with a quality string.

**Read group** — A tag (`@RG`) in SAM/BAM format that identifies which sequencing run, lane,
library, and sample a read belongs to. Required by GATK and used for tracking batch effects.

**Reference allele** — The nucleotide(s) at a given position in the reference genome assembly.
Written in the REF column of a VCF file.

**Reference genome** — A representative, consensus sequence assembly used as a coordinate
system for aligning reads and calling variants.

**RF distance** — Robinson-Foulds distance. The number of clades present in one tree but not
the other, summed for both trees. `RF = 0` means identical topologies.

---

## S

**SAM** — Sequence Alignment Map. A tab-delimited text format for storing read alignments.
The human-readable counterpart of BAM.

**scaffold** — A set of contigs joined with estimated gap sizes (indicated by `N` characters).
Scaffolding uses paired-end or long-read information to orient and order contigs.

**SNV / SNP** — Single Nucleotide Variant (or Polymorphism). A single-base change between
two DNA sequences. SNP traditionally implies population-level polymorphism; SNV is more
general (includes rare and somatic variants).

**Soft filter** — A variant filter that adds a tag to the FILTER column without removing the
record. Allows downstream tools to apply different thresholds without re-calling variants.

**Strand bias** — An artifact where the alternate allele is observed predominantly on reads
from one strand (forward or reverse). Usually indicates a technical error rather than a real
variant.

**Structural variant (SV)** — A genomic alteration > 50 bp: deletions, insertions, inversions,
duplications, or translocations. Detected from split reads, discordant read pairs, or coverage
changes.

**Swiss-Prot** — The manually reviewed (high-quality) section of UniProt. Contains well-
characterized proteins with curated functional annotations.

---

## T

**Ti/Tv ratio** — Transition/transversion ratio. Transitions (A↔G, C↔T) occur more frequently
than transversions (A↔C, A↔T, G↔C, G↔T) in most genomes. For human WGS, expected Ti/Tv ≈ 2.0–2.2.
A ratio much lower than expected indicates excess false-positive variant calls.

**Transversion** — A nucleotide substitution that changes a purine (A, G) to a pyrimidine
(C, T) or vice versa. There are 8 possible transversions (e.g., A→C, A→T, G→C, G→T and
their reverses).

**Transition** — A nucleotide substitution within the same chemical class (purine↔purine or
pyrimidine↔pyrimidine). There are 4 transitions: A↔G and C↔T.

---

## U

**UFBoot2** — Ultrafast Bootstrap Approximation version 2 (IQ-TREE). A computationally
efficient bootstrap method that reuses intermediate computations. Values ≥ 95 are
considered strongly supported.

**UTR** — Untranslated Region. The portions of an mRNA upstream (5' UTR) and downstream
(3' UTR) of the coding sequence (CDS). Important for regulating mRNA stability and
translation efficiency.

---

## V

**VCF** — Variant Call Format. The standard text format for storing genetic variants. Each
record describes one variant site with chromosome, position, reference allele, alternate
allele, quality, filter status, and per-sample genotype information.

**VEP** — Variant Effect Predictor (Ensembl). A tool that annotates each variant with its
predicted molecular consequence (missense, frameshift, etc.), gene name, HGVS notation,
population allele frequencies (gnomAD), and clinical significance (ClinVar).

**VUS** — Variant of Uncertain Significance. An ACMG classification for variants where
there is insufficient evidence to classify as benign or pathogenic.

---

## W

**WES** — Whole Exome Sequencing. Sequencing of only the protein-coding regions of the
genome (~1–2% of total, ~50 Mb). More cost-effective than WGS for finding coding variants.

**WGS** — Whole Genome Sequencing. Sequencing of the entire genome, including non-coding
regions. Higher cost than WES but captures all variant types including regulatory and
structural variants.

---

## Index by Category

**Formats:** FASTA · FASTQ · SAM · BAM · VCF · GFF3 · Newick · PHYLIP · NEXUS

**Quality metrics:** Q score · QUAL · GQ · PL · MAPQ · DP · AB · FS · Ti/Tv · N50

**Variant types:** SNV · Indel · MNV · SV · CNV

**Alignment tools:** BWA · MAFFT · Needleman-Wunsch · Smith-Waterman

**Variant callers:** FreeBayes · GATK HaplotypeCaller · Mutect2

**Annotation tools:** Prokka · DIAMOND · HMMER · VEP

**Phylogeny tools:** IQ-TREE · FastTree · MAFFT · ClipKit · gotree · MrBayes

**Databases:** gnomAD · ClinVar · Pfam · Swiss-Prot · KEGG · InterPro
