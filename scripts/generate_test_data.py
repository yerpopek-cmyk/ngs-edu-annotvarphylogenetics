#!/usr/bin/env python3
"""
scripts/generate_test_data.py — Generate minimal toy datasets for pipeline testing
===================================================================================

WHAT THIS GENERATES
───────────────────
This script creates tiny synthetic datasets so you can verify that every
pipeline step runs correctly without downloading gigabytes of real data.

Datasets created:

  1. phylogenetics/test_16S.fasta  (20 sequences × 500 bp)
     Simulates 16S rRNA gene fragments with realistic variation.
     Use with:  bash pipelines/phylogenetics/run_phylogenetics.sh -i test_16S.fasta

  2. annotation/test_genome.fasta  (1 contig × 10 000 bp)
     Synthetic bacterial genome fragment.
     Use with:  bash pipelines/annotation/run_annotation.sh --offline

  3. variants/test_ref.fasta       (1 chromosome × 50 000 bp)
     variants/test_R1.fastq.gz     (5 000 reads × 150 bp)
     variants/test_R2.fastq.gz     (paired-end reads)
     Use with:  bash pipelines/variants/run_all.sh (after editing config.sh)

HOW TO RUN
──────────
  python3 scripts/generate_test_data.py --outdir test_data/
  python3 scripts/generate_test_data.py --outdir test_data/ --seed 42

WHY SYNTHETIC DATA?
───────────────────
Real sequencing data is large (typically gigabytes per sample) and subject
to privacy restrictions when it contains human variants. Synthetic data lets
you test the pipeline logic quickly — a full run should take < 5 minutes.
"""

from __future__ import annotations

import argparse
import gzip
import random
import sys
from pathlib import Path


# =============================================================================
# CLI
# =============================================================================

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--outdir", default="test_data",
                   help="Output directory for generated files (default: test_data/)")
    p.add_argument("--seed", type=int, default=2025,
                   help="Random seed for reproducibility (default: 2025)")
    return p.parse_args()


# =============================================================================
# Sequence generation helpers
# =============================================================================

DNA_BASES = "ACGT"
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

def random_dna(length: int, rng: random.Random) -> str:
    """Generate a random DNA sequence of the given length."""
    return "".join(rng.choices(DNA_BASES, k=length))


def mutate_sequence(seq: str, mutation_rate: float, rng: random.Random) -> str:
    """
    Introduce point mutations at the given rate.
    This simulates the evolutionary divergence between sequences.
    A typical 16S rRNA divergence between different bacterial species
    is 1–10% (mutation_rate = 0.01–0.10).
    """
    bases = list(seq)
    for i, base in enumerate(bases):
        if rng.random() < mutation_rate:
            # Substitute with a different base (never the same)
            alternatives = [b for b in DNA_BASES if b != base]
            bases[i] = rng.choice(alternatives)
    return "".join(bases)


def make_phred_quality(length: int, rng: random.Random, mean_qual: int = 35) -> str:
    """
    Generate a realistic Phred quality string for a read.
    Quality is higher in the middle and lower at the ends (typical Illumina pattern).

    Phred score Q = -10 × log₁₀(P_error)
    Q=30 → P_error = 0.001  (1 error per 1000 bases)
    Q=20 → P_error = 0.01   (1 error per 100 bases)
    """
    quals = []
    for i in range(length):
        # Simulate quality drop at read ends (common in Illumina sequencing)
        position_factor = min(i, length - i - 1) / (length // 4)
        position_factor = min(1.0, position_factor)
        q = int(mean_qual * position_factor + rng.gauss(0, 3))
        q = max(2, min(40, q))   # clamp to [2, 40]
        quals.append(chr(q + 33))  # Phred+33 encoding (Illumina 1.8+)
    return "".join(quals)


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    """Write a list of (header, sequence) tuples as FASTA."""
    with open(path, "w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n")
            # Wrap sequence at 70 characters (standard FASTA)
            for i in range(0, len(seq), 70):
                fh.write(seq[i:i+70] + "\n")


def write_fastq_gz(path: Path, reads: list[tuple[str, str, str]]) -> None:
    """Write a list of (name, sequence, quality) tuples as gzipped FASTQ."""
    with gzip.open(path, "wt") as fh:
        for name, seq, qual in reads:
            fh.write(f"@{name}\n{seq}\n+\n{qual}\n")


# =============================================================================
# Dataset 1: 16S rRNA sequences for phylogenetics
# =============================================================================

# Rough conserved regions of a generic 16S rRNA (simplified).
# In reality, 16S has variable regions (V1–V9) flanked by conserved regions.
# Here we use a pseudo-conserved 500 bp sequence and introduce variation.
CONSERVED_16S_CORE = (
    "AGAGTTTGATCCTGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAGCGGTAGC"
    "ACAGAGAGCTTGCTCTCGGGTGACGAGCGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATG"
    "GAGGGGGATAACTACTGGAAACGGTAGCTAATACCGCATAACGTCGCAAGACCAAAGAGGGGGACCT"
    "TAGGGCCTTTGGGCCATCAGCTTTTGGTGGGTACGCACCGGAAATGTTGCATCAGCTTCGAGCGGGG"
    "GCAGCCGGCCATTACCGCGGCTGCTGGCACGGAGTTTAGCCCAGGGAGTCGCGCCTGTCCGAGAAG"
    "AATAACTACAAGCCCGGAGGTCATCGGAGCGGGATGTTATTAGAGCAGCAGCAGCAAGGAAAGAATA"
    "TACACCGAATCGAGAGTTCTTCGGAATCAGCTTGATCCTTCGGAGCGCATAGTCGGCGGCAGCAGGG"
    "ATTACCCGCGG"
)

def generate_16S_dataset(n_taxa: int, rng: random.Random) -> list[tuple[str, str]]:
    """
    Generate N divergent 16S rRNA sequences.
    Each sequence is a mutated version of the conserved core with variable
    insertion lengths to simulate real 16S variation.
    """
    records = []
    # Represent different "taxa" at varying divergence levels
    groups = ["Firmicutes", "Proteobacteria", "Actinobacteria", "Bacteroidetes"]

    for i in range(n_taxa):
        group = groups[i % len(groups)]
        species_num = i // len(groups) + 1
        taxon_name = f"{group}_sp{species_num:02d}"

        # Each group diverges from the core at a group-specific rate;
        # within-group variation is lower than between-group variation
        group_rate = 0.05 * (i % len(groups))    # inter-group divergence
        species_rate = rng.uniform(0.01, 0.03)   # intra-group variation

        seq = mutate_sequence(CONSERVED_16S_CORE, group_rate + species_rate, rng)

        # Append a short variable region (simulates V3-V4 hypervariable region)
        variable_len = rng.randint(20, 60)
        variable_region = random_dna(variable_len, rng)
        seq = seq[:300] + variable_region + seq[300:]

        records.append((taxon_name, seq[:500]))  # trim to 500 bp

    return records


# =============================================================================
# Dataset 2: Synthetic bacterial genome for annotation
# =============================================================================

# Simplified codon table (start codons: ATG; stop codons: TAA, TAG, TGA)
START_CODON = "ATG"
STOP_CODONS = ["TAA", "TAG", "TGA"]
CODING_CODONS = [
    "GCT", "GCC", "GCA", "GCG",   # Ala
    "CGT", "CGC", "CGA", "CGG",   # Arg
    "AAT", "AAC",                  # Asn
    "GAT", "GAC",                  # Asp
    "TGT", "TGC",                  # Cys
    "CAA", "CAG",                  # Gln
    "GAA", "GAG",                  # Glu
    "GGT", "GGC", "GGA", "GGG",   # Gly
    "CAT", "CAC",                  # His
    "ATT", "ATC", "ATA",           # Ile
    "TTA", "TTG", "CTT", "CTC",   # Leu
    "AAA", "AAG",                  # Lys
    "ATG",                         # Met
    "TTT", "TTC",                  # Phe
    "CCT", "CCC", "CCA", "CCG",   # Pro
    "TCT", "TCC", "TCA", "TCG",   # Ser
    "ACT", "ACC", "ACA", "ACG",   # Thr
    "TGG",                         # Trp
    "TAT", "TAC",                  # Tyr
    "GTT", "GTC", "GTA", "GTG",   # Val
]

def make_fake_orf(length_codons: int, rng: random.Random) -> str:
    """
    Generate a fake but plausible ORF (Open Reading Frame).
    In real genomes, ORFs start with ATG and end with a stop codon.
    Prodigal (used by Prokka) identifies these using hexamer statistics —
    the frequency of 6-base patterns differs between coding and non-coding DNA.
    """
    codons = [START_CODON]
    codons += [rng.choice(CODING_CODONS) for _ in range(length_codons - 2)]
    codons.append(rng.choice(STOP_CODONS))
    return "".join(codons)


def generate_genome(length: int, rng: random.Random) -> str:
    """
    Generate a synthetic bacterial genome fragment.
    Roughly 85% of bacterial genome is coding; we embed fake ORFs every ~900 bp
    (average bacterial gene length ≈ 900 bp) in a background of random DNA.
    """
    genome = list(random_dna(length, rng))

    # Embed ORFs at roughly regular intervals
    pos = rng.randint(50, 150)
    while pos + 900 < length:
        gene_length_codons = rng.randint(100, 400)  # 100–400 codons
        orf = make_fake_orf(gene_length_codons, rng)
        end = min(pos + len(orf), length)
        genome[pos:end] = list(orf[:end - pos])
        pos += len(orf) + rng.randint(100, 400)   # intergenic spacer

    return "".join(genome)


# =============================================================================
# Dataset 3: Reference + paired-end reads for variant calling
# =============================================================================

def generate_reads_from_ref(
    ref: str, n_reads: int, read_len: int, rng: random.Random,
    insert_size: int = 400, insert_std: int = 50,
    snv_rate: float = 0.003,
) -> tuple[list[tuple[str, str, str]], list[tuple[str, str, str]]]:
    """
    Simulate paired-end Illumina reads from a reference sequence.

    How paired-end sequencing works:
      - A DNA fragment (the "insert") is sequenced from both ends.
      - R1 reads the forward strand from the left end.
      - R2 reads the reverse-complement from the right end.
      - The insert size is the distance between the outer edges of the pair.

    We also introduce random SNVs to simulate variants in our "sample"
    relative to the reference.
    """
    # Create a "sample genome" with some variants relative to reference
    sample = list(ref)
    variant_positions = set()
    for i in range(len(sample)):
        if rng.random() < snv_rate:
            old = sample[i]
            new = rng.choice([b for b in DNA_BASES if b != old])
            sample[i] = new
            variant_positions.add(i)
    sample_str = "".join(sample)

    r1_reads: list[tuple[str, str, str]] = []
    r2_reads: list[tuple[str, str, str]] = []

    for i in range(n_reads):
        # Pick a random insert position
        frag_len = max(read_len * 2, int(rng.gauss(insert_size, insert_std)))
        frag_len = min(frag_len, len(sample_str) - 1)
        start = rng.randint(0, len(sample_str) - frag_len - 1)

        fragment = sample_str[start:start + frag_len]
        r1_seq = fragment[:read_len]
        r2_seq = fragment[-read_len:][::-1]           # reverse
        r2_seq = r2_seq.translate(str.maketrans("ACGT", "TGCA"))  # complement

        qual = make_phred_quality(read_len, rng)

        read_name = f"read{i+1:06d}:pos{start}"
        r1_reads.append((f"{read_name}/1", r1_seq, qual))
        r2_reads.append((f"{read_name}/2", r2_seq, qual[::-1]))  # R2 quality is reversed

    return r1_reads, r2_reads


# =============================================================================
# Main
# =============================================================================

def main() -> None:
    args = parse_args()
    rng = random.Random(args.seed)
    outdir = Path(args.outdir)

    # Create subdirectories
    phylo_dir = outdir / "phylogenetics"
    annot_dir = outdir / "annotation"
    var_dir   = outdir / "variants"
    for d in [phylo_dir, annot_dir, var_dir]:
        d.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("  NGS Education Hub — Test Data Generator")
    print(f"  Random seed: {args.seed}")
    print(f"  Output dir:  {outdir}")
    print("=" * 60)

    # ── Dataset 1: 16S sequences ─────────────────────────────────────────────
    print("\n[1/3] Generating 16S rRNA sequences for phylogenetics...")
    records_16S = generate_16S_dataset(n_taxa=20, rng=rng)
    fasta_16S = phylo_dir / "test_16S.fasta"
    write_fasta(fasta_16S, records_16S)
    print(f"      {len(records_16S)} sequences → {fasta_16S}")
    print(f"      Run:  bash pipelines/phylogenetics/run_phylogenetics.sh -i {fasta_16S}")

    # ── Dataset 2: Synthetic genome ──────────────────────────────────────────
    print("\n[2/3] Generating synthetic bacterial genome for annotation...")
    genome_seq = generate_genome(length=10_000, rng=rng)
    fasta_genome = annot_dir / "test_genome.fasta"
    write_fasta(fasta_genome, [("synthetic_genome_1 len=10000", genome_seq)])
    print(f"      1 contig × 10 000 bp → {fasta_genome}")
    print(f"      Edit pipelines/annotation/config.sh: ASSEMBLY={fasta_genome}")

    # ── Dataset 3: Reference + reads for variant calling ────────────────────
    print("\n[3/3] Generating reference and paired-end reads for variant calling...")
    ref_seq = generate_genome(length=50_000, rng=rng)
    fasta_ref = var_dir / "test_ref.fasta"
    write_fasta(fasta_ref, [("chr_test len=50000", ref_seq)])

    r1_reads, r2_reads = generate_reads_from_ref(
        ref=ref_seq, n_reads=5_000, read_len=150, rng=rng,
        snv_rate=0.003,  # ~150 SNVs in 50 kbp (rough WGS-like density)
    )
    r1_path = var_dir / "test_R1.fastq.gz"
    r2_path = var_dir / "test_R2.fastq.gz"
    write_fastq_gz(r1_path, r1_reads)
    write_fastq_gz(r2_path, r2_reads)

    print(f"      Reference:  {fasta_ref}  (50 000 bp)")
    print(f"      Reads R1:   {r1_path}  ({len(r1_reads)} reads × 150 bp)")
    print(f"      Reads R2:   {r2_path}")
    print(f"      Expected coverage: ~{(len(r1_reads) * 2 * 150) // 50_000}×")
    print(f"      Edit pipelines/variants/config.sh:")
    print(f"        REF_FASTA={fasta_ref}")
    print(f"        READS_R1={r1_path}")
    print(f"        READS_R2={r2_path}")
    print(f"        SAMPLE_ID=test_sample")

    print("\n" + "=" * 60)
    print("  Test data generation complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()
