# paired-kmers Architecture

## Overview
Fast, memory-efficient k-mer based method to identify conserved regions (paired kmers) in highly diverse microbial genomes for PCR primer design.

## Core Data Flow
```
Input (multi-FASTA of diverse genomes)
  → K-mer extraction per genome
    → Column-based + row-based reduction (spill to disk)
      → Paired k-mer detection (forward + reverse-complement within spacing range)
        → Output (conserved region BED/FASTA)
```

## Key Design Decisions

### 1. Two-Phase Reduction
- **Column-based**: Process genomes in batches to reduce memory
- **Row-based**: Final pairing pass with bounded memory
- Spills intermediate results to temp files when memory limit reached

### 2. Paired K-mer Definition
Two k-mers form a pair if:
- They are reverse-complements of each other
- Their distance in the genome is within `[-s, -S]` (min/max spacing)
- They appear in >= `-c` fraction of input genomes

### 3. No External Dependencies
Single Go file, no shared libraries. Self-contained for easy deployment.

## Performance Notes
- Memory usage is the primary constraint; designed to run on 16GB machines with 100+ microbial genomes
- CPU-bound during k-mer counting; benefits from fast SSD for temp files
