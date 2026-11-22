# GECKO (v2 — Rust)

A fast, memory-safe pairwise genome comparison tool for the detection of ungapped High-scoring Segment Pairs (HSPs).

GECKO (*GEnome Comparison with K-mers Out-of-core*) identifies collections of HSPs across single- and multi-FASTA genomes, from bacterial assemblies to chromosome-scale sequences. This v2 release modernizes GECKO into a single Rust binary while preserving full drop-in compatibility with GECKO-MGV (`.frags`, `.csv`, `.INF`, `.MAT`) and adding modern bioinformatics formats (**PAF**, **BLAST tabular**, **BEDPE**, **JSONL**) and transparent `.gz` FASTA support.

## Highlights vs. Legacy GECKO (v1)

- **Single unified binary (`gecko`)** replacing 19 separate C binaries and 6 Bash scripts.
- **2-bit/bp `PackedSeq` + zero-spill Reference-Slab Streaming**: packs sequences at `0.25 bytes/bp` with sparse run-length annotations for `'N'`, soft-masking, and multi-FASTA separators, cutting memory footprint and allowing chromosome-scale runs inside tight `--memory` budgets with zero disk spilling.
- **Pluggable CPU + GPU compute backends (`--backend auto|cuda|wgpu|cpu`)**:
  - **Fast AVX2/BMI2 CPU (`cpu`)**: single canonical Reference-Y Hash-CSR / Direct-CSR index, 64-wide prefetched streaming query, lock-free 2-tier diagonal cache, bit-exact 2-bit + AVX2/BMI2 prefix-sum ungapped extension, zero-spill reference-slab streaming, and out-of-core `.gidx` caching when `--cache-dir` is enabled.
  - **Cross-vendor WGPU (`wgpu`, enabled by default)**: Vulkan / Metal / DX12 WGSL compute shaders with adaptive sub-tiling and 1-bit `--freq-limit` dropped-position bitmap for AMD, Apple Silicon, Intel, and NVIDIA GPUs.
  - **Native NVIDIA CUDA + CUB (`cuda`, optional `--features cuda`)**: GPU radix-sorted canonical Hash-CSR index, 1-bit `--freq-limit` dropped-position bitmap + lock-free left/right-maximal seed filter, adaptive sub-tiling, two-tier warp-cooperative extension, and overlapped CPU query loading with GPU reference indexing.
- **High-speed I/O**: 4 MB chunked constant-line-width stride-predicted FASTA loading (~2.2 GB/s/core), concurrent query/reference loading, and 4 MB buffered output writers.
- **Any seed word length `4..=32`**: no longer restricted to multiples of 4.
- **Compressed FASTA input**: reads `.fa`, `.fasta`, `.fna`, and `.gz` seamlessly.
- **Drop-in + modern formats**: outputs GECKO-MGV `.frags` (exact 93-byte big-endian format), `.csv`, `.frags.INF`, `.frags.MAT`, plus **PAF**, **BLAST tabular (`outfmt 6`)**, **BEDPE**, **JSONL**, and **PNG** dotplots.
- **Fixes legacy bugs**: resolves the `readFragmentsv2` 104B/93B truncation bug (~10.6% dropped CSV fragments), `coverage[..][100]` OOB write on 100%-identical HSPs, `loadWordOcurrences` `realloc` UB, `reverseComplement` 2 GB fixed allocation / 8 MB stack array, `filterFrags` 100× similarity scale bug, and `workflow.sh` extension mismatch.

## Installation

### From Source (Rust 1.88+)

```bash
# Default build: Fast AVX2 CPU + cross-vendor WGPU (Vulkan / Metal / DX12)
cargo build --release --bin gecko

# With native NVIDIA CUDA + CUB GPU support (requires CUDA Toolkit / nvcc)
cargo build --release --bin gecko --features cuda
# Optional: target a specific GPU SM architecture for faster builds, e.g. T4:
# GECKO_CUDA_ARCH="-gencode arch=compute_75,code=sm_75" cargo build --release --bin gecko --features cuda

# Inspect detected backends
./target/release/gecko backends
```

### Docker

```bash
docker build -t gecko .
docker run --rm -v "$PWD:/data" gecko compare /data/query.fa.gz /data/ref.fa.gz -o /data/results
```

## Quick Start

```bash
# Automatic backend selection (CUDA -> hardware WGPU -> fast AVX2 CPU)
gecko compare query.fasta reference.fasta -l 40 -s 60 -k 12 -o results/

# Explicit backend selection (--backend auto|cuda|wgpu|cpu)
gecko compare chr1_human.fa chr1_gorilla.fa -l 100 -s 60 -k 32 --backend cuda -o results/

# Output all supported formats (frags, csv, inf, mat, paf, blast, bedpe, jsonl, png)
gecko compare query.fa.gz reference.fa.gz -l 100 -s 60 -k 16 --format all -o results/

# Bound RAM for chromosome-scale comparisons and reuse k-mer dictionaries
gecko compare chr1_human.fa.gz chr1_gorilla.fa.gz -l 100 -s 60 -k 32 -m 4G -t 16 --cache-dir .gecko_cache -o results/
```

### Migration from v1 (`workflow.sh`)

| Legacy v1 Command | GECKO v2 Equivalent |
|---|---|
| `./bin/workflow.sh q.fa r.fa 40 60 12 1` | `gecko compare q.fa r.fa -l 40 -s 60 -k 12 --format legacy -o results` |
| `./bin/frags2align.sh q-r.frags q.fa r.fa align.txt` | `gecko align q-r.frags q.fa r.fa -o align.txt` |
| `./bin/frags2borders.sh q-r.csv q.fa r.fa out.txt 50 query` | `gecko borders q-r.csv q.fa r.fa -b 50 -o out.txt --emit-query query-50.fasta` |
| `./bin/filterFrags q-r.frags 100 70 > out.csv` | `gecko filter q-r.frags -l 100 -s 70 -o out.csv` |
| `./bin/allVsAll.sh genomes/ 200 40 12 fa` | `gecko all-vs-all genomes/ -l 200 -s 40 -k 12 --ext fa -o results` |
| `./bin/fragStat q-r.frags` | `gecko stats q-r.frags` |

## Subcommands

- `gecko compare`: Pairwise genome comparison (`-l`, `-s`, `-k`, `--backend auto|cuda|wgpu|cpu`, `--freq-limit`, `--max-rep`, `--percent-length`, `--no-soft-mask`, `--forward-only`, `--format`, `-m/--memory`, `-t/--threads`, `--cache-dir`, `--json-stats`).
- `gecko backends`: Probe and display available compute backends (CUDA, WGPU, AVX2 CPU).
- `gecko index`: Precompute a compact binary `.gidx` k-mer dictionary.
- `gecko align`: Extract wrapped text pairwise alignments from `.frags` or `.csv` + original FASTAs.
- `gecko borders`: Extract flanking border alignments (`-b/--border`) and optionally emit query sequences as FASTA (`--emit-query`).
- `gecko filter`: Re-filter `.frags` or `.csv` by minimum length and similarity into CSV, `.frags`, PAF, BLAST, BEDPE, or JSONL.
- `gecko plot`: Render a 2D synteny matrix dotplot PNG directly from `.frags` or `.csv` (`--width`, `--height`, `--line-width`, `--color strand|similarity`, `--origin top-left|bottom-left`, `--dark`, `--query`/`--reference` for contig boundary dividers).
- `gecko all-vs-all`: Compare every genome pair in a directory, indexing each genome only once.
- `gecko stats`: Inspect `.frags`, `.csv`, or `.gidx` files (human-readable or `--json`).

## Output Formats

- **`frags`**: 93-byte big-endian binary format compatible with GECKO-MGV and legacy downstream tools.
- **`csv`**: GECKO-MGV CSV with metadata header and columns `Type,xStart,yStart,xEnd,yEnd,strand(f/r),block,length,score,ident,similarity,%ident,SeqX,SeqY`.
- **`inf` / `mat`**: Legacy `.frags.INF` summary metadata and `.frags.MAT` 1000×100 length-by-similarity histogram.
- **`paf`**: Minimap2 Pairwise mApping Format with 0-based per-sequence local coordinates and `AS:i:`, `si:f:`, `cg:Z:` tags.
- **`blast`**: BLAST tabular (`outfmt 6`) with 1-based inclusive coordinates, ungapped Karlin-Altschul E-value, and bitscore.
- **`bedpe`**: Pairwise 0-based BEDPE synteny intervals.
- **`jsonl`**: Structured newline-delimited JSON records.
- **`png`**: 2D synteny matrix dotplot image with axes, gridlines, legend, and optional contig boundary dividers.

## Whole-Genome Benchmark (Human `GRCh38.p14` vs. Chimpanzee `Clint_PTRv2`)

Measured across all **24 homologous chromosome pairs** (`chr1`–`chr22`, `chrX`, `chrY`, pairing Human `chr2` against Chimpanzee `chr2A + chr2B`; **3,088.28 Mb vs. 2,832.61 Mb**) on Google Colab (**2 vCPUs + NVIDIA Tesla T4**, `WL=32, L=1000, S=75%, FL=20, -m 4G`, zero disk spill) using [`notebooks/gecko_rust_gpu_benchmark.ipynb`](notebooks/gecko_rust_gpu_benchmark.ipynb):

| Backend | Total Wall Time | Total Compute Time | Peak Host RSS | Total Fragments | Speedup | Recall vs. CPU |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: |
| **CPU (`AVX2`, 2 threads)** | `670.50s` (`11.18 min`) | `625.77s` | `3,268.5 MiB` | `2,567,621` | `1.00x` | reference |
| **CUDA (`Tesla T4`, SM 7.5)** | **`111.15s` (`1.85 min`)** | **`105.89s`** | **`619.9 MiB`** | `2,567,765` | **`6.03x`** | **`99.997%`** |

<details>
<summary>Full 24-chromosome per-pair breakdown</summary>

| Chr | Human (Mb) | Chimp (Mb) | CPU Wall | CPU Compute | CPU RSS (MiB) | CPU Frags | CUDA Wall | CUDA Compute | CUDA RSS (MiB) | CUDA Frags | Speedup | Recall |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| chr1 | 248.96 | 224.24 | 55.60s | 52.47s | 3,169.9 | 182,049 | 8.65s | 8.43s | 619.1 | 182,050 | 6.43x | 99.998% |
| chr2 (2A+2B) | 242.19 | 236.78 | 57.91s | 54.51s | 3,268.5 | 198,904 | 8.75s | 8.50s | 619.9 | 198,917 | 6.62x | 99.997% |
| chr3 | 198.30 | 196.56 | 47.58s | 44.22s | 2,788.0 | 198,837 | 7.33s | 7.07s | 534.4 | 198,849 | 6.50x | 99.999% |
| chr4 | 190.21 | 189.15 | 43.72s | 41.16s | 2,692.4 | 207,091 | 7.96s | 7.72s | 514.6 | 207,104 | 5.49x | 99.999% |
| chr5 | 181.54 | 159.32 | 39.08s | 36.62s | 2,337.4 | 178,515 | 6.06s | 5.83s | 482.6 | 178,546 | 6.45x | 99.998% |
| chr6 | 170.81 | 168.37 | 39.56s | 37.27s | 2,439.5 | 161,937 | 6.24s | 6.01s | 472.8 | 161,945 | 6.34x | 99.999% |
| chr7 | 159.35 | 156.05 | 37.63s | 35.49s | 2,267.6 | 137,621 | 5.71s | 5.48s | 447.3 | 137,626 | 6.59x | 99.997% |
| chr8 | 145.14 | 143.34 | 32.53s | 29.53s | 2,124.1 | 144,662 | 5.19s | 4.97s | 422.1 | 144,665 | 6.27x | 99.995% |
| chr9 | 138.39 | 110.51 | 27.17s | 25.45s | 1,660.3 | 102,247 | 4.29s | 4.07s | 389.9 | 102,245 | 6.33x | 99.995% |
| chr10 | 133.80 | 129.81 | 30.32s | 28.50s | 1,823.4 | 110,690 | 5.05s | 4.82s | 396.8 | 110,700 | 6.01x | 99.998% |
| chr11 | 135.09 | 130.78 | 31.68s | 29.27s | 1,853.0 | 126,682 | 4.92s | 4.70s | 401.4 | 126,692 | 6.43x | 99.999% |
| chr12 | 133.28 | 131.00 | 31.62s | 28.90s | 1,849.0 | 114,395 | 4.85s | 4.63s | 393.0 | 114,404 | 6.51x | 99.999% |
| chr13 | 114.36 | 95.60 | 19.45s | 18.05s | 1,457.1 | 63,575 | 3.56s | 3.34s | 340.4 | 63,579 | 5.46x | 99.998% |
| chr14 | 107.04 | 87.72 | 18.75s | 17.43s | 1,376.0 | 73,899 | 3.48s | 3.26s | 340.4 | 73,905 | 5.39x | 100.000% |
| chr15 | 101.99 | 80.52 | 17.46s | 16.21s | 1,276.3 | 56,149 | 3.19s | 2.98s | 340.4 | 56,150 | 5.46x | 99.996% |
| chr16 | 90.34 | 75.91 | 17.27s | 16.11s | 1,153.7 | 42,774 | 3.04s | 2.83s | 340.4 | 42,774 | 5.68x | 99.998% |
| chr17 | 83.26 | 76.55 | 18.00s | 16.89s | 1,170.7 | 31,968 | 3.11s | 2.90s | 340.4 | 31,965 | 5.79x | 99.984% |
| chr18 | 80.37 | 74.77 | 15.00s | 13.92s | 1,155.5 | 54,180 | 2.79s | 2.58s | 340.4 | 54,184 | 5.38x | 100.000% |
| chr19 | 58.62 | 56.73 | 15.47s | 14.70s | 859.3 | 23,188 | 2.97s | 2.76s | 340.4 | 23,190 | 5.21x | 99.996% |
| chr20 | 64.44 | 64.04 | 14.51s | 13.66s | 950.8 | 36,305 | 2.43s | 2.23s | 340.4 | 36,303 | 5.96x | 99.992% |
| chr21 | 46.71 | 33.23 | 7.68s | 7.12s | 581.4 | 13,836 | 1.44s | 1.25s | 340.4 | 13,837 | 5.33x | 99.993% |
| chr22 | 50.82 | 33.70 | 8.37s | 7.79s | 605.4 | 12,666 | 1.55s | 1.36s | 340.4 | 12,667 | 5.40x | 100.000% |
| chrX | 156.04 | 151.58 | 37.98s | 34.84s | 2,207.3 | 266,024 | 7.23s | 7.01s | 471.6 | 266,040 | 5.25x | 99.998% |
| chrY | 57.23 | 26.35 | 6.15s | 5.65s | 504.4 | 29,427 | 1.36s | 1.16s | 398.5 | 29,428 | 4.53x | 99.997% |
| **TOTAL** | **3088.28** | **2832.61** | **670.50s** | **625.77s** | **3,268.5** | **2,567,621** | **111.15s** | **105.89s** | **619.9** | **2,567,765** | **6.03x** | **99.997%** |

</details>

## Development & Testing

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
```

## Cite

If you use GECKO, please cite:

> Torreno, Oscar, and Oswaldo Trelles. "Breaking the computational barriers of pairwise genome comparison." *BMC Bioinformatics* 16, no. 1 (2015): 250.
> https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0679-9
