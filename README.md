# GECKO (v2 — Rust)

A fast, memory-safe pairwise genome comparison tool for the detection of ungapped High-scoring Segment Pairs (HSPs).

GECKO (*GEnome Comparison with K-mers Out-of-core*) identifies collections of HSPs across single- and multi-FASTA genomes, from bacterial assemblies to chromosome-scale sequences. This v2 release modernizes GECKO into a single Rust binary while preserving full drop-in compatibility with GECKO-MGV (`.frags`, `.csv`, `.INF`, `.MAT`) and adding modern bioinformatics formats (**PAF**, **BLAST tabular**, **BEDPE**, **JSONL**) and transparent `.gz` FASTA support.

## Highlights vs. Legacy GECKO (v1)

- **Single unified binary (`gecko`)** replacing 19 separate C binaries and 6 Bash scripts.
- **Pluggable CPU + GPU compute backends (`--backend auto|cuda|wgpu|cpu`)**:
  - **Fast AVX2 CPU (`cpu`)**: single canonical Reference-Y Hash-CSR / Direct-CSR index, 64-wide prefetched streaming query, lock-free 2-tier diagonal cache, bit-exact AVX2 prefix-sum ungapped extension, and automatic out-of-core spilling when `--memory` is tight or `--cache-dir` is enabled.
  - **Cross-vendor WGPU (`wgpu`, enabled by default)**: Vulkan / Metal / DX12 WGSL compute shaders for AMD, Apple Silicon, Intel, and NVIDIA GPUs.
  - **Native NVIDIA CUDA + CUB (`cuda`, optional `--features cuda`)**: GPU radix-sorted canonical Hash-CSR index, lock-free left/right-maximal seed filter, two-tier warp-cooperative extension, and overlapped CPU query loading with GPU reference indexing.
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
