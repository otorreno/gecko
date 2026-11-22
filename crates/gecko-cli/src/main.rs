mod backend;

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{Context, Result, bail};
use backend::BackendArg;
use clap::{Parser, Subcommand, ValueEnum};
use gecko_core::formats::align::{write_alignments_path, write_borders_path};
use gecko_core::formats::csv::{
    ComparisonMetadata, read_csv_path, write_csv_path, write_inf_path, write_mat_path,
};
use gecko_core::formats::frags::{FragsFile, read_frags_path, write_frags_path};
use gecko_core::formats::modern::{
    write_bedpe_path, write_blast_tabular_path, write_jsonl_path, write_paf_path,
};
use gecko_core::formats::plot::{ColorMode, PlotOptions, PlotOrigin, write_dotplot_png_path};
use gecko_core::kmer::{DictOptions, INDEX_WORD_SIZE, KmerDict};
use gecko_core::pipeline::{CompareConfig, compare_genomes, compare_with_ref_index};
use gecko_core::ref_index::RefIndex;
use gecko_core::seq::Genome;
use gecko_core::{Fragment, Scoring};

#[derive(Parser)]
#[command(
    name = "gecko",
    version,
    about = "GECKO v2: Fast pairwise genome comparison for High-scoring Segment Pair (HSP) detection",
    long_about = None
)]
struct Cli {
    /// Number of worker threads (defaults to available logical CPUs).
    #[arg(short = 't', long, global = true)]
    threads: Option<usize>,

    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Compare a query FASTA (or .fa.gz) against a reference FASTA (or .fa.gz).
    Compare(CompareArgs),

    /// List available compute backends (CUDA, WGPU, AVX2 CPU).
    Backends(BackendsArgs),

    /// Precompute a reusable binary k-mer dictionary (.gidx) for a FASTA file.
    Index(IndexArgs),

    /// Extract human-readable pairwise alignments from a .frags or .csv file.
    Align(AlignArgs),

    /// Extract flanking border alignments (and optionally query FASTA) from .frags or .csv.
    Borders(BordersArgs),

    /// Filter a .frags or .csv file by minimum length and similarity thresholds.
    Filter(FilterArgs),

    /// Render a 2D synteny matrix dotplot PNG from a .frags or .csv file.
    Plot(PlotArgs),

    /// Compare all genome pairs in a directory (reusing dictionaries across pairs).
    AllVsAll(AllVsAllArgs),

    /// Print summary statistics for a .frags, .csv, or .gidx file.
    Stats(StatsArgs),
}

#[derive(Parser)]
struct BackendsArgs {
    /// Include software WGPU adapters (e.g. Mesa Lavapipe) when probing.
    #[arg(long)]
    allow_software_gpu: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
enum OutputFormat {
    /// Legacy GECKO 93-byte big-endian binary fragments (.frags)
    Frags,
    /// Legacy GECKO-MGV compatible CSV with metadata header (.csv)
    Csv,
    /// Legacy GECKO metadata info file (.frags.INF)
    Inf,
    /// Legacy GECKO 1000x100 length-by-similarity matrix (.frags.MAT)
    Mat,
    /// Minimap2 Pairwise mApping Format (.paf)
    Paf,
    /// BLAST tabular outfmt 6 (.tsv)
    Blast,
    /// Pairwise BEDPE synteny intervals (.bedpe)
    Bedpe,
    /// Newline-delimited JSON (.jsonl)
    Jsonl,
    /// 2D synteny matrix dotplot image (.png)
    Png,
    /// Preset: frags + csv + inf + mat (exact drop-in GECKO-MGV output set)
    Legacy,
    /// Preset: all legacy and modern output formats
    All,
}

impl std::fmt::Display for OutputFormat {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = match self {
            OutputFormat::Frags => "frags",
            OutputFormat::Csv => "csv",
            OutputFormat::Inf => "inf",
            OutputFormat::Mat => "mat",
            OutputFormat::Paf => "paf",
            OutputFormat::Blast => "blast",
            OutputFormat::Bedpe => "bedpe",
            OutputFormat::Jsonl => "jsonl",
            OutputFormat::Png => "png",
            OutputFormat::Legacy => "legacy",
            OutputFormat::All => "all",
        };
        f.write_str(s)
    }
}

impl OutputFormat {
    fn expand(formats: &[OutputFormat]) -> BTreeSet<OutputFormat> {
        let mut set = BTreeSet::new();
        for &f in formats {
            match f {
                OutputFormat::Legacy => {
                    set.insert(OutputFormat::Frags);
                    set.insert(OutputFormat::Csv);
                    set.insert(OutputFormat::Inf);
                    set.insert(OutputFormat::Mat);
                }
                OutputFormat::All => {
                    set.insert(OutputFormat::Frags);
                    set.insert(OutputFormat::Csv);
                    set.insert(OutputFormat::Inf);
                    set.insert(OutputFormat::Mat);
                    set.insert(OutputFormat::Paf);
                    set.insert(OutputFormat::Blast);
                    set.insert(OutputFormat::Bedpe);
                    set.insert(OutputFormat::Jsonl);
                    set.insert(OutputFormat::Png);
                }
                other => {
                    set.insert(other);
                }
            }
        }
        set
    }
}

#[derive(Parser)]
struct CompareArgs {
    /// Query genome in FASTA or gzipped FASTA (.gz) format.
    query: PathBuf,

    /// Reference genome in FASTA or gzipped FASTA (.gz) format.
    reference: PathBuf,

    /// Minimum HSP length in bp (or percentage of shorter sequence if --percent-length is set).
    #[arg(short = 'l', long, default_value_t = 40)]
    min_length: u64,

    /// Minimum similarity percentage threshold (0..=100).
    #[arg(short = 's', long, default_value_t = 60.0)]
    min_similarity: f32,

    /// Seed k-mer word length (4..=32; any integer supported, not limited to multiples of 4).
    #[arg(short = 'k', long, default_value_t = 12)]
    word_length: usize,

    /// Skip 32-mers occurring more than this many times in a genome.
    #[arg(long, default_value_t = 1000)]
    freq_limit: u64,

    /// Maximum sampled occurrences per word list.
    #[arg(long, default_value_t = 10_000)]
    max_rep: usize,

    /// Interpret --min-length as a percentage of the shorter sequence length (legacy fixedL=0).
    #[arg(long)]
    percent_length: bool,

    /// Include lowercase (soft-masked) bases when generating k-mer seeds.
    #[arg(long)]
    no_soft_mask: bool,

    /// Compare only the forward strand (skip reverse complement of reference).
    #[arg(long)]
    forward_only: bool,

    /// Output directory.
    #[arg(short = 'o', long, default_value = "results")]
    outdir: PathBuf,

    /// Output file prefix (defaults to `<query>-<reference>`).
    #[arg(long)]
    prefix: Option<String>,

    /// Output formats (comma-separated or repeated). Defaults to legacy GECKO set + PAF.
    #[arg(
        short = 'f',
        long = "format",
        value_delimiter = ',',
        default_values_t = [OutputFormat::Legacy, OutputFormat::Paf]
    )]
    formats: Vec<OutputFormat>,

    /// Optional directory to cache and reuse .gidx k-mer dictionaries across runs.
    #[arg(long)]
    cache_dir: Option<PathBuf>,

    /// Memory budget for in-memory sorting/buffering before spilling to disk (e.g. 512M, 2G, 8G).
    #[arg(short = 'm', long, default_value = "2G", value_parser = parse_byte_size)]
    memory: usize,

    /// Match score for ungapped extension.
    #[arg(long, default_value_t = 4)]
    match_score: i64,

    /// Mismatch penalty (positive value subtracted on mismatch).
    #[arg(long, default_value_t = 4)]
    mismatch_penalty: i64,

    /// Ambiguous 'N' base penalty (positive value subtracted on 'N').
    #[arg(long, default_value_t = 1)]
    n_penalty: i64,

    /// Compute backend (`auto` prefers CUDA -> hardware WGPU -> fast AVX2 CPU).
    #[arg(long, env = "GECKO_BACKEND", default_value = "auto")]
    backend: BackendArg,

    /// Allow software WGPU adapters (e.g. Mesa Lavapipe) when using `--backend wgpu|auto`.
    #[arg(long)]
    allow_software_gpu: bool,

    /// Print run statistics as JSON to stdout.
    #[arg(long)]
    json_stats: bool,
}

#[derive(Parser)]
struct IndexArgs {
    /// Input FASTA (or .fa.gz) file.
    fasta: PathBuf,

    /// Output .gidx dictionary path.
    #[arg(short = 'o', long)]
    output: PathBuf,

    /// Indexed k-mer size (4..=32, default 32).
    #[arg(short = 'k', long, default_value_t = INDEX_WORD_SIZE)]
    word_length: usize,

    /// Include lowercase (soft-masked) bases in k-mer seeds.
    #[arg(long)]
    no_soft_mask: bool,

    /// Memory budget for sorting (e.g. 512M, 2G).
    #[arg(short = 'm', long, default_value = "2G", value_parser = parse_byte_size)]
    memory: usize,
}

#[derive(Parser)]
struct AlignArgs {
    /// Input .frags or .csv fragment file.
    frags: PathBuf,

    /// Original query FASTA (or .fa.gz).
    query: PathBuf,

    /// Original reference FASTA (or .fa.gz).
    reference: PathBuf,

    /// Output text alignment file.
    #[arg(short = 'o', long)]
    output: PathBuf,
}

#[derive(Parser)]
struct BordersArgs {
    /// Input .frags or .csv fragment file.
    frags: PathBuf,

    /// Original query FASTA (or .fa.gz).
    query: PathBuf,

    /// Original reference FASTA (or .fa.gz).
    reference: PathBuf,

    /// Output text alignment file with borders.
    #[arg(short = 'o', long)]
    output: PathBuf,

    /// Number of flanking bases to extend on each side of every HSP.
    #[arg(short = 'b', long)]
    border: u64,

    /// Optional output path to also emit extracted query border sequences as FASTA.
    #[arg(long)]
    emit_query: Option<PathBuf>,
}

#[derive(Parser)]
struct FilterArgs {
    /// Input .frags or .csv fragment file.
    input: PathBuf,

    /// Minimum fragment length in bp (inclusive >=).
    #[arg(short = 'l', long)]
    min_length: u64,

    /// Minimum similarity percentage (0..=100, inclusive >=).
    #[arg(short = 's', long)]
    min_similarity: f32,

    /// Output file path.
    #[arg(short = 'o', long)]
    output: PathBuf,

    /// Output format (csv, frags, paf, blast, bedpe, jsonl).
    #[arg(short = 'f', long, default_value = "csv")]
    format: OutputFormat,

    /// Query FASTA (required when exporting paf, blast, bedpe, or jsonl).
    #[arg(long)]
    query: Option<PathBuf>,

    /// Reference FASTA (required when exporting paf, blast, bedpe, or jsonl).
    #[arg(long)]
    reference: Option<PathBuf>,
}

#[derive(Parser)]
struct AllVsAllArgs {
    /// Directory containing genome FASTA files.
    dir: PathBuf,

    /// Minimum HSP length in bp.
    #[arg(short = 'l', long, default_value_t = 200)]
    min_length: u64,

    /// Minimum similarity percentage threshold.
    #[arg(short = 's', long, default_value_t = 40.0)]
    min_similarity: f32,

    /// Seed k-mer word length (4..=32).
    #[arg(short = 'k', long, default_value_t = 12)]
    word_length: usize,

    /// File extensions to include (comma-separated).
    #[arg(long, value_delimiter = ',', default_values_t = ["fa".to_string(), "fasta".to_string(), "fna".to_string(), "gz".to_string()])]
    ext: Vec<String>,

    /// Output directory for all pairwise comparisons.
    #[arg(short = 'o', long, default_value = "results")]
    outdir: PathBuf,

    /// Output formats per pair.
    #[arg(
        short = 'f',
        long = "format",
        value_delimiter = ',',
        default_values_t = [OutputFormat::Legacy, OutputFormat::Paf]
    )]
    formats: Vec<OutputFormat>,

    /// Skip 32-mers occurring more than this many times.
    #[arg(long, default_value_t = 1000)]
    freq_limit: u64,

    /// Memory budget for sorting/buffering.
    #[arg(short = 'm', long, default_value = "2G", value_parser = parse_byte_size)]
    memory: usize,

    /// Skip pairs whose primary output file already exists.
    #[arg(long)]
    skip_existing: bool,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum PlotColorArg {
    Strand,
    Similarity,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum PlotOriginArg {
    TopLeft,
    BottomLeft,
}

#[derive(Parser)]
struct PlotArgs {
    /// Input .frags or .csv fragment file.
    input: PathBuf,

    /// Output PNG image path.
    #[arg(short = 'o', long)]
    output: PathBuf,

    /// Image width in pixels.
    #[arg(long, default_value_t = 1200)]
    width: u32,

    /// Image height in pixels.
    #[arg(long, default_value_t = 1200)]
    height: u32,

    /// Stroke thickness in pixels for each HSP segment.
    #[arg(long, default_value_t = 2)]
    line_width: u32,

    /// Color encoding (strand: blue fwd / red rev, or similarity heatmap).
    #[arg(long, default_value = "strand")]
    color: PlotColorArg,

    /// Matrix origin (top-left = matrix/GECKO-MGV style, bottom-left = Cartesian dotplot).
    #[arg(long, default_value = "top-left")]
    origin: PlotOriginArg,

    /// Render with a dark background theme.
    #[arg(long)]
    dark: bool,

    /// Optional minimum HSP length filter before plotting.
    #[arg(short = 'l', long)]
    min_length: Option<u64>,

    /// Optional minimum similarity filter (0..=100) before plotting.
    #[arg(short = 's', long)]
    min_similarity: Option<f32>,

    /// Optional custom plot title.
    #[arg(long)]
    title: Option<String>,

    /// Optional query FASTA (adds axis name and multi-FASTA contig boundary guides).
    #[arg(long)]
    query: Option<PathBuf>,

    /// Optional reference FASTA (adds axis name and multi-FASTA contig boundary guides).
    #[arg(long)]
    reference: Option<PathBuf>,
}

#[derive(Parser)]
struct StatsArgs {
    /// Path to a .frags, .csv, or .gidx file.
    input: PathBuf,

    /// Emit JSON instead of human-readable text.
    #[arg(long)]
    json: bool,
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    if let Some(threads) = cli.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .context("failed to configure rayon thread pool")?;
    }

    match cli.command {
        Commands::Compare(args) => cmd_compare(args),
        Commands::Backends(args) => {
            backend::print_backends(args.allow_software_gpu);
            Ok(())
        }
        Commands::Index(args) => cmd_index(args),
        Commands::Align(args) => cmd_align(args),
        Commands::Borders(args) => cmd_borders(args),
        Commands::Filter(args) => cmd_filter(args),
        Commands::Plot(args) => cmd_plot(args),
        Commands::AllVsAll(args) => cmd_all_vs_all(args),
        Commands::Stats(args) => cmd_stats(args),
    }
}

fn cmd_compare(args: CompareArgs) -> Result<()> {
    let t0 = Instant::now();
    let config = CompareConfig {
        min_length: args.min_length,
        min_similarity: args.min_similarity,
        word_length: args.word_length,
        freq_limit: args.freq_limit,
        max_rep: args.max_rep,
        fixed_length: !args.percent_length,
        soft_mask: !args.no_soft_mask,
        memory_budget: args.memory,
        compare_reverse: !args.forward_only,
        scoring: Scoring {
            match_score: args.match_score,
            mismatch_penalty: args.mismatch_penalty,
            n_penalty: args.n_penalty,
        },
    };

    let run = backend::run_compare_from_paths(
        &args.query,
        &args.reference,
        &config,
        args.cache_dir.as_deref(),
        args.backend,
        args.allow_software_gpu,
    )?;
    let query = run.query;
    let reference = run.reference;
    let res = run.result;

    let t_write_start = Instant::now();
    let formats = OutputFormat::expand(&args.formats);
    let prefix = args
        .prefix
        .unwrap_or_else(|| format!("{}-{}", query.name, reference.name));

    std::fs::create_dir_all(&args.outdir)?;
    let q_fname = args
        .query
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or(&query.name);
    let r_fname = args
        .reference
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or(&reference.name);
    let meta = res.metadata(q_fname, r_fname, &query, &reference, &config);

    write_outputs(
        &args.outdir,
        &prefix,
        &formats,
        &query,
        &reference,
        &meta,
        &res.fragments,
        Some(&res.coverage),
    )?;
    let write_secs = t_write_start.elapsed().as_secs_f64();

    let elapsed = t0.elapsed().as_secs_f64();
    if args.json_stats {
        let summary = serde_json::json!({
            "backend": run.backend_label,
            "query": query.name,
            "reference": reference.name,
            "query_bases": query.total_bases,
            "reference_bases": reference.total_bases,
            "fragments": res.fragments.len(),
            "forward_fragments": res.stats.forward_fragments,
            "reverse_fragments": res.stats.reverse_fragments,
            "word_matches": res.stats.word_matches,
            "raw_hits": res.stats.raw_hits,
            "filtered_hits": res.stats.filtered_hits,
            "hits_used": res.stats.hits_used,
            "spilled_runs": res.stats.spilled_runs,
            "load_secs": run.load_secs,
            "compute_secs": run.compute_secs,
            "write_secs": write_secs,
            "elapsed_secs": elapsed,
        });
        println!("{}", serde_json::to_string_pretty(&summary)?);
    } else {
        eprintln!(
            "GECKO [{}] done in {:.3}s (load={:.3}s compute={:.3}s write={:.3}s) | frags={} (fwd={}, rev={}) | raw_hits={} filtered={} used={} spilled_runs={}",
            run.backend_label,
            elapsed,
            run.load_secs,
            run.compute_secs,
            write_secs,
            res.fragments.len(),
            res.stats.forward_fragments,
            res.stats.reverse_fragments,
            res.stats.raw_hits,
            res.stats.filtered_hits,
            res.stats.hits_used,
            res.stats.spilled_runs
        );
    }

    Ok(())
}

fn cmd_index(args: IndexArgs) -> Result<()> {
    let t0 = Instant::now();
    let genome = Genome::from_fasta_path(&args.fasta)?;
    let dict = KmerDict::build(
        &genome,
        DictOptions {
            soft_mask: !args.no_soft_mask,
            memory_budget: args.memory,
            index_k: args.word_length,
        },
    )?;
    if let Some(parent) = args.output.parent()
        && !parent.as_os_str().is_empty()
    {
        std::fs::create_dir_all(parent)?;
    }
    dict.save_to_path(&args.output)?;
    let st = dict.stats();
    eprintln!(
        "Indexed {} in {:.3}s | k={} unique_words={} total_occ={} max_freq={} -> {}",
        genome.name,
        t0.elapsed().as_secs_f64(),
        st.index_k,
        st.unique_words,
        st.total_occurrences,
        st.max_frequency,
        args.output.display()
    );
    Ok(())
}

fn cmd_align(args: AlignArgs) -> Result<()> {
    let ff = load_frags_or_csv(&args.frags)?;
    let query = Genome::from_fasta_path(&args.query)?;
    let reference = Genome::from_fasta_path(&args.reference)?;
    write_alignments_path(&args.output, &query, &reference, &ff.fragments)?;
    eprintln!(
        "Wrote {} alignments -> {}",
        ff.fragments.len(),
        args.output.display()
    );
    Ok(())
}

fn cmd_borders(args: BordersArgs) -> Result<()> {
    let ff = load_frags_or_csv(&args.frags)?;
    let query = Genome::from_fasta_path(&args.query)?;
    let reference = Genome::from_fasta_path(&args.reference)?;
    write_borders_path(
        &args.output,
        &query,
        &reference,
        &ff.fragments,
        args.border,
        args.emit_query.as_deref(),
    )?;
    eprintln!(
        "Wrote {} border alignments (+/-{} bp) -> {}",
        ff.fragments.len(),
        args.border,
        args.output.display()
    );
    Ok(())
}

fn cmd_filter(args: FilterArgs) -> Result<()> {
    let ff = load_frags_or_csv(&args.input)?;
    let filtered: Vec<Fragment> = ff
        .fragments
        .into_iter()
        .filter(|f| f.length >= args.min_length && f.similarity >= args.min_similarity)
        .collect();

    match args.format {
        OutputFormat::Frags => {
            write_frags_path(&args.output, ff.x_total, ff.y_total, &filtered)?;
        }
        OutputFormat::Csv => {
            let meta = ComparisonMetadata {
                seq_x_filename: "Unknown".into(),
                seq_y_filename: "Unknown".into(),
                seq_x_name: "Unknown".into(),
                seq_y_name: "Unknown".into(),
                seq_x_len: ff.x_total,
                seq_y_len: ff.y_total,
                min_length: args.min_length,
                min_similarity: args.min_similarity.round() as u64,
                total_hits: 0,
                total_hits_used: 0,
                total_fragments: filtered.len() as u64,
            };
            write_csv_path(&args.output, &meta, &filtered)?;
        }
        OutputFormat::Paf | OutputFormat::Blast | OutputFormat::Bedpe | OutputFormat::Jsonl => {
            let q_path = args
                .query
                .as_ref()
                .context("--query FASTA is required for PAF/BLAST/BEDPE/JSONL output")?;
            let r_path = args
                .reference
                .as_ref()
                .context("--reference FASTA is required for PAF/BLAST/BEDPE/JSONL output")?;
            let query = Genome::from_fasta_path(q_path)?;
            let reference = Genome::from_fasta_path(r_path)?;
            match args.format {
                OutputFormat::Paf => write_paf_path(&args.output, &query, &reference, &filtered)?,
                OutputFormat::Blast => {
                    write_blast_tabular_path(&args.output, &query, &reference, &filtered)?
                }
                OutputFormat::Bedpe => {
                    write_bedpe_path(&args.output, &query, &reference, &filtered)?
                }
                OutputFormat::Jsonl => {
                    write_jsonl_path(&args.output, &query, &reference, &filtered)?
                }
                _ => unreachable!(),
            }
        }
        other => bail!("unsupported output format {other:?} for `gecko filter`"),
    }

    eprintln!(
        "Filtered to {} fragments -> {}",
        filtered.len(),
        args.output.display()
    );
    Ok(())
}

fn cmd_all_vs_all(args: AllVsAllArgs) -> Result<()> {
    let t0 = Instant::now();
    let mut files: Vec<PathBuf> = std::fs::read_dir(&args.dir)
        .with_context(|| format!("failed to read directory {}", args.dir.display()))?
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| {
            p.is_file()
                && p.extension().and_then(|e| e.to_str()).is_some_and(|ext| {
                    args.ext
                        .iter()
                        .any(|allowed| allowed.eq_ignore_ascii_case(ext))
                })
        })
        .collect();
    files.sort();

    if files.len() < 2 {
        bail!(
            "found {} matching genome files in {} (need at least 2)",
            files.len(),
            args.dir.display()
        );
    }

    std::fs::create_dir_all(&args.outdir)?;
    let formats = OutputFormat::expand(&args.formats);

    let config = CompareConfig {
        min_length: args.min_length,
        min_similarity: args.min_similarity,
        word_length: args.word_length,
        freq_limit: args.freq_limit,
        memory_budget: args.memory,
        ..Default::default()
    };
    config.validate()?;

    // Load all genomes once; index one reference genome `j` at a time (one canonical RefIndex
    // instead of 2 KmerDicts per genome) to minimize peak RAM while reusing the index across
    // all queries `0..j`.
    let mut genomes = Vec::with_capacity(files.len());
    for path in &files {
        genomes.push(Genome::from_fasta_path(path)?);
    }

    let total_pairs = files.len() * (files.len() - 1) / 2;
    let mut done = 0usize;
    for j in 1..files.len() {
        // Check if any pair `(i, j)` actually needs running before indexing `genomes[j]`.
        let needed_queries: Vec<usize> = (0..j)
            .filter(|&i| {
                let prefix = format!("{}-{}", genomes[i].name, genomes[j].name);
                !(args.skip_existing && args.outdir.join(format!("{prefix}.frags")).exists())
            })
            .collect();

        if needed_queries.is_empty() {
            for i in 0..j {
                done += 1;
                let prefix = format!("{}-{}", genomes[i].name, genomes[j].name);
                eprintln!("[{done}/{total_pairs}] Skipping existing {prefix}");
            }
            continue;
        }

        let use_ref_index = RefIndex::fits_in_budget(
            &genomes[j],
            config.word_length,
            config.soft_mask,
            config.compare_reverse,
            config.memory_budget,
        );

        let cached_ref = if use_ref_index {
            let idx = RefIndex::build(
                &genomes[j],
                config.word_length,
                config.freq_limit,
                config.soft_mask,
            )?;
            let rev = if config.compare_reverse {
                Some(genomes[j].reverse_complement())
            } else {
                None
            };
            Some((idx, rev))
        } else {
            None
        };

        for i in 0..j {
            done += 1;
            let prefix = format!("{}-{}", genomes[i].name, genomes[j].name);
            if args.skip_existing && args.outdir.join(format!("{prefix}.frags")).exists() {
                eprintln!("[{done}/{total_pairs}] Skipping existing {prefix}");
                continue;
            }
            let res = if let Some((ref idx, ref rev)) = cached_ref {
                compare_with_ref_index(&genomes[i], &genomes[j], idx, rev.as_ref(), &config)?
            } else {
                compare_genomes(&genomes[i], &genomes[j], &config, None)?
            };
            let q_fname = files[i].file_name().and_then(|s| s.to_str()).unwrap_or("");
            let r_fname = files[j].file_name().and_then(|s| s.to_str()).unwrap_or("");
            let meta = res.metadata(q_fname, r_fname, &genomes[i], &genomes[j], &config);
            write_outputs(
                &args.outdir,
                &prefix,
                &formats,
                &genomes[i],
                &genomes[j],
                &meta,
                &res.fragments,
                Some(&res.coverage),
            )?;
            eprintln!(
                "[{done}/{total_pairs}] {prefix}: {} fragments",
                res.fragments.len()
            );
        }
    }

    eprintln!(
        "All-vs-all complete ({total_pairs} pairs) in {:.3}s -> {}",
        t0.elapsed().as_secs_f64(),
        args.outdir.display()
    );
    Ok(())
}

fn cmd_stats(args: StatsArgs) -> Result<()> {
    let ext = args
        .input
        .extension()
        .and_then(|s| s.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();

    if ext == "gidx" {
        let dict = KmerDict::load_from_path(&args.input)?;
        let st = dict.stats();
        if args.json {
            println!("{}", serde_json::to_string_pretty(&st)?);
        } else {
            println!("Dictionary : {}", args.input.display());
            println!("Word size  : {}", st.index_k);
            println!("Unique k-m : {}", st.unique_words);
            println!("Occurrences: {}", st.total_occurrences);
            println!("Max freq   : {}", st.max_frequency);
        }
        return Ok(());
    }

    let ff = load_frags_or_csv(&args.input)?;
    let n = ff.fragments.len();
    let fwd = ff.fragments.iter().filter(|f| f.strand == 'f').count();
    let rev = n - fwd;
    let total_len: u64 = ff.fragments.iter().map(|f| f.length).sum();
    let max_len = ff.fragments.iter().map(|f| f.length).max().unwrap_or(0);
    let min_len = ff.fragments.iter().map(|f| f.length).min().unwrap_or(0);
    let mean_len = if n > 0 {
        total_len as f64 / n as f64
    } else {
        0.0
    };
    let mean_sim = if n > 0 {
        ff.fragments
            .iter()
            .map(|f| f.similarity as f64)
            .sum::<f64>()
            / n as f64
    } else {
        0.0
    };

    if args.json {
        let obj = serde_json::json!({
            "file": args.input,
            "x_total": ff.x_total,
            "y_total": ff.y_total,
            "fragments": n,
            "forward": fwd,
            "reverse": rev,
            "min_length": min_len,
            "max_length": max_len,
            "mean_length": mean_len,
            "mean_similarity": mean_sim,
        });
        println!("{}", serde_json::to_string_pretty(&obj)?);
    } else {
        println!("File            : {}", args.input.display());
        println!("SeqX length     : {}", ff.x_total);
        println!("SeqY length     : {}", ff.y_total);
        println!("Fragments       : {} (fwd={}, rev={})", n, fwd, rev);
        println!(
            "Length (min/avg/max): {} / {:.1} / {}",
            min_len, mean_len, max_len
        );
        println!("Mean similarity : {:.2}%", mean_sim);
    }

    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn write_outputs(
    outdir: &Path,
    prefix: &str,
    formats: &BTreeSet<OutputFormat>,
    query: &Genome,
    reference: &Genome,
    meta: &ComparisonMetadata,
    fragments: &[Fragment],
    coverage: Option<&gecko_core::CoverageMatrix>,
) -> Result<()> {
    if formats.contains(&OutputFormat::Frags) {
        write_frags_path(
            outdir.join(format!("{prefix}.frags")),
            query.total_bases,
            reference.total_bases,
            fragments,
        )?;
    }
    if formats.contains(&OutputFormat::Inf) {
        write_inf_path(outdir.join(format!("{prefix}.frags.INF")), meta)?;
    }
    if formats.contains(&OutputFormat::Mat)
        && let Some(cov) = coverage
    {
        write_mat_path(outdir.join(format!("{prefix}.frags.MAT")), cov)?;
    }
    if formats.contains(&OutputFormat::Csv) {
        write_csv_path(outdir.join(format!("{prefix}.csv")), meta, fragments)?;
    }
    if formats.contains(&OutputFormat::Paf) {
        write_paf_path(
            outdir.join(format!("{prefix}.paf")),
            query,
            reference,
            fragments,
        )?;
    }
    if formats.contains(&OutputFormat::Blast) {
        write_blast_tabular_path(
            outdir.join(format!("{prefix}.blast.tsv")),
            query,
            reference,
            fragments,
        )?;
    }
    if formats.contains(&OutputFormat::Bedpe) {
        write_bedpe_path(
            outdir.join(format!("{prefix}.bedpe")),
            query,
            reference,
            fragments,
        )?;
    }
    if formats.contains(&OutputFormat::Jsonl) {
        write_jsonl_path(
            outdir.join(format!("{prefix}.jsonl")),
            query,
            reference,
            fragments,
        )?;
    }
    if formats.contains(&OutputFormat::Png) {
        let plot_opts = PlotOptions::default().with_genomes(query, reference);
        write_dotplot_png_path(
            outdir.join(format!("{prefix}.png")),
            query.total_bases,
            reference.total_bases,
            fragments,
            &plot_opts,
        )?;
    }
    Ok(())
}

fn cmd_plot(args: PlotArgs) -> Result<()> {
    let t0 = Instant::now();
    let ff = load_frags_or_csv(&args.input)?;

    let fragments: Vec<Fragment> = ff
        .fragments
        .into_iter()
        .filter(|f| {
            args.min_length.is_none_or(|ml| f.length >= ml)
                && args.min_similarity.is_none_or(|ms| f.similarity >= ms)
        })
        .collect();

    let mut opts = PlotOptions {
        width: args.width,
        height: args.height,
        line_width: args.line_width,
        color_mode: match args.color {
            PlotColorArg::Strand => ColorMode::Strand,
            PlotColorArg::Similarity => ColorMode::Similarity,
        },
        origin: match args.origin {
            PlotOriginArg::TopLeft => PlotOrigin::TopLeft,
            PlotOriginArg::BottomLeft => PlotOrigin::BottomLeft,
        },
        dark: args.dark,
        title: args.title,
        ..Default::default()
    };

    if let (Some(qp), Some(rp)) = (args.query.as_ref(), args.reference.as_ref()) {
        let q = Genome::from_fasta_path(qp)?;
        let r = Genome::from_fasta_path(rp)?;
        opts = opts.with_genomes(&q, &r);
    } else if opts.title.is_none() {
        let stem = args
            .input
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("GECKO Dotplot");
        opts.title = Some(stem.to_string());
    }

    if let Some(parent) = args.output.parent()
        && !parent.as_os_str().is_empty()
    {
        std::fs::create_dir_all(parent)?;
    }

    write_dotplot_png_path(&args.output, ff.x_total, ff.y_total, &fragments, &opts)?;
    eprintln!(
        "Rendered {} fragments ({}x{}) in {:.3}s -> {}",
        fragments.len(),
        opts.width,
        opts.height,
        t0.elapsed().as_secs_f64(),
        args.output.display()
    );
    Ok(())
}

fn load_frags_or_csv(path: &Path) -> Result<FragsFile> {
    let is_csv = path
        .extension()
        .and_then(|s| s.to_str())
        .is_some_and(|e| e.eq_ignore_ascii_case("csv"));
    if is_csv {
        Ok(read_csv_path(path)?)
    } else {
        Ok(read_frags_path(path)?)
    }
}

fn parse_byte_size(s: &str) -> std::result::Result<usize, String> {
    let trimmed = s.trim();
    if trimmed.is_empty() {
        return Err("empty memory size".into());
    }
    let idx = trimmed
        .find(|c: char| !c.is_ascii_digit())
        .unwrap_or(trimmed.len());
    let num: usize = trimmed[..idx]
        .parse()
        .map_err(|_| format!("invalid memory size '{s}'"))?;
    let unit = trimmed[idx..].trim().to_ascii_uppercase();
    let mult: usize = match unit.as_str() {
        "" | "B" => 1,
        "K" | "KB" | "KIB" => 1024,
        "M" | "MB" | "MIB" => 1024 * 1024,
        "G" | "GB" | "GIB" => 1024 * 1024 * 1024,
        _ => return Err(format!("unknown memory unit '{unit}' in '{s}'")),
    };
    num.checked_mul(mult)
        .ok_or_else(|| format!("memory size overflow: '{s}'"))
}
