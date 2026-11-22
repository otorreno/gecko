//! Shared MetaGECKO CLI arguments and execution pipeline (used by both `gecko meta` and
//! the standalone `metagecko` binary).

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{Context, Result, bail};
use clap::{Parser, ValueEnum};
use gecko_core::Scoring;
use gecko_core::formats::csv::{
    ComparisonMetadata, write_csv_path, write_inf_path, write_mat_path,
};
use gecko_core::formats::frags::write_frags_path;
use gecko_core::meta::{
    MetaConfig, MetaReference, MetaSampleAccumulator, MetaSampleResult, MultiMappingMode,
    collect_sequence_files, write_abundance_json, write_abundance_matrix_tsv, write_abundance_tsv,
    write_read_assignments_tsv,
};
use gecko_core::pipeline::CompareConfig;
use gecko_core::seq::{ReadBatchReader, file_stem_name};

use crate::backend::{BackendArg, PreparedRefSession};

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum MultiMapArg {
    /// Split read weight 1/N across tied best-scoring genomes while tracking unique reads.
    Fractional,
    /// Keep only the single highest-scoring alignment per read.
    Best,
    /// Retain all alignments passing thresholds (exhaustive fine-grained MetaGECKO mode).
    All,
}

impl From<MultiMapArg> for MultiMappingMode {
    fn from(v: MultiMapArg) -> Self {
        match v {
            MultiMapArg::Fractional => MultiMappingMode::Fractional,
            MultiMapArg::Best => MultiMappingMode::Best,
            MultiMapArg::All => MultiMappingMode::All,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
pub enum MetaOutputFormat {
    /// Per-genome taxonomic abundance, breadth/depth, RPKM & Lander-Waterman evenness (.abundance.tsv)
    Abundance,
    /// Machine-readable JSON sample summary + genome profiles (.abundance.json)
    Json,
    /// Per-read alignment assignment table (.reads.tsv)
    Reads,
    /// GECKO-MGV compatible CSV fragments (.csv)
    Csv,
    /// Legacy GECKO 93-byte binary fragments (.frags + .INF + .MAT)
    Frags,
    /// Preset: abundance + json
    Default,
    /// Preset: all MetaGECKO and legacy GECKO-MGV outputs
    All,
}

impl std::fmt::Display for MetaOutputFormat {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = match self {
            Self::Abundance => "abundance",
            Self::Json => "json",
            Self::Reads => "reads",
            Self::Csv => "csv",
            Self::Frags => "frags",
            Self::Default => "default",
            Self::All => "all",
        };
        f.write_str(s)
    }
}

impl MetaOutputFormat {
    fn expand(formats: &[Self]) -> BTreeSet<Self> {
        let mut set = BTreeSet::new();
        for &f in formats {
            match f {
                Self::Default => {
                    set.insert(Self::Abundance);
                    set.insert(Self::Json);
                }
                Self::All => {
                    set.insert(Self::Abundance);
                    set.insert(Self::Json);
                    set.insert(Self::Reads);
                    set.insert(Self::Csv);
                    set.insert(Self::Frags);
                }
                other => {
                    set.insert(other);
                }
            }
        }
        set
    }
}

#[derive(Parser, Debug)]
pub struct MetaArgs {
    /// Metagenomic read sample file(s) or directory(ies) (FASTA, FASTQ, or .gz).
    #[arg(required = true)]
    pub reads: Vec<PathBuf>,

    /// Reference genome FASTA file(s) or directory(ies) to profile against.
    #[arg(short = 'r', long = "refs", required = true, num_args = 1..)]
    pub refs: Vec<PathBuf>,

    /// Treat each FASTA record inside a single reference file as a distinct genome/taxon.
    #[arg(long)]
    pub split_records: bool,

    /// Minimum HSP length in bp.
    #[arg(short = 'l', long, default_value_t = 40)]
    pub min_length: u64,

    /// Minimum HSP similarity percentage threshold (0..=100).
    #[arg(short = 's', long, default_value_t = 70.0)]
    pub min_similarity: f32,

    /// Minimum fraction of the read length covered by the HSP (0..=100, filters spurious short hits).
    #[arg(long, default_value_t = 0.0)]
    pub min_read_cov: f32,

    /// Minimum exact base identity percentage (0..=100).
    #[arg(long, default_value_t = 0.0)]
    pub min_identity: f32,

    /// Skip input reads shorter than this length in bp.
    #[arg(long, default_value_t = 0)]
    pub min_read_len: u64,

    /// Seed k-mer word length (4..=32; GPU backends require 8..=32).
    #[arg(short = 'k', long, default_value_t = 16)]
    pub word_length: usize,

    /// Skip reference k-mers occurring more than this many times per strand.
    #[arg(long, default_value_t = 1000)]
    pub freq_limit: u64,

    /// Enable soft-masking (skip lowercase seeds in reference/reads). Off by default for metagenomes.
    #[arg(long)]
    pub soft_mask: bool,

    /// Compare only the forward strand (skip reverse complement of reference).
    #[arg(long)]
    pub forward_only: bool,

    /// Multi-mapping read resolution mode.
    #[arg(long, value_enum, default_value_t = MultiMapArg::Fractional)]
    pub multi_map: MultiMapArg,

    /// Score margin below a read's best score still considered tied in fractional mode.
    #[arg(long, default_value_t = 0)]
    pub tie_margin: u64,

    /// Maximum base pairs per streaming read batch (e.g. 32M, 64M, 256M).
    #[arg(long, default_value = "64M", value_parser = parse_meta_byte_size)]
    pub batch_size: usize,

    /// Memory budget for CPU reference indexing before slab-streaming fallback (e.g. 2G, 8G).
    #[arg(short = 'm', long, default_value = "4G", value_parser = parse_meta_byte_size)]
    pub memory: usize,

    /// Output directory.
    #[arg(short = 'o', long, default_value = "metagecko_results")]
    pub outdir: PathBuf,

    /// Output formats (comma-separated or repeated).
    #[arg(
        short = 'f',
        long = "format",
        value_delimiter = ',',
        default_values_t = [MetaOutputFormat::Default]
    )]
    pub formats: Vec<MetaOutputFormat>,

    /// Match score for ungapped extension.
    #[arg(long, default_value_t = 4)]
    pub match_score: i64,

    /// Mismatch penalty (positive integer subtracted on mismatch).
    #[arg(long, default_value_t = 4)]
    pub mismatch_penalty: i64,

    /// Penalty for aligning an 'N' or ambiguous base.
    #[arg(long, default_value_t = 1)]
    pub n_penalty: i64,

    /// Compute backend (`auto` prefers CUDA -> WGPU -> AVX2 CPU).
    #[arg(short = 'b', long, value_enum, default_value_t = BackendArg::Auto)]
    pub backend: BackendArg,

    /// Allow software WGPU adapters (e.g. Mesa Lavapipe) when using `--backend wgpu|auto`.
    #[arg(long)]
    pub allow_software_gpu: bool,
}

pub fn run_meta(args: MetaArgs) -> Result<()> {
    let t_total = Instant::now();

    // 1. Resolve metagenome read sample files
    let mut read_files: Vec<PathBuf> = Vec::new();
    for p in &args.reads {
        if p.is_dir() {
            collect_sequence_files(p, &mut read_files)?;
        } else {
            read_files.push(p.clone());
        }
    }
    read_files.sort();
    if read_files.is_empty() {
        bail!("no metagenome read files found");
    }

    let formats = MetaOutputFormat::expand(&args.formats);
    let collect_assignments = formats.contains(&MetaOutputFormat::Reads);
    let collect_fragments =
        formats.contains(&MetaOutputFormat::Csv) || formats.contains(&MetaOutputFormat::Frags);

    let meta_cfg = MetaConfig {
        compare: CompareConfig {
            min_length: args.min_length,
            min_similarity: args.min_similarity,
            word_length: args.word_length,
            freq_limit: args.freq_limit,
            max_rep: 10_000,
            fixed_length: true,
            soft_mask: args.soft_mask,
            memory_budget: args.memory,
            compare_reverse: !args.forward_only,
            scoring: Scoring {
                match_score: args.match_score,
                mismatch_penalty: args.mismatch_penalty,
                n_penalty: args.n_penalty,
            },
        },
        min_read_coverage: args.min_read_cov,
        min_identity: args.min_identity,
        min_read_len: args.min_read_len,
        multi_mapping: args.multi_map.into(),
        tie_score_margin: args.tie_margin,
        batch_bases: args.batch_size,
        collect_assignments,
        collect_fragments,
    };
    meta_cfg.compare.validate()?;

    // 2. Load unified reference database & prepare reusable CPU/CUDA/WGPU index ONCE
    let t_ref = Instant::now();
    let meta_ref = MetaReference::load_from_paths(&args.refs, args.split_records)
        .context("failed to load reference genome database")?;
    let ref_load_s = t_ref.elapsed().as_secs_f64();

    let t_prep = Instant::now();
    let session = PreparedRefSession::prepare(
        &meta_ref.genome,
        &meta_cfg.compare,
        args.backend,
        args.allow_software_gpu,
    )?;
    let prep_s = t_prep.elapsed().as_secs_f64();

    eprintln!(
        "[metagecko] Reference DB: {} genomes ({} contigs, {} bp) loaded in {:.3}s, indexed on [{}] in {:.3}s",
        meta_ref.num_genomes(),
        meta_ref.genome.num_records(),
        meta_ref.genome.total_bases,
        ref_load_s,
        session.backend_label(),
        prep_s
    );

    std::fs::create_dir_all(&args.outdir)?;
    let mut sample_results: Vec<MetaSampleResult> = Vec::with_capacity(read_files.len());

    // 3. Stream each metagenome sample in bounded batches through the prepared reference session
    for read_path in &read_files {
        let sample_name = file_stem_name(read_path);
        let t_sample = Instant::now();
        let mut batch_reader =
            ReadBatchReader::open(read_path, meta_cfg.batch_bases, meta_cfg.min_read_len)
                .with_context(|| format!("failed to open reads '{}'", read_path.display()))?;

        let mut acc = MetaSampleAccumulator::new(&sample_name, &meta_ref);
        let mut query_secs = 0.0f64;

        while let Some(batch) = batch_reader.next_batch()? {
            let t_q = Instant::now();
            let cmp = session.query_batch(&batch.genome, &meta_cfg.compare)?;
            query_secs += t_q.elapsed().as_secs_f64();
            acc.ingest_batch(&batch, &meta_ref, cmp, &meta_cfg);
        }

        let res = acc.finish(&meta_ref);
        let sample_s = t_sample.elapsed().as_secs_f64();
        let reads_per_s = if sample_s > 0.0 {
            (res.reads_analyzed as f64) / sample_s
        } else {
            0.0
        };
        let mapped_pct = if res.reads_analyzed > 0 {
            (res.mapped_reads as f64) * 100.0 / (res.reads_analyzed as f64)
        } else {
            0.0
        };
        let detected_taxa = res
            .profiles
            .iter()
            .filter(|p| p.assigned_reads > 0.0)
            .count();

        eprintln!(
            "[metagecko] Sample '{}': {}/{} reads mapped ({:.2}%, unique={}, multi={}), {} taxa detected, {} batches | {:.3}s ({:.0} reads/s, compute={:.3}s)",
            res.sample_name,
            res.mapped_reads,
            res.reads_analyzed,
            mapped_pct,
            res.unique_mapped_reads,
            res.multi_mapped_reads,
            detected_taxa,
            res.batches_processed,
            sample_s,
            reads_per_s,
            query_secs
        );

        write_sample_outputs(
            &args.outdir,
            read_path,
            &meta_ref,
            &res,
            &meta_cfg,
            &formats,
        )?;
        sample_results.push(res);
    }

    // 4. Multi-sample comparative abundance matrices (G x M) when >= 2 samples
    if sample_results.len() >= 2 {
        let rel_path = args.outdir.join("abundance_matrix.rel.tsv");
        let cnt_path = args.outdir.join("abundance_matrix.reads.tsv");
        write_abundance_matrix_tsv(&rel_path, &meta_ref, &sample_results, true)?;
        write_abundance_matrix_tsv(&cnt_path, &meta_ref, &sample_results, false)?;
        eprintln!(
            "[metagecko] Wrote multi-sample abundance matrices -> {}, {}",
            rel_path.display(),
            cnt_path.display()
        );
    }

    eprintln!(
        "[metagecko] Done ({} sample(s) x {} genome(s)) in {:.3}s -> {}",
        sample_results.len(),
        meta_ref.num_genomes(),
        t_total.elapsed().as_secs_f64(),
        args.outdir.display()
    );

    Ok(())
}

fn write_sample_outputs(
    outdir: &Path,
    read_path: &Path,
    meta_ref: &MetaReference,
    res: &MetaSampleResult,
    cfg: &MetaConfig,
    formats: &BTreeSet<MetaOutputFormat>,
) -> Result<()> {
    let prefix = &res.sample_name;
    if formats.contains(&MetaOutputFormat::Abundance) {
        write_abundance_tsv(outdir.join(format!("{prefix}.abundance.tsv")), res)?;
    }
    if formats.contains(&MetaOutputFormat::Json) {
        write_abundance_json(outdir.join(format!("{prefix}.abundance.json")), res)?;
    }
    if formats.contains(&MetaOutputFormat::Reads) {
        write_read_assignments_tsv(outdir.join(format!("{prefix}.reads.tsv")), res)?;
    }
    if formats.contains(&MetaOutputFormat::Csv) || formats.contains(&MetaOutputFormat::Frags) {
        let meta = ComparisonMetadata {
            seq_x_filename: read_path.display().to_string(),
            seq_y_filename: meta_ref.genome.name.clone(),
            seq_x_name: res.sample_name.clone(),
            seq_y_name: meta_ref.genome.primary_id().to_string(),
            seq_x_len: res.total_read_bases,
            seq_y_len: meta_ref.genome.total_bases,
            min_length: cfg.compare.min_length,
            min_similarity: cfg.compare.min_similarity.round() as u64,
            total_hits: res.total_hits,
            total_hits_used: res.hits_used,
            total_fragments: res.fragments.len() as u64,
        };
        if formats.contains(&MetaOutputFormat::Csv) {
            write_csv_path(outdir.join(format!("{prefix}.csv")), &meta, &res.fragments)?;
        }
        if formats.contains(&MetaOutputFormat::Frags) {
            write_frags_path(
                outdir.join(format!("{prefix}.frags")),
                res.total_read_bases,
                meta_ref.genome.total_bases,
                &res.fragments,
            )?;
            write_inf_path(outdir.join(format!("{prefix}.frags.INF")), &meta)?;
            write_mat_path(outdir.join(format!("{prefix}.frags.MAT")), &res.coverage)?;
        }
    }
    Ok(())
}

fn parse_meta_byte_size(s: &str) -> std::result::Result<usize, String> {
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
