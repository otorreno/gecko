use std::path::Path;
use std::sync::OnceLock;
use std::time::Instant;

use serde::Serialize;

use crate::error::{GeckoError, Result};
use crate::extend::{CoverageMatrix, ExtendOptions, Scoring, extend_hits_parallel};
use crate::formats::csv::ComparisonMetadata;
use crate::formats::frags::Fragment;
use crate::hits::{HitOptions, generate_and_consume_hits, generate_and_consume_hits_from_index};
use crate::kmer::{DictOptions, INDEX_WORD_SIZE, KmerDict, KmerIndex};
use crate::minimizer::{MINIMIZER_MIN_K, compare_genomes_minimizer_with_coverage};
use crate::ref_index::{RefIndex, compute_partition_slabs_with_mask_and_strand};
use crate::seq::Genome;
use crate::spaced::SpacedSeed;
use crate::stream::{compare_streaming_with_coverage, dedup_fragments};

/// Return `true` when `GECKO_PROFILE` is set to a truthy value (cached once per process for zero
/// hot-path overhead when unset).
#[inline]
pub fn profile_enabled() -> bool {
    static ENABLED: OnceLock<bool> = OnceLock::new();
    *ENABLED.get_or_init(|| match std::env::var("GECKO_PROFILE") {
        Ok(v) => {
            let trimmed = v.trim();
            !trimmed.is_empty()
                && trimmed != "0"
                && !trimmed.eq_ignore_ascii_case("false")
                && !trimmed.eq_ignore_ascii_case("no")
                && !trimmed.eq_ignore_ascii_case("off")
        }
        Err(_) => false,
    })
}

/// CPU/CUDA seeding mechanism selection (`--seeding exact|minimizer|auto`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SeedingMode {
    /// Stride-1 exact canonical `RefIndex` seeding (Direct-CSR / Hash-CSR).
    Exact,
    /// Shift-invariant canonical minimizer super-seeding ($k' < k, w = k - k' + 1$) with exact-$k$
    /// verification and exact `freq_limit` resolution for contiguous $k \ge 18$ (falls back to
    /// `Exact` for $k < 18$, spaced seeds, `--two-hit`, or `.gidx` cache runs).
    Minimizer,
    /// Automatically select `Minimizer` for contiguous $k \ge 18$ (without spaced seeds,
    /// `--two-hit`, or `.gidx` caching) and `Exact` otherwise.
    #[default]
    Auto,
}

/// Full configuration for a pairwise GECKO comparison.
#[derive(Debug, Clone)]
pub struct CompareConfig {
    /// Minimum HSP length (in bp if `fixed_length`, or percentage if `!fixed_length`).
    pub min_length: u64,
    /// Minimum similarity percentage threshold (`0.0 ..= 100.0`).
    pub min_similarity: f32,
    /// Seed k-mer length (`4..=32`; no longer limited to multiples of 4).
    /// When `spaced_seed` is set, this must equal `spaced_seed.span`.
    pub word_length: usize,
    /// Optional palindromic spaced-seed pattern (`None` for contiguous seeds).
    pub spaced_seed: Option<SpacedSeed>,
    /// Skip reference words occurring more than `freq_limit` times per strand (default `1000`).
    pub freq_limit: u64,
    /// Subsampling cap per word list (default `10_000`).
    pub max_rep: usize,
    /// If true, `min_length` is in bp; if false, percentage of shorter sequence length.
    pub fixed_length: bool,
    /// If true (default), lowercase bases do not seed k-mers (soft-masking) but can extend.
    pub soft_mask: bool,
    /// Memory budget in bytes for indexing/sorting before falling back to slab streaming or disk runs.
    pub memory_budget: usize,
    /// Whether to compare the reverse complement of the reference (`true` by default).
    pub compare_reverse: bool,
    /// Ungapped alignment scoring parameters.
    pub scoring: Scoring,
    /// Seeding mechanism (`Exact`, `Minimizer`, or `Auto`).
    pub seeding_mode: SeedingMode,
    /// Optional BLAST/LASTZ two-hit diagonal window in bp (`None` = standard single-hit exact mode).
    pub two_hit_window: Option<u32>,
    /// Whether to populate the 1000x100 `.MAT` [`CoverageMatrix`] (default `true` for API callers;
    /// CLI sets `false` when `--format mat` is not requested).
    pub collect_coverage: bool,
}

impl Default for CompareConfig {
    fn default() -> Self {
        Self {
            min_length: 40,
            min_similarity: 60.0,
            word_length: 12,
            spaced_seed: None,
            freq_limit: 1000,
            max_rep: 10_000,
            fixed_length: true,
            soft_mask: true,
            memory_budget: 2 * 1024 * 1024 * 1024,
            compare_reverse: true,
            scoring: Scoring::default(),
            seeding_mode: SeedingMode::Auto,
            two_hit_window: None,
            collect_coverage: true,
        }
    }
}

impl CompareConfig {
    /// Return the active 2-bit spaced-seed mask if a gapped spaced seed is configured, else `None`.
    #[inline]
    pub fn spaced_mask_2bit(&self) -> Option<u64> {
        self.spaced_seed.as_ref().and_then(|s| s.gapped_mask_2bit())
    }

    /// Return `true` when this configuration should execute via the canonical minimizer engine.
    #[inline]
    pub fn use_minimizer_seeding(&self, cache_dir: Option<&Path>) -> bool {
        matches!(
            self.seeding_mode,
            SeedingMode::Minimizer | SeedingMode::Auto
        ) && self.word_length >= MINIMIZER_MIN_K
            && self.spaced_mask_2bit().is_none()
            && self.two_hit_window.is_none()
            && cache_dir.is_none()
    }

    pub fn validate(&self) -> Result<()> {
        if !(4..=32).contains(&self.word_length) {
            return Err(GeckoError::InvalidParameter(format!(
                "word_length must be between 4 and 32, got {}",
                self.word_length
            )));
        }
        if let Some(ref sp) = self.spaced_seed
            && sp.span != self.word_length
        {
            return Err(GeckoError::InvalidParameter(format!(
                "spaced_seed span ({}) must match word_length ({})",
                sp.span, self.word_length
            )));
        }
        if !(0.0..=100.0).contains(&self.min_similarity) {
            return Err(GeckoError::InvalidParameter(format!(
                "min_similarity must be between 0 and 100, got {}",
                self.min_similarity
            )));
        }
        if self.scoring.match_score <= 0 {
            return Err(GeckoError::InvalidParameter(
                "match_score must be positive".into(),
            ));
        }
        if self.two_hit_window == Some(0) {
            return Err(GeckoError::InvalidParameter(
                "two_hit_window must be >= 1".into(),
            ));
        }
        Ok(())
    }

    /// Choose the dictionary index k-mer length: 32 for normal sequences (exact GECKO
    /// compatibility and multi-`word_length` index reuse), or `word_length` if a sequence
    /// is shorter than 32 bp.
    pub fn effective_index_k(&self, query: &Genome, reference: &Genome) -> usize {
        let max_q = query.records.iter().map(|r| r.len).max().unwrap_or(0) as usize;
        let max_r = reference.records.iter().map(|r| r.len).max().unwrap_or(0) as usize;
        if max_q < INDEX_WORD_SIZE || max_r < INDEX_WORD_SIZE {
            self.word_length
        } else {
            INDEX_WORD_SIZE
        }
    }
}

/// Summary statistics from a pairwise comparison run.
#[derive(Debug, Clone, Default, Serialize)]
pub struct CompareStats {
    pub word_matches: u64,
    pub raw_hits: u64,
    pub filtered_hits: u64,
    pub hits_used: u64,
    pub forward_fragments: u64,
    pub reverse_fragments: u64,
    pub spilled_runs: usize,
}

/// Output of a pairwise comparison.
#[derive(Debug, Clone)]
pub struct CompareResult {
    pub fragments: Vec<Fragment>,
    pub coverage: CoverageMatrix,
    pub stats: CompareStats,
}

impl CompareResult {
    /// Build legacy GECKO `.INF`/`.csv` metadata block.
    pub fn metadata(
        &self,
        query_filename: impl Into<String>,
        ref_filename: impl Into<String>,
        query: &Genome,
        reference: &Genome,
        config: &CompareConfig,
    ) -> ComparisonMetadata {
        ComparisonMetadata {
            seq_x_filename: query_filename.into(),
            seq_y_filename: ref_filename.into(),
            seq_x_name: query.primary_id().to_string(),
            seq_y_name: reference.primary_id().to_string(),
            seq_x_len: query.total_bases,
            seq_y_len: reference.total_bases,
            min_length: config.min_length,
            min_similarity: config.min_similarity.round() as u64,
            total_hits: self.stats.filtered_hits,
            total_hits_used: self.stats.hits_used,
            total_fragments: self.fragments.len() as u64,
        }
    }
}

/// Run a full pairwise genome comparison, optionally reusing cached `.gidx` dictionaries in
/// `cache_dir`.
///
/// Strategy when `cache_dir` is `None`:
/// 1. **Single-slab streaming** (`RefIndex` over all of `reference`) when it fits in `memory_budget`.
/// 2. **Reference-slab streaming** (sequential memory-bounded slabs of `reference`, zero temp disk
///    spilling) when `memory_budget` has at least 128 MiB net headroom above 2-bit `rev_ref`.
/// 3. **Out-of-core `KmerIndex::OnDisk`** fallback for ultra-tight `--memory` budgets or `.gidx`
///    disk caching.
pub fn compare_genomes(
    query: &Genome,
    reference: &Genome,
    config: &CompareConfig,
    cache_dir: Option<&Path>,
) -> Result<CompareResult> {
    compare_genomes_with_coverage(query, reference, config, cache_dir, config.collect_coverage)
}

/// Run a full pairwise genome comparison with optional `.MAT` [`CoverageMatrix`] collection
/// (`collect_coverage = false` avoids per-worker 800 KiB matrix allocations and rejected-hit
/// FP64 divisions while producing 100% bit-identical fragments and statistics).
pub fn compare_genomes_with_coverage(
    query: &Genome,
    reference: &Genome,
    config: &CompareConfig,
    cache_dir: Option<&Path>,
    collect_coverage: bool,
) -> Result<CompareResult> {
    config.validate()?;
    let eff_cov = collect_coverage && config.collect_coverage;

    if config.use_minimizer_seeding(cache_dir) {
        return compare_genomes_minimizer_with_coverage(query, reference, config, eff_cov);
    }

    let spaced_mask = config.spaced_mask_2bit();
    let streaming_only = spaced_mask.is_some() || config.two_hit_window.is_some();

    if cache_dir.is_none() || streaming_only {
        if RefIndex::fits_in_budget_with_mask(
            reference,
            config.word_length,
            spaced_mask,
            config.soft_mask,
            config.compare_reverse,
            config.memory_budget,
        ) {
            let forward_only = !config.compare_reverse
                && config.freq_limit.saturating_mul(2) <= config.max_rep.max(1) as u64;
            // Build `RefIndex` first on Main Thread (glibc Arena 0), then compute internally
            // parallel `reverse_complement()` after transient `part_keys` are freed so `rev_ref`
            // reuses freed heap space at zero peak-RSS and zero CPU-time cost.
            let t_idx = Instant::now();
            let ref_index = RefIndex::build_with_mask_and_strand(
                reference,
                config.word_length,
                spaced_mask,
                config.freq_limit,
                config.soft_mask,
                forward_only,
            )?;
            let idx_secs = t_idx.elapsed().as_secs_f64();

            let t_rev = Instant::now();
            let rev_ref = if config.compare_reverse {
                Some(reference.reverse_complement())
            } else {
                None
            };
            let rev_secs = t_rev.elapsed().as_secs_f64();

            if profile_enabled() {
                eprintln!(
                    "[GECKO_PROFILE] index_build slabs=1 kept_entries={} heap_mib={:.2} index_ms={:.2} rev_ms={:.2}",
                    ref_index.kept_entries,
                    ref_index.heap_bytes() as f64 / (1024.0 * 1024.0),
                    idx_secs * 1e3,
                    rev_secs * 1e3,
                );
            }

            return compare_streaming_with_coverage(
                query,
                reference,
                &ref_index,
                rev_ref.as_ref(),
                config,
                eff_cov,
            );
        }

        if let Some(max_kmers_per_slab) = RefIndex::max_kmers_per_slab_for_budget_with_mask(
            reference,
            config.word_length,
            spaced_mask,
            config.compare_reverse,
            config.memory_budget,
        ) {
            return compare_genomes_slab_streaming_with_coverage(
                query,
                reference,
                config,
                max_kmers_per_slab,
                eff_cov,
            );
        }

        if streaming_only {
            // Spaced seeds and `--two-hit` use Hash-CSR/Direct-CSR slab streaming even under tight
            // memory budgets.
            return compare_genomes_slab_streaming_with_coverage(
                query, reference, config, 4_000_000, eff_cov,
            );
        }
    }

    compare_genomes_kmer_index(query, reference, config, cache_dir)
}

/// Zero-spill Hash-Partitioned Reference-Slab Streaming: partitions the 256 coarse k-mer
/// partitions (`0..256`) into disjoint contiguous ranges `[part_start .. part_end)` each holding
/// at most `max_kmers_per_slab` valid k-mers. Because every occurrence of a canonical k-mer across
/// the entire reference maps to the same coarse partition, per-strand `freq_limit` filtering is
/// **100% globally exact** (bit-identical to single-slab runs), query streaming skips `(S - 1) / S`
/// of k-mers before touching CSR tables, and full 2-bit sequences stay in RAM for exact extension.
#[cfg(test)]
pub(crate) fn compare_genomes_slab_streaming(
    query: &Genome,
    reference: &Genome,
    config: &CompareConfig,
    max_kmers_per_slab: usize,
) -> Result<CompareResult> {
    compare_genomes_slab_streaming_with_coverage(query, reference, config, max_kmers_per_slab, true)
}

pub(crate) fn compare_genomes_slab_streaming_with_coverage(
    query: &Genome,
    reference: &Genome,
    config: &CompareConfig,
    max_kmers_per_slab: usize,
    collect_coverage: bool,
) -> Result<CompareResult> {
    let spaced_mask = config.spaced_mask_2bit();
    let forward_only = !config.compare_reverse
        && config.freq_limit.saturating_mul(2) <= config.max_rep.max(1) as u64;
    let slabs = compute_partition_slabs_with_mask_and_strand(
        reference,
        config.word_length,
        spaced_mask,
        config.soft_mask,
        max_kmers_per_slab,
        forward_only,
    );
    if slabs.is_empty() {
        return Ok(CompareResult {
            fragments: Vec::new(),
            coverage: CoverageMatrix::new(),
            stats: CompareStats::default(),
        });
    }

    // Compute internally parallel `rev_ref` once before the slab loop, then build EVERY slab
    // sequentially on Main Thread (glibc Arena 0) so freed per-slab `part_keys` are reused across
    // slabs instead of stranding Slab 0 in a Rayon worker arena.
    let t_rev = Instant::now();
    let rev_ref = if config.compare_reverse {
        Some(reference.reverse_complement())
    } else {
        None
    };
    let rev_secs = t_rev.elapsed().as_secs_f64();
    let n_slabs = slabs.len();
    let build_unique_bm = n_slabs == 1;

    let mut merged = CompareResult {
        fragments: Vec::new(),
        coverage: CoverageMatrix::new(),
        stats: CompareStats::default(),
    };
    let mut total_idx_secs = 0.0f64;

    for (part_start, part_end) in slabs {
        let t_idx = Instant::now();
        let slab_index = RefIndex::build_partitions_with_options(
            reference,
            part_start,
            part_end,
            config.word_length,
            spaced_mask,
            config.freq_limit,
            config.soft_mask,
            forward_only,
            build_unique_bm,
        )?;
        total_idx_secs += t_idx.elapsed().as_secs_f64();
        let mut slab_res = compare_streaming_with_coverage(
            query,
            reference,
            &slab_index,
            rev_ref.as_ref(),
            config,
            collect_coverage,
        )?;
        drop(slab_index);

        merged.stats.word_matches += slab_res.stats.word_matches;
        merged.stats.raw_hits += slab_res.stats.raw_hits;
        merged.stats.filtered_hits += slab_res.stats.filtered_hits;
        merged.stats.hits_used += slab_res.stats.hits_used;
        if collect_coverage {
            merged.coverage.merge_from(&slab_res.coverage);
        }
        merged.fragments.append(&mut slab_res.fragments);
    }

    dedup_fragments(&mut merged.fragments);
    for f in &merged.fragments {
        if f.strand == 'f' {
            merged.stats.forward_fragments += 1;
        } else {
            merged.stats.reverse_fragments += 1;
        }
    }

    if profile_enabled() {
        eprintln!(
            "[GECKO_PROFILE] index_multislab slabs={} index_ms={:.2} rev_ms={:.2} frags={}",
            n_slabs,
            total_idx_secs * 1e3,
            rev_secs * 1e3,
            merged.fragments.len(),
        );
    }

    Ok(merged)
}

/// Out-of-core / `.gidx`-cached comparison using legacy sorted [`KmerIndex`] dictionaries.
pub(crate) fn compare_genomes_kmer_index(
    query: &Genome,
    reference: &Genome,
    config: &CompareConfig,
    cache_dir: Option<&Path>,
) -> Result<CompareResult> {
    config.validate()?;
    let index_k = config.effective_index_k(query, reference);
    let dict_opts = DictOptions {
        soft_mask: config.soft_mask,
        memory_budget: config.memory_budget,
        index_k,
    };

    let mut result = CompareResult {
        fragments: Vec::new(),
        coverage: CoverageMatrix::new(),
        stats: CompareStats::default(),
    };

    let idx_x = load_or_build_index(query, dict_opts, cache_dir)?;

    {
        let idx_y = load_or_build_index(reference, dict_opts, cache_dir)?;
        run_strand_index(query, &idx_x, reference, &idx_y, 'f', config, &mut result)?;
    }

    if config.compare_reverse {
        let rev_ref = reference.reverse_complement();
        let idx_y_rev = load_or_build_index(&rev_ref, dict_opts, cache_dir)?;
        run_strand_index(
            query,
            &idx_x,
            &rev_ref,
            &idx_y_rev,
            'r',
            config,
            &mut result,
        )?;
    }

    Ok(result)
}

/// Compare a query genome against a pre-built canonical [`RefIndex`] (ideal for all-vs-all
/// and multi-query workflows).
pub fn compare_with_ref_index(
    query: &Genome,
    reference: &Genome,
    ref_index: &RefIndex,
    rev_ref: Option<&Genome>,
    config: &CompareConfig,
) -> Result<CompareResult> {
    compare_with_ref_index_with_coverage(
        query,
        reference,
        ref_index,
        rev_ref,
        config,
        config.collect_coverage,
    )
}

/// Compare a query genome against a pre-built [`RefIndex`] with optional `.MAT` coverage
/// matrix collection.
pub fn compare_with_ref_index_with_coverage(
    query: &Genome,
    reference: &Genome,
    ref_index: &RefIndex,
    rev_ref: Option<&Genome>,
    config: &CompareConfig,
    collect_coverage: bool,
) -> Result<CompareResult> {
    compare_streaming_with_coverage(
        query,
        reference,
        ref_index,
        rev_ref,
        config,
        collect_coverage,
    )
}

/// Compare two genomes using pre-built [`KmerDict`]s (ideal for `.gidx`-cached workflows).
pub fn compare_with_dicts(
    query: &Genome,
    dict_x: &KmerDict,
    reference: &Genome,
    dict_y: &KmerDict,
    rev_ref: Option<&Genome>,
    dict_y_rev: Option<&KmerDict>,
    config: &CompareConfig,
) -> Result<CompareResult> {
    config.validate()?;

    let mut result = CompareResult {
        fragments: Vec::new(),
        coverage: CoverageMatrix::new(),
        stats: CompareStats::default(),
    };

    run_strand(query, dict_x, reference, dict_y, 'f', config, &mut result)?;

    if config.compare_reverse
        && let (Some(rr), Some(dr)) = (rev_ref, dict_y_rev)
    {
        run_strand(query, dict_x, rr, dr, 'r', config, &mut result)?;
    }

    Ok(result)
}

fn run_strand(
    query: &Genome,
    dict_x: &KmerDict,
    target: &Genome,
    dict_y: &KmerDict,
    strand: char,
    config: &CompareConfig,
    acc: &mut CompareResult,
) -> Result<()> {
    let hit_opts = HitOptions {
        word_length: config.word_length,
        freq_limit: config.freq_limit,
        max_rep: config.max_rep,
        memory_budget: config.memory_budget,
    };
    let ext_opts = ExtendOptions {
        min_length: config.min_length,
        min_similarity: config.min_similarity,
        word_length: config.word_length,
        fixed_length: config.fixed_length,
        strand,
        scoring: config.scoring,
        spaced_seed: false,
    };

    let mut strand_frags: u64 = 0;
    let mut hits_used: u64 = 0;

    let hit_stats = generate_and_consume_hits(dict_x, dict_y, query, target, hit_opts, |batch| {
        let mut ext = extend_hits_parallel(batch, query, target, ext_opts);
        strand_frags += ext.fragments.len() as u64;
        hits_used += ext.hits_used;
        acc.coverage.merge_from(&ext.coverage);
        acc.fragments.append(&mut ext.fragments);
    })?;

    record_strand_stats(acc, hit_stats, hits_used, strand, strand_frags);
    Ok(())
}

fn run_strand_index(
    query: &Genome,
    idx_x: &KmerIndex,
    target: &Genome,
    idx_y: &KmerIndex,
    strand: char,
    config: &CompareConfig,
    acc: &mut CompareResult,
) -> Result<()> {
    let hit_opts = HitOptions {
        word_length: config.word_length,
        freq_limit: config.freq_limit,
        max_rep: config.max_rep,
        memory_budget: config.memory_budget,
    };
    let ext_opts = ExtendOptions {
        min_length: config.min_length,
        min_similarity: config.min_similarity,
        word_length: config.word_length,
        fixed_length: config.fixed_length,
        strand,
        scoring: config.scoring,
        spaced_seed: false,
    };

    let mut strand_frags: u64 = 0;
    let mut hits_used: u64 = 0;

    let hit_stats =
        generate_and_consume_hits_from_index(idx_x, idx_y, query, target, hit_opts, |batch| {
            let mut ext = extend_hits_parallel(batch, query, target, ext_opts);
            strand_frags += ext.fragments.len() as u64;
            hits_used += ext.hits_used;
            acc.coverage.merge_from(&ext.coverage);
            acc.fragments.append(&mut ext.fragments);
        })?;

    record_strand_stats(acc, hit_stats, hits_used, strand, strand_frags);
    Ok(())
}

fn record_strand_stats(
    acc: &mut CompareResult,
    hit_stats: crate::hits::HitStats,
    hits_used: u64,
    strand: char,
    strand_frags: u64,
) {
    acc.stats.word_matches += hit_stats.word_matches;
    acc.stats.raw_hits += hit_stats.raw_hits;
    acc.stats.filtered_hits += hit_stats.filtered_hits;
    acc.stats.hits_used += hits_used;
    acc.stats.spilled_runs += hit_stats.spilled_runs;
    if strand == 'f' {
        acc.stats.forward_fragments += strand_frags;
    } else {
        acc.stats.reverse_fragments += strand_frags;
    }
}

fn load_or_build_index(
    genome: &Genome,
    opts: DictOptions,
    cache_dir: Option<&Path>,
) -> Result<KmerIndex> {
    if let Some(dir) = cache_dir {
        std::fs::create_dir_all(dir)?;
        let mask_tag = if opts.soft_mask { "sm" } else { "nm" };
        let idx_path = dir.join(format!("{}.k{}.{mask_tag}.gidx", genome.name, opts.index_k));
        if idx_path.exists()
            && let Ok(loaded) = KmerDict::load_from_path(&idx_path)
            && loaded.index_k == opts.index_k
        {
            return Ok(KmerIndex::InMemory(loaded));
        }
        let dict = KmerDict::build(genome, opts)?;
        dict.save_to_path(&idx_path)?;
        Ok(KmerIndex::InMemory(dict))
    } else {
        KmerIndex::build(genome, opts)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn detects_forward_and_inverted_hsps() {
        let block_a = "ACGTCCGGTTAAACGTTTGGCCAAAACCCGGGTTTACGTACGGTTAACCGGAATTCCGGAA";
        let block_b = "GCTAGCTAGGCCTTAAACCGGTTTGGCCAATTCCGGATGCTAGCTAGGCCTTAAACCGGTT";
        let rc_b: String = block_b
            .bytes()
            .rev()
            .map(crate::seq::complement_base)
            .map(|b| b as char)
            .collect();

        let query = Genome::from_records(
            "q",
            [(
                "q0",
                format!("NNNNNNNNNN{block_a}NNNNNNNNNN{block_b}NNNNNNNNNN"),
            )],
        );
        let reference = Genome::from_records(
            "r",
            [(
                "r0",
                format!("TTTTTTTTTT{block_a}CCCCCCCCCC{rc_b}AAAAAAAAAA"),
            )],
        );

        let cfg = CompareConfig {
            min_length: 40,
            min_similarity: 60.0,
            word_length: 12,
            ..Default::default()
        };

        let res = compare_genomes(&query, &reference, &cfg, None).unwrap();
        assert_eq!(res.stats.forward_fragments, 1);
        assert_eq!(res.stats.reverse_fragments, 1);
        assert_eq!(res.fragments.len(), 2);

        let fwd = res.fragments.iter().find(|f| f.strand == 'f').unwrap();
        assert_eq!(fwd.length, block_a.len() as u64);
        assert_eq!(fwd.x_start, 10);
        assert_eq!(fwd.y_start, 10);

        let rev = res.fragments.iter().find(|f| f.strand == 'r').unwrap();
        assert_eq!(rev.length, block_b.len() as u64);
        assert_eq!(rev.x_start, (10 + block_a.len() + 10) as u64);
    }

    #[test]
    fn on_disk_index_matches_in_memory_exactly() {
        let mut s_x = String::new();
        let mut s_y = String::new();
        let mut seed = 12345u64;
        let bases = *b"ACGT";
        let mut next_dna = |n: usize| {
            let mut out = Vec::with_capacity(n);
            for _ in 0..n {
                seed = seed
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                out.push(bases[((seed >> 32) & 3) as usize]);
            }
            String::from_utf8(out).unwrap()
        };

        let fwd_block = next_dna(180);
        let rev_block = next_dna(160);
        let rev_block_rc: String = rev_block
            .bytes()
            .rev()
            .map(crate::seq::complement_base)
            .map(|b| b as char)
            .collect();

        // Make sequences > 1500 bp so a 32 KiB budget forces KmerIndex::OnDisk spilling!
        s_x.push_str(&next_dna(700));
        s_x.push_str(&fwd_block);
        s_x.push_str(&next_dna(700));
        s_x.push_str(&rev_block);
        s_x.push_str(&next_dna(700));

        s_y.push_str(&next_dna(650));
        s_y.push_str(&fwd_block);
        s_y.push_str(&next_dna(650));
        s_y.push_str(&rev_block_rc);
        s_y.push_str(&next_dna(650));

        let q = Genome::from_records("q", [("q0", s_x.as_bytes())]);
        let r = Genome::from_records("r", [("r0", s_y.as_bytes())]);

        let cfg_mem = CompareConfig {
            min_length: 80,
            min_similarity: 70.0,
            word_length: 12,
            memory_budget: 64 * 1024 * 1024,
            ..Default::default()
        };
        let cfg_disk = CompareConfig {
            memory_budget: 32 * 1024,
            ..cfg_mem.clone()
        };

        let res_mem = compare_genomes_kmer_index(&q, &r, &cfg_mem, None).unwrap();
        let res_disk = compare_genomes(&q, &r, &cfg_disk, None).unwrap();
        let res_stream = compare_genomes(&q, &r, &cfg_mem, None).unwrap();
        assert_eq!(res_mem.fragments.len(), 2);
        assert_eq!(res_mem.fragments, res_disk.fragments);
        assert_eq!(res_stream.fragments, res_disk.fragments);
        assert_eq!(res_mem.stats.raw_hits, res_disk.stats.raw_hits);
        assert_eq!(res_mem.stats.filtered_hits, res_disk.stats.filtered_hits);
    }

    #[test]
    fn streaming_ref_index_matches_legacy_multi_fasta() {
        let mut seed = 0x5EED_CAFE_1234_5678u64;
        let bases = *b"ACGT";
        let mut next_dna = |n: usize| {
            let mut out = Vec::with_capacity(n);
            for _ in 0..n {
                seed = seed
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                out.push(bases[((seed >> 32) & 3) as usize]);
            }
            String::from_utf8(out).unwrap()
        };

        let fwd1 = next_dna(150);
        let fwd2 = next_dna(120);
        let rev1 = next_dna(140);
        let rev1_rc: String = rev1
            .bytes()
            .rev()
            .map(crate::seq::complement_base)
            .map(|b| b as char)
            .collect();
        // Palindromic 64 bp HSP (identical to its reverse complement).
        let half_pal = next_dna(32);
        let pal: String = format!(
            "{half_pal}{}",
            half_pal
                .bytes()
                .rev()
                .map(crate::seq::complement_base)
                .map(|b| b as char)
                .collect::<String>()
        );
        // Soft-masked lowercase flank around an uppercase seed.
        let masked_hsp = format!(
            "{}{}{}",
            fwd2[..30].to_ascii_lowercase(),
            &fwd2[30..90],
            fwd2[90..].to_ascii_lowercase()
        );

        let q = Genome::from_records(
            "q",
            [
                (
                    "q0",
                    format!(
                        "{}{fwd1}{}{rev1}{}",
                        next_dna(200),
                        next_dna(200),
                        next_dna(200)
                    ),
                ),
                (
                    "q1",
                    format!(
                        "{}{masked_hsp}{}{pal}{}",
                        next_dna(180),
                        next_dna(180),
                        next_dna(180)
                    ),
                ),
            ],
        );
        let r = Genome::from_records(
            "r",
            [
                (
                    "r0",
                    format!(
                        "{}{rev1_rc}{}{pal}{}",
                        next_dna(190),
                        next_dna(190),
                        next_dna(190)
                    ),
                ),
                (
                    "r1",
                    format!(
                        "{}{fwd1}{}{masked_hsp}{}",
                        next_dna(210),
                        next_dna(210),
                        next_dna(210)
                    ),
                ),
            ],
        );

        // Test both Direct CSR (wl = 8) and Hash-CSR (wl = 16, wl = 32), plus forward-only.
        for &(wl, compare_reverse) in &[(8, true), (16, true), (32, true), (16, false)] {
            let cfg = CompareConfig {
                min_length: 50,
                min_similarity: 75.0,
                word_length: wl,
                compare_reverse,
                ..Default::default()
            };

            let res_stream = compare_genomes(&q, &r, &cfg, None).unwrap();
            let res_slab = compare_genomes_slab_streaming(&q, &r, &cfg, 300).unwrap();
            let res_legacy = compare_genomes_kmer_index(&q, &r, &cfg, None).unwrap();

            let mut f_stream = res_stream.fragments.clone();
            let mut f_slab = res_slab.fragments.clone();
            let mut f_legacy = res_legacy.fragments.clone();
            let key = |f: &Fragment| (f.strand, f.seq_x, f.seq_y, f.diag, f.x_start, f.x_end);
            f_stream.sort_by_key(key);
            f_slab.sort_by_key(key);
            f_legacy.sort_by_key(key);

            assert_eq!(
                f_stream, f_legacy,
                "fragment mismatch for wl={wl}, compare_reverse={compare_reverse}"
            );
            assert_eq!(
                f_slab, f_legacy,
                "slab streaming mismatch for wl={wl}, compare_reverse={compare_reverse}"
            );
            assert!(
                !f_stream.is_empty(),
                "expected non-empty fragments for wl={wl}"
            );
        }
    }

    #[test]
    fn partition_slab_streaming_preserves_exact_global_freq_limit() {
        // Place 5 copies of `rep_block` spread far apart across `r`, plus 1 unique `uniq_block`.
        // With `freq_limit = 2`, ALL 5 copies of `rep_block` must be masked globally — even when
        // `r` is split into many partition slabs (`max_kmers_per_slab = 250`).
        let mut seed = 0xfeed_face_cafe_beefu64;
        let mut next_dna = |n: usize| -> String {
            let bases = *b"ACGT";
            let mut s = String::with_capacity(n);
            for _ in 0..n {
                seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                s.push(bases[((seed >> 33) & 3) as usize] as char);
            }
            s
        };

        let rep_block = next_dna(120);
        let uniq_block = next_dna(140);
        let q = Genome::from_records(
            "q",
            [(
                "q0",
                format!(
                    "{}NNNNNNNN{rep_block}NNNNNNNN{}{uniq_block}{}",
                    next_dna(200),
                    next_dna(200),
                    next_dna(200)
                ),
            )],
        );
        let r = Genome::from_records(
            "r",
            [(
                "r0",
                format!(
                    "{}NNNNNNNN{rep_block}NNNNNNNN{}NNNNNNNN{rep_block}NNNNNNNN{}{uniq_block}{}NNNNNNNN{rep_block}NNNNNNNN{}NNNNNNNN{rep_block}NNNNNNNN{}NNNNNNNN{rep_block}NNNNNNNN{}",
                    next_dna(300),
                    next_dna(300),
                    next_dna(300),
                    next_dna(300),
                    next_dna(300),
                    next_dna(300),
                    next_dna(300)
                ),
            )],
        );

        for wl in [8usize, 16, 32] {
            let cfg = CompareConfig {
                min_length: 80,
                min_similarity: 85.0,
                word_length: wl,
                freq_limit: 2,
                compare_reverse: true,
                ..Default::default()
            };
            let single = compare_genomes(&q, &r, &cfg, None).unwrap();
            let multi = compare_genomes_slab_streaming(&q, &r, &cfg, 250).unwrap();

            let key = |f: &Fragment| (f.strand, f.seq_x, f.seq_y, f.diag, f.x_start, f.x_end);
            let mut f_single = single.fragments.clone();
            let mut f_multi = multi.fragments.clone();
            f_single.sort_by_key(key);
            f_multi.sort_by_key(key);

            assert_eq!(
                f_multi, f_single,
                "multi-slab violated global freq_limit invariance for wl={wl}"
            );
            assert_eq!(
                f_multi.len(),
                1,
                "expected only uniq_block HSP (rep_block masked by FL=2) for wl={wl}"
            );
        }
    }

    #[test]
    fn spaced_seed_discovers_divergent_homologs_on_both_strands() {
        let mut rng = 0x5A17_C0DE_BEEF_1234u64;
        let bases = *b"ACGT";
        let mut next_dna = |len: usize| -> String {
            let mut s = Vec::with_capacity(len);
            for _ in 0..len {
                rng = rng
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                s.push(bases[((rng >> 33) & 3) as usize]);
            }
            String::from_utf8(s).unwrap()
        };
        let mutate_every_8 = |src: &str| -> String {
            src.bytes()
                .enumerate()
                .map(|(i, b)| {
                    if i % 8 == 3 {
                        match b {
                            b'A' => b'C',
                            b'C' => b'G',
                            b'G' => b'T',
                            _ => b'A',
                        }
                    } else {
                        b
                    }
                })
                .map(|b| b as char)
                .collect()
        };
        let revcomp = |src: &str| -> String {
            src.bytes()
                .rev()
                .map(|b| match b {
                    b'A' => 'T',
                    b'C' => 'G',
                    b'G' => 'C',
                    _ => 'A',
                })
                .collect()
        };

        let fwd_orig = next_dna(180);
        let rev_orig = next_dna(180);
        // Every 8th base is mutated (87.5% identity, max contiguous run is 7 bp!):
        // contiguous k=12 finds 0 hits, while spaced seed `1110111110111` (w=11/s=13) or `w12s18`
        // finds BOTH forward and reverse homologs!
        let fwd_mut = mutate_every_8(&fwd_orig);
        let rev_mut_rc = revcomp(&mutate_every_8(&rev_orig));

        let q = Genome::from_records(
            "q",
            [(
                "q0",
                format!(
                    "{}{}{}{}{}",
                    next_dna(200),
                    fwd_orig,
                    next_dna(200),
                    rev_orig,
                    next_dna(200)
                ),
            )],
        );
        let r = Genome::from_records(
            "r",
            [(
                "r0",
                format!(
                    "{}{}{}{}{}",
                    next_dna(200),
                    fwd_mut,
                    next_dna(200),
                    rev_mut_rc,
                    next_dna(200)
                ),
            )],
        );

        let cfg_contig = CompareConfig {
            min_length: 100,
            min_similarity: 70.0,
            word_length: 12,
            spaced_seed: None,
            ..Default::default()
        };
        let res_contig = compare_genomes(&q, &r, &cfg_contig, None).unwrap();
        assert_eq!(
            res_contig.fragments.len(),
            0,
            "contiguous k=12 must miss 7-bp-capped matches"
        );

        // Pattern `111011110111`? Check period-8 alignment: positions i%8==3 are mutated (distance 8).
        // `111011111110111` has don't-care at indices 3 and 11 (distance 8, palindromic, w=13/s=15)!
        let sp = SpacedSeed::parse("111011111110111").unwrap();
        let cfg_spaced = CompareConfig {
            min_length: 100,
            min_similarity: 70.0,
            word_length: sp.span,
            spaced_seed: Some(sp),
            ..Default::default()
        };
        let res_spaced = compare_genomes(&q, &r, &cfg_spaced, None).unwrap();
        assert_eq!(
            res_spaced.fragments.len(),
            2,
            "spaced seed must recover both fwd & rev HSPs"
        );
        assert!(res_spaced.fragments.iter().any(|f| f.strand == 'f'));
        assert!(res_spaced.fragments.iter().any(|f| f.strand == 'r'));
        for f in &res_spaced.fragments {
            assert!(
                f.length >= 160 && f.ident >= 150 && f.similarity > 72.0,
                "unexpected frag: {f:?}"
            );
        }

        // Verify multi-slab partition streaming matches single-slab bit-for-bit with spaced seeds.
        let res_slab = compare_genomes_slab_streaming(&q, &r, &cfg_spaced, 150).unwrap();
        let key = |f: &Fragment| (f.strand, f.seq_x, f.seq_y, f.diag, f.x_start, f.x_end);
        let mut a = res_spaced.fragments.clone();
        let mut b = res_slab.fragments.clone();
        a.sort_by_key(key);
        b.sort_by_key(key);
        assert_eq!(a, b);
    }

    #[test]
    fn lazy_coverage_false_matches_fragments_and_stats_bit_exact() {
        let alphabet = *b"ACGT";
        let mut rng = 0xCAFE_1234_ABCD_EF01u64;
        let mut make_dna = |n: usize| -> Vec<u8> {
            (0..n)
                .map(|_| {
                    rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
                    alphabet[((rng >> 32) & 3) as usize]
                })
                .collect()
        };
        let mut sx = make_dna(6_000);
        let mut sy = make_dna(6_000);
        sy[500..900].copy_from_slice(&sx[200..600]);
        let rc: Vec<u8> = sx[1_500..1_900]
            .iter()
            .rev()
            .map(|&b| match b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                _ => b'A',
            })
            .collect();
        sy[3_000..3_400].copy_from_slice(&rc);
        sx[350] = b'T';

        let q = Genome::from_records("q", [("q0", sx.as_slice())]);
        let r = Genome::from_records("r", [("r0", sy.as_slice())]);

        for wl in [8usize, 16] {
            let cfg = CompareConfig {
                word_length: wl,
                min_length: 40,
                min_similarity: 60.0,
                freq_limit: 100,
                ..Default::default()
            };
            let res_cov = compare_genomes_with_coverage(&q, &r, &cfg, None, true).unwrap();
            let res_nocov = compare_genomes_with_coverage(&q, &r, &cfg, None, false).unwrap();

            assert_eq!(res_cov.fragments, res_nocov.fragments);
            assert_eq!(res_cov.stats.word_matches, res_nocov.stats.word_matches);
            assert_eq!(res_cov.stats.raw_hits, res_nocov.stats.raw_hits);
            assert_eq!(res_cov.stats.hits_used, res_nocov.stats.hits_used);
            assert_eq!(
                res_cov.stats.forward_fragments,
                res_nocov.stats.forward_fragments
            );
            assert_eq!(
                res_cov.stats.reverse_fragments,
                res_nocov.stats.reverse_fragments
            );
            assert_eq!(res_nocov.coverage.bins.len(), 1000);
            assert!(
                res_nocov
                    .coverage
                    .bins
                    .iter()
                    .all(|row| row.iter().all(|&c| c == 0))
            );
            assert!(
                res_cov
                    .coverage
                    .bins
                    .iter()
                    .any(|row| row.iter().any(|&c| c > 0))
            );
        }
    }

    #[test]
    fn forward_only_specialization_is_bit_exact_across_all_modes_and_edge_cases() {
        let alphabet = *b"ACGT";
        let mut rng = 0x9876_5432_10FE_DCBAu64;
        let mut make_dna = |n: usize| -> Vec<u8> {
            (0..n)
                .map(|_| {
                    rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
                    alphabet[((rng >> 32) & 3) as usize]
                })
                .collect()
        };
        let mut sx = make_dna(7_000);
        let mut sy = make_dna(7_000);
        // Forward homolog with palindrome block and soft-masked interior.
        let pal = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        sx[300..332].copy_from_slice(pal);
        sy[800..1_200].copy_from_slice(&sx[200..600]);
        for b in &mut sx[360..390] {
            *b = b.to_ascii_lowercase();
        }
        // Reverse-only homolog (must NOT appear in forward-only results).
        let rc: Vec<u8> = sx[2_500..2_900]
            .iter()
            .rev()
            .map(|&b| match b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                _ => b'A',
            })
            .collect();
        sy[4_000..4_400].copy_from_slice(&rc);
        // Repeat a 16-mer 5 times on forward strand and 1 time on reverse strand to test freq_limit.
        let rep = b"AACCTTGGCCAATTGG";
        for k in 0..5 {
            sy[5_000 + k * 40..5_000 + k * 40 + 16].copy_from_slice(rep);
        }

        let q = Genome::from_records("q", [("q0", sx.as_slice())]);
        let r = Genome::from_records("r", [("r0", sy.as_slice())]);

        let cases: Vec<(usize, Option<SpacedSeed>, u64, usize, bool)> = vec![
            (8, None, 100, 10_000, true),   // DirectCsr, soft_mask=true
            (12, None, 100, 10_000, false), // HashCsr (small ref), soft_mask=false
            (16, None, 3, 10_000, true),    // HashCsr + tight freq_limit
            (32, None, 100, 10_000, true),  // HashCsr WL32
            (
                20,
                Some(SpacedSeed::parse("w14s20").unwrap()),
                100,
                10_000,
                true,
            ), // Spaced seed
            (16, None, 50, 20, true),       // Adversarial 2*freq_limit > max_rep fallback
        ];

        for (wl, spaced, fl, max_rep, soft_mask) in cases {
            let cfg_dual = CompareConfig {
                word_length: wl,
                spaced_seed: spaced.clone(),
                min_length: 40,
                min_similarity: 60.0,
                freq_limit: fl,
                max_rep,
                soft_mask,
                compare_reverse: true,
                ..Default::default()
            };
            let cfg_fwd = CompareConfig {
                compare_reverse: false,
                ..cfg_dual.clone()
            };

            let res_dual = compare_genomes(&q, &r, &cfg_dual, None).unwrap();
            let res_fwd = compare_genomes(&q, &r, &cfg_fwd, None).unwrap();
            let res_fwd_slab = compare_genomes_slab_streaming(&q, &r, &cfg_fwd, 1_500).unwrap();

            // Canonical RefIndex queried with `compare_reverse = false` (all-vs-all reuse path).
            let canon_idx =
                RefIndex::build_with_mask(&r, wl, cfg_fwd.spaced_mask_2bit(), fl, soft_mask)
                    .unwrap();
            assert!(!canon_idx.is_forward_only());
            let res_canon_fwd = compare_with_ref_index(&q, &r, &canon_idx, None, &cfg_fwd).unwrap();

            let expected_fwd: Vec<Fragment> = res_dual
                .fragments
                .iter()
                .filter(|f| f.strand == 'f')
                .cloned()
                .collect();

            assert_eq!(
                res_fwd.fragments,
                expected_fwd,
                "fwd-only mismatch for wl={wl}, spaced={}, fl={fl}, max_rep={max_rep}",
                spaced.is_some()
            );
            assert_eq!(res_fwd_slab.fragments, expected_fwd);
            assert_eq!(res_canon_fwd.fragments, expected_fwd);
            assert!(res_fwd.fragments.iter().all(|f| f.strand == 'f'));
            assert_eq!(res_fwd.stats.reverse_fragments, 0);

            // Querying a forward-only RefIndex with compare_reverse=true must return InvalidParameter.
            let fwd_idx = RefIndex::build_with_mask_and_strand(
                &r,
                wl,
                cfg_fwd.spaced_mask_2bit(),
                fl,
                soft_mask,
                true,
            )
            .unwrap();
            assert!(fwd_idx.is_forward_only());
            let rev_r = r.reverse_complement();
            assert!(compare_with_ref_index(&q, &r, &fwd_idx, Some(&rev_r), &cfg_dual).is_err());
        }
    }

    #[test]
    fn multislab_without_unique_bm_is_bit_exact_with_single_slab() {
        let alphabet = *b"ACGT";
        let mut rng = 0x7777_1234_ABCD_EF01u64;
        let mut sx: Vec<u8> = (0..12_000)
            .map(|_| {
                rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
                alphabet[((rng >> 32) & 3) as usize]
            })
            .collect();
        let mut sy: Vec<u8> = (0..12_000)
            .map(|_| {
                rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
                alphabet[((rng >> 32) & 3) as usize]
            })
            .collect();

        // Embed conserved forward and inverted HSP blocks plus a repeated motif.
        let fwd_block = sx[1_000..1_350].to_vec();
        sy[2_000..2_350].copy_from_slice(&fwd_block);
        let rev_block: Vec<u8> = sx[4_000..4_320]
            .iter()
            .rev()
            .map(|&b| match b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                _ => b'A',
            })
            .collect();
        sy[6_000..6_320].copy_from_slice(&rev_block);
        // Inject a low-complexity duplicate block so duplicate-marking paths are also exercised.
        let dup = b"ACGTACGTCCGGTTAAACGTTTGGCCAAAACC";
        for k in 0..4 {
            sx[8_000 + k * 60..8_000 + k * 60 + 32].copy_from_slice(dup);
            sy[8_500 + k * 60..8_500 + k * 60 + 32].copy_from_slice(dup);
        }

        let q = Genome::from_records("q", [("q0", sx.as_slice())]);
        let r = Genome::from_records("r", [("r0", sy.as_slice())]);
        let rev_r = r.reverse_complement();

        for wl in [8usize, 16, 32] {
            for compare_reverse in [true, false] {
                let cfg = CompareConfig {
                    word_length: wl,
                    min_length: 40,
                    min_similarity: 60.0,
                    freq_limit: 20,
                    compare_reverse,
                    seeding_mode: SeedingMode::Exact,
                    ..Default::default()
                };
                let single = compare_genomes(&q, &r, &cfg, None).unwrap();
                assert!(!single.fragments.is_empty());

                // Multi-slab path automatically sets `build_unique_bm = false` when `slabs.len() > 1`.
                let multi = compare_genomes_slab_streaming(&q, &r, &cfg, 2_500).unwrap();
                assert_eq!(
                    multi.fragments, single.fragments,
                    "multi-slab fragments mismatch for wl={wl}, compare_reverse={compare_reverse}"
                );
                assert_eq!(multi.stats.word_matches, single.stats.word_matches);
                assert_eq!(multi.stats.raw_hits, single.stats.raw_hits);

                // Also directly verify `RefIndex::build_partitions_with_options(..., build_unique_bm=false)`
                // leaves `unique_bm.is_empty()` and produces 100% bit-exact fragments, stats, AND
                // `coverage.bins` against `build_unique_bm = true` (both single-slab and partial-slab).
                for (ps, pe) in [(0u16, 256u16), (0, 128), (128, 256)] {
                    let idx_with_bm = RefIndex::build_partitions_with_options(
                        &r,
                        ps,
                        pe,
                        wl,
                        None,
                        cfg.freq_limit,
                        cfg.soft_mask,
                        !compare_reverse,
                        true,
                    )
                    .unwrap();
                    let idx_no_bm = RefIndex::build_partitions_with_options(
                        &r,
                        ps,
                        pe,
                        wl,
                        None,
                        cfg.freq_limit,
                        cfg.soft_mask,
                        !compare_reverse,
                        false,
                    )
                    .unwrap();
                    assert!(idx_no_bm.unique_bm.is_empty());
                    let res_with_bm = compare_with_ref_index(
                        &q,
                        &r,
                        &idx_with_bm,
                        compare_reverse.then_some(&rev_r),
                        &cfg,
                    )
                    .unwrap();
                    let res_no_bm = compare_with_ref_index(
                        &q,
                        &r,
                        &idx_no_bm,
                        compare_reverse.then_some(&rev_r),
                        &cfg,
                    )
                    .unwrap();
                    assert_eq!(res_no_bm.fragments, res_with_bm.fragments);
                    assert_eq!(res_no_bm.stats.word_matches, res_with_bm.stats.word_matches);
                    assert_eq!(res_no_bm.stats.raw_hits, res_with_bm.stats.raw_hits);
                    assert_eq!(res_no_bm.stats.hits_used, res_with_bm.stats.hits_used);
                    assert_eq!(res_no_bm.coverage.bins, res_with_bm.coverage.bins);
                    if ps == 0 && pe == 256 {
                        assert_eq!(res_no_bm.fragments, single.fragments);
                        assert_eq!(res_no_bm.coverage.bins, single.coverage.bins);
                    }
                }
            }
        }
    }
}
