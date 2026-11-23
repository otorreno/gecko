use std::path::Path;

use serde::Serialize;

use crate::error::{GeckoError, Result};
use crate::extend::{CoverageMatrix, ExtendOptions, Scoring, extend_hits_parallel};
use crate::formats::csv::ComparisonMetadata;
use crate::formats::frags::Fragment;
use crate::hits::{HitOptions, generate_and_consume_hits, generate_and_consume_hits_from_index};
use crate::kmer::{DictOptions, INDEX_WORD_SIZE, KmerDict, KmerIndex};
use crate::ref_index::{RefIndex, compute_partition_slabs};
use crate::seq::Genome;
use crate::stream::{compare_streaming, dedup_fragments};

/// Full configuration for a pairwise GECKO comparison.
#[derive(Debug, Clone)]
pub struct CompareConfig {
    /// Minimum HSP length (in bp if `fixed_length`, or percentage if `!fixed_length`).
    pub min_length: u64,
    /// Minimum similarity percentage threshold (`0.0 ..= 100.0`).
    pub min_similarity: f32,
    /// Seed k-mer length (`4..=32`; no longer limited to multiples of 4).
    pub word_length: usize,
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
}

impl Default for CompareConfig {
    fn default() -> Self {
        Self {
            min_length: 40,
            min_similarity: 60.0,
            word_length: 12,
            freq_limit: 1000,
            max_rep: 10_000,
            fixed_length: true,
            soft_mask: true,
            memory_budget: 2 * 1024 * 1024 * 1024,
            compare_reverse: true,
            scoring: Scoring::default(),
        }
    }
}

impl CompareConfig {
    pub fn validate(&self) -> Result<()> {
        if !(4..=32).contains(&self.word_length) {
            return Err(GeckoError::InvalidParameter(format!(
                "word_length must be between 4 and 32, got {}",
                self.word_length
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
    config.validate()?;

    if cache_dir.is_none() {
        if RefIndex::fits_in_budget(
            reference,
            config.word_length,
            config.soft_mask,
            config.compare_reverse,
            config.memory_budget,
        ) {
            let ref_index = RefIndex::build(
                reference,
                config.word_length,
                config.freq_limit,
                config.soft_mask,
            )?;
            let rev_ref = if config.compare_reverse {
                Some(reference.reverse_complement())
            } else {
                None
            };
            return compare_streaming(query, reference, &ref_index, rev_ref.as_ref(), config);
        }

        if let Some(max_kmers_per_slab) = RefIndex::max_kmers_per_slab_for_budget(
            reference,
            config.word_length,
            config.compare_reverse,
            config.memory_budget,
        ) {
            return compare_genomes_slab_streaming(query, reference, config, max_kmers_per_slab);
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
pub(crate) fn compare_genomes_slab_streaming(
    query: &Genome,
    reference: &Genome,
    config: &CompareConfig,
    max_kmers_per_slab: usize,
) -> Result<CompareResult> {
    let slabs = compute_partition_slabs(
        reference,
        config.word_length,
        config.soft_mask,
        max_kmers_per_slab,
    );
    if slabs.is_empty() {
        return Ok(CompareResult {
            fragments: Vec::new(),
            coverage: CoverageMatrix::new(),
            stats: CompareStats::default(),
        });
    }

    let rev_ref = if config.compare_reverse {
        Some(reference.reverse_complement())
    } else {
        None
    };

    let mut merged = CompareResult {
        fragments: Vec::new(),
        coverage: CoverageMatrix::new(),
        stats: CompareStats::default(),
    };

    for (part_start, part_end) in slabs {
        let slab_index = RefIndex::build_partitions(
            reference,
            part_start,
            part_end,
            config.word_length,
            config.freq_limit,
            config.soft_mask,
        )?;
        let mut slab_res =
            compare_streaming(query, reference, &slab_index, rev_ref.as_ref(), config)?;
        drop(slab_index);

        merged.stats.word_matches += slab_res.stats.word_matches;
        merged.stats.raw_hits += slab_res.stats.raw_hits;
        merged.stats.filtered_hits += slab_res.stats.filtered_hits;
        merged.stats.hits_used += slab_res.stats.hits_used;
        merged.coverage.merge_from(&slab_res.coverage);
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
    compare_streaming(query, reference, ref_index, rev_ref, config)
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
            let bases = [b'A', b'C', b'G', b'T'];
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
}
