use std::path::Path;

use serde::Serialize;

use crate::error::{GeckoError, Result};
use crate::extend::{CoverageMatrix, ExtendOptions, Scoring, extend_hits_parallel};
use crate::formats::csv::ComparisonMetadata;
use crate::formats::frags::Fragment;
use crate::hits::{HitOptions, generate_and_consume_hits};
use crate::kmer::{DictOptions, INDEX_WORD_SIZE, KmerDict};
use crate::seq::Genome;

/// Full configuration for a pairwise GECKO comparison.
#[derive(Debug, Clone)]
pub struct CompareConfig {
    /// Minimum HSP length (in bp if `fixed_length`, or percentage if `!fixed_length`).
    pub min_length: u64,
    /// Minimum similarity percentage threshold (`0.0 ..= 100.0`).
    pub min_similarity: f32,
    /// Seed k-mer length (`4..=32`; no longer limited to multiples of 4).
    pub word_length: usize,
    /// Skip dictionary 32-mers occurring more than `freq_limit` times (default `1000`).
    pub freq_limit: u64,
    /// Subsampling cap per word list (default `10_000`).
    pub max_rep: usize,
    /// If true, `min_length` is in bp; if false, percentage of shorter sequence length.
    pub fixed_length: bool,
    /// If true (default), lowercase bases do not seed k-mers (soft-masking) but can extend.
    pub soft_mask: bool,
    /// Memory budget in bytes for sorting/buffering before spilling to temporary disk runs.
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
pub fn compare_genomes(
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

    let dict_x = load_or_build_dict(query, dict_opts, cache_dir)?;

    {
        let dict_y = load_or_build_dict(reference, dict_opts, cache_dir)?;
        run_strand(query, &dict_x, reference, &dict_y, 'f', config, &mut result)?;
    }

    if config.compare_reverse {
        let rev_ref = reference.reverse_complement();
        let dict_y_rev = load_or_build_dict(&rev_ref, dict_opts, cache_dir)?;
        run_strand(
            query,
            &dict_x,
            &rev_ref,
            &dict_y_rev,
            'r',
            config,
            &mut result,
        )?;
    }

    Ok(result)
}

/// Compare two genomes using pre-built [`KmerDict`]s (ideal for all-vs-all workflows).
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

    Ok(())
}

fn load_or_build_dict(
    genome: &Genome,
    opts: DictOptions,
    cache_dir: Option<&Path>,
) -> Result<KmerDict> {
    if let Some(dir) = cache_dir {
        std::fs::create_dir_all(dir)?;
        let mask_tag = if opts.soft_mask { "sm" } else { "nm" };
        let idx_path = dir.join(format!("{}.k{}.{mask_tag}.gidx", genome.name, opts.index_k));
        if idx_path.exists()
            && let Ok(loaded) = KmerDict::load_from_path(&idx_path)
            && loaded.index_k == opts.index_k
        {
            return Ok(loaded);
        }
        let dict = KmerDict::build(genome, opts)?;
        dict.save_to_path(&idx_path)?;
        Ok(dict)
    } else {
        KmerDict::build(genome, opts)
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
}
