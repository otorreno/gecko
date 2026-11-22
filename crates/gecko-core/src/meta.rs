//! MetaGECKO metagenomic read-to-genome profiling, spurious-match filtering,
//! multi-mapping resolution, and multi-sample differential abundance analysis.
//!
//! Built directly on top of GECKO v2's 2-bit [`PackedSeq`], canonical [`RefIndex`],
//! AVX2/BMI2 extension, and CUDA/WGPU GPU backends so all backends share the same
//! prepare-once, query-many metagenomics pipeline.

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use serde::Serialize;

use crate::error::{GeckoError, Result};
use crate::extend::CoverageMatrix;
use crate::formats::frags::Fragment;
use crate::pipeline::{CompareConfig, CompareResult};
use crate::seq::{Genome, PackedSeqBuilder, ReadBatch, SeqRecord, file_stem_name};

/// Metadata for a single reference genome (taxon/assembly) inside a [`MetaReference`].
#[derive(Debug, Clone, Serialize)]
pub struct RefGenomeInfo {
    pub genome_id: usize,
    pub name: String,
    pub path: String,
    pub first_record: usize,
    pub num_records: usize,
    pub total_bases: u64,
    /// Start coordinate in concatenated base-only (`start_no_sep`) reference space.
    pub start_no_sep: u64,
}

/// Unified multi-genome reference database for MetaGECKO.
///
/// Concatenates one or more genome FASTA files (or splits records of a single multi-FASTA
/// into individual taxa when `split_records = true`) into a single 2-bit [`Genome`] so the
/// CPU, CUDA, or WGPU reference index is built **once** and queried across all read batches
/// and all metagenomic samples.
#[derive(Debug, Clone)]
pub struct MetaReference {
    pub genome: Genome,
    pub genomes: Vec<RefGenomeInfo>,
    /// O(1) lookup from `Fragment::seq_y` (record index) -> `genome_id`.
    pub record_to_genome: Vec<usize>,
}

impl MetaReference {
    /// Load one or more reference genome files (or directories of genome FASTAs) into a
    /// unified [`MetaReference`].
    ///
    /// - If `split_records` is `true` and a single file is supplied, each FASTA record is
    ///   treated as a distinct taxon/genome (useful for concatenated reference databases).
    /// - Otherwise, each input file represents one genome (potentially multi-contig).
    pub fn load_from_paths(paths: &[PathBuf], split_records: bool) -> Result<Self> {
        let mut files: Vec<PathBuf> = Vec::new();
        for p in paths {
            if p.is_dir() {
                collect_sequence_files(p, &mut files)?;
            } else {
                files.push(p.clone());
            }
        }
        files.sort();
        if files.is_empty() {
            return Err(GeckoError::InvalidParameter(
                "no reference genome FASTA files found".into(),
            ));
        }

        if files.len() == 1 {
            let g = Genome::from_fasta_path(&files[0])?;
            return Ok(Self::from_single_genome(
                g,
                files[0].display().to_string(),
                split_records,
            ));
        }

        // Load each genome file in parallel, then pack into one unified Genome.
        use rayon::prelude::*;
        let loaded: Vec<(PathBuf, Genome)> = files
            .par_iter()
            .map(|p| Genome::from_fasta_path(p).map(|g| (p.clone(), g)))
            .collect::<Result<Vec<_>>>()?;

        let est_bases: usize = loaded.iter().map(|(_, g)| g.total_bases as usize).sum();
        let mut builder = PackedSeqBuilder::with_capacity_bases(est_bases);
        let mut all_records: Vec<SeqRecord> = Vec::new();
        let mut genomes: Vec<RefGenomeInfo> = Vec::with_capacity(loaded.len());
        let mut record_to_genome: Vec<usize> = Vec::new();
        let mut total_bases = 0u64;

        for (genome_id, (path, g)) in loaded.into_iter().enumerate() {
            let first_record = all_records.len();
            let genome_start_no_sep = total_bases;
            let raw_ascii = g.seq.to_ascii_vec();

            for rec in &g.records {
                if !all_records.is_empty() {
                    builder.push_sep();
                }
                let start_with_sep = builder.len() as u64;
                let s = rec.start_with_sep as usize;
                let e = s + rec.len as usize;
                let added = builder.push_fasta_bytes(&raw_ascii[s..e]);
                all_records.push(SeqRecord {
                    id: rec.id.clone(),
                    header: rec.header.clone(),
                    start_with_sep,
                    start_no_sep: total_bases,
                    len: added,
                });
                record_to_genome.push(genome_id);
                total_bases += added;
            }

            genomes.push(RefGenomeInfo {
                genome_id,
                name: file_stem_name(&path),
                path: path.display().to_string(),
                first_record,
                num_records: g.records.len(),
                total_bases: g.total_bases,
                start_no_sep: genome_start_no_sep,
            });
        }

        Ok(Self {
            genome: Genome {
                name: format!("meta_ref_{}_genomes", genomes.len()),
                seq: builder.finish(),
                records: all_records,
                total_bases,
            },
            genomes,
            record_to_genome,
        })
    }

    /// Wrap an already-loaded [`Genome`] into a [`MetaReference`].
    pub fn from_single_genome(genome: Genome, path: String, split_records: bool) -> Self {
        if split_records {
            let mut genomes = Vec::with_capacity(genome.records.len());
            let mut record_to_genome = Vec::with_capacity(genome.records.len());
            for (idx, rec) in genome.records.iter().enumerate() {
                genomes.push(RefGenomeInfo {
                    genome_id: idx,
                    name: rec.id.clone(),
                    path: path.clone(),
                    first_record: idx,
                    num_records: 1,
                    total_bases: rec.len,
                    start_no_sep: rec.start_no_sep,
                });
                record_to_genome.push(idx);
            }
            Self {
                genome,
                genomes,
                record_to_genome,
            }
        } else {
            let nrecs = genome.records.len();
            let info = RefGenomeInfo {
                genome_id: 0,
                name: genome.name.clone(),
                path,
                first_record: 0,
                num_records: nrecs,
                total_bases: genome.total_bases,
                start_no_sep: 0,
            };
            Self {
                genome,
                genomes: vec![info],
                record_to_genome: vec![0; nrecs],
            }
        }
    }

    /// O(1) lookup from reference record index (`Fragment::seq_y`) to genome ID.
    #[inline]
    pub fn genome_for_record(&self, record_idx: usize) -> usize {
        self.record_to_genome.get(record_idx).copied().unwrap_or(0)
    }

    #[inline]
    pub fn num_genomes(&self) -> usize {
        self.genomes.len()
    }
}

/// Strategy for resolving reads that align to multiple reference genomes/loci.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub enum MultiMappingMode {
    /// Keep all alignments within `tie_score_margin` of the read's best score, splitting
    /// abundance weight `1.0 / n_tied_genomes` across tied genomes while tracking unique reads.
    Fractional,
    /// Keep only the single highest-scoring alignment per read (deterministic tie-break).
    Best,
    /// Keep all alignments passing quality/coverage thresholds (exhaustive MetaGECKO mode).
    All,
}

/// Configuration for a MetaGECKO metagenomic profiling run.
#[derive(Debug, Clone)]
pub struct MetaConfig {
    /// Underlying GECKO v2 seed-and-extend comparison parameters.
    pub compare: CompareConfig,
    /// Minimum fraction of the read length covered by the HSP (`0.0 ..= 100.0`, default `0.0`).
    pub min_read_coverage: f32,
    /// Minimum exact sequence identity percentage (`0.0 ..= 100.0`, default `0.0`).
    pub min_identity: f32,
    /// Skip reads shorter than this length in bp (default `0`).
    pub min_read_len: u64,
    /// How to resolve reads hitting multiple reference genomes.
    pub multi_mapping: MultiMappingMode,
    /// Score margin below a read's best score still considered tied in `Fractional` mode.
    pub tie_score_margin: u64,
    /// Maximum base pairs per streaming read batch (default 64 Mbp; keeps RAM & GPU VRAM bounded).
    pub batch_bases: usize,
    /// Whether to collect per-read assignment records (`reads.tsv`).
    pub collect_assignments: bool,
    /// Whether to retain individual `Fragment` rows for `.csv` / `.frags` / `.paf` / `.jsonl`.
    pub collect_fragments: bool,
}

impl Default for MetaConfig {
    fn default() -> Self {
        Self {
            compare: CompareConfig {
                min_length: 40,
                min_similarity: 70.0,
                word_length: 16,
                freq_limit: 1000,
                soft_mask: false, // Metagenomic reads are typically unmasked
                ..CompareConfig::default()
            },
            min_read_coverage: 0.0,
            min_identity: 0.0,
            min_read_len: 0,
            multi_mapping: MultiMappingMode::Fractional,
            tie_score_margin: 0,
            batch_bases: 64 * 1024 * 1024,
            collect_assignments: false,
            collect_fragments: false,
        }
    }
}

/// Per-read alignment assignment emitted in `<sample>.reads.tsv`.
#[derive(Debug, Clone, Serialize)]
pub struct ReadAssignment {
    pub read_id: String,
    pub read_index: u64,
    pub read_len: u64,
    pub genome_id: usize,
    pub genome_name: String,
    pub contig_id: String,
    pub read_start: u64,
    pub read_end: u64,
    pub ref_start: u64,
    pub ref_end: u64,
    pub strand: char,
    pub length: u64,
    pub score: u64,
    pub identity_pct: f32,
    pub similarity_pct: f32,
    pub read_coverage_pct: f32,
    pub weight: f64,
    pub is_unique_genome: bool,
}

/// Aggregated taxonomic abundance, coverage breadth/depth, and low-abundance confidence
/// metrics for one reference genome in a metagenomic sample.
#[derive(Debug, Clone, Serialize)]
pub struct GenomeProfile {
    pub genome_id: usize,
    pub genome_name: String,
    pub genome_length: u64,
    /// Reads that mapped exclusively to this genome (key signal for low-abundance taxa).
    pub unique_reads: u64,
    /// Reads that also tied with at least one other genome.
    pub multi_reads: u64,
    /// Fractional-weighted read count (`unique_reads + sum(1 / n_tied_genomes)`).
    pub assigned_reads: f64,
    /// Total HSPs passing all filters on this genome.
    pub total_hsps: u64,
    /// Total aligned read bases on this genome (sum of HSP lengths).
    pub mapped_bases: u64,
    /// Distinct reference bases covered by >= 1 passing HSP.
    pub covered_bases: u64,
    /// `100 * covered_bases / genome_length`.
    pub breadth_pct: f64,
    /// `mapped_bases / genome_length`.
    pub mean_depth: f64,
    /// Expected Lander-Waterman breadth `100 * (1 - exp(-mean_depth))`.
    pub expected_breadth_pct: f64,
    /// `breadth_pct / expected_breadth_pct` (~1.0 = uniform coverage; <<0.3 = localized repeat pileup).
    pub evenness_ratio: f64,
    /// Percentage of all mapped reads (`100 * assigned_reads / total_mapped_reads`).
    pub read_abundance_pct: f64,
    /// Genome-length-normalized relative abundance (`RPK_i / sum(RPK)` in percent).
    pub relative_abundance_pct: f64,
    /// Reads Per Kilobase of genome per Million mapped reads.
    pub rpkm: f64,
    pub mean_identity_pct: f64,
    pub mean_similarity_pct: f64,
    pub mean_hsp_length: f64,
}

/// Complete profiling result for one metagenomic sample against a [`MetaReference`].
#[derive(Debug, Clone, Serialize)]
pub struct MetaSampleResult {
    pub sample_name: String,
    pub total_raw_reads: u64,
    pub short_reads_skipped: u64,
    pub reads_analyzed: u64,
    pub total_read_bases: u64,
    pub mapped_reads: u64,
    pub unique_mapped_reads: u64,
    pub multi_mapped_reads: u64,
    pub total_raw_hsps: u64,
    pub total_kept_hsps: u64,
    pub total_hits: u64,
    pub hits_used: u64,
    pub batches_processed: usize,
    pub profiles: Vec<GenomeProfile>,
    #[serde(skip)]
    pub assignments: Vec<ReadAssignment>,
    #[serde(skip)]
    pub fragments: Vec<Fragment>,
    #[serde(skip)]
    pub coverage: CoverageMatrix,
}

#[derive(Default)]
struct GenomeAccum {
    unique_reads: u64,
    multi_reads: u64,
    assigned_reads: f64,
    total_hsps: u64,
    mapped_bases: u64,
    sum_identity: f64,
    sum_similarity: f64,
    sum_length: u64,
    /// Half-open `[start, end)` intervals in reference `start_no_sep` space for exact breadth.
    intervals: Vec<(u64, u64)>,
}

/// Streaming per-sample accumulator consuming batch [`CompareResult`] outputs from
/// CPU, CUDA, or WGPU backends in bounded memory.
pub struct MetaSampleAccumulator {
    sample_name: String,
    num_genomes: usize,
    total_raw_reads: u64,
    short_reads_skipped: u64,
    reads_analyzed: u64,
    total_read_bases: u64,
    mapped_reads: u64,
    unique_mapped_reads: u64,
    multi_mapped_reads: u64,
    total_raw_hsps: u64,
    total_kept_hsps: u64,
    total_hits: u64,
    hits_used: u64,
    batches_processed: usize,
    genomes: Vec<GenomeAccum>,
    assignments: Vec<ReadAssignment>,
    fragments: Vec<Fragment>,
    coverage: CoverageMatrix,
}

impl MetaSampleAccumulator {
    pub fn new(sample_name: impl Into<String>, reference: &MetaReference) -> Self {
        let num_genomes = reference.num_genomes();
        let genomes = (0..num_genomes).map(|_| GenomeAccum::default()).collect();
        Self {
            sample_name: sample_name.into(),
            num_genomes,
            total_raw_reads: 0,
            short_reads_skipped: 0,
            reads_analyzed: 0,
            total_read_bases: 0,
            mapped_reads: 0,
            unique_mapped_reads: 0,
            multi_mapped_reads: 0,
            total_raw_hsps: 0,
            total_kept_hsps: 0,
            total_hits: 0,
            hits_used: 0,
            batches_processed: 0,
            genomes,
            assignments: Vec::new(),
            fragments: Vec::new(),
            coverage: CoverageMatrix::new(),
        }
    }

    /// Ingest one read batch's comparison results from any backend (CPU, CUDA, or WGPU).
    pub fn ingest_batch(
        &mut self,
        batch: &ReadBatch,
        reference: &MetaReference,
        result: CompareResult,
        config: &MetaConfig,
    ) {
        self.batches_processed += 1;
        self.total_raw_reads += batch.raw_reads_in_batch;
        self.short_reads_skipped += batch.short_reads_skipped;
        self.reads_analyzed += batch.genome.records.len() as u64;
        self.total_read_bases += batch.genome.total_bases;
        self.total_hits += result.stats.filtered_hits;
        self.hits_used += result.stats.hits_used;
        self.total_raw_hsps += result.fragments.len() as u64;

        let y_total = reference.genome.total_bases;

        // 1. Filter batch fragments by read coverage & identity, then group by batch-local seq_x.
        let mut by_read: BTreeMap<usize, Vec<Fragment>> = BTreeMap::new();
        for f in result.fragments {
            let local_read_idx = f.seq_x as usize;
            if local_read_idx >= batch.genome.records.len() {
                continue;
            }
            let read_len = batch.genome.records[local_read_idx].len.max(1);
            let read_cov = (f.length as f32) * 100.0 / (read_len as f32);
            let ident_pct = (f.ident as f32) * 100.0 / (f.length.max(1) as f32);
            if read_cov < config.min_read_coverage || ident_pct < config.min_identity {
                continue;
            }
            by_read.entry(local_read_idx).or_default().push(f);
        }

        // 2. Resolve multi-mapping per read and accumulate exact per-genome metrics.
        for (local_read_idx, mut hits) in by_read {
            let read_rec = &batch.genome.records[local_read_idx];
            let read_len = read_rec.len.max(1);
            let global_read_idx = batch.first_read_id + local_read_idx as u64;

            // Deterministic ordering: highest score, then length, then identity, then genome_id, then y_start.
            hits.sort_unstable_by(|a, b| {
                b.score
                    .cmp(&a.score)
                    .then_with(|| b.length.cmp(&a.length))
                    .then_with(|| b.ident.cmp(&a.ident))
                    .then_with(|| a.seq_y.cmp(&b.seq_y))
                    .then_with(|| a.y_start.cmp(&b.y_start))
            });

            let best_score = hits[0].score;
            let selected: Vec<Fragment> = match config.multi_mapping {
                MultiMappingMode::Best => {
                    hits.truncate(1);
                    hits
                }
                MultiMappingMode::Fractional => {
                    let cutoff = best_score.saturating_sub(config.tie_score_margin);
                    hits.into_iter().filter(|h| h.score >= cutoff).collect()
                }
                MultiMappingMode::All => hits,
            };

            if selected.is_empty() {
                continue;
            }

            // Collect distinct genomes hit by this read's selected alignments.
            let mut hit_genomes: Vec<usize> = selected
                .iter()
                .map(|f| reference.genome_for_record(f.seq_y as usize))
                .collect();
            hit_genomes.sort_unstable();
            hit_genomes.dedup();

            let n_genomes_hit = hit_genomes.len();
            let is_unique = n_genomes_hit == 1;
            let weight = 1.0 / (n_genomes_hit as f64);

            self.mapped_reads += 1;
            if is_unique {
                self.unique_mapped_reads += 1;
            } else {
                self.multi_mapped_reads += 1;
            }

            for &gid in &hit_genomes {
                let gacc = &mut self.genomes[gid];
                gacc.assigned_reads += weight;
                if is_unique {
                    gacc.unique_reads += 1;
                } else {
                    gacc.multi_reads += 1;
                }
            }

            // Record HSP-level coverage, intervals, assignments, and optional fragments.
            for mut f in selected {
                let gid = reference.genome_for_record(f.seq_y as usize);
                let gacc = &mut self.genomes[gid];
                let ident_pct = (f.ident as f64) * 100.0 / (f.length.max(1) as f64);
                let read_cov_pct = (f.length as f32) * 100.0 / (read_len as f32);

                gacc.total_hsps += 1;
                gacc.mapped_bases += f.length;
                gacc.sum_identity += ident_pct;
                gacc.sum_similarity += f.similarity as f64;
                gacc.sum_length += f.length;

                // Convert to forward reference [lo, hi) interval in `start_no_sep` coordinates.
                let (ref_lo, ref_hi_incl) = if f.strand == 'r' {
                    let y_hi = y_total.saturating_sub(f.y_start + 1);
                    let y_lo = y_total.saturating_sub(f.y_end + 1);
                    (y_lo.min(y_hi), y_lo.max(y_hi))
                } else {
                    (f.y_start.min(f.y_end), f.y_start.max(f.y_end))
                };
                gacc.intervals.push((ref_lo, ref_hi_incl + 1));

                self.coverage.record(f.length, f.similarity);
                self.total_kept_hsps += 1;

                if config.collect_assignments {
                    let ref_rec = &reference.genome.records[f.seq_y as usize];
                    let local_ref_start = ref_lo.saturating_sub(ref_rec.start_no_sep);
                    let local_ref_end = ref_hi_incl.saturating_sub(ref_rec.start_no_sep);
                    let local_read_start = f.x_start.saturating_sub(read_rec.start_no_sep);
                    let local_read_end = f.x_end.saturating_sub(read_rec.start_no_sep);
                    self.assignments.push(ReadAssignment {
                        read_id: read_rec.id.clone(),
                        read_index: global_read_idx,
                        read_len,
                        genome_id: gid,
                        genome_name: reference.genomes[gid].name.clone(),
                        contig_id: ref_rec.id.clone(),
                        read_start: local_read_start,
                        read_end: local_read_end,
                        ref_start: local_ref_start,
                        ref_end: local_ref_end,
                        strand: f.strand,
                        length: f.length,
                        score: f.score,
                        identity_pct: ident_pct as f32,
                        similarity_pct: f.similarity,
                        read_coverage_pct: read_cov_pct,
                        weight,
                        is_unique_genome: is_unique,
                    });
                }

                if config.collect_fragments {
                    f.seq_x = global_read_idx;
                    self.fragments.push(f);
                }
            }
        }
    }

    /// Finalize exact coverage breadth, Lander-Waterman evenness, RPKM, and relative abundance.
    pub fn finish(self, reference: &MetaReference) -> MetaSampleResult {
        let total_mapped = self.mapped_reads as f64;

        // First compute RPK per genome for length-normalized relative abundance
        let rpk: Vec<f64> = self
            .genomes
            .iter()
            .enumerate()
            .map(|(gid, acc)| {
                let kb = (reference.genomes[gid].total_bases.max(1) as f64) / 1000.0;
                acc.assigned_reads / kb
            })
            .collect();
        let total_rpk: f64 = rpk.iter().sum();

        let mut profiles: Vec<GenomeProfile> = self
            .genomes
            .into_iter()
            .enumerate()
            .map(|(gid, mut acc)| {
                let info = &reference.genomes[gid];
                let glen = info.total_bases.max(1);
                let covered_bases = union_interval_len(&mut acc.intervals).min(glen);
                let breadth_pct = (covered_bases as f64) * 100.0 / (glen as f64);
                let mean_depth = (acc.mapped_bases as f64) / (glen as f64);
                let expected_breadth_pct = (1.0 - (-mean_depth).exp()) * 100.0;
                let evenness_ratio = if expected_breadth_pct > 1e-9 {
                    (breadth_pct / expected_breadth_pct).min(1.5)
                } else {
                    0.0
                };
                let read_abundance_pct = if total_mapped > 0.0 {
                    acc.assigned_reads * 100.0 / total_mapped
                } else {
                    0.0
                };
                let relative_abundance_pct = if total_rpk > 0.0 {
                    rpk[gid] * 100.0 / total_rpk
                } else {
                    0.0
                };
                let rpkm = if total_mapped > 0.0 {
                    rpk[gid] * 1_000_000.0 / total_mapped
                } else {
                    0.0
                };
                let denom = acc.total_hsps.max(1) as f64;
                GenomeProfile {
                    genome_id: gid,
                    genome_name: info.name.clone(),
                    genome_length: info.total_bases,
                    unique_reads: acc.unique_reads,
                    multi_reads: acc.multi_reads,
                    assigned_reads: acc.assigned_reads,
                    total_hsps: acc.total_hsps,
                    mapped_bases: acc.mapped_bases,
                    covered_bases,
                    breadth_pct,
                    mean_depth,
                    expected_breadth_pct,
                    evenness_ratio,
                    read_abundance_pct,
                    relative_abundance_pct,
                    rpkm,
                    mean_identity_pct: if acc.total_hsps > 0 {
                        acc.sum_identity / denom
                    } else {
                        0.0
                    },
                    mean_similarity_pct: if acc.total_hsps > 0 {
                        acc.sum_similarity / denom
                    } else {
                        0.0
                    },
                    mean_hsp_length: if acc.total_hsps > 0 {
                        (acc.sum_length as f64) / denom
                    } else {
                        0.0
                    },
                }
            })
            .collect();

        // Sort genomes by assigned reads descending (most abundant taxa first)
        profiles.sort_by(|a, b| {
            b.assigned_reads
                .total_cmp(&a.assigned_reads)
                .then_with(|| a.genome_id.cmp(&b.genome_id))
        });

        let _ = self.num_genomes;
        MetaSampleResult {
            sample_name: self.sample_name,
            total_raw_reads: self.total_raw_reads,
            short_reads_skipped: self.short_reads_skipped,
            reads_analyzed: self.reads_analyzed,
            total_read_bases: self.total_read_bases,
            mapped_reads: self.mapped_reads,
            unique_mapped_reads: self.unique_mapped_reads,
            multi_mapped_reads: self.multi_mapped_reads,
            total_raw_hsps: self.total_raw_hsps,
            total_kept_hsps: self.total_kept_hsps,
            total_hits: self.total_hits,
            hits_used: self.hits_used,
            batches_processed: self.batches_processed,
            profiles,
            assignments: self.assignments,
            fragments: self.fragments,
            coverage: self.coverage,
        }
    }
}

/// Compute exact total base length of the union of half-open `[start, end)` intervals.
fn union_interval_len(intervals: &mut [(u64, u64)]) -> u64 {
    if intervals.is_empty() {
        return 0;
    }
    intervals.sort_unstable_by_key(|&(s, _)| s);
    let mut total = 0u64;
    let mut cur_s = intervals[0].0;
    let mut cur_e = intervals[0].1;
    for &(s, e) in &intervals[1..] {
        if s <= cur_e {
            if e > cur_e {
                cur_e = e;
            }
        } else {
            total += cur_e.saturating_sub(cur_s);
            cur_s = s;
            cur_e = e;
        }
    }
    total + cur_e.saturating_sub(cur_s)
}

/// Collect sequence files (FASTA/FASTQ/gzipped) from a directory.
pub fn collect_sequence_files(dir: &Path, out: &mut Vec<PathBuf>) -> Result<()> {
    for entry in std::fs::read_dir(dir)? {
        let entry = entry?;
        let path = entry.path();
        if path.is_file() && is_sequence_filename(&path) {
            out.push(path);
        }
    }
    Ok(())
}

/// Check whether a path looks like a FASTA/FASTQ sequence file (plain or `.gz`).
pub fn is_sequence_filename(path: &Path) -> bool {
    let fname = path
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();
    let base = fname.strip_suffix(".gz").unwrap_or(&fname);
    matches!(
        base.rsplit_once('.').map(|(_, ext)| ext),
        Some("fa" | "fna" | "fasta" | "fas" | "ffn" | "fq" | "fastq")
    )
}

/// Write per-genome taxonomic abundance & coverage profile as TSV.
pub fn write_abundance_tsv(path: impl AsRef<Path>, res: &MetaSampleResult) -> Result<()> {
    let file = File::create(path)?;
    let mut w = BufWriter::with_capacity(4 << 20, file);
    writeln!(
        w,
        "genome\tgenome_len\tunique_reads\tmulti_reads\tassigned_reads\ttotal_hsps\t\
         mapped_bases\tcovered_bases\tbreadth_pct\tmean_depth\texpected_breadth_pct\t\
         evenness_ratio\tread_abundance_pct\trel_abundance_pct\trpkm\t\
         mean_identity_pct\tmean_similarity_pct\tmean_hsp_len"
    )?;
    for p in &res.profiles {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{:.2}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.2}\t{:.2}\t{:.2}\t{:.1}",
            p.genome_name,
            p.genome_length,
            p.unique_reads,
            p.multi_reads,
            p.assigned_reads,
            p.total_hsps,
            p.mapped_bases,
            p.covered_bases,
            p.breadth_pct,
            p.mean_depth,
            p.expected_breadth_pct,
            p.evenness_ratio,
            p.read_abundance_pct,
            p.relative_abundance_pct,
            p.rpkm,
            p.mean_identity_pct,
            p.mean_similarity_pct,
            p.mean_hsp_length
        )?;
    }
    w.flush()?;
    Ok(())
}

/// Write full sample summary + per-genome profiles as pretty JSON.
pub fn write_abundance_json(path: impl AsRef<Path>, res: &MetaSampleResult) -> Result<()> {
    let file = File::create(path)?;
    let mut w = BufWriter::new(file);
    serde_json::to_writer_pretty(&mut w, res)
        .map_err(|e| GeckoError::InvalidParameter(format!("JSON serialization failed: {e}")))?;
    writeln!(w)?;
    w.flush()?;
    Ok(())
}

/// Write per-read alignment assignments as TSV (`<sample>.reads.tsv`).
pub fn write_read_assignments_tsv(path: impl AsRef<Path>, res: &MetaSampleResult) -> Result<()> {
    let file = File::create(path)?;
    let mut w = BufWriter::with_capacity(4 << 20, file);
    writeln!(
        w,
        "read_id\tread_index\tread_len\tgenome\tcontig\tread_start\tread_end\t\
         ref_start\tref_end\tstrand\tlength\tscore\tidentity_pct\tsimilarity_pct\t\
         read_coverage_pct\tweight\tunique_genome"
    )?;
    for a in &res.assignments {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{:.2}\t{:.2}\t{:.4}\t{}",
            a.read_id,
            a.read_index,
            a.read_len,
            a.genome_name,
            a.contig_id,
            a.read_start,
            a.read_end,
            a.ref_start,
            a.ref_end,
            a.strand,
            a.length,
            a.score,
            a.identity_pct,
            a.similarity_pct,
            a.read_coverage_pct,
            a.weight,
            a.is_unique_genome
        )?;
    }
    w.flush()?;
    Ok(())
}

/// Write a multi-sample `genomes × samples` relative-abundance & assigned-read matrix TSV
/// for differential abundance comparisons across metagenomic conditions.
pub fn write_abundance_matrix_tsv(
    path: impl AsRef<Path>,
    reference: &MetaReference,
    samples: &[MetaSampleResult],
    use_relative_pct: bool,
) -> Result<()> {
    let file = File::create(path)?;
    let mut w = BufWriter::with_capacity(4 << 20, file);
    write!(w, "genome\tgenome_len")?;
    for s in samples {
        write!(w, "\t{}", s.sample_name)?;
    }
    writeln!(w)?;

    // Build lookup maps (sample_idx -> genome_id -> value)
    let maps: Vec<BTreeMap<usize, f64>> = samples
        .iter()
        .map(|s| {
            s.profiles
                .iter()
                .map(|p| {
                    let v = if use_relative_pct {
                        p.relative_abundance_pct
                    } else {
                        p.assigned_reads
                    };
                    (p.genome_id, v)
                })
                .collect()
        })
        .collect();

    for info in &reference.genomes {
        write!(w, "{}\t{}", info.name, info.total_bases)?;
        for m in &maps {
            let v = m.get(&info.genome_id).copied().unwrap_or(0.0);
            write!(w, "\t{:.4}", v)?;
        }
        writeln!(w)?;
    }
    w.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pipeline::compare_genomes;
    use crate::seq::complement_base;

    fn make_dna(seed: u64, len: usize) -> String {
        const ACGT: [u8; 4] = *b"ACGT";
        let mut x = seed;
        let mut s = String::with_capacity(len);
        for _ in 0..len {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            s.push(ACGT[((x >> 62) & 3) as usize] as char);
        }
        s
    }

    fn rc_str(s: &str) -> String {
        s.bytes()
            .rev()
            .map(complement_base)
            .map(|b| b as char)
            .collect()
    }

    #[test]
    fn metagecko_profiles_unique_multi_and_low_abundance_breadth() {
        // Two distinct 4000 bp bacterial genomes sharing one 200 bp conserved 16S-like region
        let shared = make_dna(777, 200);
        let mut g1 = make_dna(101, 4000);
        let mut g2 = make_dna(202, 4000);
        g1.replace_range(1000..1200, &shared);
        g2.replace_range(2000..2200, &shared);

        let ref_genome = Genome::from_records("refs", [("taxon_A", &g1), ("taxon_B", &g2)]);
        let meta_ref = MetaReference::from_single_genome(ref_genome, "refs.fa".into(), true);

        // Reads:
        // - r0, r1: unique forward & reverse reads from taxon_A
        // - r2: unique reverse read from taxon_B
        // - r3: ambiguous read from the 200 bp conserved region (hits both taxon_A and taxon_B equally)
        // - r4: unmapped noise read
        let r0 = g1[100..250].to_string();
        let r1 = rc_str(&g1[2500..2650]);
        let r2 = rc_str(&g2[500..650]);
        let r3 = shared[10..160].to_string();
        let r4 = make_dna(999, 150);

        let batch_genome = Genome::from_records(
            "sample1",
            [
                ("r0", r0.as_str()),
                ("r1", r1.as_str()),
                ("r2", r2.as_str()),
                ("r3", r3.as_str()),
                ("r4", r4.as_str()),
            ],
        );
        let batch = ReadBatch {
            first_read_id: 0,
            raw_reads_in_batch: 5,
            short_reads_skipped: 0,
            genome: batch_genome,
        };

        let cfg = MetaConfig {
            compare: CompareConfig {
                min_length: 80,
                min_similarity: 80.0,
                word_length: 16,
                soft_mask: false,
                ..CompareConfig::default()
            },
            min_read_coverage: 70.0,
            multi_mapping: MultiMappingMode::Fractional,
            collect_assignments: true,
            collect_fragments: true,
            ..MetaConfig::default()
        };

        let cmp = compare_genomes(&batch.genome, &meta_ref.genome, &cfg.compare, None).unwrap();
        let mut acc = MetaSampleAccumulator::new("sample1", &meta_ref);
        acc.ingest_batch(&batch, &meta_ref, cmp, &cfg);
        let res = acc.finish(&meta_ref);

        assert_eq!(res.mapped_reads, 4);
        assert_eq!(res.unique_mapped_reads, 3);
        assert_eq!(res.multi_mapped_reads, 1);

        let prof_a = res
            .profiles
            .iter()
            .find(|p| p.genome_name == "taxon_A")
            .unwrap();
        let prof_b = res
            .profiles
            .iter()
            .find(|p| p.genome_name == "taxon_B")
            .unwrap();

        assert_eq!(prof_a.unique_reads, 2);
        assert_eq!(prof_a.multi_reads, 1);
        assert!((prof_a.assigned_reads - 2.5).abs() < 1e-6);
        assert_eq!(prof_a.covered_bases, 450); // 150 + 150 + 150 non-overlapping bases

        assert_eq!(prof_b.unique_reads, 1);
        assert_eq!(prof_b.multi_reads, 1);
        assert!((prof_b.assigned_reads - 1.5).abs() < 1e-6);
        assert_eq!(prof_b.covered_bases, 300); // 150 + 150 non-overlapping bases
    }
}
