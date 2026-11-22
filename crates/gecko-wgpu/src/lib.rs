//! Cross-vendor WGPU (Vulkan / Metal / DX12) GPU backend for GECKO.

mod gpu;
mod index;

use gecko_core::error::{GeckoError, Result};
use gecko_core::extend::CoverageMatrix;
use gecko_core::formats::frags::Fragment;
use gecko_core::pipeline::{CompareConfig, CompareResult, CompareStats};
use gecko_core::seq::Genome;
use gecko_core::stream::dedup_fragments;

use self::gpu::{WgpuContext, WgpuPreparedRef};
use self::index::WgpuRefIndex;

/// Summary of an available WGPU adapter.
#[derive(Debug, Clone)]
pub struct WgpuAdapterInfo {
    pub name: String,
    pub backend: String,
    pub is_software: bool,
}

/// Probe whether a WGPU GPU adapter is available.
pub fn probe_wgpu(allow_software: bool) -> std::result::Result<WgpuAdapterInfo, String> {
    let ctx = WgpuContext::init(allow_software)?;
    Ok(WgpuAdapterInfo {
        name: ctx.adapter_name,
        backend: ctx.adapter_backend,
        is_software: ctx.is_software,
    })
}

/// Reusable WGPU reference session (prepares `WgpuContext`, uploads the reference Hash-CSR
/// index + 2-bit sequences to GPU VRAM, and compiles WGSL pipelines once, then queries
/// unlimited metagenomic read batches or samples without re-indexing or re-uploading).
pub struct WgpuRefSession<'a> {
    ctx: WgpuContext,
    info: WgpuAdapterInfo,
    reference: &'a Genome,
    ref_rc: Option<Genome>,
    prep: WgpuPreparedRef,
    word_length: usize,
    soft_mask: bool,
    compare_reverse: bool,
}

impl<'a> WgpuRefSession<'a> {
    /// Initialize WGPU adapter, build the canonical Hash-CSR reference index, and upload
    /// persistent GPU buffers + compiled WGSL compute pipelines once.
    pub fn prepare(
        reference: &'a Genome,
        config: &CompareConfig,
        allow_software: bool,
    ) -> Result<Self> {
        config.validate()?;
        if !(8..=32).contains(&config.word_length) {
            return Err(GeckoError::InvalidParameter(format!(
                "WGPU backend requires word_length in 8..=32 (got {})",
                config.word_length
            )));
        }
        if reference.seq.len >= (1usize << 31) {
            return Err(GeckoError::InvalidParameter(
                "WGPU backend requires reference < 2.14 Gb".into(),
            ));
        }

        let ctx = WgpuContext::init(allow_software).map_err(GeckoError::InvalidParameter)?;
        let info = WgpuAdapterInfo {
            name: ctx.adapter_name.clone(),
            backend: ctx.adapter_backend.clone(),
            is_software: ctx.is_software,
        };

        let (idx, ref_rc) = if config.compare_reverse {
            rayon::join(
                || {
                    WgpuRefIndex::build(
                        reference,
                        config.word_length,
                        config.freq_limit as u32,
                        config.soft_mask,
                    )
                },
                || Some(reference.reverse_complement()),
            )
        } else {
            (
                WgpuRefIndex::build(
                    reference,
                    config.word_length,
                    config.freq_limit as u32,
                    config.soft_mask,
                ),
                None,
            )
        };

        let prep = ctx
            .prepare_ref(reference, ref_rc.as_ref(), &idx, config.soft_mask)
            .map_err(GeckoError::InvalidParameter)?;

        Ok(Self {
            ctx,
            info,
            reference,
            ref_rc,
            prep,
            word_length: config.word_length,
            soft_mask: config.soft_mask,
            compare_reverse: config.compare_reverse,
        })
    }

    pub fn info(&self) -> &WgpuAdapterInfo {
        &self.info
    }

    /// Query a single genome or metagenomic read batch against the GPU-resident reference.
    pub fn query(&self, query: &Genome, config: &CompareConfig) -> Result<CompareResult> {
        if query.seq.len >= (1usize << 31) {
            return Err(GeckoError::InvalidParameter(
                "WGPU backend requires query batch < 2.14 Gb".into(),
            ));
        }

        let min_len = if config.fixed_length {
            config.min_length
        } else {
            let shorter = query.total_bases.min(self.reference.total_bases);
            ((shorter as f64) * (config.min_length as f64 / 100.0)) as u64
        };

        let gpu_stats = self
            .ctx
            .query_prepared_ref(
                query,
                &self.prep,
                self.word_length as u32,
                min_len as u32,
                config.min_similarity,
                self.soft_mask,
                self.compare_reverse,
                config.scoring,
            )
            .map_err(GeckoError::InvalidParameter)?;

        let mut frags_f = Vec::new();
        let mut frags_r = Vec::new();
        let mut coverage = CoverageMatrix::new();
        let n_seqs_y = self.reference.num_records() as u64;
        let match_sc = config.scoring.match_score as f64;

        for gf in &gpu_stats.raw_frags {
            let is_rev = gf.strand != 0;
            let seq_x = query.seq_index_for_data_pos(gf.x_start as u64) as u64;
            let seq_y_buf = if is_rev {
                self.ref_rc
                    .as_ref()
                    .map_or(0, |rc| rc.seq_index_for_data_pos(gf.y_start as u64) as u64)
            } else {
                self.reference.seq_index_for_data_pos(gf.y_start as u64) as u64
            };
            let length = (gf.x_end - gf.x_start + 1) as u64;
            let sim = ((gf.score as f64) * 100.0 / ((length as f64) * match_sc)) as f32;

            coverage.record(length, sim);

            let frag = Fragment {
                diag: (gf.x_start as i64) - (gf.y_start as i64),
                x_start: (gf.x_start as u64) - seq_x,
                y_start: (gf.y_start as u64) - seq_y_buf,
                x_end: (gf.x_end as u64) - seq_x,
                y_end: (gf.y_end as u64) - seq_y_buf,
                length,
                ident: gf.ident as u64,
                score: gf.score as u64,
                similarity: sim,
                seq_x,
                seq_y: if is_rev {
                    n_seqs_y - seq_y_buf - 1
                } else {
                    seq_y_buf
                },
                block: 0,
                strand: if is_rev { 'r' } else { 'f' },
            };
            if is_rev {
                frags_r.push(frag);
            } else {
                frags_f.push(frag);
            }
        }

        dedup_fragments(&mut frags_f);
        dedup_fragments(&mut frags_r);

        let forward_fragments = frags_f.len() as u64;
        let reverse_fragments = frags_r.len() as u64;
        let mut fragments = frags_f;
        fragments.extend(frags_r);

        let stats = CompareStats {
            word_matches: gpu_stats.leftmax_hits,
            raw_hits: gpu_stats.leftmax_hits,
            filtered_hits: gpu_stats.leftmax_hits,
            hits_used: gpu_stats.extended_hits,
            forward_fragments,
            reverse_fragments,
            spilled_runs: 0,
        };

        Ok(CompareResult {
            fragments,
            coverage,
            stats,
        })
    }
}

/// Run a full pairwise genome comparison on the WGPU compute backend.
pub fn compare_wgpu(
    query: &Genome,
    reference: &Genome,
    config: &CompareConfig,
    allow_software: bool,
) -> Result<(CompareResult, WgpuAdapterInfo)> {
    let session = WgpuRefSession::prepare(reference, config, allow_software)?;
    let info = session.info().clone();
    let result = session.query(query, config)?;
    Ok((result, info))
}

#[cfg(test)]
mod tests {
    use super::*;
    use gecko_core::pipeline::compare_genomes;
    use gecko_core::seq::complement_base;

    #[test]
    fn freq_bounded_interior_seeds_match_cpu_both_strands() {
        // Construct a 120 bp block: 55 bp left flank + 10 bp unique core + 55 bp right flank.
        // Repeat the left and right flanks 4 times in the reference (> freq_limit = 2) so all
        // boundary k-mers of the exact match are dropped by freq_limit, leaving ONLY interior
        // k-mers whose left & right flanking bases still match.
        let left_flank = "ACGTCCGGTTAAACGTTTGGCCAAAACCCGGGTTTACGTACGGTTAACCGGAATT";
        let right_flank = "GCTAGCTAGGCCTTAAACCGGTTTGGCCAATTCCGGATGCTAGCTAGGCCTTAAA";
        let unique_fwd = "TGCATGCAAT";
        let unique_rev = "CGATCGATCG";

        let rc_str = |s: &str| -> String {
            s.bytes()
                .rev()
                .map(complement_base)
                .map(|b| b as char)
                .collect()
        };
        let left_flank_rc = rc_str(left_flank);
        let right_flank_rc = rc_str(right_flank);
        let hsp_fwd = format!("{left_flank}{unique_fwd}{right_flank}");
        let hsp_rev = format!("{left_flank}{unique_rev}{right_flank}");
        let hsp_rev_rc = rc_str(&hsp_rev);

        let query = Genome::from_records(
            "q",
            [(
                "q0",
                format!("NNNNNNNNNN{hsp_fwd}NNNNNNNNNN{hsp_rev}NNNNNNNNNN"),
            )],
        );
        let reference = Genome::from_records(
            "r",
            [(
                "r0",
                format!(
                    "NNNN{left_flank}NNNN{left_flank}NNNN{left_flank}NNNN\
                     {right_flank}NNNN{right_flank}NNNN{right_flank}NNNN\
                     {left_flank_rc}NNNN{left_flank_rc}NNNN{left_flank_rc}NNNN\
                     {right_flank_rc}NNNN{right_flank_rc}NNNN{right_flank_rc}NNNN\
                     {hsp_fwd}NNNN{hsp_rev_rc}NNNN"
                ),
            )],
        );

        let cfg = CompareConfig {
            min_length: 80,
            min_similarity: 75.0,
            word_length: 16,
            freq_limit: 2,
            ..Default::default()
        };

        let cpu_res = compare_genomes(&query, &reference, &cfg, None).unwrap();
        assert_eq!(cpu_res.stats.forward_fragments, 1);
        assert_eq!(cpu_res.stats.reverse_fragments, 1);

        // Run WGPU if any adapter (including software Lavapipe) is present on this host.
        if probe_wgpu(true).is_ok() {
            let (gpu_res, _) = compare_wgpu(&query, &reference, &cfg, true).unwrap();
            assert_eq!(gpu_res.fragments, cpu_res.fragments);
        }
    }

    #[test]
    fn multi_fasta_softmask_n_and_iupac_match_cpu() {
        let block_a = "acgtacgtacgtACGTCCGGTTAAACGTTTGGCCAAAACCCGGGTTTACGTACGGTTAACCGGAATTNNNNTTGGCCAACCGGTTAAccggttaa";
        let block_b =
            "GCTAGCTAGGCCTTAAACCGGTTTGGCCAATTCCGGATRYGCTAGCTAGGCCTTAAACCGGTTTGGCCAATTCCGGAT";
        let block_b_rc: String = block_b
            .bytes()
            .rev()
            .map(complement_base)
            .map(|b| b as char)
            .collect();

        let query = Genome::from_records(
            "q",
            [
                ("q0", format!("TTTTTTTTTT{block_a}GGGGGGGGGG")),
                ("q1", format!("CCCCCCCCCC{block_b}AAAAAAAAAA")),
            ],
        );
        let reference = Genome::from_records(
            "r",
            [
                ("r0", format!("AAAAAAAAAA{block_b_rc}TTTTTTTTTT")),
                ("r1", format!("GGGGGGGGGG{block_a}CCCCCCCCCC")),
            ],
        );

        for soft_mask in [true, false] {
            let cfg = CompareConfig {
                min_length: 40,
                min_similarity: 60.0,
                word_length: 12,
                soft_mask,
                ..Default::default()
            };
            let cpu_res = compare_genomes(&query, &reference, &cfg, None).unwrap();
            assert_eq!(cpu_res.stats.forward_fragments, 1);
            assert_eq!(cpu_res.stats.reverse_fragments, 1);

            if probe_wgpu(true).is_ok() {
                let (gpu_res, _) = compare_wgpu(&query, &reference, &cfg, true).unwrap();
                let mut f_cpu = cpu_res.fragments.clone();
                let mut f_gpu = gpu_res.fragments.clone();
                let key = |f: &Fragment| (f.strand, f.seq_x, f.seq_y, f.diag, f.x_start, f.x_end);
                f_cpu.sort_by_key(key);
                f_gpu.sort_by_key(key);
                assert_eq!(
                    f_gpu, f_cpu,
                    "2-bit WGPU mismatch for soft_mask={soft_mask}"
                );
            }
        }
    }
}
