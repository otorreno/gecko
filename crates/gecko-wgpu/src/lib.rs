//! Cross-vendor WGPU (Vulkan / Metal / DX12) GPU backend for GECKO.

mod gpu;
mod index;

use gecko_core::error::{GeckoError, Result};
use gecko_core::extend::CoverageMatrix;
use gecko_core::formats::frags::Fragment;
use gecko_core::pipeline::{CompareConfig, CompareResult, CompareStats};
use gecko_core::seq::Genome;
use gecko_core::stream::dedup_fragments;

use self::gpu::WgpuContext;
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

/// Run a full pairwise genome comparison on the WGPU compute backend.
pub fn compare_wgpu(
    query: &Genome,
    reference: &Genome,
    config: &CompareConfig,
    allow_software: bool,
) -> Result<(CompareResult, WgpuAdapterInfo)> {
    config.validate()?;
    if !(8..=32).contains(&config.word_length) {
        return Err(GeckoError::InvalidParameter(format!(
            "WGPU backend requires word_length in 8..=32 (got {})",
            config.word_length
        )));
    }
    if query.data.len() >= (1usize << 31) || reference.data.len() >= (1usize << 31) {
        return Err(GeckoError::InvalidParameter(
            "WGPU backend requires sequences < 2.14 Gb".into(),
        ));
    }

    let ctx = WgpuContext::init(allow_software).map_err(GeckoError::InvalidParameter)?;
    let info = WgpuAdapterInfo {
        name: ctx.adapter_name.clone(),
        backend: ctx.adapter_backend.clone(),
        is_software: ctx.is_software,
    };

    let min_len = if config.fixed_length {
        config.min_length
    } else {
        let shorter = query.total_bases.min(reference.total_bases);
        ((shorter as f64) * (config.min_length as f64 / 100.0)) as u64
    };

    // Build canonical Hash-CSR index and reverse complement concurrently
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

    let gpu_stats = ctx
        .run_query_and_extend(
            query,
            reference,
            ref_rc.as_ref(),
            &idx,
            config.word_length as u32,
            min_len as u32,
            config.min_similarity,
            config.soft_mask,
            config.compare_reverse,
            config.scoring,
        )
        .map_err(GeckoError::InvalidParameter)?;

    let mut frags_f = Vec::new();
    let mut frags_r = Vec::new();
    let mut coverage = CoverageMatrix::new();
    let n_seqs_y = reference.num_records() as u64;
    let match_sc = config.scoring.match_score as f64;

    for gf in &gpu_stats.raw_frags {
        let is_rev = gf.strand != 0;
        let seq_x = query.seq_index_for_data_pos(gf.x_start as u64) as u64;
        let seq_y_buf = if is_rev {
            ref_rc
                .as_ref()
                .map_or(0, |rc| rc.seq_index_for_data_pos(gf.y_start as u64) as u64)
        } else {
            reference.seq_index_for_data_pos(gf.y_start as u64) as u64
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

    Ok((
        CompareResult {
            fragments,
            coverage,
            stats,
        },
        info,
    ))
}
