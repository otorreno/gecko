//! Safe Rust wrapper around the optional native NVIDIA CUDA + CUB GECKO backend.
//!
//! Compile with `--features cuda-native` (or `gecko-cli --features cuda`) on machines
//! with the NVIDIA CUDA Toolkit (`nvcc`) installed. Without that feature, this crate
//! compiles as a lightweight stub so pure-Rust/WGPU builds succeed anywhere.

use std::path::Path;

use gecko_core::error::Result;
use gecko_core::pipeline::{CompareConfig, CompareResult};
use gecko_core::seq::Genome;

/// Hardware summary for a detected NVIDIA CUDA GPU.
#[derive(Debug, Clone)]
pub struct CudaDeviceInfo {
    pub name: String,
    pub sm_major: i32,
    pub sm_minor: i32,
    pub vram_bytes: u64,
}

/// Telemetry counters for the process-wide persistent [`CudaScratchPool`].
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct CudaScratchPoolStats {
    pub hits: u64,
    pub misses: u64,
    pub recycles: u64,
    pub evictions: u64,
    pub idle_count: usize,
    pub cached_bytes: u64,
}

/// RAII guard that releases all idle pooled CUDA scratch buffers upon [`Drop`].
pub struct CudaScratchPoolGuard;

impl Drop for CudaScratchPoolGuard {
    fn drop(&mut self) {
        CudaScratchPool::release();
    }
}

/// Persistent process-wide CUDA scratch workspace pool (`W3-F15a`).
///
/// Automatically checked out in [`CudaPreparedRef::build`] and recycled in [`Drop`] after
/// draining device streams and verifying CUDA health, eliminating repeated ~608 MiB scratch
/// allocations across `all-vs-all` reference transitions.
pub struct CudaScratchPool;

#[cfg(not(has_cuda_kernels))]
static STUB_POOL_ENABLED: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(true);

impl CudaScratchPool {
    /// Return cumulative hit/miss/recycle/eviction and idle residency counters.
    pub fn stats() -> CudaScratchPoolStats {
        #[cfg(has_cuda_kernels)]
        {
            ffi::scratch_pool_stats()
        }
        #[cfg(not(has_cuda_kernels))]
        {
            CudaScratchPoolStats::default()
        }
    }

    /// Release all idle scratch workspaces currently cached in the pool.
    /// Safe no-op when the pool is empty or no CUDA GPU is attached.
    pub fn release() {
        #[cfg(has_cuda_kernels)]
        {
            ffi::scratch_pool_release();
        }
    }

    /// Returns `true` if scratch pooling is currently enabled.
    pub fn is_enabled() -> bool {
        #[cfg(has_cuda_kernels)]
        {
            ffi::scratch_pool_is_enabled()
        }
        #[cfg(not(has_cuda_kernels))]
        {
            STUB_POOL_ENABLED.load(std::sync::atomic::Ordering::Relaxed)
        }
    }

    /// Enable or disable scratch pooling. Disabling immediately releases all idle workspaces.
    pub fn set_enabled(enabled: bool) {
        #[cfg(has_cuda_kernels)]
        {
            ffi::scratch_pool_set_enabled(enabled);
        }
        #[cfg(not(has_cuda_kernels))]
        {
            STUB_POOL_ENABLED.store(enabled, std::sync::atomic::Ordering::Relaxed);
        }
    }

    /// Configure the maximum number of idle scratch workspaces retained in the pool.
    pub fn set_max_idle(max_idle: usize) {
        #[cfg(has_cuda_kernels)]
        {
            ffi::scratch_pool_set_max_idle(max_idle);
        }
        #[cfg(not(has_cuda_kernels))]
        {
            let _ = max_idle;
        }
    }

    /// Return an RAII guard that calls [`CudaScratchPool::release`] when dropped.
    pub fn scoped_release() -> CudaScratchPoolGuard {
        CudaScratchPoolGuard
    }
}

/// Returns `true` if native CUDA kernels were compiled into this binary.
pub const fn is_compiled_in() -> bool {
    cfg!(has_cuda_kernels)
}

/// 8-byte half-open `[start, end)` run descriptor shared with `cuda_common.cuh::GeckoRun`.
#[cfg(any(has_cuda_kernels, test))]
#[repr(C)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct GeckoRun {
    pub start: u32,
    pub end: u32,
}

#[cfg(any(has_cuda_kernels, test))]
const _: () = assert!(std::mem::size_of::<GeckoRun>() == 8);

/// 512 u32 bitmap words = 16,384 bp per chunk (32-bp aligned splits avoid inter-chunk atomic contention).
#[cfg(any(has_cuda_kernels, test))]
pub(crate) const MAX_GPU_RUN_CHUNK_BP: u64 = 16_384;

/// Convert sorted half-open `(u64, u64)` runs to `GeckoRun`, splitting giant runs (e.g. centromeric
/// N-runs) on exact 32-bp bitmap-word boundaries so GPU SMs stay saturated with zero boundary contention.
#[cfg(any(has_cuda_kernels, test))]
pub(crate) fn pack_runs_for_gpu(runs: &[(u64, u64)], seq_len: u64) -> Vec<GeckoRun> {
    let mut out = Vec::with_capacity(runs.len());
    for &(raw_s, raw_e) in runs {
        let mut s = raw_s.min(seq_len);
        let e = raw_e.min(seq_len);
        if s >= e {
            continue;
        }
        while e - s > MAX_GPU_RUN_CHUNK_BP {
            let cut = (((s >> 5) + (MAX_GPU_RUN_CHUNK_BP >> 5)) << 5).min(e);
            out.push(GeckoRun {
                start: s as u32,
                end: cut as u32,
            });
            s = cut;
        }
        if s < e {
            out.push(GeckoRun {
                start: s as u32,
                end: e as u32,
            });
        }
    }
    out
}

#[cfg(has_cuda_kernels)]
mod ffi {
    use std::ffi::{CStr, c_char, c_int};
    use std::time::Instant;

    use gecko_core::error::{GeckoError, Result};
    use gecko_core::extend::CoverageMatrix;
    use gecko_core::formats::frags::Fragment;
    use gecko_core::pipeline::{CompareConfig, CompareResult, CompareStats, SeedingMode};
    use gecko_core::seq::Genome;
    use gecko_core::stream::dedup_fragments;

    use super::{CudaDeviceInfo, CudaScratchPoolStats, GeckoRun, pack_runs_for_gpu};

    #[inline]
    fn seeding_mode_to_u8(mode: SeedingMode) -> u8 {
        match mode {
            SeedingMode::Exact => 0,
            SeedingMode::Minimizer => 1,
            SeedingMode::Auto => 2,
        }
    }

    #[repr(C)]
    pub struct GeckoCudaRefHandle {
        _private: [u8; 0],
    }

    #[repr(C)]
    pub struct GeckoCudaUploadedQuery {
        _private: [u8; 0],
    }

    #[repr(C)]
    #[derive(Clone, Copy)]
    pub struct GpuFrag {
        pub x_start: u32,
        pub y_start: u32,
        pub length: u32,
        pub ident: u32,
        pub score: u32,
        pub flags: u32,
    }
    const _: () = assert!(std::mem::size_of::<GpuFrag>() == 24);

    #[repr(C)]
    #[derive(Clone, Copy)]
    pub struct GeckoCudaConfig {
        pub seed_mask: u64,
        pub freq_limit: u32,
        pub min_len: u32,
        pub two_hit_window: u32,
        pub sim_th: f32,
        pub wl: i32,
        pub match_score: i32,
        pub mismatch_penalty: i32,
        pub n_penalty: i32,
        pub soft_mask: u8,
        pub compare_reverse: u8,
        pub seeding_mode: u8,
        pub _pad0: u8,
        pub _pad1: u32,
    }
    const _: () = assert!(std::mem::size_of::<GeckoCudaConfig>() == 48);

    impl GeckoCudaConfig {
        fn from_compare_config(config: &CompareConfig, min_len: u64) -> Self {
            Self {
                seed_mask: config.spaced_mask_2bit().unwrap_or(0),
                freq_limit: config.freq_limit.min(u32::MAX as u64) as u32,
                min_len: min_len.min(u32::MAX as u64) as u32,
                two_hit_window: config.two_hit_window.unwrap_or(0),
                sim_th: config.min_similarity,
                wl: config.word_length as i32,
                match_score: config.scoring.match_score as i32,
                mismatch_penalty: config.scoring.mismatch_penalty as i32,
                n_penalty: config.scoring.n_penalty as i32,
                soft_mask: u8::from(config.soft_mask),
                compare_reverse: u8::from(config.compare_reverse),
                seeding_mode: seeding_mode_to_u8(config.seeding_mode),
                _pad0: 0,
                _pad1: 0,
            }
        }
    }

    #[repr(C)]
    #[derive(Default)]
    pub struct GeckoCudaScratchPoolStats {
        pub hits: u64,
        pub misses: u64,
        pub recycles: u64,
        pub evictions: u64,
        pub idle_count: u64,
        pub cached_bytes: u64,
    }
    const _: () = assert!(std::mem::size_of::<GeckoCudaScratchPoolStats>() == 48);

    #[repr(C)]
    #[derive(Default)]
    pub struct GeckoCudaResult {
        pub frags: *mut GpuFrag,
        pub n_frags: u64,
        pub leftmax_hits: u64,
        pub filtered_hits: u64,
        pub extended_hits: u64,
        pub index_kept: u32,
        pub index_valid: u32,
        pub bucket_bits: i32,
        pub _pad: u32,
    }
    const _: () = assert!(std::mem::size_of::<GeckoCudaResult>() == 56);

    #[repr(C)]
    pub struct GeckoPackedSeqHost {
        pub bits: *const u64,
        pub bits_words: u64,
        pub special_runs: *const GeckoRun,
        pub soft_runs: *const GeckoRun,
        pub sep_pos: *const u32,
        pub iupac_pos: *const u32,
        pub iupac_byte: *const u8,
        pub len: u64,
        pub n_special_runs: u32,
        pub n_soft_runs: u32,
        pub n_seps: u32,
        pub n_iupac: u32,
        pub soft_mask: u32,
        pub _pad: u32,
    }
    const _: () = assert!(std::mem::size_of::<GeckoPackedSeqHost>() == 88);

    /// Owns compact run descriptors & u32 special index slices for the duration of a CUDA upload
    /// call while borrowing `PackedSeq::bits` zero-copy (zero 60-120 MiB host dense bitmap allocs!).
    struct HostPackedView<'a> {
        bits: &'a [u64],
        special_runs: Vec<GeckoRun>,
        soft_runs: Vec<GeckoRun>,
        sep_pos: Vec<u32>,
        iupac_pos: Vec<u32>,
        iupac_byte: Vec<u8>,
        len: u64,
        soft_mask: bool,
    }

    impl<'a> HostPackedView<'a> {
        fn new(genome: &'a Genome, soft_mask: bool) -> Self {
            let seq = &genome.seq;
            let len = seq.len as u64;
            let special_runs = pack_runs_for_gpu(&seq.special_runs, len);
            let soft_runs = if soft_mask && !seq.soft_runs.is_empty() {
                pack_runs_for_gpu(&seq.soft_runs, len)
            } else {
                Vec::new()
            };
            let sep_pos: Vec<u32> = seq.sep_pos.iter().map(|&p| p as u32).collect();
            let mut iupac_pos = Vec::with_capacity(seq.iupac.len());
            let mut iupac_byte = Vec::with_capacity(seq.iupac.len());
            for &(p, b) in &seq.iupac {
                iupac_pos.push(p as u32);
                iupac_byte.push(b);
            }
            Self {
                bits: &seq.bits,
                special_runs,
                soft_runs,
                sep_pos,
                iupac_pos,
                iupac_byte,
                len,
                soft_mask,
            }
        }

        fn as_ffi(&self) -> GeckoPackedSeqHost {
            GeckoPackedSeqHost {
                bits: self.bits.as_ptr(),
                bits_words: self.bits.len() as u64,
                special_runs: if self.special_runs.is_empty() {
                    std::ptr::null()
                } else {
                    self.special_runs.as_ptr()
                },
                soft_runs: if self.soft_runs.is_empty() {
                    std::ptr::null()
                } else {
                    self.soft_runs.as_ptr()
                },
                sep_pos: if self.sep_pos.is_empty() {
                    std::ptr::null()
                } else {
                    self.sep_pos.as_ptr()
                },
                iupac_pos: if self.iupac_pos.is_empty() {
                    std::ptr::null()
                } else {
                    self.iupac_pos.as_ptr()
                },
                iupac_byte: if self.iupac_byte.is_empty() {
                    std::ptr::null()
                } else {
                    self.iupac_byte.as_ptr()
                },
                len: self.len,
                n_special_runs: self.special_runs.len() as u32,
                n_soft_runs: self.soft_runs.len() as u32,
                n_seps: self.sep_pos.len() as u32,
                n_iupac: self.iupac_pos.len() as u32,
                soft_mask: u32::from(self.soft_mask),
                _pad: 0,
            }
        }
    }

    unsafe extern "C" {
        fn gecko_cuda_probe(
            name_buf: *mut c_char,
            name_cap: usize,
            sm_major: *mut c_int,
            sm_minor: *mut c_int,
            vram_bytes: *mut u64,
        ) -> c_int;

        fn gecko_cuda_warmup() -> c_int;

        fn gecko_cuda_upload_query_async(
            query_seq: *const GeckoPackedSeqHost,
        ) -> *mut GeckoCudaUploadedQuery;

        fn gecko_cuda_free_uploaded_query(uq: *mut GeckoCudaUploadedQuery);

        fn gecko_cuda_prepare_ref(
            ref_seq: *const GeckoPackedSeqHost,
            cfg: *const GeckoCudaConfig,
        ) -> *mut GeckoCudaRefHandle;

        fn gecko_cuda_ref_is_minimizer(handle: *const GeckoCudaRefHandle) -> c_int;

        fn gecko_cuda_query_and_extend(
            handle: *mut GeckoCudaRefHandle,
            query_seq: *const GeckoPackedSeqHost,
            cfg: *const GeckoCudaConfig,
            out: *mut GeckoCudaResult,
        ) -> c_int;

        fn gecko_cuda_query_and_extend_uploaded(
            handle: *mut GeckoCudaRefHandle,
            uq: *mut GeckoCudaUploadedQuery,
            cfg: *const GeckoCudaConfig,
            out: *mut GeckoCudaResult,
        ) -> c_int;

        fn gecko_cuda_free_ref(handle: *mut GeckoCudaRefHandle);
        fn gecko_cuda_free_frags(frags: *mut GpuFrag);

        fn gecko_cuda_scratch_pool_release();
        fn gecko_cuda_scratch_pool_stats(out: *mut GeckoCudaScratchPoolStats);
        fn gecko_cuda_scratch_pool_set_enabled(enabled: c_int);
        fn gecko_cuda_scratch_pool_is_enabled() -> c_int;
        fn gecko_cuda_scratch_pool_set_max_idle(max_idle: usize);
    }

    pub fn scratch_pool_release() {
        unsafe { gecko_cuda_scratch_pool_release() };
    }

    pub fn scratch_pool_stats() -> CudaScratchPoolStats {
        let mut raw = GeckoCudaScratchPoolStats::default();
        unsafe { gecko_cuda_scratch_pool_stats(&mut raw) };
        CudaScratchPoolStats {
            hits: raw.hits,
            misses: raw.misses,
            recycles: raw.recycles,
            evictions: raw.evictions,
            idle_count: raw.idle_count as usize,
            cached_bytes: raw.cached_bytes,
        }
    }

    pub fn scratch_pool_set_enabled(enabled: bool) {
        unsafe { gecko_cuda_scratch_pool_set_enabled(c_int::from(enabled)) };
    }

    pub fn scratch_pool_is_enabled() -> bool {
        unsafe { gecko_cuda_scratch_pool_is_enabled() != 0 }
    }

    pub fn scratch_pool_set_max_idle(max_idle: usize) {
        unsafe { gecko_cuda_scratch_pool_set_max_idle(max_idle) };
    }

    pub struct RefGuard(*mut GeckoCudaRefHandle);
    impl Drop for RefGuard {
        fn drop(&mut self) {
            if !self.0.is_null() {
                unsafe { gecko_cuda_free_ref(self.0) };
            }
        }
    }

    pub struct UploadedQueryGuard(*mut GeckoCudaUploadedQuery);
    // Safety: CUDA Runtime handles (streams, events, device/pinned pointers on device 0's
    // primary context) are valid across all host threads; ownership is exclusive.
    unsafe impl Send for UploadedQueryGuard {}

    impl Drop for UploadedQueryGuard {
        fn drop(&mut self) {
            if !self.0.is_null() {
                unsafe { gecko_cuda_free_uploaded_query(self.0) };
            }
        }
    }

    pub fn probe_cuda() -> std::result::Result<CudaDeviceInfo, String> {
        let mut name_buf = [0 as c_char; 256];
        let mut sm_major: c_int = 0;
        let mut sm_minor: c_int = 0;
        let mut vram_bytes: u64 = 0;
        let rc = unsafe {
            gecko_cuda_probe(
                name_buf.as_mut_ptr(),
                name_buf.len(),
                &mut sm_major,
                &mut sm_minor,
                &mut vram_bytes,
            )
        };
        if rc != 0 {
            return Err("No CUDA-capable GPU detected".to_string());
        }
        let name = unsafe { CStr::from_ptr(name_buf.as_ptr()) }
            .to_string_lossy()
            .into_owned();
        Ok(CudaDeviceInfo {
            name,
            sm_major,
            sm_minor,
            vram_bytes,
        })
    }

    /// O(n_seqs) reverse-complement sequence start offsets for resolving `seq_y` IDs
    /// without building `Y_rc` on the host CPU.
    fn build_rc_starts(reference: &Genome) -> Vec<u64> {
        let n_seqs = reference.records.len();
        let buf_len = reference.seq.len as u64;
        let mut starts = Vec::with_capacity(n_seqs);
        starts.push(0);
        for k in (1..n_seqs).rev() {
            let orig_start = reference.records[k].start_with_sep;
            starts.push(buf_len - orig_start + 1);
        }
        starts
    }

    #[inline]
    fn find_rc_seq_id(rc_starts: &[u64], pos: u64) -> u64 {
        if rc_starts.len() <= 1 {
            return 0;
        }
        let idx = rc_starts.partition_point(|&s| s <= pos);
        idx.saturating_sub(1) as u64
    }

    fn validate_cuda_config(config: &CompareConfig) -> Result<()> {
        config.validate()?;
        if !(8..=32).contains(&config.word_length) {
            return Err(GeckoError::InvalidParameter(format!(
                "CUDA backend requires word_length in 8..=32 (got {})",
                config.word_length
            )));
        }
        Ok(())
    }

    fn finalize_cuda_result(
        query: &Genome,
        reference: &Genome,
        config: &CompareConfig,
        raw: &GeckoCudaResult,
    ) -> CompareResult {
        let rc_starts = build_rc_starts(reference);
        let n_seqs_y = reference.num_records() as u64;
        let match_sc = config.scoring.match_score as f64;

        let mut frags_f = Vec::new();
        let mut frags_r = Vec::new();
        let mut coverage = CoverageMatrix::new();

        let slice = if !raw.frags.is_null() && raw.n_frags > 0 {
            unsafe { std::slice::from_raw_parts(raw.frags, raw.n_frags as usize) }
        } else {
            &[]
        };

        for gf in slice {
            let is_rev = (gf.flags & 1) != 0;
            let x_start_raw = gf.x_start as u64;
            let y_start_raw = gf.y_start as u64;
            let length = gf.length as u64;
            let x_end_raw = x_start_raw + length - 1;
            let y_end_raw = y_start_raw + length - 1;
            let diag = (x_start_raw as i64) - (y_start_raw as i64);

            let seq_x = query.seq_index_for_data_pos(x_start_raw) as u64;
            let seq_y_buf = if is_rev {
                find_rc_seq_id(&rc_starts, y_start_raw)
            } else {
                reference.seq_index_for_data_pos(y_start_raw) as u64
            };
            let sim = ((gf.score as f64) * 100.0 / ((length as f64) * match_sc)) as f32;

            let frag = Fragment {
                diag,
                x_start: x_start_raw - seq_x,
                y_start: y_start_raw - seq_y_buf,
                x_end: x_end_raw - seq_x,
                y_end: y_end_raw - seq_y_buf,
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
        for f in &fragments {
            coverage.record(f.length, f.similarity);
        }

        let stats = CompareStats {
            word_matches: raw.leftmax_hits,
            raw_hits: raw.leftmax_hits,
            filtered_hits: raw.filtered_hits,
            hits_used: raw.extended_hits,
            forward_fragments,
            reverse_fragments,
            spilled_runs: 0,
        };

        CompareResult {
            fragments,
            coverage,
            stats,
        }
    }

    pub struct PreparedRefInner<'a> {
        handle: RefGuard,
        reference: &'a Genome,
        soft_mask: bool,
        is_minimizer: bool,
    }

    impl<'a> PreparedRefInner<'a> {
        pub fn build(
            reference: &'a Genome,
            config: &CompareConfig,
        ) -> Result<(Self, CudaDeviceInfo)> {
            let prof = gecko_core::profile_enabled();
            validate_cuda_config(config)?;
            if reference.seq.len >= (1usize << 31) {
                return Err(GeckoError::InvalidParameter(
                    "CUDA backend requires reference < 2.14 Gb".into(),
                ));
            }
            let t_probe = prof.then(Instant::now);
            let info = probe_cuda().map_err(GeckoError::InvalidParameter)?;
            let probe_ms = t_probe.map_or(0.0, |t| t.elapsed().as_secs_f64() * 1000.0);

            let t_view = prof.then(Instant::now);
            let ref_view = HostPackedView::new(reference, config.soft_mask);
            let ref_ffi = ref_view.as_ffi();
            let view_ref_ms = t_view.map_or(0.0, |t| t.elapsed().as_secs_f64() * 1000.0);

            let cuda_cfg = GeckoCudaConfig::from_compare_config(config, 0);
            let t_prep = prof.then(Instant::now);
            let handle = RefGuard(unsafe { gecko_cuda_prepare_ref(&ref_ffi, &cuda_cfg) });
            let prepare_ref_ms = t_prep.map_or(0.0, |t| t.elapsed().as_secs_f64() * 1000.0);
            if prof {
                eprintln!(
                    "[gecko-cuda-profile] host_prepare_ref probe_ms={probe_ms:.2} view_ref_ms={view_ref_ms:.2} prepare_ref_ms={prepare_ref_ms:.2}"
                );
            }
            if handle.0.is_null() {
                return Err(GeckoError::InvalidParameter(
                    "CUDA reference index preparation failed".into(),
                ));
            }
            let is_minimizer = unsafe { gecko_cuda_ref_is_minimizer(handle.0) != 0 };
            Ok((
                Self {
                    handle,
                    reference,
                    soft_mask: config.soft_mask,
                    is_minimizer,
                },
                info,
            ))
        }

        pub fn query(&mut self, query: &Genome, config: &CompareConfig) -> Result<CompareResult> {
            let prof = gecko_core::profile_enabled();
            validate_cuda_config(config)?;
            if self.is_minimizer && config.two_hit_window.is_some() {
                return Err(GeckoError::InvalidParameter(
                    "cannot query a minimizer CudaPreparedRef with --two-hit".into(),
                ));
            }
            if query.seq.len >= (1usize << 31) {
                return Err(GeckoError::InvalidParameter(
                    "CUDA backend requires query < 2.14 Gb".into(),
                ));
            }
            let min_len = if config.fixed_length {
                config.min_length
            } else {
                let shorter = query.total_bases.min(self.reference.total_bases);
                ((shorter as f64) * (config.min_length as f64 / 100.0)) as u64
            };
            let cuda_cfg = GeckoCudaConfig::from_compare_config(config, min_len);
            let mut raw = GeckoCudaResult::default();
            let t_vq = prof.then(Instant::now);
            let q_view = HostPackedView::new(query, self.soft_mask);
            let q_ffi = q_view.as_ffi();
            let view_query_ms = t_vq.map_or(0.0, |t| t.elapsed().as_secs_f64() * 1000.0);

            let t_qe = prof.then(Instant::now);
            let rc =
                unsafe { gecko_cuda_query_and_extend(self.handle.0, &q_ffi, &cuda_cfg, &mut raw) };
            let query_extend_ms = t_qe.map_or(0.0, |t| t.elapsed().as_secs_f64() * 1000.0);
            if rc != 0 {
                if !raw.frags.is_null() {
                    unsafe { gecko_cuda_free_frags(raw.frags) };
                }
                return Err(GeckoError::InvalidParameter(format!(
                    "CUDA query/extend execution failed (rc={rc})"
                )));
            }
            let t_fin = prof.then(Instant::now);
            let result = finalize_cuda_result(query, self.reference, config, &raw);
            if !raw.frags.is_null() {
                unsafe { gecko_cuda_free_frags(raw.frags) };
            }
            let finalize_ms = t_fin.map_or(0.0, |t| t.elapsed().as_secs_f64() * 1000.0);
            if prof {
                eprintln!(
                    "[gecko-cuda-profile] host_query view_query_ms={view_query_ms:.2} query_extend_ms={query_extend_ms:.2} finalize_ms={finalize_ms:.2}"
                );
            }
            Ok(result)
        }
    }

    pub fn compare_cuda(
        query: &Genome,
        reference: &Genome,
        config: &CompareConfig,
    ) -> Result<(CompareResult, CudaDeviceInfo)> {
        let (mut prepared, info) = PreparedRefInner::build(reference, config)?;
        let result = prepared.query(query, config)?;
        Ok((result, info))
    }

    /// 3-Phase Overlapped I/O + GPU pipeline (`W3-F12a`):
    /// - Phase 1: Uncontended `ref_path` FASTA load on all CPU cores || background `gecko_cuda_warmup()`.
    /// - Phase 2: Main thread builds ref index on hot GPU || background thread loads `query_path`
    ///   and launches async `seqX` H2D upload + GPU `fill_runs_bitmap_kernel` on a non-blocking stream.
    /// - Phase 3: Query & extend using pre-uploaded `seqX` (with transparent OOM fallback).
    pub fn compare_cuda_paths_overlapped(
        query_path: &std::path::Path,
        ref_path: &std::path::Path,
        config: &CompareConfig,
    ) -> Result<(Genome, Genome, CompareResult, CudaDeviceInfo)> {
        let prof = gecko_core::profile_enabled();
        let t_total = prof.then(Instant::now);
        validate_cuda_config(config)?;

        let t_probe = prof.then(Instant::now);
        let info = probe_cuda().map_err(GeckoError::InvalidParameter)?;
        let probe_ms = t_probe.map_or(0.0, |t| t.elapsed().as_secs_f64() * 1000.0);

        // Phase 1: Uncontended ref FASTA load || background CUDA primary context warmup (`cudaFree(0)`).
        let (reference, warmup_rc, load_ref_ms, warmup_bg_ms, warmup_wait_ms) =
            std::thread::scope(|s| -> Result<_> {
                let warmup_thread = s.spawn(|| {
                    let tw = prof.then(Instant::now);
                    let rc = unsafe { gecko_cuda_warmup() };
                    let ms = tw.map_or(0.0, |t| t.elapsed().as_secs_f64() * 1000.0);
                    (rc, ms)
                });
                let tr = prof.then(Instant::now);
                let ref_res = Genome::from_fasta_path(ref_path);
                let load_ref_ms = tr.map_or(0.0, |t| t.elapsed().as_secs_f64() * 1000.0);

                let tw_wait = prof.then(Instant::now);
                let (warmup_rc, warmup_bg_ms) =
                    warmup_thread.join().expect("CUDA warmup thread panicked");
                let warmup_wait_ms = tw_wait.map_or(0.0, |t| t.elapsed().as_secs_f64() * 1000.0);
                let reference = ref_res?;
                Ok((
                    reference,
                    warmup_rc,
                    load_ref_ms,
                    warmup_bg_ms,
                    warmup_wait_ms,
                ))
            })?;

        if warmup_rc != 0 {
            return Err(GeckoError::InvalidParameter(
                "CUDA primary context warmup failed".into(),
            ));
        }
        if reference.seq.len >= (1usize << 31) {
            return Err(GeckoError::InvalidParameter(
                "CUDA backend requires reference < 2.14 Gb".into(),
            ));
        }

        // Phase 2: Main thread builds ref index || bg thread loads query & uploads seqX async.
        let soft_mask = config.soft_mask;
        let (
            query,
            uploaded_q,
            handle,
            view_ref_ms,
            prepare_ref_ms,
            load_q_bg_ms,
            upload_q_bg_ms,
            q_wait_ms,
        ) = std::thread::scope(|s| -> Result<_> {
            let q_thread = s.spawn(
                || -> Result<(Genome, Option<UploadedQueryGuard>, f64, f64)> {
                    let tq = prof.then(Instant::now);
                    let q = Genome::from_fasta_path(query_path)?;
                    let load_q_ms = tq.map_or(0.0, |t| t.elapsed().as_secs_f64() * 1000.0);
                    if q.seq.len >= (1usize << 31) {
                        return Err(GeckoError::InvalidParameter(
                            "CUDA backend requires query < 2.14 Gb".into(),
                        ));
                    }
                    let tu = prof.then(Instant::now);
                    let q_view = HostPackedView::new(&q, soft_mask);
                    let q_ffi = q_view.as_ffi();
                    let ptr = unsafe { gecko_cuda_upload_query_async(&q_ffi) };
                    let upload_q_ms = tu.map_or(0.0, |t| t.elapsed().as_secs_f64() * 1000.0);
                    let guard = if ptr.is_null() {
                        None
                    } else {
                        Some(UploadedQueryGuard(ptr))
                    };
                    // SAFETY: `upload_h2d_staged_stream` guarantees on async streams that every byte
                    // of `q_view` (including <64 KiB arrays) is copied into pinned staging or synced
                    // before returning, so dropping `q_view` here is 100% DMA-safe.
                    Ok((q, guard, load_q_ms, upload_q_ms))
                },
            );

            let tvr = prof.then(Instant::now);
            let ref_view = HostPackedView::new(&reference, config.soft_mask);
            let ref_ffi = ref_view.as_ffi();
            let view_ref_ms = tvr.map_or(0.0, |t| t.elapsed().as_secs_f64() * 1000.0);

            let prep_cfg = GeckoCudaConfig::from_compare_config(config, 0);
            let tp = prof.then(Instant::now);
            let handle = RefGuard(unsafe { gecko_cuda_prepare_ref(&ref_ffi, &prep_cfg) });
            let prepare_ref_ms = tp.map_or(0.0, |t| t.elapsed().as_secs_f64() * 1000.0);

            let twq = prof.then(Instant::now);
            let (query, uploaded_q, load_q_bg_ms, upload_q_bg_ms) =
                q_thread.join().expect("FASTA query loader panicked")?;
            let q_wait_ms = twq.map_or(0.0, |t| t.elapsed().as_secs_f64() * 1000.0);

            if handle.0.is_null() {
                return Err(GeckoError::InvalidParameter(
                    "CUDA reference index preparation failed".into(),
                ));
            }
            Ok((
                query,
                uploaded_q,
                handle,
                view_ref_ms,
                prepare_ref_ms,
                load_q_bg_ms,
                upload_q_bg_ms,
                q_wait_ms,
            ))
        })?;

        let min_len = if config.fixed_length {
            config.min_length
        } else {
            let shorter = query.total_bases.min(reference.total_bases);
            ((shorter as f64) * (config.min_length as f64 / 100.0)) as u64
        };
        let query_cfg = GeckoCudaConfig::from_compare_config(config, min_len);

        // Phase 3: Query & extend using pre-uploaded seqX (or fallback if async upload OOM'd).
        let mut raw = GeckoCudaResult::default();
        let (view_query_ms, query_extend_ms, rc) = if let Some(uq) = uploaded_q.as_ref() {
            let t_qe = prof.then(Instant::now);
            let rc = unsafe {
                gecko_cuda_query_and_extend_uploaded(handle.0, uq.0, &query_cfg, &mut raw)
            };
            let qe_ms = t_qe.map_or(0.0, |t| t.elapsed().as_secs_f64() * 1000.0);
            (0.0, qe_ms, rc)
        } else {
            let t_vq = prof.then(Instant::now);
            let q_view = HostPackedView::new(&query, config.soft_mask);
            let q_ffi = q_view.as_ffi();
            let vq_ms = t_vq.map_or(0.0, |t| t.elapsed().as_secs_f64() * 1000.0);
            let t_qe = prof.then(Instant::now);
            let rc = unsafe { gecko_cuda_query_and_extend(handle.0, &q_ffi, &query_cfg, &mut raw) };
            let qe_ms = t_qe.map_or(0.0, |t| t.elapsed().as_secs_f64() * 1000.0);
            (vq_ms, qe_ms, rc)
        };
        drop(uploaded_q);

        if rc != 0 {
            if !raw.frags.is_null() {
                unsafe { gecko_cuda_free_frags(raw.frags) };
            }
            return Err(GeckoError::InvalidParameter(format!(
                "CUDA query/extend execution failed (rc={rc})"
            )));
        }

        let t_fin = prof.then(Instant::now);
        let result = finalize_cuda_result(&query, &reference, config, &raw);
        if !raw.frags.is_null() {
            unsafe { gecko_cuda_free_frags(raw.frags) };
        }
        let finalize_ms = t_fin.map_or(0.0, |t| t.elapsed().as_secs_f64() * 1000.0);

        if prof {
            let total_ms = t_total.map_or(0.0, |t| t.elapsed().as_secs_f64() * 1000.0);
            eprintln!(
                "[gecko-cuda-profile] host_overlapped total_ms={total_ms:.2} probe_ms={probe_ms:.2} \
                 load_ref_ms={load_ref_ms:.2} warmup_bg_ms={warmup_bg_ms:.2} warmup_wait_ms={warmup_wait_ms:.2} \
                 view_ref_ms={view_ref_ms:.2} prepare_ref_ms={prepare_ref_ms:.2} \
                 load_query_bg_ms={load_q_bg_ms:.2} upload_x_bg_ms={upload_q_bg_ms:.2} \
                 query_wait_ms={q_wait_ms:.2} view_query_ms={view_query_ms:.2} \
                 query_extend_ms={query_extend_ms:.2} finalize_ms={finalize_ms:.2}"
            );
        }

        Ok((query, reference, result, info))
    }
}

/// VRAM-resident prepared CUDA reference index + pooled tile scratch workspace.
///
/// Build once per reference genome and call [`CudaPreparedRef::query`] across arbitrarily many
/// query genomes (e.g., all-vs-all matrices) without rebuilding the GPU index or reallocating
/// tile scratch buffers.
pub struct CudaPreparedRef<'a> {
    #[cfg(has_cuda_kernels)]
    inner: ffi::PreparedRefInner<'a>,
    #[cfg(not(has_cuda_kernels))]
    _marker: std::marker::PhantomData<&'a Genome>,
}

impl<'a> CudaPreparedRef<'a> {
    /// Upload `reference` to VRAM and build its canonical GPU Hash-CSR index once.
    pub fn build(reference: &'a Genome, config: &CompareConfig) -> Result<(Self, CudaDeviceInfo)> {
        #[cfg(has_cuda_kernels)]
        {
            let (inner, info) = ffi::PreparedRefInner::build(reference, config)?;
            Ok((Self { inner }, info))
        }
        #[cfg(not(has_cuda_kernels))]
        {
            let _ = (reference, config);
            Err(gecko_core::error::GeckoError::InvalidParameter(
                "CUDA backend not compiled in (rebuild with `--features cuda`)".into(),
            ))
        }
    }

    /// Query `query` against the VRAM-resident reference index, reusing pooled scratch buffers.
    pub fn query(&mut self, query: &Genome, config: &CompareConfig) -> Result<CompareResult> {
        #[cfg(has_cuda_kernels)]
        {
            self.inner.query(query, config)
        }
        #[cfg(not(has_cuda_kernels))]
        {
            let _ = (query, config);
            Err(gecko_core::error::GeckoError::InvalidParameter(
                "CUDA backend not compiled in (rebuild with `--features cuda`)".into(),
            ))
        }
    }
}

/// Probe whether a native CUDA GPU is available in this build.
pub fn probe_cuda() -> std::result::Result<CudaDeviceInfo, String> {
    #[cfg(has_cuda_kernels)]
    {
        ffi::probe_cuda()
    }
    #[cfg(not(has_cuda_kernels))]
    {
        Err("CUDA backend not compiled in (rebuild with `--features cuda`)".to_string())
    }
}

/// Run a pairwise comparison on the native CUDA + CUB backend.
pub fn compare_cuda(
    query: &Genome,
    reference: &Genome,
    config: &CompareConfig,
) -> Result<(CompareResult, CudaDeviceInfo)> {
    #[cfg(has_cuda_kernels)]
    {
        ffi::compare_cuda(query, reference, config)
    }
    #[cfg(not(has_cuda_kernels))]
    {
        let _ = (query, reference, config);
        Err(gecko_core::error::GeckoError::InvalidParameter(
            "CUDA backend not compiled in (rebuild with `--features cuda`)".into(),
        ))
    }
}

/// Load FASTA files and run CUDA comparison with overlapped CPU query loading & GPU ref indexing.
pub fn compare_cuda_paths_overlapped(
    query_path: &Path,
    ref_path: &Path,
    config: &CompareConfig,
) -> Result<(Genome, Genome, CompareResult, CudaDeviceInfo)> {
    #[cfg(has_cuda_kernels)]
    {
        ffi::compare_cuda_paths_overlapped(query_path, ref_path, config)
    }
    #[cfg(not(has_cuda_kernels))]
    {
        let _ = (query_path, ref_path, config);
        Err(gecko_core::error::GeckoError::InvalidParameter(
            "CUDA backend not compiled in (rebuild with `--features cuda`)".into(),
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::{CudaScratchPool, GeckoRun, MAX_GPU_RUN_CHUNK_BP, pack_runs_for_gpu};
    use gecko_core::seq::PackedSeq;

    #[cfg(has_cuda_kernels)]
    unsafe extern "C" {
        fn gecko_cuda_fill_runs_bitmap_host_ref(
            runs: *const GeckoRun,
            n_runs: u32,
            len: u32,
            out_bm: *mut u32,
        );
        fn gecko_cuda_parse_profile_env_str(s: *const std::ffi::c_char) -> std::ffi::c_int;
        fn gecko_cuda_mix48_host(canon_kp: u64) -> u64;
        fn gecko_cuda_mix48_inv_host(h48: u64) -> u64;
        fn gecko_cuda_minimizer_params_host(
            wl: std::ffi::c_int,
            out_kp: *mut std::ffi::c_int,
            out_w: *mut std::ffi::c_int,
        ) -> std::ffi::c_int;
        fn gecko_cuda_two_hit_scan_host(
            sorted_hits: *const u64,
            n_hits: u32,
            tile_start: u64,
            tile_end: u64,
            len_y: u64,
            wl: std::ffi::c_int,
            compare_reverse: std::ffi::c_int,
            two_hit_window: u32,
            out_flags: *mut u8,
        );
    }

    fn build_packed_seq(raw: &[u8]) -> PackedSeq {
        gecko_core::seq::Genome::from_records(
            "test",
            raw.split(|&c| c == 0)
                .enumerate()
                .map(|(i, chunk)| (format!("s{i}"), chunk)),
        )
        .seq
    }

    fn parse_profile_flag(val: Option<&str>) -> bool {
        let Some(raw) = val else {
            return false;
        };
        let t = raw.trim();
        !t.is_empty()
            && t != "0"
            && !t.eq_ignore_ascii_case("false")
            && !t.eq_ignore_ascii_case("no")
            && !t.eq_ignore_ascii_case("off")
    }

    fn sim_fill_runs_bitmap(runs: &[GeckoRun], len: u32) -> Vec<u32> {
        let bm_words = (len as usize).div_ceil(32) + 2;
        let mut bm = vec![0u32; bm_words];
        for r in runs {
            let s = r.start.min(len);
            let e = r.end.min(len);
            if s >= e {
                continue;
            }
            let w_s = (s >> 5) as usize;
            let b_s = s & 31;
            let last = e - 1;
            let w_e = (last >> 5) as usize;
            let b_e = last & 31;
            if w_s == w_e {
                let mask = (0xFFFF_FFFFu32 << b_s) & (0xFFFF_FFFFu32 >> (31 - b_e));
                bm[w_s] |= mask;
            } else {
                bm[w_s] |= 0xFFFF_FFFFu32 << b_s;
                for w in &mut bm[(w_s + 1)..w_e] {
                    *w = 0xFFFF_FFFFu32;
                }
                bm[w_e] |= 0xFFFF_FFFFu32 >> (31 - b_e);
            }
        }
        bm
    }

    fn sim_materialize_pair_bitmaps(
        seq: &PackedSeq,
        soft_mask: bool,
    ) -> (Option<Vec<u32>>, Option<Vec<u32>>, bool) {
        let len = seq.len as u32;
        let special_runs = pack_runs_for_gpu(&seq.special_runs, seq.len as u64);
        let soft_runs = if soft_mask && !seq.soft_runs.is_empty() {
            pack_runs_for_gpu(&seq.soft_runs, seq.len as u64)
        } else {
            Vec::new()
        };
        let has_special = len > 0 && !special_runs.is_empty();
        let has_soft = len > 0 && soft_mask && !soft_runs.is_empty();

        let fill = |runs: &[GeckoRun], dst: &mut [u32]| {
            let sim_tmp = sim_fill_runs_bitmap(runs, len);
            #[cfg(has_cuda_kernels)]
            {
                let mut ffi_tmp = vec![0u32; dst.len()];
                unsafe {
                    gecko_cuda_fill_runs_bitmap_host_ref(
                        runs.as_ptr(),
                        runs.len() as u32,
                        len,
                        ffi_tmp.as_mut_ptr(),
                    );
                }
                assert_eq!(
                    ffi_tmp, sim_tmp,
                    "CUDA C++ host helper vs Rust sim mismatch"
                );
            }
            for (d, &t) in dst.iter_mut().zip(sim_tmp.iter()) {
                *d |= t;
            }
        };

        let bm_words = seq.len.div_ceil(32) + 2;
        match (has_special, has_soft) {
            (false, false) => (None, None, false),
            (true, false) => {
                let mut sp = vec![0u32; bm_words];
                fill(&special_runs, &mut sp);
                let inv = sp.clone();
                (Some(sp), Some(inv), true) // aliased
            }
            (false, true) => {
                let mut inv = vec![0u32; bm_words];
                fill(&soft_runs, &mut inv);
                (None, Some(inv), false)
            }
            (true, true) => {
                let mut sp = vec![0u32; bm_words];
                fill(&special_runs, &mut sp);
                let mut inv = vec![0u32; bm_words];
                fill(&special_runs, &mut inv);
                fill(&soft_runs, &mut inv);
                (Some(sp), Some(inv), false)
            }
        }
    }

    // ---- Test Group 1: Run-Bitmap Fill Parity (W3-F12b) ----

    #[test]
    fn test_pack_runs_for_gpu_alignment_and_clamping() {
        let runs = vec![
            (0u64, 0u64),             // empty
            (50u64, 40u64),           // inverted
            (5u64, 200_005u64),       // giant run spanning many chunks
            (250_000u64, 300_000u64), // clamped to seq_len = 260_000
        ];
        let packed = pack_runs_for_gpu(&runs, 260_000);
        assert!(packed.len() > 10);
        // Verify strict disjointness, ascending order, <= MAX_GPU_RUN_CHUNK_BP, and 32-bp aligned splits.
        for win in packed.windows(2) {
            assert!(win[0].start < win[0].end);
            assert!(win[0].end <= win[1].start);
            assert!((win[0].end - win[0].start) as u64 <= MAX_GPU_RUN_CHUNK_BP);
            if win[0].end == win[1].start {
                assert_eq!(win[0].end & 31, 0, "split boundary must be 32-bp aligned");
            }
        }
        // Check exact coverage of [5, 200_005) and [250_000, 260_000)
        let mut covered = 0u64;
        for r in &packed {
            assert!(r.end <= 260_000);
            covered += (r.end - r.start) as u64;
        }
        assert_eq!(covered, 200_000 + 10_000);
    }

    #[test]
    fn test_gpu_run_bitmap_fill_parity_vs_packed_seq() {
        let mut cases: Vec<Vec<u8>> = vec![
            b"".to_vec(),
            b"A".to_vec(),
            b"N".to_vec(),
            b"a".to_vec(),
            b"ACGTACGTACGTACGTACGTACGTACGTACGT".to_vec(), // 32 clean
            b"ACGTACGTACGTACGTACGTACGTACGTACG".to_vec(),  // 31 clean
            b"ACGTACGTACGTACGTACGTACGTACGTACGTA".to_vec(), // 33 clean
            // Single-base runs at 0, 1, 30, 31, 32, 33, and end
            {
                let mut s = vec![b'A'; 65];
                for &p in &[0usize, 1, 30, 31, 32, 33, 64] {
                    s[p] = b'N';
                }
                s
            },
            // Multiple short runs sharing same 32-bp word + IUPAC + SEQ_SEP + soft-masked
            {
                let mut s = Vec::new();
                s.extend_from_slice(b"acgNNNtgcaRYSWKM");
                s.push(0u8); // SEQ_SEP
                s.extend_from_slice(b"ACGTacgtNNN");
                s.extend_from_slice(&[b't'; 75]);
                s.extend_from_slice(b"GCTAN");
                s.push(0u8);
                s.extend_from_slice(b"acgT");
                s
            },
            // Giant 100,000-bp N-run + flanking soft-masked regions
            {
                let mut s = vec![b'A'; 500];
                s.extend_from_slice(&[b'c'; 37]);
                s.extend_from_slice(&vec![b'N'; 100_000]);
                s.extend_from_slice(&[b'g'; 41]);
                s.extend_from_slice(&vec![b'T'; 513]);
                s
            },
        ];
        // Also test exact word-aligned lengths 64 & 96 with straddling runs [28..36) & [5..90)
        {
            let mut s = vec![b'C'; 64];
            for b in &mut s[28..36] {
                *b = b'N';
            }
            cases.push(s);
        }
        {
            let mut s = vec![b'G'; 96];
            for b in &mut s[5..90] {
                *b = b'a';
            }
            cases.push(s);
        }

        for raw in &cases {
            let seq = build_packed_seq(raw);
            for soft_mask in [false, true] {
                let ref_sp = seq.special_bitmap();
                let ref_inv = seq.invalid_seed_bitmap(soft_mask);
                let (gpu_sp, gpu_inv, aliased) = sim_materialize_pair_bitmaps(&seq, soft_mask);

                assert_eq!(
                    gpu_sp, ref_sp,
                    "special_bm mismatch (len={}, soft_mask={soft_mask})",
                    seq.len
                );
                assert_eq!(
                    gpu_inv, ref_inv,
                    "invalid_bm mismatch (len={}, soft_mask={soft_mask})",
                    seq.len
                );
                if ref_sp.is_some() {
                    assert_eq!(
                        aliased,
                        seq.invalid_bitmap_aliases_special(soft_mask),
                        "aliasing mismatch"
                    );
                }
                // Trailing 2 sentinel words MUST be strictly 0u32
                for bm in [gpu_sp.as_ref(), gpu_inv.as_ref()].into_iter().flatten() {
                    assert!(bm.len() >= 2);
                    assert_eq!(bm[bm.len() - 1], 0u32);
                    assert_eq!(bm[bm.len() - 2], 0u32);
                }
            }
        }
    }

    // ---- Shared bit-exact helpers mirroring cuda_common.cuh ----

    fn fmix64(mut x: u64) -> u64 {
        x ^= x >> 33;
        x = x.wrapping_mul(0xff51afd7ed558ccd);
        x ^= x >> 33;
        x = x.wrapping_mul(0xc4ceb9fe1a85ec53);
        x ^= x >> 33;
        x
    }

    fn fmix64_inv(mut x: u64) -> u64 {
        x ^= x >> 33;
        x = x.wrapping_mul(0x9cb4b2f8129337db);
        x ^= x >> 33;
        x = x.wrapping_mul(0x4f74430c22a54005);
        x ^= x >> 33;
        x
    }

    fn reverse_32_2bit(x: u64) -> u64 {
        let mut y = x.swap_bytes();
        y = ((y & 0x0F0F_0F0F_0F0F_0F0F) << 4) | ((y >> 4) & 0x0F0F_0F0F_0F0F_0F0F);
        y = ((y & 0x3333_3333_3333_3333) << 2) | ((y >> 2) & 0x3333_3333_3333_3333);
        y
    }

    fn is_palindromic_canon(canon: u64, wl: usize, eff_mask: u64) -> bool {
        let rev = ((!reverse_32_2bit(canon)) >> (64 - 2 * wl)) & eff_mask;
        rev == canon
    }

    fn is_bit_set(bm: Option<&[u32]>, pos: usize) -> bool {
        match bm {
            None => false,
            Some(w) => ((w[pos >> 5] >> (pos & 31)) & 1) != 0,
        }
    }

    fn set_bit(bm: &mut [u32], pos: usize) {
        bm[pos >> 5] |= 1u32 << (pos & 31);
    }

    fn load_32_2bit(bits: &[u64], pos: usize) -> u64 {
        let w = pos >> 5;
        let shift = (pos & 31) << 1;
        let lo = bits[w];
        if shift == 0 {
            lo
        } else {
            let hi = bits[w + 1];
            (lo >> shift) | (hi << (64 - shift))
        }
    }

    fn base_2bit(bits: &[u64], pos: usize) -> u32 {
        ((bits[pos >> 5] >> ((pos & 31) << 1)) & 3) as u32
    }

    fn make_revcomp_bytes(s: &[u8]) -> Vec<u8> {
        s.iter()
            .rev()
            .map(|&c| match c {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                b'a' => b't',
                b'c' => b'g',
                b'g' => b'c',
                b't' => b'a',
                other => other,
            })
            .collect()
    }

    // ---- Test Group 2: Inverted d_ref_non_unique_bm + Fused d_dropped_bitmap Parity (W3-F13a) ----

    #[test]
    fn test_inverted_non_unique_and_fused_dropped_bm_parity() {
        // Construct deterministic reference covering all 6 critical k-mer classes.
        let mut ref_bytes = Vec::new();
        let alphabet = *b"ACGT";
        // Deterministic pseudo-random unique background (512 bp)
        let mut rng = 0x1234_5678_9abc_def0u64;
        for _ in 0..512 {
            rng = rng
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ref_bytes.push(alphabet[((rng >> 32) & 3) as usize]);
        }
        // Inject a 32-bp palindrome (16 bp + revcomp(16 bp)) surrounded by N separators
        ref_bytes.push(b'N');
        let half16 = b"ACGTACGTTAGCGCTA";
        let pal32: Vec<u8> = half16
            .iter()
            .copied()
            .chain(make_revcomp_bytes(half16))
            .collect();
        ref_bytes.extend_from_slice(&pal32);
        ref_bytes.push(b'N');
        // Inject a 20-bp palindrome and 12-bp palindrome
        let half10 = b"ACGTCTAGCA";
        let pal20: Vec<u8> = half10
            .iter()
            .copied()
            .chain(make_revcomp_bytes(half10))
            .collect();
        ref_bytes.extend_from_slice(&pal20);
        ref_bytes.push(b'N');
        let half6 = b"ACGTCG";
        let pal12: Vec<u8> = half6
            .iter()
            .copied()
            .chain(make_revcomp_bytes(half6))
            .collect();
        ref_bytes.extend_from_slice(&pal12);
        ref_bytes.push(b'N');

        // Low-frequency repeat (2 fwd + 1 rev <= fl=3)
        let rep_low = b"GGTAACCGTTAACCGGTTAACCGGTTAACCGG";
        let rep_low_rc = make_revcomp_bytes(rep_low);
        for _ in 0..2 {
            ref_bytes.extend_from_slice(rep_low);
            ref_bytes.push(b'N');
        }
        ref_bytes.extend_from_slice(&rep_low_rc);
        ref_bytes.push(b'N');

        // One-strand over-frequent repeat (fl=2: 2 fwd kept, 4 rev dropped)
        let rep_one = b"CCATTGGAACCATTGGAACCATTGGAACCATT";
        let rep_one_rc = make_revcomp_bytes(rep_one);
        for _ in 0..2 {
            ref_bytes.extend_from_slice(rep_one);
            ref_bytes.push(b'N');
        }
        for _ in 0..4 {
            ref_bytes.extend_from_slice(&rep_one_rc);
            ref_bytes.push(b'N');
        }

        // Both-strands over-frequent repeat (fl=2: 4 fwd + 4 rev => tests early-break trap!)
        let rep_both = b"TGCATGCAGTACGTACTGCAACTGGTCAGTCA";
        let rep_both_rc = make_revcomp_bytes(rep_both);
        for _ in 0..4 {
            ref_bytes.extend_from_slice(rep_both);
            ref_bytes.push(b'N');
            ref_bytes.extend_from_slice(&rep_both_rc);
            ref_bytes.push(b'N');
        }

        let seq = build_packed_seq(&ref_bytes);
        let inv_bm = seq.invalid_seed_bitmap(true);
        let bm_words = seq.len.div_ceil(32) + 2;

        for (wl, is_direct) in [(12usize, true), (16, false), (20, false), (32, false)] {
            let span_mask = if wl == 32 {
                !0u64
            } else {
                (1u64 << (2 * wl)) - 1
            };
            let rev_rshift = 64 - 2 * wl;

            // Extract all valid k-mers
            let mut pairs: Vec<(u64, u32)> = Vec::new();
            let mut valid_positions = Vec::new();
            for pos in 0..=seq.len.saturating_sub(wl) {
                let mut ok = true;
                for d in 0..wl {
                    if is_bit_set(inv_bm.as_deref(), pos + d) {
                        ok = false;
                        break;
                    }
                }
                if !ok {
                    continue;
                }
                valid_positions.push(pos);
                let w_lsb = load_32_2bit(&seq.bits, pos);
                let kf = reverse_32_2bit(w_lsb) >> rev_rshift;
                let kr = (!w_lsb) & span_mask;
                let strand = if kf <= kr { 0u32 } else { 1u32 };
                let canon = if strand == 0 { kf } else { kr };
                let key = if is_direct { canon } else { fmix64(canon) };
                pairs.push((key, ((pos as u32) << 1) | strand));
            }
            pairs.sort_by_key(|&(k, _)| k);
            let n = pairs.len();

            for &freq_limit in &[2u32, 3u32, n as u32 + 10] {
                let no_freq_drops = (freq_limit as usize) >= n;

                // OLD reference model (positive ref_unique_bm + separate dropped_bm)
                let mut old_keep = vec![0u8; n];
                let mut old_unique_bm = vec![0u32; bm_words];
                let mut old_dropped_bm = vec![0u32; bm_words];
                if no_freq_drops {
                    old_keep.fill(1);
                    for i in 0..n {
                        let k = pairs[i].0;
                        let unique =
                            (i == 0 || pairs[i - 1].0 != k) && (i + 1 >= n || pairs[i + 1].0 != k);
                        if unique {
                            let v = pairs[i].1;
                            let is_pal = if (v & 1) == 0 {
                                let canon = if is_direct { k } else { fmix64_inv(k) };
                                is_palindromic_canon(canon, wl, span_mask)
                            } else {
                                false
                            };
                            if !is_pal {
                                set_bit(&mut old_unique_bm, (v >> 1) as usize);
                            }
                        }
                    }
                } else {
                    let mut i = 0;
                    while i < n {
                        let k = pairs[i].0;
                        let mut j = i + 1;
                        while j < n && pairs[j].0 == k {
                            j += 1;
                        }
                        let run_len = j - i;
                        if run_len == 1 {
                            old_keep[i] = 1;
                            let v = pairs[i].1;
                            let is_pal = if (v & 1) == 0 {
                                let canon = if is_direct { k } else { fmix64_inv(k) };
                                is_palindromic_canon(canon, wl, span_mask)
                            } else {
                                false
                            };
                            if !is_pal {
                                set_bit(&mut old_unique_bm, (v >> 1) as usize);
                            }
                        } else if run_len <= freq_limit as usize {
                            for f in &mut old_keep[i..j] {
                                *f = 1;
                            }
                        } else {
                            let mut c0 = 0usize;
                            let mut c1 = 0usize;
                            for &(_, v) in &pairs[i..j] {
                                if (v & 1) != 0 {
                                    c1 += 1;
                                } else {
                                    c0 += 1;
                                }
                            }
                            let k0 = u8::from(c0 <= freq_limit as usize);
                            let k1 = u8::from(c1 <= freq_limit as usize);
                            for m in i..j {
                                old_keep[m] = if (pairs[m].1 & 1) != 0 { k1 } else { k0 };
                            }
                        }
                        i = j;
                    }
                    for m in 0..n {
                        if old_keep[m] == 0 {
                            set_bit(&mut old_dropped_bm, (pairs[m].1 >> 1) as usize);
                        }
                    }
                }

                // NEW W3-F13a model (inverted ref_non_unique_bm + fused dropped_bm with break not return)
                let mut new_keep = vec![0u8; n];
                let mut new_non_unique_bm = vec![0u32; bm_words];
                let mut new_dropped_bm = vec![0u32; bm_words];
                let mut any_dropped = false;
                if no_freq_drops {
                    new_keep.fill(1);
                    for i in 0..n {
                        let k = pairs[i].0;
                        let mut non_unique =
                            (i > 0 && pairs[i - 1].0 == k) || (i + 1 < n && pairs[i + 1].0 == k);
                        if !non_unique && (wl & 1) == 0 {
                            let canon = if is_direct { k } else { fmix64_inv(k) };
                            non_unique = is_palindromic_canon(canon, wl, span_mask);
                        }
                        if non_unique {
                            set_bit(&mut new_non_unique_bm, (pairs[i].1 >> 1) as usize);
                        }
                    }
                } else {
                    for i in 0..n {
                        let k = pairs[i].0;
                        let same_prev = i > 0 && pairs[i - 1].0 == k;
                        if same_prev {
                            set_bit(&mut new_non_unique_bm, (pairs[i].1 >> 1) as usize);
                            continue;
                        }
                        let same_next = i + 1 < n && pairs[i + 1].0 == k;
                        if !same_next {
                            new_keep[i] = 1;
                            if (wl & 1) == 0 {
                                let canon = if is_direct { k } else { fmix64_inv(k) };
                                if is_palindromic_canon(canon, wl, span_mask) {
                                    set_bit(&mut new_non_unique_bm, (pairs[i].1 >> 1) as usize);
                                }
                            }
                            continue;
                        }
                        set_bit(&mut new_non_unique_bm, (pairs[i].1 >> 1) as usize);
                        let mut run_end = i + 1;
                        while run_end < n && pairs[run_end].0 == k {
                            run_end += 1;
                        }
                        let run_len = run_end - i;
                        if run_len <= freq_limit as usize {
                            for f in &mut new_keep[i..run_end] {
                                *f = 1;
                            }
                            continue;
                        }
                        let mut c0 = 0usize;
                        let mut c1 = 0usize;
                        for &(_, v) in &pairs[i..run_end] {
                            if (v & 1) != 0 {
                                c1 += 1;
                            } else {
                                c0 += 1;
                            }
                            if c0 > freq_limit as usize && c1 > freq_limit as usize {
                                break;
                            }
                        }
                        let k0 = u8::from(c0 <= freq_limit as usize);
                        let k1 = u8::from(c1 <= freq_limit as usize);
                        if k0 == 0 || k1 == 0 {
                            any_dropped = true;
                            for m in i..run_end {
                                let v = pairs[m].1;
                                let keep = if (v & 1) != 0 { k1 } else { k0 };
                                if keep != 0 {
                                    new_keep[m] = 1;
                                } else {
                                    set_bit(&mut new_dropped_bm, (v >> 1) as usize);
                                }
                            }
                        } else {
                            new_keep[i..run_end].fill(1);
                        }
                    }
                }

                assert_eq!(
                    new_keep, old_keep,
                    "keep flags mismatch (wl={wl}, fl={freq_limit})"
                );
                assert_eq!(
                    new_dropped_bm, old_dropped_bm,
                    "fused dropped_bm mismatch (wl={wl}, fl={freq_limit})"
                );
                let has_any_drop_ref = old_keep.contains(&0);
                assert_eq!(any_dropped, has_any_drop_ref);

                // De Morgan parity on EVERY valid ref k-mer position!
                for &pos in &valid_positions {
                    assert_eq!(
                        !is_bit_set(Some(&new_non_unique_bm), pos),
                        is_bit_set(Some(&old_unique_bm), pos),
                        "De Morgan parity failed at pos={pos}, wl={wl}, fl={freq_limit}"
                    );
                }

                // Multi-slab (2 slabs: [0..128) and [128..256)) global dropped_bm parity for Hash-CSR
                if !is_direct && !no_freq_drops {
                    let mut slab_dropped_bm = vec![0u32; bm_words];
                    for (ps, pe) in [(0u32, 128u32), (128u32, 256u32)] {
                        let sub: Vec<(u64, u32)> = pairs
                            .iter()
                            .copied()
                            .filter(|&(k, _)| {
                                let p = (k >> 56) as u32;
                                p >= ps && p < pe
                            })
                            .collect();
                        let mut i = 0;
                        while i < sub.len() {
                            let k = sub[i].0;
                            let mut j = i + 1;
                            while j < sub.len() && sub[j].0 == k {
                                j += 1;
                            }
                            if j - i > freq_limit as usize {
                                let mut c0 = 0usize;
                                let mut c1 = 0usize;
                                for &(_, v) in &sub[i..j] {
                                    if (v & 1) != 0 {
                                        c1 += 1;
                                    } else {
                                        c0 += 1;
                                    }
                                    if c0 > freq_limit as usize && c1 > freq_limit as usize {
                                        break;
                                    }
                                }
                                let k0 = c0 <= freq_limit as usize;
                                let k1 = c1 <= freq_limit as usize;
                                for &(_, v) in &sub[i..j] {
                                    let keep = if (v & 1) != 0 { k1 } else { k0 };
                                    if !keep {
                                        set_bit(&mut slab_dropped_bm, (v >> 1) as usize);
                                    }
                                }
                            }
                            i = j;
                        }
                    }
                    assert_eq!(slab_dropped_bm, old_dropped_bm);
                }
            }
        }
    }

    // ---- Test Group 3: Hybrid Pre-Bypass + CTA-Privatized L2 Scatter Parity (W3-F13b) ----

    #[test]
    fn test_hybrid_prebypass_and_privatized_l2_scatter_parity() {
        let wl = 16usize;
        let mask = (1u64 << (2 * wl)) - 1;
        let rev_rshift = 64 - 2 * wl;
        let freq_limit = 2usize;
        let alphabet = *b"ACGT";

        // Build deterministic reference (4,096 bp) with unique segments, palindrome, repeats, N, soft-mask
        let mut rng = 0xCAFE_BABE_DEAD_BEEFu64;
        let mut next_base = || {
            rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
            alphabet[((rng >> 33) & 3) as usize]
        };

        let mut ref_bytes: Vec<u8> = (0..1500).map(|_| next_base()).collect();
        let unique_fwd_block: Vec<u8> = (0..160).map(|_| next_base()).collect();
        let unique_rev_block: Vec<u8> = (0..160).map(|_| next_base()).collect();
        let pal_half = b"ACGTCTAG";
        let pal16: Vec<u8> = pal_half
            .iter()
            .copied()
            .chain(make_revcomp_bytes(pal_half))
            .collect();
        let rep_keep: Vec<u8> = (0..24).map(|_| next_base()).collect();
        let rep_drop: Vec<u8> = (0..24).map(|_| next_base()).collect();

        ref_bytes.extend_from_slice(&unique_fwd_block);
        ref_bytes.extend_from_slice(&pal16);
        ref_bytes.extend_from_slice(&unique_fwd_block[..80]); // makes prefix non-unique!
        ref_bytes.extend_from_slice(&unique_rev_block);
        // Low-freq repeat (2x <= fl=2)
        for _ in 0..2 {
            ref_bytes.extend_from_slice(&rep_keep);
            ref_bytes.push(b'G');
        }
        // Over-freq repeat (4x > fl=2 -> dropped_bitmap == 1)
        for _ in 0..4 {
            ref_bytes.extend_from_slice(&rep_drop);
            ref_bytes.push(b'T');
        }
        ref_bytes.extend_from_slice(b"NNNNacgtacgtNNNN");
        while ref_bytes.len() < 4096 {
            ref_bytes.push(next_base());
        }

        // Build query (non-multiple-of-256 length = 1,373 bp) containing:
        //  - long unique fwd match (triggers Case 2 bypass)
        //  - SNP inside unique match (resets act_valid)
        //  - palindrome on anchor and inside diagonal (Case 3 / Subcase 4c)
        //  - reverse-complement long match
        //  - repeated kept & dropped k-mers
        //  - N / soft-masked bases (Case 1)
        let mut q_bytes = Vec::new();
        q_bytes.extend_from_slice(&unique_fwd_block);
        q_bytes.extend_from_slice(&pal16);
        // Copy with 1 SNP in the middle
        let mut snp_block = unique_fwd_block.clone();
        snp_block[64] = if snp_block[64] == b'A' { b'C' } else { b'A' };
        q_bytes.extend_from_slice(&snp_block);
        q_bytes.extend_from_slice(&make_revcomp_bytes(&unique_rev_block));
        q_bytes.extend_from_slice(&rep_keep);
        q_bytes.extend_from_slice(&rep_drop);
        q_bytes.extend_from_slice(b"NNNacgtNNN");
        q_bytes.extend_from_slice(&pal16);
        while q_bytes.len() < 1373 {
            q_bytes.push(next_base());
        }

        // Also test a 100%-bypassed tile query (pure unique forward exact match slice!)
        let pure_unique_query = unique_fwd_block[80..160].to_vec();

        let seq_y = build_packed_seq(&ref_bytes);
        let inv_y = seq_y.invalid_seed_bitmap(true);
        let bm_words_y = seq_y.len.div_ceil(32) + 2;

        // Build Hash-CSR index + inverted ref_non_unique_bm + dropped_bitmap
        let mut kvs: Vec<(u64, u32)> = Vec::new();
        for pos in 0..=seq_y.len - wl {
            if (0..wl).any(|d| is_bit_set(inv_y.as_deref(), pos + d)) {
                continue;
            }
            let w_lsb = load_32_2bit(&seq_y.bits, pos);
            let kf = reverse_32_2bit(w_lsb) >> rev_rshift;
            let kr = (!w_lsb) & mask;
            let strand = if kf <= kr { 0u32 } else { 1u32 };
            let canon = if strand == 0 { kf } else { kr };
            kvs.push((fmix64(canon), ((pos as u32) << 1) | strand));
        }
        kvs.sort_by_key(|&(h, _)| h);

        let mut ref_non_unique_bm = vec![0u32; bm_words_y];
        let mut dropped_bm = vec![0u32; bm_words_y];
        let mut kept: Vec<(u64, u32)> = Vec::new();
        {
            let mut i = 0;
            while i < kvs.len() {
                let h = kvs[i].0;
                let mut j = i + 1;
                while j < kvs.len() && kvs[j].0 == h {
                    j += 1;
                }
                if j - i == 1 {
                    let v = kvs[i].1;
                    if is_palindromic_canon(fmix64_inv(h), wl, mask) {
                        set_bit(&mut ref_non_unique_bm, (v >> 1) as usize);
                    }
                    kept.push(kvs[i]);
                } else {
                    for &(_, v) in &kvs[i..j] {
                        set_bit(&mut ref_non_unique_bm, (v >> 1) as usize);
                    }
                    let c0 = kvs[i..j].iter().filter(|&&(_, v)| (v & 1) == 0).count();
                    let c1 = (j - i) - c0;
                    for &(kh, v) in &kvs[i..j] {
                        let ok = if (v & 1) != 0 {
                            c1 <= freq_limit
                        } else {
                            c0 <= freq_limit
                        };
                        if ok {
                            kept.push((kh, v));
                        } else {
                            set_bit(&mut dropped_bm, (v >> 1) as usize);
                        }
                    }
                }
                i = j;
            }
        }

        let bucket_bits = 10usize;
        let num_buckets = 1usize << bucket_bits;
        let shift_bits = 64 - bucket_bits;
        let tag_shift = shift_bits - 32;
        let mut bucket_off = vec![0usize; num_buckets + 1];
        let mut idx_tags = Vec::with_capacity(kept.len());
        let mut idx_vals = Vec::with_capacity(kept.len());
        for &(h, v) in &kept {
            let b = (h >> shift_bits) as usize;
            bucket_off[b + 1] += 1;
            idx_tags.push((h >> tag_shift) as u32);
            idx_vals.push(v);
        }
        for b in 0..num_buckets {
            bucket_off[b + 1] += bucket_off[b];
        }

        let check_lr_fwd =
            |seq_x: &PackedSeq, inv_x: Option<&[u32]>, pos_x: usize, pos_y: usize| -> bool {
                let prev_x_valid = pos_x > 0 && !is_bit_set(inv_x, pos_x - 1);
                let left_ext = prev_x_valid
                    && pos_y > 0
                    && !is_bit_set(inv_y.as_deref(), pos_y - 1)
                    && base_2bit(&seq_y.bits, pos_y - 1) == base_2bit(&seq_x.bits, pos_x - 1)
                    && !is_bit_set(Some(&dropped_bm), pos_y - 1);
                if !left_ext {
                    return true;
                }
                let next_x_valid = pos_x + wl < seq_x.len && !is_bit_set(inv_x, pos_x + wl);
                let right_ext = next_x_valid
                    && pos_y + wl < seq_y.len
                    && !is_bit_set(inv_y.as_deref(), pos_y + wl)
                    && base_2bit(&seq_y.bits, pos_y + wl) == base_2bit(&seq_x.bits, pos_x + wl)
                    && !is_bit_set(Some(&dropped_bm), pos_y + 1);
                !right_ext
            };

        let check_lr_rev =
            |seq_x: &PackedSeq, inv_x: Option<&[u32]>, pos_x: usize, pos_y: usize| -> bool {
                let pos_y_rc = seq_y.len - pos_y - wl;
                let prev_x_valid = pos_x > 0 && !is_bit_set(inv_x, pos_x - 1);
                let left_ext = prev_x_valid
                    && pos_y_rc > 0
                    && !is_bit_set(inv_y.as_deref(), pos_y + wl)
                    && (base_2bit(&seq_y.bits, pos_y + wl) ^ 3)
                        == base_2bit(&seq_x.bits, pos_x - 1)
                    && !is_bit_set(Some(&dropped_bm), pos_y + 1);
                if !left_ext {
                    return true;
                }
                let next_x_valid = pos_x + wl < seq_x.len && !is_bit_set(inv_x, pos_x + wl);
                let right_ext = next_x_valid
                    && pos_y > 0
                    && !is_bit_set(inv_y.as_deref(), pos_y - 1)
                    && (base_2bit(&seq_y.bits, pos_y - 1) ^ 3)
                        == base_2bit(&seq_x.bits, pos_x + wl)
                    && !is_bit_set(Some(&dropped_bm), pos_y - 1);
                !right_ext
            };

        for (q_label, raw_q) in [
            ("mixed", &q_bytes),
            ("pure_unique_100pct", &pure_unique_query),
        ] {
            let seq_x = build_packed_seq(raw_q);
            let inv_x = seq_x.invalid_seed_bitmap(true);
            let num_pos = seq_x.len - (wl - 1);
            let bp_words = (num_pos + 31) >> 5;

            for compare_reverse in [true, false] {
                // 1. Monolithic Oracle (no filter_sorted_hits min_gap heuristic!)
                let mut oracle_hits: Vec<(u8, usize, usize)> = Vec::new();
                for pos_x in 0..num_pos {
                    if (0..wl).any(|d| is_bit_set(inv_x.as_deref(), pos_x + d)) {
                        continue;
                    }
                    let w_lsb = load_32_2bit(&seq_x.bits, pos_x);
                    let kf = reverse_32_2bit(w_lsb) >> rev_rshift;
                    let kr = (!w_lsb) & mask;
                    let x_strand = u32::from(kf > kr);
                    let canon = if x_strand == 0 { kf } else { kr };
                    let is_pal = kf == kr;
                    let h = fmix64(canon);
                    let bkt = (h >> shift_bits) as usize;
                    let tag = (h >> tag_shift) as u32;
                    let x_seed_lsb = w_lsb & mask;
                    for j in bucket_off[bkt]..bucket_off[bkt + 1] {
                        if idx_tags[j] != tag {
                            continue;
                        }
                        let ps = idx_vals[j];
                        let pos_y = (ps >> 1) as usize;
                        let y_strand = ps & 1;
                        let y_word = load_32_2bit(&seq_y.bits, pos_y);
                        if (is_pal || x_strand == y_strand)
                            && ((y_word ^ x_seed_lsb) & mask) == 0
                            && check_lr_fwd(&seq_x, inv_x.as_deref(), pos_x, pos_y)
                        {
                            oracle_hits.push((0, pos_x, pos_y));
                        }
                        if compare_reverse && (is_pal || x_strand != y_strand) {
                            let y_rc_lsb = (!reverse_32_2bit(y_word)) >> rev_rshift;
                            if ((y_rc_lsb ^ x_seed_lsb) & mask) == 0
                                && check_lr_rev(&seq_x, inv_x.as_deref(), pos_x, pos_y)
                            {
                                let pos_y_rc = seq_y.len - pos_y - wl;
                                oracle_hits.push((1, pos_x, pos_y_rc));
                            }
                        }
                    }
                }
                oracle_hits.sort_unstable();

                // 2. Hybrid Pre-Bypass + CTA-Privatized L2 Scatter + L2 Probe across P in [2, 8, 64]
                for p_parts in [2usize, 8, 64] {
                    let p_bits = p_parts.trailing_zeros() as usize;
                    let part_shift = bucket_bits - p_bits;

                    let mut hybrid_hits: Vec<(u8, usize, usize)> = Vec::new();
                    let mut d_q_bypass_bm = vec![0u32; bp_words];
                    let mut total_valid_q = 0usize;
                    let mut total_bypassed_q = 0usize;

                    // Simulate l2_pre_bypass_kernel (1 thread = 4 consecutive 32-bp words)
                    let n_threads_bp = (bp_words + 3) >> 2;
                    for gtid in 0..n_threads_bp {
                        let mut act_valid = false;
                        let mut act_strand = 0u32;
                        let mut act_pos_y_prev = 0usize;
                        for wave in 0..4usize {
                            let w_out = gtid * 4 + wave;
                            if w_out >= bp_words {
                                break;
                            }
                            let w_rel0 = w_out << 5;
                            let w_cnt = (num_pos - w_rel0).min(32);
                            let mut bp_word = if w_cnt < 32 { !0u32 << w_cnt } else { 0u32 };

                            for s in 0..w_cnt {
                                let pos_x = w_rel0 + s;
                                if (0..wl).any(|d| is_bit_set(inv_x.as_deref(), pos_x + d)) {
                                    act_valid = false;
                                    bp_word |= 1u32 << s;
                                    continue;
                                }
                                total_valid_q += 1;
                                if act_valid {
                                    let b = base_2bit(&seq_x.bits, pos_x + wl - 1);
                                    let (ext_ok, cand_y) = if act_strand == 0 {
                                        let cy = act_pos_y_prev + 1;
                                        let new_y = cy + wl - 1;
                                        let ok = new_y < seq_y.len
                                            && !is_bit_set(inv_y.as_deref(), new_y)
                                            && base_2bit(&seq_y.bits, new_y) == b;
                                        (ok, cy)
                                    } else if act_pos_y_prev > 0 {
                                        let cy = act_pos_y_prev - 1;
                                        let ok = !is_bit_set(inv_y.as_deref(), cy)
                                            && (base_2bit(&seq_y.bits, cy) ^ 3) == b;
                                        (ok, cy)
                                    } else {
                                        (false, 0)
                                    };
                                    if ext_ok {
                                        act_pos_y_prev = cand_y;
                                        if !is_bit_set(Some(&ref_non_unique_bm), cand_y) {
                                            if act_strand == 0 {
                                                let left_ext =
                                                    !is_bit_set(Some(&dropped_bm), cand_y - 1);
                                                let emit = if !left_ext {
                                                    true
                                                } else {
                                                    let next_x_ok = pos_x + wl < seq_x.len
                                                        && !is_bit_set(
                                                            inv_x.as_deref(),
                                                            pos_x + wl,
                                                        );
                                                    let right_ext = next_x_ok
                                                        && cand_y + wl < seq_y.len
                                                        && !is_bit_set(
                                                            inv_y.as_deref(),
                                                            cand_y + wl,
                                                        )
                                                        && base_2bit(&seq_y.bits, cand_y + wl)
                                                            == base_2bit(&seq_x.bits, pos_x + wl)
                                                        && !is_bit_set(
                                                            Some(&dropped_bm),
                                                            cand_y + 1,
                                                        );
                                                    !right_ext
                                                };
                                                if emit {
                                                    hybrid_hits.push((0, pos_x, cand_y));
                                                }
                                            } else {
                                                let left_ext =
                                                    !is_bit_set(Some(&dropped_bm), cand_y + 1);
                                                let emit = if !left_ext {
                                                    true
                                                } else {
                                                    let next_x_ok = pos_x + wl < seq_x.len
                                                        && !is_bit_set(
                                                            inv_x.as_deref(),
                                                            pos_x + wl,
                                                        );
                                                    let right_ext = next_x_ok
                                                        && cand_y > 0
                                                        && !is_bit_set(
                                                            inv_y.as_deref(),
                                                            cand_y - 1,
                                                        )
                                                        && (base_2bit(&seq_y.bits, cand_y - 1) ^ 3)
                                                            == base_2bit(&seq_x.bits, pos_x + wl)
                                                        && !is_bit_set(
                                                            Some(&dropped_bm),
                                                            cand_y - 1,
                                                        );
                                                    !right_ext
                                                };
                                                if emit {
                                                    hybrid_hits.push((
                                                        1,
                                                        pos_x,
                                                        seq_y.len - cand_y - wl,
                                                    ));
                                                }
                                            }
                                            bp_word |= 1u32 << s;
                                            total_bypassed_q += 1;
                                        }
                                        continue;
                                    }
                                    act_valid = false;
                                }

                                if (s & 15) == 0 {
                                    let w_lsb = load_32_2bit(&seq_x.bits, pos_x);
                                    let kf = reverse_32_2bit(w_lsb) >> rev_rshift;
                                    let kr = (!w_lsb) & mask;
                                    let x_strand = u32::from(kf > kr);
                                    let canon = if x_strand == 0 { kf } else { kr };
                                    let is_pal = kf == kr;
                                    let h = fmix64(canon);
                                    let bkt = (h >> shift_bits) as usize;
                                    let tag = (h >> tag_shift) as u32;
                                    let x_seed_lsb = w_lsb & mask;

                                    let mut found = None;
                                    for j in bucket_off[bkt]..bucket_off[bkt + 1] {
                                        if idx_tags[j] != tag {
                                            continue;
                                        }
                                        let ps = idx_vals[j];
                                        let pos_y = (ps >> 1) as usize;
                                        let y_strand = ps & 1;
                                        let y_word = load_32_2bit(&seq_y.bits, pos_y);
                                        if (is_pal || x_strand == y_strand)
                                            && ((y_word ^ x_seed_lsb) & mask) == 0
                                        {
                                            found = Some((pos_y, 0u32));
                                            break;
                                        }
                                        if compare_reverse && (is_pal || x_strand != y_strand) {
                                            let y_rc_lsb = (!reverse_32_2bit(y_word)) >> rev_rshift;
                                            if ((y_rc_lsb ^ x_seed_lsb) & mask) == 0 {
                                                found = Some((pos_y, 1u32));
                                                break;
                                            }
                                        }
                                    }
                                    match found {
                                        None => bp_word |= 1u32 << s, // Subcase 4a: anchor miss
                                        Some((m_pos_y, m_strand)) => {
                                            act_valid = true;
                                            act_strand = m_strand;
                                            act_pos_y_prev = m_pos_y;
                                            if !is_pal
                                                && !is_bit_set(Some(&ref_non_unique_bm), m_pos_y)
                                            {
                                                if m_strand == 0 {
                                                    if check_lr_fwd(
                                                        &seq_x,
                                                        inv_x.as_deref(),
                                                        pos_x,
                                                        m_pos_y,
                                                    ) {
                                                        hybrid_hits.push((0, pos_x, m_pos_y));
                                                    }
                                                } else if check_lr_rev(
                                                    &seq_x,
                                                    inv_x.as_deref(),
                                                    pos_x,
                                                    m_pos_y,
                                                ) {
                                                    hybrid_hits.push((
                                                        1,
                                                        pos_x,
                                                        seq_y.len - m_pos_y - wl,
                                                    ));
                                                }
                                                bp_word |= 1u32 << s;
                                                total_bypassed_q += 1;
                                            }
                                        }
                                    }
                                }
                            }
                            d_q_bypass_bm[w_out] = bp_word;
                        }
                    }

                    // Simulate l2_count<true> + prefix sum + CTA-privatized l2_scatter<true>
                    let mut d_part_counts = vec![0usize; p_parts];
                    for (w, &bp_w) in d_q_bypass_bm.iter().enumerate() {
                        let mut todo = !bp_w;
                        while todo != 0 {
                            let s = todo.trailing_zeros() as usize;
                            todo &= todo - 1;
                            let pos_x = (w << 5) + s;
                            let w_lsb = load_32_2bit(&seq_x.bits, pos_x);
                            let kf = reverse_32_2bit(w_lsb) >> rev_rshift;
                            let kr = (!w_lsb) & mask;
                            let canon = kf.min(kr);
                            let h = fmix64(canon);
                            let part = ((h >> shift_bits) as usize) >> part_shift;
                            d_part_counts[part] += 1;
                        }
                    }
                    let mut d_part_offsets = vec![0usize; p_parts + 1];
                    for p in 0..p_parts {
                        d_part_offsets[p + 1] = d_part_offsets[p] + d_part_counts[p];
                    }
                    let total_q = d_part_offsets[p_parts];
                    if q_label == "pure_unique_100pct" {
                        // 100%-bypassed tile: verify total_q == 0 still preserves pre-bypass hits!
                        assert_eq!(total_q, 0);
                        assert_eq!(total_bypassed_q, total_valid_q);
                    } else {
                        assert!(total_bypassed_q > 0, "expected active-diagonal bypass hits");
                    }

                    // Two-phase CTA-privatized scatter (simulate 256-thread CTAs, 2 words/thread)
                    let mut d_part_cursor = d_part_offsets[..p_parts].to_vec();
                    let mut d_q_hash = vec![0u64; total_q];
                    let mut d_q_meta = vec![u32::MAX; total_q];
                    let n_threads_p1 = (bp_words + 1) >> 1;
                    let n_ctas = n_threads_p1.div_ceil(256);
                    for cta in 0..n_ctas {
                        let mut s_cnt = vec![0usize; p_parts];
                        // Phase 1: CTA-local histogram
                        for tid in 0..256 {
                            let w0 = (cta * 256 + tid) * 2;
                            for u in 0..2 {
                                let w = w0 + u;
                                if w >= bp_words {
                                    break;
                                }
                                let mut todo = !d_q_bypass_bm[w];
                                while todo != 0 {
                                    let s = todo.trailing_zeros() as usize;
                                    todo &= todo - 1;
                                    let pos_x = (w << 5) + s;
                                    let w_lsb = load_32_2bit(&seq_x.bits, pos_x);
                                    let kf = reverse_32_2bit(w_lsb) >> rev_rshift;
                                    let kr = (!w_lsb) & mask;
                                    let h = fmix64(kf.min(kr));
                                    let part = ((h >> shift_bits) as usize) >> part_shift;
                                    s_cnt[part] += 1;
                                }
                            }
                        }
                        // Phase 1.5: Reserve 1 global slice per non-empty partition per CTA
                        for p in 0..p_parts {
                            let c = s_cnt[p];
                            if c > 0 {
                                s_cnt[p] = d_part_cursor[p];
                                d_part_cursor[p] += c;
                            }
                        }
                        // Phase 2: CTA-local scatter
                        for tid in 0..256 {
                            let w0 = (cta * 256 + tid) * 2;
                            for u in 0..2 {
                                let w = w0 + u;
                                if w >= bp_words {
                                    break;
                                }
                                let mut todo = !d_q_bypass_bm[w];
                                while todo != 0 {
                                    let s = todo.trailing_zeros() as usize;
                                    todo &= todo - 1;
                                    let rel_x = (w << 5) + s;
                                    let w_lsb = load_32_2bit(&seq_x.bits, rel_x);
                                    let kf = reverse_32_2bit(w_lsb) >> rev_rshift;
                                    let kr = (!w_lsb) & mask;
                                    let x_strand = u32::from(kf > kr);
                                    let is_pal = kf == kr;
                                    let h = fmix64(if x_strand == 0 { kf } else { kr });
                                    let part = ((h >> shift_bits) as usize) >> part_shift;
                                    let slot = s_cnt[part];
                                    s_cnt[part] += 1;
                                    d_q_hash[slot] = h;
                                    d_q_meta[slot] =
                                        ((rel_x as u32) << 2) | (u32::from(is_pal) << 1) | x_strand;
                                }
                            }
                        }
                    }

                    // Verify partition bounds & zero collisions/overwrites
                    for p in 0..p_parts {
                        assert_eq!(d_part_cursor[p], d_part_offsets[p + 1]);
                        for q_idx in d_part_offsets[p]..d_part_offsets[p + 1] {
                            assert_ne!(d_q_meta[q_idx], u32::MAX);
                            let h = d_q_hash[q_idx];
                            assert_eq!(((h >> shift_bits) as usize) >> part_shift, p);
                            let meta = d_q_meta[q_idx];
                            let pos_x = (meta >> 2) as usize;
                            let is_pal = ((meta >> 1) & 1) != 0;
                            let x_strand = meta & 1;
                            let bkt = (h >> shift_bits) as usize;
                            let tag = (h >> tag_shift) as u32;
                            let x_seed_lsb = load_32_2bit(&seq_x.bits, pos_x) & mask;
                            for j in bucket_off[bkt]..bucket_off[bkt + 1] {
                                if idx_tags[j] != tag {
                                    continue;
                                }
                                let ps = idx_vals[j];
                                let pos_y = (ps >> 1) as usize;
                                let y_strand = ps & 1;
                                let y_word = load_32_2bit(&seq_y.bits, pos_y);
                                if (is_pal || x_strand == y_strand)
                                    && ((y_word ^ x_seed_lsb) & mask) == 0
                                    && check_lr_fwd(&seq_x, inv_x.as_deref(), pos_x, pos_y)
                                {
                                    hybrid_hits.push((0, pos_x, pos_y));
                                }
                                if compare_reverse && (is_pal || x_strand != y_strand) {
                                    let y_rc_lsb = (!reverse_32_2bit(y_word)) >> rev_rshift;
                                    if ((y_rc_lsb ^ x_seed_lsb) & mask) == 0
                                        && check_lr_rev(&seq_x, inv_x.as_deref(), pos_x, pos_y)
                                    {
                                        hybrid_hits.push((1, pos_x, seq_y.len - pos_y - wl));
                                    }
                                }
                            }
                        }
                    }

                    hybrid_hits.sort_unstable();
                    assert_eq!(
                        hybrid_hits, oracle_hits,
                        "Hybrid L2-Part mismatch vs Monolithic Oracle (q={q_label}, rev={compare_reverse}, P={p_parts})"
                    );
                }
            }
        }
    }

    // ---- Test Group 4: GECKO_PROFILE Env Gating + g++ Header Cross-Check (W3-F11a / W3-F12b) ----

    #[test]
    fn test_gecko_profile_env_gating_and_cpp_header_parity() {
        let cases: &[(&str, bool)] = &[
            ("", false),
            ("   ", false),
            ("0", false),
            (" 0\n", false),
            ("false", false),
            ("FALSE", false),
            ("False", false),
            (" false ", false),
            ("no", false),
            ("NO", false),
            ("off", false),
            ("OFF", false),
            ("1", true),
            (" 1\n", true),
            ("true", true),
            ("TRUE", true),
            ("yes", true),
            ("on", true),
            ("2", true),
        ];

        assert!(!parse_profile_flag(None));
        for &(raw, expected) in cases {
            assert_eq!(
                parse_profile_flag(Some(raw)),
                expected,
                "Rust mismatch on {raw:?}"
            );
            #[cfg(has_cuda_kernels)]
            {
                let cstr = std::ffi::CString::new(raw).unwrap();
                let got = unsafe { gecko_cuda_parse_profile_env_str(cstr.as_ptr()) } != 0;
                assert_eq!(got, expected, "CUDA FFI mismatch on {raw:?}");
            }
        }
        #[cfg(has_cuda_kernels)]
        {
            assert_eq!(
                unsafe { gecko_cuda_parse_profile_env_str(std::ptr::null()) },
                0
            );
        }

        // Cross-compile cuda_common.cuh directly with host `g++ -std=c++17` if available
        if std::process::Command::new("g++")
            .arg("--version")
            .output()
            .is_ok_and(|o| o.status.success())
        {
            let tmp_dir =
                std::env::temp_dir().join(format!("gecko_cuda_hdr_test_{}", std::process::id()));
            let _ = std::fs::create_dir_all(&tmp_dir);
            struct Cleanup(std::path::PathBuf);
            impl Drop for Cleanup {
                fn drop(&mut self) {
                    let _ = std::fs::remove_dir_all(&self.0);
                }
            }
            let _guard = Cleanup(tmp_dir.clone());

            let cpp_path = tmp_dir.join("check_hdr.cpp");
            let bin_path = tmp_dir.join("check_hdr");
            let cpp_src = r#"
#include <cassert>
#include <cstdint>
#include <vector>
#include "cuda_common.cuh"

int main() {
    assert(!gecko_parse_profile_str(nullptr));
    assert(!gecko_parse_profile_str(""));
    assert(!gecko_parse_profile_str("   "));
    assert(!gecko_parse_profile_str("0"));
    assert(!gecko_parse_profile_str(" 0\n"));
    assert(!gecko_parse_profile_str("FALSE"));
    assert(!gecko_parse_profile_str("no"));
    assert(!gecko_parse_profile_str("Off"));
    assert(gecko_parse_profile_str("1"));
    assert(gecko_parse_profile_str(" true "));
    assert(gecko_parse_profile_str("YES"));

    for (uint64_t x : {0ULL, 1ULL, 0x123456789ABCDEF0ULL, ~0ULL}) {
        assert(gecko_fmix64_inv(gecko_fmix64(x)) == x);
    }
    // 4-bp palindrome ACGT: A=0,C=1,G=2,T=3 -> 0b00_01_10_11 = 0x1B
    assert(gecko_is_palindromic_canon(0x1BULL, 4, 0xFFULL));
    // Non-palindrome AAAA = 0
    assert(!gecko_is_palindromic_canon(0ULL, 4, 0xFFULL));

    GeckoRun runs[2] = {{3U, 10U}, {30U, 98U}};
    std::vector<uint32_t> bm(((100 + 31) / 32) + 2, 0U);
    gecko_fill_runs_bitmap_host(runs, 2, 100, bm.data());
    assert(bm[0] == (((0xFFFFFFFFU << 3) & (0xFFFFFFFFU >> (31 - 9))) | (0xFFFFFFFFU << 30)));
    assert(bm[1] == 0xFFFFFFFFU);
    assert(bm[2] == 0xFFFFFFFFU);
    assert(bm[3] == (0xFFFFFFFFU >> (31 - 1)));
    assert(bm[4] == 0U);
    assert(bm[5] == 0U);

    // W3-M4 C-ABI layout & minimizer param parity checks
    static_assert(sizeof(GeckoCudaConfig) == 48, "48B config");
    static_assert(sizeof(GeckoCudaResult) == 56, "56B result");
    static_assert(sizeof(GeckoCudaScratchPoolStats) == 48, "48B pool stats");
    for (int wl = 4; wl <= 34; wl++) {
        int kp = 0, w = 0;
        bool ok = gecko_minimizer_params(wl, kp, w);
        if (wl < 18 || wl > 32) {
            assert(!ok);
        } else {
            assert(ok);
            assert(kp >= 14 && kp <= 21);
            assert(w >= 5 && w <= 12);
            assert(kp + w - 1 == wl);
        }
    }
    assert(!gecko_cuda_should_use_minimizer(GECKO_SEEDING_EXACT, 32, false, 0, 20));
    assert(gecko_cuda_should_use_minimizer(GECKO_SEEDING_MINIMIZER, 32, false, 0, 20));
    assert(gecko_cuda_should_use_minimizer(GECKO_SEEDING_AUTO, 18, false, 0, 20));
    assert(!gecko_cuda_should_use_minimizer(GECKO_SEEDING_AUTO, 17, false, 0, 20));
    assert(!gecko_cuda_should_use_minimizer(GECKO_SEEDING_AUTO, 32, true, 0, 20));
    assert(!gecko_cuda_should_use_minimizer(GECKO_SEEDING_AUTO, 32, false, 40, 20));
    assert(!gecko_cuda_should_use_minimizer(GECKO_SEEDING_AUTO, 32, false, 0, 0));

    const uint64_t mix_xs[] = {0ULL, 1ULL, 0x3FFFFFFFFFFULL, 0x0000ABCDEF012345ULL, GECKO_MASK48};
    for (uint64_t x : mix_xs) {
        uint64_t h = gecko_mix48(x);
        assert((h & ~GECKO_MASK48) == 0ULL);
        assert(gecko_mix48_inv(h) == x);
    }
    return 0;
}
"#;
            std::fs::write(&cpp_path, cpp_src).expect("write check_hdr.cpp");
            let inc_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("cuda");
            let compile = std::process::Command::new("g++")
                .args([
                    "-O2",
                    "-std=c++17",
                    "-fsanitize=address,undefined",
                    "-Wall",
                    "-Wextra",
                    "-Werror",
                ])
                .arg(format!("-I{}", inc_dir.display()))
                .arg(&cpp_path)
                .arg("-o")
                .arg(&bin_path)
                .output()
                .expect("run g++");
            assert!(
                compile.status.success(),
                "g++ failed to compile cuda_common.cuh:\n{}",
                String::from_utf8_lossy(&compile.stderr)
            );
            let run = std::process::Command::new(&bin_path)
                .output()
                .expect("run check_hdr");
            assert!(
                run.status.success(),
                "check_hdr failed:\n{}",
                String::from_utf8_lossy(&run.stderr)
            );
        }

        // Cross-check compiled CUDA host C-ABI wrappers against gecko-core
        #[cfg(has_cuda_kernels)]
        {
            for wl in 4..=34 {
                let mut kp = 0;
                let mut w = 0;
                let ok = unsafe { gecko_cuda_minimizer_params_host(wl, &mut kp, &mut w) } != 0;
                let core = gecko_core::minimizer::minimizer_params(wl as usize);
                match core {
                    Some((ckp, cw)) => {
                        assert!(ok);
                        assert_eq!(kp as usize, ckp);
                        assert_eq!(w as usize, cw);
                    }
                    None => assert!(!ok),
                }
            }
            for x in [0u64, 1, 0x3FF_FFFF_FFFF, 0x1234_5678_9ABC, 0xFFFF_FFFF_FFFF] {
                let h = unsafe { gecko_cuda_mix48_host(x) };
                assert!(h <= 0x0000_FFFF_FFFF_FFFF);
                assert_eq!(unsafe { gecko_cuda_mix48_inv_host(h) }, x);
            }
            // Smoke test two-hit host scan FFI on a single-diag hit stream (wl=12, W=30)
            let tile_end = 1000u64;
            let len_y = 1000u64;
            let hits = [100u64, 111u64, 112u64, 155u64]; // same diag=0, rel_x=100,111,112,155
            let mut flags = [0u8; 4];
            unsafe {
                gecko_cuda_two_hit_scan_host(
                    hits.as_ptr(),
                    hits.len() as u32,
                    0,
                    tile_end,
                    len_y,
                    12,
                    0,
                    30,
                    flags.as_mut_ptr(),
                );
            }
            // 100: first -> 0; 111: dist 11 < 12 -> 0; 112: dist 12 to 100 -> 1; 155: dist 43 > 42 -> 0
            assert_eq!(flags, [0, 0, 1, 0]);
        }
    }

    // ---- Test Group 5: Persistent CudaScratchPool API, State Machine & Concurrency (W3-F15a) ----

    #[test]
    fn test_cuda_scratch_pool_api_no_gpu_safe() {
        let orig = CudaScratchPool::is_enabled();
        let st = CudaScratchPool::stats();
        assert_eq!(st.idle_count, 0);
        CudaScratchPool::set_max_idle(2);
        CudaScratchPool::set_enabled(false);
        assert!(!CudaScratchPool::is_enabled());
        CudaScratchPool::set_enabled(true);
        assert!(CudaScratchPool::is_enabled());
        CudaScratchPool::release();
        {
            let _guard = CudaScratchPool::scoped_release();
        }
        CudaScratchPool::set_max_idle(1);
        CudaScratchPool::set_enabled(orig);
    }

    #[test]
    fn test_cuda_scratch_pool_state_machine_cpp_and_concurrency() {
        if !std::process::Command::new("g++")
            .arg("--version")
            .output()
            .is_ok_and(|o| o.status.success())
        {
            return;
        }
        let tmp_dir =
            std::env::temp_dir().join(format!("gecko_cuda_pool_test_{}", std::process::id()));
        let _ = std::fs::create_dir_all(&tmp_dir);
        struct Cleanup(std::path::PathBuf);
        impl Drop for Cleanup {
            fn drop(&mut self) {
                let _ = std::fs::remove_dir_all(&self.0);
            }
        }
        let _guard = Cleanup(tmp_dir.clone());

        let cpp_path = tmp_dir.join("test_pool.cpp");
        let bin_path = tmp_dir.join("test_pool");
        let cpp_src = r#"
#include <atomic>
#include <cassert>
#include <cstdint>
#include <thread>
#include <vector>
#include "cuda_common.cuh"

struct MockScratch {
    int id = 0;
    uint64_t bytes = 0;
};

int main() {
    std::atomic<uint64_t> freed{0};
    auto bytes_fn = [](const MockScratch &s) -> uint64_t { return s.bytes; };
    auto free_fn  = [&](MockScratch &s) {
        if (s.id != 0) freed.fetch_add(1, std::memory_order_relaxed);
        s = MockScratch{};
    };

    // 1. State machine: miss -> healthy checkin -> LIFO hit -> poisoned eviction -> cap eviction
    GeckoScratchPoolCore<MockScratch, decltype(bytes_fn), decltype(free_fn)> pool(bytes_fn, free_fn, true, 1);
    MockScratch item{};
    assert(!pool.checkout(item));
    assert(pool.stats().misses == 1 && pool.stats().hits == 0);

    item = MockScratch{1, 1024};
    assert(pool.checkin(std::move(item), /*healthy=*/true));
    assert(pool.stats().recycles == 1 && pool.stats().idle_count == 1 && pool.stats().cached_bytes == 1024);

    MockScratch got{};
    assert(pool.checkout(got));
    assert(got.id == 1 && got.bytes == 1024);
    assert(pool.stats().hits == 1 && pool.stats().idle_count == 0 && pool.stats().cached_bytes == 0);

    // Poisoned / unhealthy checkin must immediately free and NEVER recycle
    assert(!pool.checkin(std::move(got), /*healthy=*/false));
    assert(freed.load() == 1 && pool.stats().evictions == 1 && pool.stats().idle_count == 0);
    assert(!pool.checkout(got));

    // Capacity bound max_idle=1: second healthy checkin evicts cleanly outside lock
    assert(pool.checkin(MockScratch{2, 2048}, true));
    assert(!pool.checkin(MockScratch{3, 4096}, true));
    assert(freed.load() == 2 && pool.stats().idle_count == 1 && pool.stats().cached_bytes == 2048);

    // Increase cap to 2, push LIFO, verify LIFO pop & set_max_idle trim
    pool.set_max_idle(2);
    assert(pool.checkin(MockScratch{4, 4096}, true));
    assert(pool.stats().idle_count == 2 && pool.stats().cached_bytes == 6144);
    pool.set_max_idle(1);
    assert(freed.load() == 3 && pool.stats().idle_count == 1 && pool.stats().cached_bytes == 2048);

    // set_enabled(false) drains all idle entries and rejects caching
    pool.set_enabled(false);
    assert(!pool.is_enabled());
    assert(freed.load() == 4 && pool.stats().idle_count == 0);
    assert(!pool.checkin(MockScratch{5, 1024}, true));
    assert(freed.load() == 5);

    // Re-enable and verify release_all()
    pool.set_enabled(true);
    pool.set_max_idle(4);
    assert(pool.checkin(MockScratch{6, 512}, true));
    assert(pool.checkin(MockScratch{7, 512}, true));
    assert(pool.release_all() == 2);
    assert(freed.load() == 7 && pool.stats().idle_count == 0);

    // 2. Multi-threaded concurrency stress: 8 threads x 500 cycles
    std::atomic<uint64_t> created{0};
    std::atomic<uint64_t> destroyed{0};
    auto stress_free = [&](MockScratch &s) {
        if (s.id != 0) destroyed.fetch_add(1, std::memory_order_relaxed);
        s = MockScratch{};
    };
    GeckoScratchPoolCore<MockScratch, decltype(bytes_fn), decltype(stress_free)> spool(
        bytes_fn, stress_free, true, 2);

    std::vector<std::thread> threads;
    for (int t = 0; t < 8; t++) {
        threads.emplace_back([&, t]() {
            for (int iter = 0; iter < 500; iter++) {
                MockScratch s{};
                if (!spool.checkout(s)) {
                    uint64_t id = created.fetch_add(1, std::memory_order_relaxed) + 1;
                    s = MockScratch{(int)id, 608ULL << 20};
                }
                bool healthy = ((iter + t) % 13) != 0;
                spool.checkin(std::move(s), healthy);
                if ((iter % 97) == 0) (void)spool.stats();
                if ((iter % 199) == 0 && t == 0) spool.release_all();
            }
        });
    }
    for (auto &th : threads) th.join();
    uint64_t idle_before = spool.stats().idle_count;
    assert(created.load() == destroyed.load() + idle_before);
    spool.release_all();
    assert(created.load() == destroyed.load());
    assert(spool.stats().idle_count == 0 && spool.stats().cached_bytes == 0);
    return 0;
}
"#;
        std::fs::write(&cpp_path, cpp_src).expect("write test_pool.cpp");
        let inc_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("cuda");
        let compile = std::process::Command::new("g++")
            .args([
                "-O2",
                "-std=c++17",
                "-pthread",
                "-fsanitize=address,undefined",
                "-Wall",
                "-Wextra",
                "-Werror",
            ])
            .arg(format!("-I{}", inc_dir.display()))
            .arg(&cpp_path)
            .arg("-o")
            .arg(&bin_path)
            .output()
            .expect("run g++");
        assert!(
            compile.status.success(),
            "g++ failed to compile test_pool.cpp:\n{}",
            String::from_utf8_lossy(&compile.stderr)
        );
        let run = std::process::Command::new(&bin_path)
            .output()
            .expect("run test_pool");
        assert!(
            run.status.success(),
            "test_pool failed:\nstdout:\n{}\nstderr:\n{}",
            String::from_utf8_lossy(&run.stdout),
            String::from_utf8_lossy(&run.stderr)
        );
    }

    // ---- Test Group 6: Host C++ ASan/UBSan Parity Suites P1–P7 (W3-F14 + W3-F15b) ----

    fn read_fasta_dna_slice(path: &std::path::Path, start: usize, len: usize) -> Option<Vec<u8>> {
        use std::io::{BufRead, BufReader};
        let file = std::fs::File::open(path).ok()?;
        let reader = BufReader::new(file);
        let mut out = Vec::with_capacity(len);
        let mut seen = 0usize;
        let end = start + len;
        for line in reader.lines() {
            let line = line.ok()?;
            if line.starts_with('>') {
                continue;
            }
            let bytes = line.as_bytes();
            let next_seen = seen + bytes.len();
            if next_seen > start {
                let lo = start.saturating_sub(seen);
                let hi = (end - seen).min(bytes.len());
                out.extend_from_slice(&bytes[lo..hi]);
                if out.len() >= len {
                    break;
                }
            }
            seen = next_seen;
        }
        if out.len() < 1024 {
            return None;
        }
        Some(out)
    }

    #[test]
    fn test_w3_m4_cpp_asan_ubsan_parity_suites_p1_to_p7() {
        if !std::process::Command::new("g++")
            .arg("--version")
            .output()
            .is_ok_and(|o| o.status.success())
        {
            return;
        }
        let tmp_dir =
            std::env::temp_dir().join(format!("gecko_cuda_w3_m4_parity_{}", std::process::id()));
        let _ = std::fs::create_dir_all(&tmp_dir);
        struct Cleanup(std::path::PathBuf);
        impl Drop for Cleanup {
            fn drop(&mut self) {
                let _ = std::fs::remove_dir_all(&self.0);
            }
        }
        let _guard = Cleanup(tmp_dir.clone());

        // Prepare optional real chimp_chr21 x gorilla_chr21 slices from repo root (skip gracefully if absent)
        let repo_root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .parent()
            .unwrap();
        let chimp_path = repo_root.join("chimp_chr21.fa");
        let gorilla_path = repo_root.join("gorilla_chr21.fa");
        let mut real_args: Vec<String> = Vec::new();
        if chimp_path.exists() && gorilla_path.exists() {
            let slices = [
                ("slice200k", 25_000_000usize, 29_650_000usize, 200_000usize),
                ("slice64k", 10_000_000usize, 14_634_652usize, 65_536usize),
            ];
            for (tag, c_off, g_off, slen) in slices {
                if let (Some(c_dna), Some(g_dna)) = (
                    read_fasta_dna_slice(&chimp_path, c_off, slen),
                    read_fasta_dna_slice(&gorilla_path, g_off, slen),
                ) {
                    let cp = tmp_dir.join(format!("{tag}_chimp.raw"));
                    let gp = tmp_dir.join(format!("{tag}_gorilla.raw"));
                    std::fs::write(&cp, &c_dna).expect("write chimp slice");
                    std::fs::write(&gp, &g_dna).expect("write gorilla slice");
                    real_args.push(cp.to_string_lossy().into_owned());
                    real_args.push(gp.to_string_lossy().into_owned());
                }
            }
        }

        let cpp_path = tmp_dir.join("parity_p1_p7.cpp");
        let bin_path = tmp_dir.join("parity_p1_p7");
        let cpp_src = r##"
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <deque>
#include <fstream>
#include <map>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include "cuda_common.cuh"

// ---- Host 2-bit PackedSeq representation mirroring gecko_core::seq::PackedSeq ----
struct HostSeq {
    std::vector<uint64_t> bits;
    std::vector<uint32_t> inv_bm; // always allocated for simplicity
    uint64_t len = 0;

    static HostSeq from_ascii(const std::string &raw, bool soft_mask) {
        HostSeq s;
        s.len = raw.size();
        s.bits.assign((s.len + 31ULL) / 32ULL + 2ULL, 0ULL);
        s.inv_bm.assign((s.len + 31ULL) / 32ULL + 2ULL, 0U);
        for (size_t i = 0; i < raw.size(); i++) {
            unsigned char c = (unsigned char)raw[i];
            bool is_lower = (c >= 'a' && c <= 'z');
            unsigned char u = is_lower ? (unsigned char)(c - 'a' + 'A') : c;
            uint64_t b = 0;
            bool valid = true;
            if (u == 'A') b = 0;
            else if (u == 'C') b = 1;
            else if (u == 'G') b = 2;
            else if (u == 'T') b = 3;
            else valid = false;
            if (valid) {
                s.bits[i >> 5] |= (b << ((i & 31U) << 1));
            }
            if (!valid || (soft_mask && is_lower)) {
                s.inv_bm[i >> 5] |= (1U << (i & 31U));
            }
        }
        return s;
    }

    HostSeq revcomp() const {
        HostSeq rc;
        rc.len = len;
        rc.bits.assign((len + 31ULL) / 32ULL + 2ULL, 0ULL);
        rc.inv_bm.assign((len + 31ULL) / 32ULL + 2ULL, 0U);
        for (uint64_t i = 0; i < len; i++) {
            uint64_t src = len - 1ULL - i;
            uint64_t b = ((bits[src >> 5] >> ((src & 31U) << 1)) & 3ULL) ^ 3ULL;
            rc.bits[i >> 5] |= (b << ((i & 31U) << 1));
            if (gecko_is_bit_set(inv_bm.data(), (uint32_t)src)) {
                rc.inv_bm[i >> 5] |= (1U << (i & 31U));
            }
        }
        return rc;
    }
};

static inline uint32_t host_base_2bit(const HostSeq &s, uint32_t pos) {
    return (uint32_t)((s.bits[pos >> 5] >> ((pos & 31U) << 1)) & 3ULL);
}

static inline uint64_t host_kp_score(const HostSeq &s, uint32_t pos, int kp) {
    if ((uint64_t)pos + (uint64_t)kp > s.len) return ~0ULL;
    uint32_t vmask = (1U << kp) - 1U;
    if ((gecko_load_32_1bit(s.inv_bm.data(), pos) & vmask) != 0U) return ~0ULL;
    uint64_t w_lsb = gecko_load_32_2bit(s.bits.data(), pos);
    uint64_t kf = gecko_reverse_32_2bit(w_lsb) >> (64 - 2 * kp);
    uint64_t kr = (~w_lsb) & ((1ULL << (2 * kp)) - 1ULL);
    uint64_t canon = (kf <= kr) ? kf : kr;
    return gecko_fmix64(canon);
}

static inline void host_kp_canon(const HostSeq &s, uint32_t pos, int kp,
                                 uint64_t &canon, uint32_t &strand, uint32_t &is_pal) {
    uint64_t w_lsb = gecko_load_32_2bit(s.bits.data(), pos);
    uint64_t kf = gecko_reverse_32_2bit(w_lsb) >> (64 - 2 * kp);
    uint64_t kr = (~w_lsb) & ((1ULL << (2 * kp)) - 1ULL);
    strand = (kf <= kr) ? 0U : 1U;
    is_pal = (kf == kr) ? 1U : 0U;
    canon = strand ? kr : kf;
}

static inline void host_kmer_canon(const HostSeq &s, uint32_t pos, int wl, uint64_t eff_mask,
                                   uint64_t &canon, uint32_t &strand, bool &is_pal) {
    uint64_t w_lsb = gecko_load_32_2bit(s.bits.data(), pos);
    uint64_t kf = (gecko_reverse_32_2bit(w_lsb) >> (64 - 2 * wl)) & eff_mask;
    uint64_t kr = (~w_lsb) & eff_mask;
    strand = (kf <= kr) ? 0U : 1U;
    is_pal = (kf == kr);
    canon = strand ? kr : kf;
}

static std::string make_revcomp_str(const std::string &s) {
    std::string out(s.size(), 'N');
    for (size_t i = 0; i < s.size(); i++) {
        char c = s[s.size() - 1 - i];
        switch (c) {
            case 'A': out[i] = 'T'; break; case 'C': out[i] = 'G'; break;
            case 'G': out[i] = 'C'; break; case 'T': out[i] = 'A'; break;
            case 'a': out[i] = 't'; break; case 'c': out[i] = 'g'; break;
            case 'g': out[i] = 'c'; break; case 't': out[i] = 'a'; break;
            default: out[i] = c; break;
        }
    }
    return out;
}

static std::string rand_dna(size_t n, uint64_t seed) {
    static const char A[4] = {'A', 'C', 'G', 'T'};
    std::string out(n, 'A');
    uint64_t rng = seed;
    for (size_t i = 0; i < n; i++) {
        rng = rng * 6364136223846793005ULL + 1442695040888963407ULL;
        out[i] = A[(rng >> 32) & 3ULL];
    }
    return out;
}

// ---- Exact Stride-1 Reference Oracle (dropped_bitmap + L/R-maximal exact hits) ----
struct OracleResult {
    std::vector<uint32_t> dropped_bm;
    std::vector<std::tuple<uint8_t, uint32_t, uint32_t>> hits;
};

static OracleResult run_exact_stride1_oracle(
    const HostSeq &seqX, const HostSeq &seqY, int wl, uint32_t freq_limit, bool compare_reverse
) {
    OracleResult res;
    uint64_t span_mask = (wl == 32) ? ~0ULL : ((1ULL << (2 * wl)) - 1ULL);
    uint32_t vmask = (wl == 32) ? 0xFFFFFFFFU : ((1U << wl) - 1U);
    res.dropped_bm.assign((seqY.len + 31ULL) / 32ULL + 2ULL, 0U);

    // Build ref map: canon -> vector of (posY, strand)
    std::unordered_map<uint64_t, std::vector<std::pair<uint32_t, uint32_t>>> idx;
    if (seqY.len >= (uint64_t)wl) {
        for (uint32_t p = 0; p + (uint32_t)wl <= (uint32_t)seqY.len; p++) {
            if ((gecko_load_32_1bit(seqY.inv_bm.data(), p) & vmask) != 0U) continue;
            uint64_t canon; uint32_t st; bool pal;
            host_kmer_canon(seqY, p, wl, span_mask, canon, st, pal);
            idx[canon].push_back({p, st});
        }
    }
    for (auto &kv : idx) {
        uint32_t c0 = 0, c1 = 0;
        for (auto &ps : kv.second) { if (ps.second) c1++; else c0++; }
        std::vector<std::pair<uint32_t, uint32_t>> kept;
        for (auto &ps : kv.second) {
            bool ok = ps.second ? (c1 <= freq_limit) : (c0 <= freq_limit);
            if (ok) kept.push_back(ps);
            else res.dropped_bm[ps.first >> 5] |= (1U << (ps.first & 31U));
        }
        kv.second.swap(kept);
    }

    auto check_lr_fwd = [&](uint32_t px, uint32_t py) -> bool {
        bool left = (px > 0 && py > 0 &&
                     !gecko_is_bit_set(seqX.inv_bm.data(), px - 1) &&
                     !gecko_is_bit_set(seqY.inv_bm.data(), py - 1) &&
                     host_base_2bit(seqX, px - 1) == host_base_2bit(seqY, py - 1) &&
                     !gecko_is_bit_set(res.dropped_bm.data(), py - 1));
        if (!left) return true;
        bool right = (px + (uint32_t)wl < seqX.len && py + (uint32_t)wl < seqY.len &&
                      !gecko_is_bit_set(seqX.inv_bm.data(), px + wl) &&
                      !gecko_is_bit_set(seqY.inv_bm.data(), py + wl) &&
                      host_base_2bit(seqX, px + wl) == host_base_2bit(seqY, py + wl) &&
                      !gecko_is_bit_set(res.dropped_bm.data(), py + 1));
        return !right;
    };
    auto check_lr_rev = [&](uint32_t px, uint32_t py) -> bool {
        uint32_t py_rc = (uint32_t)(seqY.len - py - wl);
        bool left = (px > 0 && py_rc > 0 &&
                     !gecko_is_bit_set(seqX.inv_bm.data(), px - 1) &&
                     !gecko_is_bit_set(seqY.inv_bm.data(), py + wl) &&
                     (host_base_2bit(seqY, py + wl) ^ 3U) == host_base_2bit(seqX, px - 1) &&
                     !gecko_is_bit_set(res.dropped_bm.data(), py + 1));
        if (!left) return true;
        bool right = (px + (uint32_t)wl < seqX.len && py > 0 &&
                      !gecko_is_bit_set(seqX.inv_bm.data(), px + wl) &&
                      !gecko_is_bit_set(seqY.inv_bm.data(), py - 1) &&
                      (host_base_2bit(seqY, py - 1) ^ 3U) == host_base_2bit(seqX, px + wl) &&
                      !gecko_is_bit_set(res.dropped_bm.data(), py - 1));
        return !right;
    };

    if (seqX.len >= (uint64_t)wl) {
        for (uint32_t px = 0; px + (uint32_t)wl <= (uint32_t)seqX.len; px++) {
            if ((gecko_load_32_1bit(seqX.inv_bm.data(), px) & vmask) != 0U) continue;
            uint64_t canon; uint32_t xst; bool pal;
            host_kmer_canon(seqX, px, wl, span_mask, canon, xst, pal);
            auto it = idx.find(canon);
            if (it == idx.end()) continue;
            for (auto &ps : it->second) {
                uint32_t py = ps.first, yst = ps.second;
                if ((pal || xst == yst) && check_lr_fwd(px, py)) {
                    res.hits.push_back({0, px, py});
                }
                if (compare_reverse && (pal || xst != yst) && check_lr_rev(px, py)) {
                    res.hits.push_back({1, px, (uint32_t)(seqY.len - py - wl)});
                }
            }
        }
    }
    std::sort(res.hits.begin(), res.hits.end());
    return res;
}

// ---- Full Host Simulation of GPU Minimizer Pipeline (using EXACT cuda_common.cuh helpers!) ----
struct GpuMinSimResult {
    std::vector<uint32_t> dropped_bm;
    std::vector<std::tuple<uint8_t, uint32_t, uint32_t>> hits;
    uint32_t ovf_routed_keys = 0;
    uint32_t ovf_kmers = 0;
    uint64_t exact_k_accept = 0;
    uint64_t exact_k_reject = 0;
};

static GpuMinSimResult run_gpu_min_pipeline_sim(
    const HostSeq &seqX, const HostSeq &seqY, int wl, uint32_t freq_limit,
    bool compare_reverse, uint64_t tile_bp
) {
    GpuMinSimResult res;
    int kp = 0, w = 0;
    assert(gecko_minimizer_params(wl, kp, w));
    uint32_t m_cap = gecko_minimizer_m_cap(freq_limit);
    uint64_t span_mask = (wl == 32) ? ~0ULL : ((1ULL << (2 * wl)) - 1ULL);
    int rev_rshift = 64 - 2 * wl;
    uint32_t H = (uint32_t)(w - 1);
    res.dropped_bm.assign((seqY.len + 31ULL) / 32ULL + 2ULL, 0U);

    // Step 1: CTA wave-by-wave extraction (512 bp/wave + halo H) with EMIT_TIES = true
    std::vector<std::pair<uint64_t, uint32_t>> raw_min;
    uint64_t total_kp_y = (seqY.len >= (uint64_t)kp) ? (seqY.len - (uint64_t)kp + 1ULL) : 0ULL;
    for (uint64_t w_start = 0; w_start < total_kp_y; w_start += 512ULL) {
        uint32_t wave_len = (uint32_t)std::min<uint64_t>(512ULL, total_kp_y - w_start);
        std::vector<uint64_t> s_scores(wave_len + 2U * H, ~0ULL);
        for (uint32_t i = 0; i < wave_len + 2U * H; i++) {
            int64_t gp = (int64_t)w_start - (int64_t)H + (int64_t)i;
            if (gp >= 0 && (uint64_t)gp < total_kp_y) {
                s_scores[i] = host_kp_score(seqY, (uint32_t)gp, kp);
            }
        }
        for (uint32_t p = 0; p < wave_len; p++) {
            int L = 0, R = 0;
            if (gecko_min_window_eval(s_scores.data(), (int)(p + H), w, true, L, R)) {
                uint32_t gpos = (uint32_t)(w_start + p);
                uint64_t canon; uint32_t st, pal;
                host_kp_canon(seqY, gpos, kp, canon, st, pal);
                raw_min.push_back({gecko_mix48(canon), (gpos << 1) | st});
            }
        }
    }

    // Step 2: 48-bit sort
    std::sort(raw_min.begin(), raw_min.end(), [](const auto &a, const auto &b) {
        return a.first < b.first;
    });
    uint32_t n_raw = (uint32_t)raw_min.size();
    std::vector<uint8_t> min_keep(n_raw, 1U);

    // Step 3: Exact HF freq_limit + M_CAP overflow CSR
    struct OvfCsr {
        std::vector<uint32_t> off, tags, vals;
        int bbits = 10, shift = 54, tshift = 22;
    } ovf;

    auto scan_y = [&](uint32_t p_r, bool emit_ties, uint32_t &oL, uint32_t &oR) -> bool {
        uint64_t s = host_kp_score(seqY, p_r, kp);
        if (s == ~0ULL) return false;
        uint32_t L = 0, R = 0;
        for (uint32_t d = 1; d <= H; d++) {
            if (d > p_r) break;
            uint64_t ls = host_kp_score(seqY, p_r - d, kp);
            if (ls == ~0ULL || (emit_ties ? (ls < s) : (ls <= s))) break;
            L = d;
        }
        for (uint32_t d = 1; d <= H; d++) {
            uint64_t rs = host_kp_score(seqY, p_r + d, kp);
            if (rs == ~0ULL || rs < s) break;
            R = d;
        }
        oL = L; oR = R;
        return (L + R >= H);
    };

    if (freq_limit < n_raw) {
        // Step 3a: find HF indices (run_len > freq_limit)
        std::vector<uint32_t> hf_idx;
        for (uint32_t i = 0; i < n_raw; ) {
            uint32_t j = i + 1;
            while (j < n_raw && raw_min[j].first == raw_min[i].first) j++;
            if (j - i > freq_limit) {
                for (uint32_t m = i; m < j; m++) hf_idx.push_back(m);
            }
            i = j;
        }
        if (!hf_idx.empty()) {
            // Step 3b: exact-once L_strict / R_ge candidate k-mers
            std::vector<uint32_t> hf_r_lo(hf_idx.size(), 0U);
            std::vector<uint8_t> hf_r_cnt(hf_idx.size(), 0U);
            std::vector<std::pair<uint64_t, uint32_t>> cand;
            for (size_t t = 0; t < hf_idx.size(); t++) {
                uint32_t p_r = raw_min[hf_idx[t]].second >> 1;
                uint32_t L = 0, R = 0;
                if (scan_y(p_r, false, L, R)) {
                    uint32_t r_lo = p_r - L;
                    uint32_t cnt = L + R - H + 1U;
                    hf_r_lo[t] = r_lo;
                    hf_r_cnt[t] = (uint8_t)cnt;
                    for (uint32_t d = 0; d < cnt; d++) {
                        uint32_t r = r_lo + d;
                        uint64_t canon; uint32_t st; bool pal;
                        host_kmer_canon(seqY, r, wl, span_mask, canon, st, pal);
                        cand.push_back({gecko_fmix64(canon), (r << 1) | st});
                    }
                }
            }
            // Step 3c: exact per-strand freq_limit -> dropped_bm
            std::sort(cand.begin(), cand.end(), [](const auto &a, const auto &b) {
                return a.first < b.first;
            });
            bool any_drop = false;
            for (size_t i = 0; i < cand.size(); ) {
                size_t j = i + 1;
                while (j < cand.size() && cand[j].first == cand[i].first) j++;
                if (j - i > freq_limit) {
                    uint32_t c0 = 0, c1 = 0;
                    for (size_t m = i; m < j; m++) { if (cand[m].second & 1U) c1++; else c0++; }
                    bool k0 = (c0 <= freq_limit), k1 = (c1 <= freq_limit);
                    if (!k0 || !k1) {
                        any_drop = true;
                        for (size_t m = i; m < j; m++) {
                            bool ok = (cand[m].second & 1U) ? k1 : k0;
                            if (!ok) {
                                uint32_t r = cand[m].second >> 1;
                                res.dropped_bm[r >> 5] |= (1U << (r & 31U));
                            }
                        }
                    }
                }
                i = j;
            }
            // Step 3d: L_ge/R_ge survival check + M_CAP overflow classification
            if (any_drop) {
                for (size_t t = 0; t < hf_idx.size(); t++) {
                    uint32_t idx = hf_idx[t];
                    uint32_t p_r = raw_min[idx].second >> 1;
                    uint32_t L = 0, R = 0;
                    if (!scan_y(p_r, true, L, R)) { min_keep[idx] = 0U; continue; }
                    uint32_t r_lo = p_r - L;
                    uint32_t cnt = L + R - H + 1U;
                    uint32_t db = gecko_load_32_1bit(res.dropped_bm.data(), r_lo);
                    uint32_t mask = (1U << cnt) - 1U;
                    min_keep[idx] = ((db & mask) != mask) ? 1U : 0U;
                }
            }
            std::vector<uint32_t> keep_pref(n_raw + 1U, 0U);
            for (uint32_t i = 0; i < n_raw; i++) keep_pref[i + 1U] = keep_pref[i] + min_keep[i];

            std::vector<uint8_t> hf_is_ovf(hf_idx.size(), 0U);
            for (size_t t = 0; t < hf_idx.size(); t++) {
                uint32_t idx = hf_idx[t];
                uint64_t k = raw_min[idx].first;
                uint32_t rs = idx;
                while (rs > 0 && raw_min[rs - 1].first == k) rs--;
                uint32_t re = idx + 1;
                while (re < n_raw && raw_min[re].first == k) re++;
                uint32_t kept_cnt = keep_pref[re] - keep_pref[rs];
                if (kept_cnt > m_cap) {
                    hf_is_ovf[t] = 1U;
                    if (idx == rs) { min_keep[idx] = 2U; res.ovf_routed_keys++; }
                    else min_keep[idx] = 0U;
                }
            }

            // Step 3e: Build exact-k Overflow CSR
            if (res.ovf_routed_keys > 0U) {
                std::vector<std::pair<uint64_t, uint32_t>> ovf_kv;
                for (size_t t = 0; t < hf_idx.size(); t++) {
                    if (!hf_is_ovf[t]) continue;
                    for (uint32_t d = 0; d < hf_r_cnt[t]; d++) {
                        uint32_t r = hf_r_lo[t] + d;
                        if (gecko_is_bit_set(res.dropped_bm.data(), r)) continue;
                        uint64_t canon; uint32_t st; bool pal;
                        host_kmer_canon(seqY, r, wl, span_mask, canon, st, pal);
                        ovf_kv.push_back({gecko_fmix64(canon), (r << 1) | st});
                    }
                }
                res.ovf_kmers = (uint32_t)ovf_kv.size();
                if (!ovf_kv.empty()) {
                    std::sort(ovf_kv.begin(), ovf_kv.end(), [](const auto &a, const auto &b) {
                        return a.first < b.first;
                    });
                    ovf.bbits = 10;
                    while (ovf.bbits < 24 && (1ULL << ovf.bbits) * 4ULL < ovf_kv.size()) ovf.bbits++;
                    ovf.shift = 64 - ovf.bbits;
                    ovf.tshift = ovf.shift - 32;
                    uint32_t ob_cnt = 1U << ovf.bbits;
                    ovf.off.assign(ob_cnt + 1U, 0U);
                    ovf.tags.resize(ovf_kv.size());
                    ovf.vals.resize(ovf_kv.size());
                    for (size_t i = 0; i < ovf_kv.size(); i++) {
                        uint64_t h = ovf_kv[i].first;
                        uint32_t b = (uint32_t)(h >> ovf.shift);
                        ovf.off[b + 1U]++;
                        ovf.tags[i] = (uint32_t)(h >> ovf.tshift);
                        ovf.vals[i] = ovf_kv[i].second;
                    }
                    for (uint32_t b = 0; b < ob_cnt; b++) ovf.off[b + 1U] += ovf.off[b];
                }
            }
        }
    }

    // Step 4: Main 48-bit CSR with bucket_bits in [17 .. 26] + bit-31 sentinel
    uint32_t n_kept = 0;
    for (uint8_t a : min_keep) if (a != 0U) n_kept++;
    int bbits = 17;
    while (bbits < 26 && (1ULL << bbits) * 4ULL < (uint64_t)n_kept) bbits++;
    uint32_t nb = 1U << bbits;
    int shift_bits = 48 - bbits;
    uint32_t tag_mask = (1U << shift_bits) - 1U;
    assert((tag_mask & GECKO_MIN_OVF_SENTINEL) == 0U);

    std::vector<uint32_t> bucket_off(nb + 1U, 0U), tags(n_kept), vals(n_kept);
    uint32_t dst = 0;
    for (uint32_t i = 0; i < n_raw; i++) {
        uint8_t act = min_keep[i];
        if (act == 0U) continue;
        uint64_t h48 = raw_min[i].first;
        uint32_t b = (uint32_t)(h48 >> shift_bits);
        uint32_t tag = ((uint32_t)h48) & tag_mask;
        bucket_off[b + 1U]++;
        if (act == 2U) {
            tags[dst] = tag | GECKO_MIN_OVF_SENTINEL;
            vals[dst] = GECKO_MIN_OVF_VAL_SENTINEL;
        } else {
            tags[dst] = tag;
            vals[dst] = raw_min[i].second;
        }
        dst++;
    }
    for (uint32_t b = 0; b < nb; b++) bucket_off[b + 1U] += bucket_off[b];

    // Query tiled with exact tile boundaries & CTA wave window eval
    auto expand_emit = [&](bool is_rev, uint32_t p_q, uint32_t p_r, uint32_t r_lo, uint32_t r_hi) {
        uint32_t need_l = p_q - r_lo;
        uint32_t need_r = (r_hi + (uint32_t)wl) - (p_q + (uint32_t)kp);
        uint32_t ext_l = 0, ext_r = 0;
        if (!is_rev) {
            gecko_min_expand_extents<false>(
                seqX.bits.data(), seqX.inv_bm.data(), seqX.len,
                seqY.bits.data(), seqY.inv_bm.data(), seqY.len,
                kp, p_q, p_r, need_l, need_r, ext_l, ext_r);
        } else {
            gecko_min_expand_extents<true>(
                seqX.bits.data(), seqX.inv_bm.data(), seqX.len,
                seqY.bits.data(), seqY.inv_bm.data(), seqY.len,
                kp, p_q, p_r, need_l, need_r, ext_l, ext_r);
        }
        if (ext_l + (uint32_t)kp + ext_r < (uint32_t)wl) {
            res.exact_k_reject++;
            return;
        }
        uint32_t m_lo = p_q - ext_l;
        uint32_t m_hi = p_q + (uint32_t)kp + ext_r - (uint32_t)wl;
        uint32_t e_lo = std::max(m_lo, r_lo);
        uint32_t e_hi = std::min(m_hi, r_hi);
        if (e_lo > e_hi) {
            res.exact_k_reject++;
            return;
        }
        res.exact_k_accept++;

        uint32_t y_min = (!is_rev) ? (p_r - (p_q - e_lo)) : (p_r + (p_q - e_hi) - H);
        uint32_t y_max = y_min + (e_hi - e_lo);
        uint32_t base_y = (y_min > 0U) ? (y_min - 1U) : 0U;
        uint32_t dbits = gecko_load_32_1bit(res.dropped_bm.data(), base_y);
        uint32_t dmask = (1U << (y_max + 2U - base_y)) - 1U;
        if ((dbits & dmask) == 0U) {
            if (e_lo == m_lo) {
                uint32_t ry = (!is_rev) ? (p_r - (p_q - e_lo)) : (p_r + (p_q - e_lo) - H);
                uint32_t oy = (!is_rev) ? ry : (uint32_t)(seqY.len - ry - wl);
                res.hits.push_back({is_rev ? 1 : 0, e_lo, oy});
            }
            if (e_hi == m_hi && (e_hi != e_lo || e_lo != m_lo)) {
                uint32_t ry = (!is_rev) ? (p_r - (p_q - e_hi)) : (p_r + (p_q - e_hi) - H);
                uint32_t oy = (!is_rev) ? ry : (uint32_t)(seqY.len - ry - wl);
                res.hits.push_back({is_rev ? 1 : 0, e_hi, oy});
            }
        } else {
            for (uint32_t px = e_lo; px <= e_hi; px++) {
                uint32_t ry = (!is_rev) ? (p_r - (p_q - px)) : (p_r + (p_q - px) - H);
                if (((dbits >> (ry - base_y)) & 1U) != 0U) continue;
                bool left_ext = (px > m_lo);
                if (left_ext) {
                    uint32_t ly = (!is_rev) ? (ry - 1U) : (ry + 1U);
                    left_ext = (((dbits >> (ly - base_y)) & 1U) == 0U);
                }
                bool right_ext = (px < m_hi);
                if (right_ext) {
                    uint32_t rry = (!is_rev) ? (ry + 1U) : (ry - 1U);
                    right_ext = (((dbits >> (rry - base_y)) & 1U) == 0U);
                }
                if (!left_ext || !right_ext) {
                    uint32_t oy = (!is_rev) ? ry : (uint32_t)(seqY.len - ry - wl);
                    res.hits.push_back({is_rev ? 1 : 0, px, oy});
                }
            }
        }
    };

    auto probe_ovf = [&](uint32_t r_lo, uint32_t r_hi) {
        if (ovf.off.empty()) return;
        for (uint32_t px = r_lo; px <= r_hi; px++) {
            uint64_t w_lsb = gecko_load_32_2bit(seqX.bits.data(), px);
            uint64_t canon; uint32_t xst; bool is_pal;
            host_kmer_canon(seqX, px, wl, span_mask, canon, xst, is_pal);
            uint64_t h = gecko_fmix64(canon);
            uint32_t ob = (uint32_t)(h >> ovf.shift);
            uint32_t otag = (uint32_t)(h >> ovf.tshift);
            uint32_t lo = ovf.off[ob], hi = ovf.off[ob + 1U];
            uint64_t x_seed = w_lsb & span_mask;
            for (uint32_t j = lo; j < hi; j++) {
                if (ovf.tags[j] > otag) break;
                if (ovf.tags[j] < otag) continue;
                uint32_t ps = ovf.vals[j];
                uint32_t py = ps >> 1, yst = ps & 1U;
                uint64_t yw = gecko_load_32_2bit(seqY.bits.data(), py);
                if ((is_pal || xst == yst) && (((yw ^ x_seed) & span_mask) == 0ULL)) {
                    bool left = (px > 0 && py > 0 &&
                                 !gecko_is_bit_set(seqX.inv_bm.data(), px - 1) &&
                                 !gecko_is_bit_set(seqY.inv_bm.data(), py - 1) &&
                                 host_base_2bit(seqX, px - 1) == host_base_2bit(seqY, py - 1) &&
                                 !gecko_is_bit_set(res.dropped_bm.data(), py - 1));
                    bool right = left && (px + (uint32_t)wl < seqX.len && py + (uint32_t)wl < seqY.len &&
                                  !gecko_is_bit_set(seqX.inv_bm.data(), px + wl) &&
                                  !gecko_is_bit_set(seqY.inv_bm.data(), py + wl) &&
                                  host_base_2bit(seqX, px + wl) == host_base_2bit(seqY, py + wl) &&
                                  !gecko_is_bit_set(res.dropped_bm.data(), py + 1));
                    if (!left || !right) res.hits.push_back({0, px, py});
                }
                if (compare_reverse && (is_pal || xst != yst)) {
                    uint64_t y_rc = (~gecko_reverse_32_2bit(yw)) >> rev_rshift;
                    if (((y_rc ^ x_seed) & span_mask) == 0ULL) {
                        uint32_t py_rc = (uint32_t)(seqY.len - py - wl);
                        bool left = (px > 0 && py_rc > 0 &&
                                     !gecko_is_bit_set(seqX.inv_bm.data(), px - 1) &&
                                     !gecko_is_bit_set(seqY.inv_bm.data(), py + wl) &&
                                     (host_base_2bit(seqY, py + wl) ^ 3U) == host_base_2bit(seqX, px - 1) &&
                                     !gecko_is_bit_set(res.dropped_bm.data(), py + 1));
                        bool right = left && (px + (uint32_t)wl < seqX.len && py > 0 &&
                                      !gecko_is_bit_set(seqX.inv_bm.data(), px + wl) &&
                                      !gecko_is_bit_set(seqY.inv_bm.data(), py - 1) &&
                                      (host_base_2bit(seqY, py - 1) ^ 3U) == host_base_2bit(seqX, px + wl) &&
                                      !gecko_is_bit_set(res.dropped_bm.data(), py - 1));
                        if (!left || !right) res.hits.push_back({1, px, py_rc});
                    }
                }
            }
        }
    };

    uint64_t total_kp_x = (seqX.len >= (uint64_t)kp) ? (seqX.len - (uint64_t)kp + 1ULL) : 0ULL;
    for (uint64_t ts = 0; ts < seqX.len; ts += tile_bp) {
        uint64_t te = std::min<uint64_t>(ts + tile_bp, seqX.len);
        if (te < (uint64_t)wl) continue;
        uint32_t tmin = (ts >= (uint64_t)(wl - 1)) ? (uint32_t)(ts - (uint64_t)(wl - 1)) : 0U;
        uint32_t tmax = (uint32_t)(te - (uint64_t)wl);
        if (tmin > tmax) continue;
        uint32_t p_begin = tmin, p_end = tmax + (uint32_t)w;
        for (uint32_t ws = p_begin; ws < p_end; ws += 512U) {
            uint32_t wlen = std::min<uint32_t>(512U, p_end - ws);
            std::vector<uint64_t> sc(wlen + 2U * H, ~0ULL);
            for (uint32_t i = 0; i < wlen + 2U * H; i++) {
                int64_t gp = (int64_t)ws - (int64_t)H + (int64_t)i;
                if (gp >= 0 && (uint64_t)gp < total_kp_x) sc[i] = host_kp_score(seqX, (uint32_t)gp, kp);
            }
            for (uint32_t li = 0; li < wlen; li++) {
                int L = 0, R = 0;
                if (!gecko_min_window_eval(sc.data(), (int)(li + H), w, false, L, R)) continue;
                uint32_t p_q = ws + li;
                uint32_t r_lo = std::max(tmin, p_q - (uint32_t)L);
                uint32_t r_hi = std::min(tmax, p_q + (uint32_t)R - H);
                if (r_lo > r_hi) continue;

                uint64_t canon; uint32_t qst, qpal;
                host_kp_canon(seqX, p_q, kp, canon, qst, qpal);
                uint64_t h48 = gecko_mix48(canon);
                uint32_t bkt = (uint32_t)(h48 >> shift_bits);
                uint32_t tag = ((uint32_t)h48) & tag_mask;
                uint32_t lo = bucket_off[bkt], hi = bucket_off[bkt + 1U];
                if (hi - lo > 8U) {
                    uint32_t l = lo, r = hi;
                    while (l < r) {
                        uint32_t mid = l + ((r - l) >> 1);
                        if ((tags[mid] & GECKO_MIN_TAG_MASK) < tag) l = mid + 1U;
                        else r = mid;
                    }
                    lo = l;
                }
                for (uint32_t j = lo; j < hi; j++) {
                    uint32_t raw_t = tags[j];
                    uint32_t t = raw_t & GECKO_MIN_TAG_MASK;
                    if (t > tag) break;
                    if (t < tag) continue;
                    if ((raw_t & GECKO_MIN_OVF_SENTINEL) != 0U) {
                        assert(vals[j] == GECKO_MIN_OVF_VAL_SENTINEL);
                        probe_ovf(r_lo, r_hi);
                        break;
                    }
                    uint32_t p_r = vals[j] >> 1, rst = vals[j] & 1U;
                    if (qpal || qst == rst) expand_emit(false, p_q, p_r, r_lo, r_hi);
                    if (compare_reverse && (qpal || qst != rst)) expand_emit(true, p_q, p_r, r_lo, r_hi);
                }
            }
        }
    }
    std::sort(res.hits.begin(), res.hits.end());
    return res;
}

// ---- P1: Sliding-Window Canonical Minimizer vs Exact Monotone Deque ----
static void test_p1_sliding_window() {
    std::string seq_str = rand_dna(5200, 0x1111ULL);
    seq_str.replace(100, 120, std::string(120, 'A'));
    seq_str.replace(300, 120, std::string(120, 'T'));
    for (int i = 0; i < 60; i++) { seq_str[500 + 2*i] = 'A'; seq_str[501 + 2*i] = 'T'; }
    for (int i = 0; i < 30; i++) {
        seq_str[700 + 4*i] = 'A'; seq_str[701 + 4*i] = 'C';
        seq_str[702 + 4*i] = 'G'; seq_str[703 + 4*i] = 'T';
    }
    // CTA wave/chunk boundary straddles at 512, 1024, 4096
    seq_str.replace(505, 16, "GGGGGGGGGGGGGGGG");
    seq_str.replace(1018, 14, "ACACACACACACAC");
    seq_str.replace(4088, 20, "TTTTTTTTTTTTTTTTTTTT");
    // N-runs, soft-masked lowercase, SEQ_SEP, short valid runs (len = wl-1, len = wl)
    seq_str.replace(1500, 10, "NNNNNNNNNN");
    for (size_t i = 1600; i < 1640; i++) seq_str[i] = (char)(seq_str[i] - 'A' + 'a');
    seq_str[1800] = '\0';
    seq_str[1800 + 18] = 'N'; // 17-bp run (< 18)
    seq_str[1850] = 'N';
    seq_str[1850 + 33] = 'N'; // 32-bp exact run

    for (bool soft_mask : {false, true}) {
        HostSeq seq = HostSeq::from_ascii(seq_str, soft_mask);
        for (int wl : {18, 20, 24, 28, 32}) {
            int kp = 0, w = 0;
            assert(gecko_minimizer_params(wl, kp, w));
            uint32_t H = (uint32_t)(w - 1);
            uint64_t total_kp = seq.len - (uint64_t)kp + 1ULL;

            for (bool emit_ties : {true, false}) {
                // Reference monotone deque (matching gecko_core::minimizer::for_each_minimizer)
                std::vector<uint32_t> deque_pos;
                size_t rs = 0;
                while (rs < seq.len) {
                    if (gecko_is_bit_set(seq.inv_bm.data(), (uint32_t)rs)) { rs++; continue; }
                    size_t re = rs;
                    while (re < seq.len && !gecko_is_bit_set(seq.inv_bm.data(), (uint32_t)re)) re++;
                    if (re - rs >= (size_t)wl) {
                        std::deque<std::pair<uint64_t, uint32_t>> dq;
                        uint32_t last_emit = 0xFFFFFFFFU;
                        for (uint32_t pos = (uint32_t)rs; pos + (uint32_t)kp <= (uint32_t)re; pos++) {
                            uint64_t sc = host_kp_score(seq, pos, kp);
                            while (!dq.empty() && dq.back().first > sc) dq.pop_back();
                            dq.push_back({sc, pos});
                            if (pos >= (uint32_t)(rs + w - 1)) {
                                uint32_t win_s = pos + 1U - (uint32_t)w;
                                while (!dq.empty() && dq.front().second < win_s) dq.pop_front();
                                if (emit_ties) {
                                    uint64_t ms = dq.front().first;
                                    for (auto &e : dq) {
                                        if (e.first != ms) break;
                                        if (last_emit == 0xFFFFFFFFU || e.second > last_emit) {
                                            deque_pos.push_back(e.second);
                                            last_emit = e.second;
                                        }
                                    }
                                } else if (dq.front().second != last_emit) {
                                    deque_pos.push_back(dq.front().second);
                                    last_emit = dq.front().second;
                                }
                            }
                        }
                    }
                    rs = re;
                }

                // CTA 512-bp wave + halo H using gecko_min_window_eval
                std::vector<uint32_t> cta_pos;
                for (uint64_t ws = 0; ws < total_kp; ws += 512ULL) {
                    uint32_t wlen = (uint32_t)std::min<uint64_t>(512ULL, total_kp - ws);
                    std::vector<uint64_t> sc(wlen + 2U * H, ~0ULL);
                    for (uint32_t i = 0; i < wlen + 2U * H; i++) {
                        int64_t gp = (int64_t)ws - (int64_t)H + (int64_t)i;
                        if (gp >= 0 && (uint64_t)gp < total_kp) sc[i] = host_kp_score(seq, (uint32_t)gp, kp);
                    }
                    for (uint32_t p = 0; p < wlen; p++) {
                        int L = 0, R = 0;
                        if (gecko_min_window_eval(sc.data(), (int)(p + H), w, emit_ties, L, R)) {
                            cta_pos.push_back((uint32_t)(ws + p));
                        }
                    }
                }
                assert(cta_pos == deque_pos);
            }
        }
    }
}

// ---- P2: mix48 Bijection, 0% Tag Collisions & Bucket Uniformity ----
static void test_p2_mix48_bijection() {
    uint64_t rng = 0xDEADBEEFCAFEBABEULL;
    for (int i = 0; i < 100000; i++) {
        rng = rng * 6364136223846793005ULL + 1442695040888963407ULL;
        uint64_t x = rng & GECKO_MASK48;
        uint64_t h = gecko_mix48(x);
        assert((h & ~GECKO_MASK48) == 0ULL);
        assert(gecko_mix48_inv(h) == x);
    }
    // Verify bucket_bits in [17..26] gives bit 31 == 0 and 0.000% (bkt, tag) collisions
    for (int bb = 17; bb <= 26; bb++) {
        int sb = 48 - bb;
        uint32_t tmask = (1U << sb) - 1U;
        assert((tmask & GECKO_MIN_OVF_SENTINEL) == 0U);
        std::unordered_map<uint64_t, uint64_t> seen;
        for (uint64_t i = 0; i < 5000; i++) {
            uint64_t c = (i * 0x1234567ULL) & 0x3FFFFFFFFFFULL;
            uint64_t h48 = gecko_mix48(c);
            uint32_t bkt = (uint32_t)(h48 >> sb);
            uint32_t tag = ((uint32_t)h48) & tmask;
            assert((tag & GECKO_MIN_OVF_SENTINEL) == 0U);
            uint64_t recon = ((uint64_t)bkt << sb) | (uint64_t)tag;
            assert(recon == h48);
            auto it = seen.find(recon);
            if (it != seen.end()) assert(it->second == c);
            seen[recon] = c;
        }
    }
    // Bucket uniformity on minimum-skewed real minimizer keys (256 buckets, max/mean < 2.0)
    HostSeq seq = HostSeq::from_ascii(rand_dna(120000, 0x9999ULL), false);
    std::vector<uint32_t> bcnt(256, 0U);
    uint64_t total_kp = seq.len - 20ULL;
    uint32_t n_emitted = 0;
    for (uint64_t ws = 0; ws < total_kp; ws += 512ULL) {
        uint32_t wlen = (uint32_t)std::min<uint64_t>(512ULL, total_kp - ws);
        std::vector<uint64_t> sc(wlen + 22U, ~0ULL);
        for (uint32_t i = 0; i < wlen + 22U; i++) {
            int64_t gp = (int64_t)ws - 11 + (int64_t)i;
            if (gp >= 0 && (uint64_t)gp < total_kp) sc[i] = host_kp_score(seq, (uint32_t)gp, 21);
        }
        for (uint32_t p = 0; p < wlen; p++) {
            int L = 0, R = 0;
            if (gecko_min_window_eval(sc.data(), (int)(p + 11U), 12, false, L, R)) {
                uint64_t c; uint32_t st, pal;
                host_kp_canon(seq, (uint32_t)(ws + p), 21, c, st, pal);
                bcnt[(gecko_mix48(c) >> 40) & 255U]++;
                n_emitted++;
            }
        }
    }
    uint32_t mx = *std::max_element(bcnt.begin(), bcnt.end());
    double mean = (double)n_emitted / 256.0;
    assert(mean > 20.0 && (double)mx / mean < 2.0);
}

// ---- P3: Exact Per-Strand freq_limit Resolution (palindromes + asymmetric strands + hash ties) ----
static void test_p3_exact_freq_limit() {
    std::string ref = rand_dna(600, 0x3333ULL);
    std::string k_asym = rand_dna(32, 0x3344ULL);
    std::string k_asym_rc = make_revcomp_str(k_asym);
    // 2 fwd + 3 rev -> at FL=2: fwd kept, rev dropped!
    for (int i = 0; i < 2; i++) ref += "N" + k_asym;
    for (int i = 0; i < 3; i++) ref += "N" + k_asym_rc;
    // Palindromic 32-mer (wl=32, kp=21: mirrored tie inside!) & palindromic 14-mer (wl=18, kp=14)
    std::string h16 = "ACGTACGTTAGCGCTA";
    std::string pal32 = h16 + make_revcomp_str(h16);
    std::string h7 = "ACGTACG";
    std::string pal18 = "GG" + h7 + make_revcomp_str(h7) + "CC";
    for (int i = 0; i < 3; i++) ref += "N" + pal32 + "N" + pal18;
    // Periodic tie block
    std::string per;
    for (int i = 0; i < 16; i++) per += "ACGT";
    ref += "N" + per + "N" + per + "N" + per;

    std::string qry = k_asym + "N" + k_asym_rc + "N" + pal32 + "N" + pal18 + "N" + per +
                      "N" + make_revcomp_str(ref.substr(50, 120));

    HostSeq sx = HostSeq::from_ascii(qry, true);
    HostSeq sy = HostSeq::from_ascii(ref, true);
    for (int wl : {18, 24, 32}) {
        for (uint32_t fl : {1U, 2U, 20U, 100000U}) {
            for (bool rev : {true, false}) {
                auto o = run_exact_stride1_oracle(sx, sy, wl, fl, rev);
                auto m = run_gpu_min_pipeline_sim(sx, sy, wl, fl, rev, 128ULL);
                assert(m.dropped_bm == o.dropped_bm);
                assert(m.hits == o.hits);
            }
        }
    }
}

// ---- P4: High-Freq Repeat Minimizer Exact-k Overflow Routing (M_CAP) ----
static void test_p4_hf_overflow_routing() {
    // 500x copies of a shared 21-bp core (fwd + rev) with unique flanks + 20x K_keep + 21x K_drop
    // Choose core with minimal fmix64 score so it is guaranteed the minimizer!
    uint64_t best_sc = ~0ULL;
    std::string core(21, 'A');
    for (int t = 0; t < 2000; t++) {
        std::string cand = rand_dna(21, 0x7000ULL + (uint64_t)t);
        HostSeq hs = HostSeq::from_ascii(cand, false);
        uint64_t sc = host_kp_score(hs, 0, 21);
        if (sc < best_sc) { best_sc = sc; core = cand; }
    }
    std::string core_rc = make_revcomp_str(core);

    std::string ref;
    std::vector<std::string> unique_32;
    for (int i = 0; i < 250; i++) {
        std::string lf = rand_dna(6, 0xA000ULL + (uint64_t)i * 2ULL);
        std::string rf = rand_dna(5, 0xA001ULL + (uint64_t)i * 2ULL);
        std::string k32 = lf + core + rf;
        unique_32.push_back(k32);
        ref += "N" + k32;
    }
    for (int i = 0; i < 250; i++) {
        std::string lf = rand_dna(6, 0xB000ULL + (uint64_t)i * 2ULL);
        std::string rf = rand_dna(5, 0xB001ULL + (uint64_t)i * 2ULL);
        std::string k32 = lf + core_rc + rf;
        unique_32.push_back(k32);
        ref += "N" + k32;
    }
    std::string k_keep = "ACGTAC" + core + "TGCAT"; // 20x == FL -> KEPT
    std::string k_drop = "TGCAAC" + core + "ACGTA"; // 21x == FL+1 -> DROPPED
    for (int i = 0; i < 20; i++) ref += "N" + k_keep;
    for (int i = 0; i < 21; i++) ref += "N" + k_drop;

    // Query tests all 500 unique copies (fwd & rev), K_keep (fwd & rev), and K_drop (fwd & rev)
    std::string qry = k_keep + "N" + make_revcomp_str(k_keep) + "N" +
                      k_drop + "N" + make_revcomp_str(k_drop);
    for (int i = 0; i < 500; i += 5) {
        qry += "N" + unique_32[i] + "N" + make_revcomp_str(unique_32[i]);
    }

    HostSeq sx = HostSeq::from_ascii(qry, false);
    HostSeq sy = HostSeq::from_ascii(ref, false);
    auto o = run_exact_stride1_oracle(sx, sy, 32, 20U, true);
    auto m = run_gpu_min_pipeline_sim(sx, sy, 32, 20U, true, 1024ULL);

    assert(m.ovf_routed_keys >= 1U);     // HF minimizer bucket MUST route via sentinel!
    assert(m.ovf_kmers >= 500U);         // Non-dropped exact k-mers routed to overflow CSR
    assert(m.exact_k_reject < 5000ULL);  // Zero 500x500 = 250,000 quadratic FP explosion!
    assert(m.dropped_bm == o.dropped_bm);
    assert(m.hits == o.hits);
    assert(!o.hits.empty());
}

// ---- P5: Exact-k 2-Bit Expand + L/R-Maximal Filter + Reverse-Strand Dropped-Gap Path ----
static void test_p5_expand_lr_max_and_reverse_dropped_gaps() {
    std::string ref = rand_dna(300, 0x5500ULL);
    std::string qry = rand_dna(100, 0x5501ULL);

    // Exact MEMs of lengths 31, 32, 33, 40, 44, 48, 49, 64, 200 on BOTH forward & reverse strands
    for (int L : {31, 32, 33, 40, 44, 48, 49, 64, 200}) {
        std::string blk = rand_dna((size_t)L, 0x5600ULL + (uint64_t)L);
        ref += "A" + blk + "C" + "N" + "G" + make_revcomp_str(blk) + "T" + "N";
        qry += "T" + blk + "G" + "N" + "C" + blk + "A" + "N";
    }

    // Reverse-strand & Forward-strand MEMs interrupted by INTERIOR dropped 32-mers (FL = 2)!
    // Specifically exercises min_expand_and_emit<true> dropped-gap path (ly = ry + 1, rry = ry - 1)!
    std::string gap_mem = rand_dna(120, 0x5777ULL);
    std::string dropped_mid1 = gap_mem.substr(40, 32);
    std::string dropped_mid2 = gap_mem.substr(55, 32);
    ref += "N" + gap_mem + "N" + make_revcomp_str(gap_mem) + "N";
    // Add 3 extra copies of dropped_mid1 & dropped_mid2 (so count > FL=2 -> dropped_bitmap=1!)
    for (int i = 0; i < 3; i++) {
        ref += "N" + dropped_mid1 + "N" + make_revcomp_str(dropped_mid1);
        ref += "N" + dropped_mid2 + "N" + make_revcomp_str(dropped_mid2);
    }
    qry += "N" + gap_mem + "N" + make_revcomp_str(gap_mem) + "N";

    // Soft-masked flanks & SEQ_SEP boundary directly adjacent to MEM
    std::string sm_mem = rand_dna(60, 0x5888ULL);
    ref += "acgt" + sm_mem + "tgca";
    ref.push_back('\0');
    ref += sm_mem;
    qry += "acgt" + sm_mem + "tgca";
    qry.push_back('\0');
    qry += sm_mem;

    for (bool soft_mask : {true, false}) {
        HostSeq sx = HostSeq::from_ascii(qry, soft_mask);
        HostSeq sy = HostSeq::from_ascii(ref, soft_mask);
        auto o = run_exact_stride1_oracle(sx, sy, 32, 2U, true);
        // Verify small tile sizes (64, 100, 1024) have ZERO boundary misses/duplicates!
        for (uint64_t tile_bp : {64ULL, 100ULL, 1024ULL}) {
            auto m = run_gpu_min_pipeline_sim(sx, sy, 32, 2U, true, tile_bp);
            assert(m.dropped_bm == o.dropped_bm);
            assert(m.hits == o.hits);
        }
    }
}

// ---- P7: two_hit_filter Parity vs CPU TwoHitTable (contiguous + spaced + long MEMs + tile overlap) ----
static std::vector<uint64_t> cpu_two_hit_reference(
    const std::vector<uint64_t> &sorted_hits, HitKeyLayout lay, uint64_t wl, uint64_t window
) {
    // Mirrors gecko_core::stream::DiagSlot32::two_hit_check exactly
    std::unordered_map<uint64_t, uint64_t> last_x;
    std::vector<uint64_t> kept;
    for (uint64_t key : sorted_hits) {
        uint64_t sd = key >> lay.x_bits;
        uint64_t rx = key & (uint64_t)lay.x_mask;
        auto it = last_x.find(sd);
        if (it == last_x.end()) {
            last_x[sd] = rx;
            continue;
        }
        uint64_t prev = it->second;
        if (rx <= prev) continue;
        uint64_t dist = rx - prev;
        if (dist < wl) continue; // overlapping -> keep anchor unchanged!
        it->second = rx;
        if (dist - wl <= window) kept.push_back(key);
    }
    return kept;
}

static void test_p7_two_hit_parity(const HostSeq *real_x, const HostSeq *real_y) {
    // 1. Exact synthetic boundary matrix: wl in {10, 12, 16, 32} and spaced w14s20 (wl=20)
    for (int wl : {10, 12, 16, 20, 32}) {
        for (uint32_t W : {1U, 30U, 31U, 32U, 40U}) {
            HitKeyLayout lay = make_hit_key_layout(0, 5000, 5000, wl, true);
            // (a) Isolated single hit (n_hits == 1) -> MUST drop!
            uint64_t single = gecko_pack_hit_key_tile(lay, 0, 500, 500);
            assert(gecko_two_hit_scan_backward(&single, 0, lay, (uint64_t)wl, (uint64_t)wl + W) == 0U);

            // (b) Exact boundary tests on both strands:
            for (uint32_t rev : {0U, 1U}) {
                // dist == wl - 1 (overlapping) -> DROP
                uint64_t h_ov[2] = {
                    gecko_pack_hit_key_tile(lay, rev, 200, 200),
                    gecko_pack_hit_key_tile(lay, rev, 200 + (uint32_t)wl - 1U, 200 + (uint32_t)wl - 1U)
                };
                assert(gecko_two_hit_scan_backward(h_ov, 1, lay, wl, (uint64_t)wl + W) == 0U);

                // dist == wl (gap 0) -> KEEP
                uint64_t h_eq[2] = {
                    gecko_pack_hit_key_tile(lay, rev, 200, 200),
                    gecko_pack_hit_key_tile(lay, rev, 200 + (uint32_t)wl, 200 + (uint32_t)wl)
                };
                assert(gecko_two_hit_scan_backward(h_eq, 1, lay, wl, (uint64_t)wl + W) == 1U);

                // dist == wl + W (exact inclusive boundary) -> KEEP
                uint64_t h_max[2] = {
                    gecko_pack_hit_key_tile(lay, rev, 200, 200),
                    gecko_pack_hit_key_tile(lay, rev, 200 + (uint32_t)wl + W, 200 + (uint32_t)wl + W)
                };
                assert(gecko_two_hit_scan_backward(h_max, 1, lay, wl, (uint64_t)wl + W) == 1U);

                // dist == wl + W + 1 (1 bp past boundary) -> DROP
                uint64_t h_past[2] = {
                    gecko_pack_hit_key_tile(lay, rev, 200, 200),
                    gecko_pack_hit_key_tile(lay, rev, 200 + (uint32_t)wl + W + 1U, 200 + (uint32_t)wl + W + 1U)
                };
                assert(gecko_two_hit_scan_backward(h_past, 1, lay, wl, (uint64_t)wl + W) == 0U);

                // (c) Overlapping triple [x1, x2, x3] where x2 - x1 < wl, wl <= x3 - x1 <= wl + W -> x3 KEPT
                uint64_t h_tri[3] = {
                    gecko_pack_hit_key_tile(lay, rev, 200, 200),
                    gecko_pack_hit_key_tile(lay, rev, 200 + (uint32_t)wl / 2U, 200 + (uint32_t)wl / 2U),
                    gecko_pack_hit_key_tile(lay, rev, 200 + (uint32_t)wl + (W / 2U), 200 + (uint32_t)wl + (W / 2U))
                };
                assert(gecko_two_hit_scan_backward(h_tri, 0, lay, wl, (uint64_t)wl + W) == 0U);
                assert(gecko_two_hit_scan_backward(h_tri, 1, lay, wl, (uint64_t)wl + W) == 0U);
                assert(gecko_two_hit_scan_backward(h_tri, 2, lay, wl, (uint64_t)wl + W) == 1U);
            }
        }
    }

    // 2. 100-bp contiguous exact match (wl=12, W=30) on BOTH forward & reverse strands!
    // L/R-maximal emits ONLY endpoints [1000, 1088] (dist=88 > 42); gecko_two_hit_self_prev MUST keep 1088!
    {
        std::string mem100 = rand_dna(100, 0x7777ULL);
        std::string ref = rand_dna(500, 0x7701ULL) + mem100 + rand_dna(500, 0x7702ULL) +
                          make_revcomp_str(mem100) + rand_dna(500, 0x7703ULL);
        std::string qry = rand_dna(1000, 0x7704ULL) + mem100 + rand_dna(500, 0x7705ULL);
        HostSeq sx = HostSeq::from_ascii(qry, false);
        HostSeq sy = HostSeq::from_ascii(ref, false);
        HostSeq sy_rc = sy.revcomp();
        HitKeyLayout lay = make_hit_key_layout(0, sx.len, sy.len, 12, true);
        uint32_t rev_y_rc_start = (uint32_t)(sy.len - (1100ULL + 100ULL));

        // Forward endpoints: (1000, 500) and (1088, 588)
        uint64_t fwd_hits[2] = {
            gecko_pack_hit_key_tile(lay, 0, 1000, 500),
            gecko_pack_hit_key_tile(lay, 0, 1088, 588)
        };
        assert(gecko_two_hit_scan_backward(fwd_hits, 1, lay, 12, 42) == 0U);
        assert(gecko_two_hit_self_prev(fwd_hits[0], lay, 12, sx.bits.data(), sx.inv_bm.data(),
               sy.bits.data(), sy.inv_bm.data(), sy_rc.bits.data(), nullptr, sy.len) == 0U);
        assert(gecko_two_hit_self_prev(fwd_hits[1], lay, 12, sx.bits.data(), sx.inv_bm.data(),
               sy.bits.data(), sy.inv_bm.data(), sy_rc.bits.data(), nullptr, sy.len) == 1U);

        // Reverse endpoints: (1000, rev_y_rc_start) and (1088, rev_y_rc_start + 88)
        uint64_t rev_hits[2] = {
            gecko_pack_hit_key_tile(lay, 1, 1000, rev_y_rc_start),
            gecko_pack_hit_key_tile(lay, 1, 1088, rev_y_rc_start + 88U)
        };
        assert(gecko_two_hit_scan_backward(rev_hits, 1, lay, 12, 42) == 0U);
        assert(gecko_two_hit_self_prev(rev_hits[0], lay, 12, sx.bits.data(), sx.inv_bm.data(),
               sy.bits.data(), sy.inv_bm.data(), sy_rc.bits.data(), nullptr, sy.len) == 0U);
        assert(gecko_two_hit_self_prev(rev_hits[1], lay, 12, sx.bits.data(), sx.inv_bm.data(),
               sy.bits.data(), sy.inv_bm.data(), sy_rc.bits.data(), nullptr, sy.len) == 1U);
    }

    // 3. Tile-boundary straddling gapped two-hit pair (tile_end - 25, tile_end + 10) KEPT via overlap!
    {
        uint64_t lenX = 2048, init_tile = 1024, wl = 12, W = 40;
        uint64_t overlap = std::min<uint64_t>(wl + W, init_tile / 4ULL); // 52
        uint64_t step = init_tile - overlap; // 972
        // Pair at x1 = 1024 - 25 = 999, x2 = 1024 + 10 = 1034 (dist = 35 in [12, 52])
        uint32_t x1 = 999, x2 = 1034;
        bool kept_in_any_tile = false;
        for (uint64_t ts = 0; ts < lenX; ) {
            uint64_t te = std::min<uint64_t>(ts + init_tile, lenX);
            HitKeyLayout lay = make_hit_key_layout(ts, te, 2048, (int)wl, false);
            std::vector<uint64_t> th;
            uint32_t tmin = (ts >= wl - 1) ? (uint32_t)(ts - (wl - 1)) : 0U;
            uint32_t tmax = (uint32_t)(te - wl);
            if (x1 >= tmin && x1 <= tmax) th.push_back(gecko_pack_hit_key_tile(lay, 0, x1, x1));
            if (x2 >= tmin && x2 <= tmax) th.push_back(gecko_pack_hit_key_tile(lay, 0, x2, x2));
            for (uint32_t i = 0; i < th.size(); i++) {
                if (gecko_two_hit_scan_backward(th.data(), i, lay, wl, wl + W)) kept_in_any_tile = true;
            }
            if (te == lenX) break;
            ts += step;
        }
        assert(kept_in_any_tile);
    }

    // 4. Real chimp_chr21 x gorilla_chr21 k=12 & spaced w14s20 Keep(CPU TwoHitTable) subseteq Keep(GPU)
    if (real_x && real_y) {
        struct SeedSpec { int wl; uint64_t mask; };
        // w14s20 pattern 11101101011010110111 -> 2-bit mask
        const char *pat = "11101101011010110111";
        uint64_t w14s20_mask = 0ULL;
        for (int i = 0; i < 20; i++) {
            if (pat[i] == '1') w14s20_mask |= (3ULL << (2 * i));
        }
        SeedSpec specs[2] = {{12, (1ULL << 24) - 1ULL}, {20, w14s20_mask}};
        for (auto sp : specs) {
            uint64_t slice_len = std::min<uint64_t>(real_x->len, 32768ULL);
            uint32_t vmask = (1U << sp.wl) - 1U;
            std::unordered_map<uint64_t, std::vector<std::pair<uint32_t, uint32_t>>> ref_idx;
            for (uint32_t py = 0; py + (uint32_t)sp.wl <= (uint32_t)slice_len; py++) {
                if ((gecko_load_32_1bit(real_y->inv_bm.data(), py) & vmask) != 0U) continue;
                uint64_t c; uint32_t st; bool pal;
                host_kmer_canon(*real_y, py, sp.wl, sp.mask, c, st, pal);
                ref_idx[c].push_back({py, st});
            }
            HitKeyLayout lay = make_hit_key_layout(0, slice_len, slice_len, sp.wl, true);
            std::vector<uint64_t> all_hits;
            for (uint32_t px = 0; px + (uint32_t)sp.wl <= (uint32_t)slice_len; px++) {
                if ((gecko_load_32_1bit(real_x->inv_bm.data(), px) & vmask) != 0U) continue;
                uint64_t c; uint32_t xst; bool pal;
                host_kmer_canon(*real_x, px, sp.wl, sp.mask, c, xst, pal);
                auto it = ref_idx.find(c);
                if (it == ref_idx.end() || it->second.size() > 20U) continue;
                for (auto &ps : it->second) {
                    if (pal || xst == ps.second) all_hits.push_back(gecko_pack_hit_key_tile(lay, 0, px, ps.first));
                    if (pal || xst != ps.second) {
                        uint32_t py_rc = (uint32_t)(slice_len - ps.first - sp.wl);
                        all_hits.push_back(gecko_pack_hit_key_tile(lay, 1, px, py_rc));
                    }
                }
            }
            std::sort(all_hits.begin(), all_hits.end());
            auto cpu_kept = cpu_two_hit_reference(all_hits, lay, (uint64_t)sp.wl, 40ULL);
            std::vector<uint64_t> gpu_kept;
            for (uint32_t i = 0; i < all_hits.size(); i++) {
                if (gecko_two_hit_scan_backward(all_hits.data(), i, lay, (uint64_t)sp.wl, (uint64_t)sp.wl + 40ULL)) {
                    gpu_kept.push_back(all_hits[i]);
                }
            }
            // Assert Keep(CPU TwoHitTable) subseteq Keep(GPU two_hit_filter) (0 false negatives!)
            for (uint64_t k : cpu_kept) {
                assert(std::binary_search(gpu_kept.begin(), gpu_kept.end(), k));
            }
            assert(gpu_kept.size() < all_hits.size()); // confirms substantial single-hit filtering!
        }
    }
}

static std::string read_file_raw(const char *path) {
    std::ifstream ifs(path, std::ios::binary);
    assert(ifs.good());
    return std::string((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
}

int main(int argc, char **argv) {
    test_p1_sliding_window();
    test_p2_mix48_bijection();
    test_p3_exact_freq_limit();
    test_p4_hf_overflow_routing();
    test_p5_expand_lr_max_and_reverse_dropped_gaps();

    HostSeq real_x_for_p7, real_y_for_p7;
    bool have_real = false;

    // P6: Real chimp_chr21 x gorilla_chr21 slices (BOTH soft_mask=true AND soft_mask=false!)
    for (int a = 1; a + 1 < argc; a += 2) {
        std::string c_raw = read_file_raw(argv[a]);
        std::string g_raw = read_file_raw(argv[a + 1]);
        for (bool soft_mask : {true, false}) {
            HostSeq sx = HostSeq::from_ascii(c_raw, soft_mask);
            HostSeq sy = HostSeq::from_ascii(g_raw, soft_mask);
            if (!have_real && !soft_mask) {
                real_x_for_p7 = sx;
                real_y_for_p7 = sy;
                have_real = true;
            }
            auto o = run_exact_stride1_oracle(sx, sy, 32, 20U, true);
            auto m = run_gpu_min_pipeline_sim(sx, sy, 32, 20U, true, 65536ULL);
            assert(m.dropped_bm == o.dropped_bm);
            assert(m.hits == o.hits);
            // Bounded exact_k_reject under --no-soft-mask thanks to M_CAP overflow routing
            assert(m.exact_k_reject < (uint64_t)sx.len * 4ULL);
            fprintf(stderr, "[P6] slice=%d soft_mask=%d hits=%zu ovf_keys=%u ovf_kmers=%u acc=%llu rej=%llu\n",
                    (a - 1) / 2, soft_mask ? 1 : 0, o.hits.size(), m.ovf_routed_keys, m.ovf_kmers,
                    (unsigned long long)m.exact_k_accept, (unsigned long long)m.exact_k_reject);
        }
    }

    test_p7_two_hit_parity(have_real ? &real_x_for_p7 : nullptr, have_real ? &real_y_for_p7 : nullptr);
    return 0;
}
"##;
        std::fs::write(&cpp_path, cpp_src).expect("write parity_p1_p7.cpp");
        let inc_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("cuda");
        let compile = std::process::Command::new("g++")
            .args([
                "-O2",
                "-std=c++17",
                "-fsanitize=address,undefined",
                "-Wall",
                "-Wextra",
                "-Werror",
            ])
            .arg(format!("-I{}", inc_dir.display()))
            .arg(&cpp_path)
            .arg("-o")
            .arg(&bin_path)
            .output()
            .expect("run g++");
        assert!(
            compile.status.success(),
            "g++ failed to compile parity_p1_p7.cpp:\n{}",
            String::from_utf8_lossy(&compile.stderr)
        );
        let run = std::process::Command::new(&bin_path)
            .args(&real_args)
            .output()
            .expect("run parity_p1_p7");
        assert!(
            run.status.success(),
            "parity_p1_p7 failed:\nstdout:\n{}\nstderr:\n{}",
            String::from_utf8_lossy(&run.stdout),
            String::from_utf8_lossy(&run.stderr)
        );
    }
}
