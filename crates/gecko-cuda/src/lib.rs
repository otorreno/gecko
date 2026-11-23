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

/// Returns `true` if native CUDA kernels were compiled into this binary.
pub const fn is_compiled_in() -> bool {
    cfg!(has_cuda_kernels)
}

#[cfg(has_cuda_kernels)]
mod ffi {
    use std::ffi::{CStr, c_char, c_int};

    use gecko_core::error::{GeckoError, Result};
    use gecko_core::extend::CoverageMatrix;
    use gecko_core::formats::frags::Fragment;
    use gecko_core::pipeline::{CompareConfig, CompareResult, CompareStats};
    use gecko_core::seq::Genome;
    use gecko_core::stream::dedup_fragments;

    use super::CudaDeviceInfo;

    #[repr(C)]
    pub struct GeckoCudaRefHandle {
        _private: [u8; 0],
    }

    #[repr(C)]
    #[derive(Clone, Copy)]
    pub struct GpuFrag {
        pub diag: i64,
        pub x_start: u32,
        pub y_start: u32,
        pub x_end: u32,
        pub y_end: u32,
        pub ident: u32,
        pub score: u32,
        pub strand: u8,
        pub pad: [u8; 7],
    }

    #[repr(C)]
    #[derive(Default)]
    pub struct GeckoCudaResult {
        pub frags: *mut GpuFrag,
        pub n_frags: u64,
        pub leftmax_hits: u64,
        pub extended_hits: u64,
        pub index_kept: u32,
        pub index_valid: u32,
        pub bucket_bits: i32,
    }

    #[repr(C)]
    pub struct GeckoPackedSeqHost {
        pub bits: *const u64,
        pub bits_words: u64,
        pub special_bm: *const u32,
        pub invalid_bm: *const u32,
        pub bm_words: u64,
        pub sep_pos: *const u32,
        pub n_seps: u32,
        pub iupac_pos: *const u32,
        pub iupac_byte: *const u8,
        pub n_iupac: u32,
        pub len: u64,
    }

    /// Owns transient sparse bitmaps & u32 index slices for the duration of a CUDA upload call
    /// while borrowing `PackedSeq::bits` zero-copy.
    struct HostPackedView<'a> {
        bits: &'a [u64],
        special_bm: Option<Vec<u32>>,
        invalid_bm: Option<Vec<u32>>,
        alias_invalid_to_special: bool,
        sep_pos: Vec<u32>,
        iupac_pos: Vec<u32>,
        iupac_byte: Vec<u8>,
        len: u64,
    }

    impl<'a> HostPackedView<'a> {
        fn new(genome: &'a Genome, soft_mask: bool) -> Self {
            let seq = &genome.seq;
            let special_bm = seq.special_bitmap();
            let alias = seq.invalid_bitmap_aliases_special(soft_mask);
            let invalid_bm = if alias {
                None
            } else {
                seq.invalid_seed_bitmap(soft_mask)
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
                special_bm,
                invalid_bm,
                alias_invalid_to_special: alias,
                sep_pos,
                iupac_pos,
                iupac_byte,
                len: seq.len as u64,
            }
        }

        fn as_ffi(&self) -> GeckoPackedSeqHost {
            let special_ptr = self
                .special_bm
                .as_ref()
                .map_or(std::ptr::null(), Vec::as_ptr);
            let invalid_ptr = if self.alias_invalid_to_special {
                special_ptr
            } else {
                self.invalid_bm
                    .as_ref()
                    .map_or(std::ptr::null(), Vec::as_ptr)
            };
            let bm_words = self
                .special_bm
                .as_ref()
                .or(self.invalid_bm.as_ref())
                .map_or(0u64, |v| v.len() as u64);
            GeckoPackedSeqHost {
                bits: self.bits.as_ptr(),
                bits_words: self.bits.len() as u64,
                special_bm: special_ptr,
                invalid_bm: invalid_ptr,
                bm_words,
                sep_pos: if self.sep_pos.is_empty() {
                    std::ptr::null()
                } else {
                    self.sep_pos.as_ptr()
                },
                n_seps: self.sep_pos.len() as u32,
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
                n_iupac: self.iupac_pos.len() as u32,
                len: self.len,
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

        fn gecko_cuda_prepare_ref(
            ref_seq: *const GeckoPackedSeqHost,
            wl: c_int,
            freq_limit: u32,
            soft_mask: c_int,
            compare_reverse: c_int,
        ) -> *mut GeckoCudaRefHandle;

        fn gecko_cuda_query_and_extend(
            handle: *mut GeckoCudaRefHandle,
            query_seq: *const GeckoPackedSeqHost,
            min_len: u32,
            sim_th: f32,
            match_score: c_int,
            mismatch_penalty: c_int,
            n_penalty: c_int,
            out: *mut GeckoCudaResult,
        ) -> c_int;

        fn gecko_cuda_free_ref(handle: *mut GeckoCudaRefHandle);
        fn gecko_cuda_free_frags(frags: *mut GpuFrag);
    }

    pub struct RefGuard(*mut GeckoCudaRefHandle);
    impl Drop for RefGuard {
        fn drop(&mut self) {
            if !self.0.is_null() {
                unsafe { gecko_cuda_free_ref(self.0) };
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
            let is_rev = gf.strand != 0;
            let seq_x = query.seq_index_for_data_pos(gf.x_start as u64) as u64;
            let seq_y_buf = if is_rev {
                find_rc_seq_id(&rc_starts, gf.y_start as u64)
            } else {
                reference.seq_index_for_data_pos(gf.y_start as u64) as u64
            };
            let length = (gf.x_end - gf.x_start + 1) as u64;
            let sim = ((gf.score as f64) * 100.0 / ((length as f64) * match_sc)) as f32;
            coverage.record(length, sim);

            let frag = Fragment {
                diag: gf.diag,
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
            word_matches: raw.leftmax_hits,
            raw_hits: raw.leftmax_hits,
            filtered_hits: raw.leftmax_hits,
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
    }

    impl<'a> PreparedRefInner<'a> {
        pub fn build(
            reference: &'a Genome,
            config: &CompareConfig,
        ) -> Result<(Self, CudaDeviceInfo)> {
            validate_cuda_config(config)?;
            if reference.seq.len >= (1usize << 31) {
                return Err(GeckoError::InvalidParameter(
                    "CUDA backend requires reference < 2.14 Gb".into(),
                ));
            }
            let info = probe_cuda().map_err(GeckoError::InvalidParameter)?;
            let handle = {
                let ref_view = HostPackedView::new(reference, config.soft_mask);
                let ref_ffi = ref_view.as_ffi();
                RefGuard(unsafe {
                    gecko_cuda_prepare_ref(
                        &ref_ffi,
                        config.word_length as c_int,
                        config.freq_limit as u32,
                        c_int::from(config.soft_mask),
                        c_int::from(config.compare_reverse),
                    )
                })
            };
            if handle.0.is_null() {
                return Err(GeckoError::InvalidParameter(
                    "CUDA reference index preparation failed".into(),
                ));
            }
            Ok((
                Self {
                    handle,
                    reference,
                    soft_mask: config.soft_mask,
                },
                info,
            ))
        }

        pub fn query(&mut self, query: &Genome, config: &CompareConfig) -> Result<CompareResult> {
            validate_cuda_config(config)?;
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
            let mut raw = GeckoCudaResult::default();
            let rc = {
                let q_view = HostPackedView::new(query, self.soft_mask);
                let q_ffi = q_view.as_ffi();
                unsafe {
                    gecko_cuda_query_and_extend(
                        self.handle.0,
                        &q_ffi,
                        min_len as u32,
                        config.min_similarity,
                        config.scoring.match_score as c_int,
                        config.scoring.mismatch_penalty as c_int,
                        config.scoring.n_penalty as c_int,
                        &mut raw,
                    )
                }
            };
            if rc != 0 {
                if !raw.frags.is_null() {
                    unsafe { gecko_cuda_free_frags(raw.frags) };
                }
                return Err(GeckoError::InvalidParameter(format!(
                    "CUDA query/extend execution failed (rc={rc})"
                )));
            }
            let result = finalize_cuda_result(query, self.reference, config, &raw);
            if !raw.frags.is_null() {
                unsafe { gecko_cuda_free_frags(raw.frags) };
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

    /// Overlapped I/O + GPU pipeline: loads `query_path` on a background thread while
    /// the main thread loads `ref_path`, zero-copy uploads 2-bit `PackedSeq` to the GPU,
    /// builds 2-bit `Y_rc`, and constructs the CUB Hash-CSR reference index.
    pub fn compare_cuda_paths_overlapped(
        query_path: &std::path::Path,
        ref_path: &std::path::Path,
        config: &CompareConfig,
    ) -> Result<(Genome, Genome, CompareResult, CudaDeviceInfo)> {
        validate_cuda_config(config)?;
        let info = probe_cuda().map_err(GeckoError::InvalidParameter)?;

        let (query, reference, handle) = std::thread::scope(|s| -> Result<_> {
            let q_thread = s.spawn(|| Genome::from_fasta_path(query_path));
            let reference = Genome::from_fasta_path(ref_path)?;
            if reference.seq.len >= (1usize << 31) {
                let _ = q_thread.join();
                return Err(GeckoError::InvalidParameter(
                    "CUDA backend requires reference < 2.14 Gb".into(),
                ));
            }
            let handle = {
                let ref_view = HostPackedView::new(&reference, config.soft_mask);
                let ref_ffi = ref_view.as_ffi();
                RefGuard(unsafe {
                    gecko_cuda_prepare_ref(
                        &ref_ffi,
                        config.word_length as c_int,
                        config.freq_limit as u32,
                        c_int::from(config.soft_mask),
                        c_int::from(config.compare_reverse),
                    )
                })
            };
            let query = q_thread.join().expect("FASTA query loader panicked")?;
            if handle.0.is_null() {
                return Err(GeckoError::InvalidParameter(
                    "CUDA reference index preparation failed".into(),
                ));
            }
            Ok((query, reference, handle))
        })?;

        if query.seq.len >= (1usize << 31) {
            return Err(GeckoError::InvalidParameter(
                "CUDA backend requires query < 2.14 Gb".into(),
            ));
        }

        let min_len = if config.fixed_length {
            config.min_length
        } else {
            let shorter = query.total_bases.min(reference.total_bases);
            ((shorter as f64) * (config.min_length as f64 / 100.0)) as u64
        };

        let mut raw = GeckoCudaResult::default();
        let rc = {
            let q_view = HostPackedView::new(&query, config.soft_mask);
            let q_ffi = q_view.as_ffi();
            unsafe {
                gecko_cuda_query_and_extend(
                    handle.0,
                    &q_ffi,
                    min_len as u32,
                    config.min_similarity,
                    config.scoring.match_score as c_int,
                    config.scoring.mismatch_penalty as c_int,
                    config.scoring.n_penalty as c_int,
                    &mut raw,
                )
            }
        };
        if rc != 0 {
            if !raw.frags.is_null() {
                unsafe { gecko_cuda_free_frags(raw.frags) };
            }
            return Err(GeckoError::InvalidParameter(format!(
                "CUDA query/extend execution failed (rc={rc})"
            )));
        }

        let result = finalize_cuda_result(&query, &reference, config, &raw);
        if !raw.frags.is_null() {
            unsafe { gecko_cuda_free_frags(raw.frags) };
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
