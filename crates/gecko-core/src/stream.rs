use std::cmp::Reverse;
use std::time::Instant;

use rayon::prelude::*;

use crate::error::{GeckoError, Result};
use crate::extend::{
    CoverageMatrix, ExtendOptions, SkipScratch, exact_rejected_skip,
    extend_single_hit_with_seed_portable, has_avx2, has_bmi2, has_popcnt, pack_even_bits_swar,
    seed_init_stats_from_diff,
};
#[cfg(target_arch = "x86_64")]
use crate::extend::{extend_single_hit_with_seed_avx2_bmi2, extend_single_hit_with_seed_avx2_swar};
use crate::formats::frags::Fragment;
use crate::hits::Hit;
use crate::pipeline::{CompareConfig, CompareResult, CompareStats, profile_enabled};
use crate::ref_index::{
    BlockedBloom, PosStrandVec, RefIndex, RefIndexKind, fmix64, prefetch_t1, prefetch_t2,
    probe_bucket_binary_u32, probe_sorted_tags_scalar,
};
use crate::seq::{Genome, PackedSeq};

const QUERY_CHUNK_BYTES: usize = 1 << 20; // ~1 Mb dynamic Rayon chunks
const BATCH_SIZE: usize = 64;
const MISS_BUF_CAP: usize = 16; // 16 * 16 B = 256 B (matches Haswell/Skylake LFB depth)

const HSP_CACHE_BITS: u32 = 17;
const HSP_CACHE_LEN: usize = 1 << HSP_CACHE_BITS; // 131,072 slots (1 MiB @ 8 B, 2 MiB @ 16 B)
const NOISE_CACHE_BITS: u32 = 16;
const NOISE_CACHE_LEN: usize = 1 << NOISE_CACHE_BITS; // 65,536 slots (512 KiB @ 8 B, 1 MiB @ 16 B)
const FIB_MULT: u64 = 0x9E37_79B9_7F4A_7C15;

/// Upper sequence length bound for compact 8-byte `DiagSlot32` (128+ headroom from `i32::MIN/MAX`).
const COMPACT_DIAG_MAX_LEN: usize = (i32::MAX as usize) - 128;

const ACTIVE_RING_LEN: usize = 8;
const STRAND_BIT_MASK: u64 = 1u64 << 63;

// ---------------------------------------------------------------------------
// Compact 8-byte (`DiagSlot32`) and 16-byte (`DiagSlot64`) 2-tier diagonal cache
// ---------------------------------------------------------------------------

trait DiagSlotOps: Copy + Send {
    const EMPTY: Self;
    fn covers(self, diag: i64, pos_x: u64) -> bool;
    fn update(&mut self, diag: i64, max_x: u64);
    /// BLAST/LASTZ non-overlapping two-hit check on `diag`:
    /// - If `diag` has no prior anchor, records `pos_x` and returns `false`.
    /// - If `pos_x - last_x < wl` (overlapping seed), keeps `last_x` unchanged and returns `false`.
    /// - If `wl <= pos_x - last_x <= wl + window`, updates anchor to `pos_x` and returns `true`.
    /// - Otherwise (`pos_x - last_x > wl + window`), updates anchor to `pos_x` and returns `false`.
    fn two_hit_check(&mut self, diag: i64, pos_x: u64, wl: u64, window: u64) -> bool;
}

#[derive(Clone, Copy)]
#[repr(C, align(8))]
struct DiagSlot32 {
    diag: i32,
    max_x: u32,
}
const _: () = assert!(std::mem::size_of::<DiagSlot32>() == 8);

impl DiagSlotOps for DiagSlot32 {
    const EMPTY: Self = Self {
        diag: i32::MIN,
        max_x: 0,
    };

    #[inline(always)]
    fn covers(self, diag: i64, pos_x: u64) -> bool {
        self.diag == diag as i32 && (pos_x as u32) <= self.max_x
    }

    #[inline(always)]
    fn update(&mut self, diag: i64, max_x: u64) {
        let d = diag as i32;
        let mx = max_x as u32;
        if self.diag == d {
            if mx > self.max_x {
                self.max_x = mx;
            }
        } else {
            *self = Self { diag: d, max_x: mx };
        }
    }

    #[inline(always)]
    fn two_hit_check(&mut self, diag: i64, pos_x: u64, wl: u64, window: u64) -> bool {
        let d = diag as i32;
        let px = pos_x as u32;
        if self.diag == d {
            let last_x = self.max_x as u64;
            if pos_x <= last_x {
                return false;
            }
            let dist = pos_x - last_x;
            if dist < wl {
                return false;
            }
            self.max_x = px;
            dist - wl <= window
        } else {
            *self = Self { diag: d, max_x: px };
            false
        }
    }
}

#[derive(Clone, Copy)]
#[repr(C)]
struct DiagSlot64 {
    diag: i64,
    max_x: u64,
}
const _: () = assert!(std::mem::size_of::<DiagSlot64>() == 16);

impl DiagSlotOps for DiagSlot64 {
    const EMPTY: Self = Self {
        diag: i64::MIN,
        max_x: 0,
    };

    #[inline(always)]
    fn covers(self, diag: i64, pos_x: u64) -> bool {
        self.diag == diag && pos_x <= self.max_x
    }

    #[inline(always)]
    fn update(&mut self, diag: i64, max_x: u64) {
        if self.diag == diag {
            if max_x > self.max_x {
                self.max_x = max_x;
            }
        } else {
            *self = Self { diag, max_x };
        }
    }

    #[inline(always)]
    fn two_hit_check(&mut self, diag: i64, pos_x: u64, wl: u64, window: u64) -> bool {
        if self.diag == diag {
            if pos_x <= self.max_x {
                return false;
            }
            let dist = pos_x - self.max_x;
            if dist < wl {
                return false;
            }
            self.max_x = pos_x;
            dist - wl <= window
        } else {
            *self = Self { diag, max_x: pos_x };
            false
        }
    }
}

/// Lock-free 2-tier per-thread, per-strand Fibonacci-hashed diagonal cache.
///
/// Separates accepted HSPs / exact reject-skip bounds (`hsp`, 1 MiB @ 8 B) from short rejected
/// hits (`noise`, 512 KiB @ 8 B) while preserving the exact 64-bit `FIB_MULT` slot mapping.
struct StrandDiagCache<S: DiagSlotOps> {
    hsp: Vec<S>,
    noise: Vec<S>,
}

impl<S: DiagSlotOps> StrandDiagCache<S> {
    fn new() -> Self {
        Self {
            hsp: vec![S::EMPTY; HSP_CACHE_LEN],
            noise: vec![S::EMPTY; NOISE_CACHE_LEN],
        }
    }

    fn clear(&mut self) {
        self.hsp.fill(S::EMPTY);
        self.noise.fill(S::EMPTY);
    }

    #[inline(always)]
    fn is_covered(&self, diag: i64, pos_x: u64) -> bool {
        let h = (diag as u64).wrapping_mul(FIB_MULT);
        let hi = (h >> (64 - HSP_CACHE_BITS)) as usize;
        // Safety: `hi < HSP_CACHE_LEN` by right-shift of 64-bit value.
        let hslot = unsafe { *self.hsp.get_unchecked(hi) };
        if hslot.covers(diag, pos_x) {
            return true;
        }
        let ni = (h >> (64 - NOISE_CACHE_BITS)) as usize;
        // Safety: `ni < NOISE_CACHE_LEN` by right-shift of 64-bit value.
        let nslot = unsafe { *self.noise.get_unchecked(ni) };
        nslot.covers(diag, pos_x)
    }

    #[inline(always)]
    fn update_hsp(&mut self, diag: i64, max_x: u64) {
        let h = (diag as u64).wrapping_mul(FIB_MULT);
        let hi = (h >> (64 - HSP_CACHE_BITS)) as usize;
        // Safety: `hi < HSP_CACHE_LEN`.
        unsafe { self.hsp.get_unchecked_mut(hi) }.update(diag, max_x);
    }

    #[inline(always)]
    fn update_noise(&mut self, diag: i64, max_x: u64) {
        let h = (diag as u64).wrapping_mul(FIB_MULT);
        let ni = (h >> (64 - NOISE_CACHE_BITS)) as usize;
        // Safety: `ni < NOISE_CACHE_LEN`.
        unsafe { self.noise.get_unchecked_mut(ni) }.update(diag, max_x);
    }
}

/// Opt-in per-thread, per-strand Fibonacci-hashed diagonal anchor table for `--two-hit` (`F10`).
/// Only allocated when `CompareConfig::two_hit_window.is_some()`.
struct TwoHitTable<S: DiagSlotOps> {
    slots: Vec<S>,
}

impl<S: DiagSlotOps> TwoHitTable<S> {
    fn new() -> Self {
        Self {
            slots: vec![S::EMPTY; NOISE_CACHE_LEN],
        }
    }

    fn clear(&mut self) {
        self.slots.fill(S::EMPTY);
    }

    #[inline(always)]
    fn accept(&mut self, diag: i64, pos_x: u64, wl: u64, window: u64) -> bool {
        let h = (diag as u64).wrapping_mul(FIB_MULT);
        let idx = (h >> (64 - NOISE_CACHE_BITS)) as usize;
        // Safety: `idx < NOISE_CACHE_LEN` by right-shift of 64-bit value.
        unsafe { self.slots.get_unchecked_mut(idx) }.two_hit_check(diag, pos_x, wl, window)
    }
}

// ---------------------------------------------------------------------------
// 128 B L1d Active-Diagonal Ring for 1-bit `unique_bm` exact-probe bypass (F4)
// ---------------------------------------------------------------------------

/// 16-byte active-diagonal slot (8 slots = 128 B = 2 L1d cache lines).
#[derive(Clone, Copy, Default)]
#[repr(C)]
struct ActiveDiag {
    diag: i64,
    /// Bits `[62..0]`: exclusive upper bound `max_x + 1` (`0` = empty slot).
    /// Bit `63`: strand (`0` = forward, `1` = reverse).
    end_x_strand: u64,
}
const _: () = assert!(std::mem::size_of::<ActiveDiag>() == 16);

#[repr(C, align(64))]
struct ActiveDiagRing {
    entries: [ActiveDiag; ACTIVE_RING_LEN],
    /// Tight max of `end_x` across all 8 entries -> 1-cycle fast reject `pos_x >= max_end_x`
    /// when outside all active HSP / `exact_rejected_skip` spans.
    max_end_x: u64,
}

impl ActiveDiagRing {
    fn new() -> Self {
        Self {
            entries: [ActiveDiag::default(); ACTIVE_RING_LEN],
            max_end_x: 0,
        }
    }

    fn clear(&mut self) {
        self.entries = [ActiveDiag::default(); ACTIVE_RING_LEN];
        self.max_end_x = 0;
    }

    /// Record an accepted HSP (`outcome.xfil_max`) or `exact_rejected_skip` (`skip_end`) span.
    /// Never call for short `update_noise` bounds.
    #[inline(always)]
    fn record(&mut self, diag: i64, max_x: u64, strand: u8) {
        let end_x = max_x.saturating_add(1) & !STRAND_BIT_MASK;
        let packed = end_x | ((strand as u64) << 63);
        for i in 0..ACTIVE_RING_LEN {
            let e = self.entries[i];
            let e_end = e.end_x_strand & !STRAND_BIT_MASK;
            if e_end > 0 && e.diag == diag && ((e.end_x_strand >> 63) as u8) == strand {
                let new_packed = if end_x > e_end {
                    packed
                } else {
                    e.end_x_strand
                };
                if i > 0 {
                    self.entries.copy_within(0..i, 1);
                }
                self.entries[0] = ActiveDiag {
                    diag,
                    end_x_strand: new_packed,
                };
                if end_x > self.max_end_x {
                    self.max_end_x = end_x;
                }
                return;
            }
        }
        self.entries.copy_within(0..ACTIVE_RING_LEN - 1, 1);
        self.entries[0] = ActiveDiag {
            diag,
            end_x_strand: packed,
        };
        let mut m = 0u64;
        for e in &self.entries {
            m = m.max(e.end_x_strand & !STRAND_BIT_MASK);
        }
        self.max_end_x = m;
    }
}

/// Check if `pos_x` hits a unique, non-palindromic, in-slab reference k-mer on an already-covered
/// active diagonal. If so, increments `word_matches` and `raw_hits` by 1 and returns `true`,
/// bypassing `fmix64` + Bloom + CSR DRAM lookups completely with 100% bit-exactness.
#[inline(always)]
#[allow(clippy::too_many_arguments)]
unsafe fn try_active_diag_bypass<const INVERT_UNIQUE: bool, S: DiagSlotOps>(
    pos_x: u64,
    seq_x: &PackedSeq,
    seq_y: &PackedSeq,
    rev_seq_y: Option<&PackedSeq>,
    unique_bm: &[u64],
    max_y_start: u64,
    seed_mask: u64,
    state: &mut WorkerState<S>,
) -> bool {
    if pos_x >= state.active_ring.max_end_x {
        return false;
    }
    if unique_bm.is_empty() {
        return false;
    }
    let q_word = seq_x.load_32(pos_x as usize);
    for i in 0..ACTIVE_RING_LEN {
        let e = unsafe { *state.active_ring.entries.get_unchecked(i) };
        let end_x = e.end_x_strand & !STRAND_BIT_MASK;
        if pos_x >= end_x {
            continue;
        }
        let target_y = (pos_x as i64).wrapping_sub(e.diag) as u64;
        if target_y > max_y_start {
            continue;
        }
        let is_rev = (e.end_x_strand & STRAND_BIT_MASK) != 0;
        let fwd_pos_y = if !is_rev {
            target_y as usize
        } else {
            (max_y_start - target_y) as usize
        };
        // Check `unique_bm` BEFORE loading `seq_y.bits` so non-unique positions never touch `bits`.
        let raw_w = unsafe { *unique_bm.get_unchecked(fwd_pos_y >> 6) };
        let w = if INVERT_UNIQUE { !raw_w } else { raw_w };
        if ((w >> (fwd_pos_y & 63)) & 1) == 0 {
            continue;
        }
        let ref_word = if !is_rev {
            seq_y.load_32(target_y as usize)
        } else if let Some(rs) = rev_seq_y {
            rs.load_32(target_y as usize)
        } else {
            continue;
        };
        if ((q_word ^ ref_word) & seed_mask) == 0 {
            if i > 0 {
                state.active_ring.entries.swap(0, i);
            }
            state.word_matches += 1;
            state.raw_hits += 1;
            state.bypass_hits += 1;
            return true;
        }
    }
    false
}

// ---------------------------------------------------------------------------
// Pure-scalar 8-base LUT short-death fast-reject engine (F5d, zero AVX2 / zero `vzeroupper`)
// ---------------------------------------------------------------------------

#[derive(Clone, Copy)]
#[repr(C, align(4))]
struct ByteStep {
    delta: i8,   // net unit-score change over 8 bases (-8..=+8)
    min_pre: i8, // min inclusive prefix sum over 1..=8 bases (-8..=+1)
    max_pre: i8, // max inclusive prefix sum over 1..=8 bases (-1..=+8)
    max_idx: u8, // FIRST 0..=7 index achieving `max_pre` (strict `>` tie-break)
}

const fn build_byte_step_lut(rev: bool) -> [ByteStep; 256] {
    let mut lut = [ByteStep {
        delta: 0,
        min_pre: 0,
        max_pre: 0,
        max_idx: 0,
    }; 256];
    let mut b = 0usize;
    while b < 256 {
        let mut sum = 0i8;
        let mut min_pre = 127i8;
        let mut max_pre = -128i8;
        let mut max_idx = 0u8;
        let mut i = 0u32;
        while i < 8 {
            let bit = if rev {
                (b >> (7 - i)) & 1
            } else {
                (b >> i) & 1
            };
            sum += if bit != 0 { -1 } else { 1 };
            if sum < min_pre {
                min_pre = sum;
            }
            if sum > max_pre {
                max_pre = sum;
                max_idx = i as u8;
            }
            i += 1;
        }
        lut[b] = ByteStep {
            delta: sum,
            min_pre,
            max_pre,
            max_idx,
        };
        b += 1;
    }
    lut
}

const BYTE_STEP_LUT_FWD: [ByteStep; 256] = build_byte_step_lut(false);
const BYTE_STEP_LUT_REV: [ByteStep; 256] = build_byte_step_lut(true);

#[inline(always)]
fn step32_lut<const REV: bool>(
    mm32: u32,
    u: &mut i32,
    u_max: &mut i32,
    best_off: &mut u32,
    base_off: u32,
) -> bool {
    let lut = if REV {
        &BYTE_STEP_LUT_REV
    } else {
        &BYTE_STEP_LUT_FWD
    };
    let mut cur = *u;
    for k in 0..4u32 {
        let b = if REV {
            ((mm32 >> (24 - 8 * k)) & 0xFF) as usize
        } else {
            ((mm32 >> (8 * k)) & 0xFF) as usize
        };
        let e = unsafe { *lut.get_unchecked(b) };
        if cur + (e.min_pre as i32) < 0 {
            if cur + (e.max_pre as i32) > *u_max {
                for i in 0..8u32 {
                    let bit = if REV {
                        (b >> (7 - i)) & 1
                    } else {
                        (b >> i) & 1
                    };
                    cur += if bit != 0 { -1 } else { 1 };
                    if cur > *u_max {
                        *u_max = cur;
                        *best_off = base_off + 8 * k + i + 1;
                    }
                    if cur < 0 {
                        break;
                    }
                }
            }
            return true;
        }
        if cur + (e.max_pre as i32) > *u_max {
            *u_max = cur + (e.max_pre as i32);
            *best_off = base_off + 8 * k + e.max_idx as u32 + 1;
        }
        cur += e.delta as i32;
    }
    *u = cur;
    false
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "bmi2")]
#[inline]
unsafe fn pack_mm32_bmi2(mm_even: u64) -> u32 {
    core::arch::x86_64::_pext_u64(mm_even, 0x5555_5555_5555_5555) as u32
}

#[inline(always)]
unsafe fn pack_mm32_swar(mm_even: u64) -> u32 {
    pack_even_bits_swar(mm_even)
}

/// Pure-scalar 64-bp (2x32-bp) bidirectional drop-to-zero short-death fast-reject helper.
///
/// Returns `Some((length, u_max))` iff BOTH forward and backward passes strictly terminate
/// (`score < 0`) within 64 bp of the seed AND `length <= min_len`, guaranteeing 100% bit-exact
/// rejection without touching AVX2, `special_runs` binary search, `seq_index_for_data_pos`,
/// FP64 division, or `$miss_slow`.
#[inline(always)]
#[allow(clippy::too_many_arguments)]
unsafe fn try_fast_short_reject_64<F: FnMut(u64) -> u32>(
    seq_x: &PackedSeq,
    seq_y: &PackedSeq,
    pos_x: u64,
    pos_y: u64,
    wl: u64,
    seed_score: i64,
    min_len: u64,
    mut pack_mm32: F,
) -> Option<(usize, u32)> {
    let mut u = (seed_score >> 2) as i32;
    let mut u_max = u;
    let mut fwd_best = 0u32;

    // Forward chunk 1: `[pos + wl .. pos + wl + 32)`
    let fx1 = (pos_x + wl) as usize;
    let fy1 = (pos_y + wl) as usize;
    let d_f1 = seq_x.load_32(fx1) ^ seq_y.load_32(fy1);
    let mm_f1 = pack_mm32((d_f1 | (d_f1 >> 1)) & 0x5555_5555_5555_5555);
    if !step32_lut::<false>(mm_f1, &mut u, &mut u_max, &mut fwd_best, 0) {
        // Forward chunk 2: `[pos + wl + 32 .. pos + wl + 64)`
        let d_f2 = seq_x.load_32(fx1 + 32) ^ seq_y.load_32(fy1 + 32);
        let mm_f2 = pack_mm32((d_f2 | (d_f2 >> 1)) & 0x5555_5555_5555_5555);
        if !step32_lut::<false>(mm_f2, &mut u, &mut u_max, &mut fwd_best, 32) {
            return None;
        }
    }

    // Backward chunk 1: `[pos - 32 .. pos)` (lane 31 = `pos - 1`, lane 0 = `pos - 32`)
    u = u_max;
    let mut bwd_best = 0u32;
    let bx1 = (pos_x - 32) as usize;
    let by1 = (pos_y - 32) as usize;
    let d_b1 = seq_x.load_32(bx1) ^ seq_y.load_32(by1);
    let mm_b1 = pack_mm32((d_b1 | (d_b1 >> 1)) & 0x5555_5555_5555_5555);
    if !step32_lut::<true>(mm_b1, &mut u, &mut u_max, &mut bwd_best, 0) {
        // Backward chunk 2: `[pos - 64 .. pos - 32)`
        let d_b2 = seq_x.load_32(bx1 - 32) ^ seq_y.load_32(by1 - 32);
        let mm_b2 = pack_mm32((d_b2 | (d_b2 >> 1)) & 0x5555_5555_5555_5555);
        if !step32_lut::<true>(mm_b2, &mut u, &mut u_max, &mut bwd_best, 32) {
            return None;
        }
    }

    let length = wl + fwd_best as u64 + bwd_best as u64;
    if length > min_len {
        return None;
    }
    Some((length as usize, u_max as u32))
}

// ---------------------------------------------------------------------------
// Typed `pos_strand` readers, batch items, Stage-3c miss items, and worker state
// ---------------------------------------------------------------------------

/// Typed `pos_strand` slice abstraction so `flush_hash` / `flush_direct` monomorphize once over
/// `u32` vs `u64` instead of matching on [`PosStrandVec`] per raw hit.
trait PsRead: Copy {
    unsafe fn unpack_at(self, idx: usize) -> (u64, u8);
    unsafe fn prefetch_at(self, idx: usize);
}

impl PsRead for *const u32 {
    #[inline(always)]
    unsafe fn unpack_at(self, idx: usize) -> (u64, u8) {
        let raw = unsafe { *self.add(idx) };
        ((raw >> 1) as u64, (raw & 1) as u8)
    }

    #[inline(always)]
    unsafe fn prefetch_at(self, idx: usize) {
        prefetch_t2(unsafe { self.add(idx) });
    }
}

impl PsRead for *const u64 {
    #[inline(always)]
    unsafe fn unpack_at(self, idx: usize) -> (u64, u8) {
        let raw = unsafe { *self.add(idx) };
        (raw >> 1, (raw & 1) as u8)
    }

    #[inline(always)]
    unsafe fn prefetch_at(self, idx: usize) {
        prefetch_t2(unsafe { self.add(idx) });
    }
}

trait PsSlice: Copy + Send + Sync {
    type Ptr: PsRead;
    fn as_ps_ptr(self) -> Self::Ptr;
}

impl PsSlice for &[u32] {
    type Ptr = *const u32;
    #[inline(always)]
    fn as_ps_ptr(self) -> *const u32 {
        self.as_ptr()
    }
}

impl PsSlice for &[u64] {
    type Ptr = *const u64;
    #[inline(always)]
    fn as_ps_ptr(self) -> *const u64 {
        self.as_ptr()
    }
}

/// Compact 32-byte batch item storing query metadata, 32-bit hash tag, Bloom block index, and
/// in-place bucket / hit range `[start, end)` across all prefetch stages.
#[derive(Clone, Copy, Default)]
#[repr(C)]
struct QueryItem {
    pos_x: u64,
    b_idx: u32,
    tag: u32,
    start: u32,
    end: u32,
    bloom_blk: u32,
    /// `0` = forward canonical (`kf < kr`), `1` = revcomp canonical (`kr < kf`), `2` = palindrome (`kf == kr`).
    q_kind: u8,
}
const _: () = assert!(std::mem::size_of::<QueryItem>() == 32);

/// Compact 16-byte uncovered miss queued in Stage 3b and drained in Stage 3c.
/// Bit 63 of `pos_y_rev` encodes `is_rev` (`0` = forward `reference`, `1` = `rev_ref`).
#[derive(Clone, Copy, Default)]
#[repr(C)]
struct MissItem {
    pos_x: u64,
    pos_y_rev: u64,
}
const _: () = assert!(std::mem::size_of::<MissItem>() == 16);

struct WorkerState<S: DiagSlotOps> {
    fwd_cache: StrandDiagCache<S>,
    rev_cache: Option<StrandDiagCache<S>>,
    fwd_two_hit: Option<TwoHitTable<S>>,
    rev_two_hit: Option<TwoHitTable<S>>,
    two_hit_window: u32,
    active_ring: ActiveDiagRing,
    has_unique_bm: bool,
    skip_scratch: SkipScratch,
    last_chunk_end: usize,
    sp_x_cursor: usize,
    fragments: Vec<Fragment>,
    coverage: Option<CoverageMatrix>,
    word_matches: u64,
    raw_hits: u64,
    hits_used: u64,
    bypass_hits: u64,
    bloom_rejected: u64,
    fast_rejects: u64,
    two_hit_dropped: u64,
}

impl<S: DiagSlotOps> WorkerState<S> {
    fn new(
        compare_reverse: bool,
        collect_coverage: bool,
        has_unique_bm: bool,
        two_hit_window: Option<u32>,
    ) -> Self {
        let th_win = two_hit_window.unwrap_or(0);
        Self {
            fwd_cache: StrandDiagCache::new(),
            rev_cache: if compare_reverse {
                Some(StrandDiagCache::new())
            } else {
                None
            },
            fwd_two_hit: if th_win > 0 {
                Some(TwoHitTable::new())
            } else {
                None
            },
            rev_two_hit: if th_win > 0 && compare_reverse {
                Some(TwoHitTable::new())
            } else {
                None
            },
            two_hit_window: th_win,
            active_ring: ActiveDiagRing::new(),
            has_unique_bm,
            skip_scratch: SkipScratch::default(),
            last_chunk_end: 0,
            sp_x_cursor: 0,
            fragments: Vec::new(),
            coverage: if collect_coverage {
                Some(CoverageMatrix::new())
            } else {
                None
            },
            word_matches: 0,
            raw_hits: 0,
            hits_used: 0,
            bypass_hits: 0,
            bloom_rejected: 0,
            fast_rejects: 0,
            two_hit_dropped: 0,
        }
    }

    /// Advance the monotonic `query.seq.special_runs` cursor so `sp_x[..cursor].end <= pos_x`.
    #[inline(always)]
    fn sp_x_cursor_for(&mut self, sp_x: &[(u64, u64)], pos_x: u64) -> usize {
        let mut c = self.sp_x_cursor;
        while c < sp_x.len() && unsafe { sp_x.get_unchecked(c).1 } <= pos_x {
            c += 1;
        }
        self.sp_x_cursor = c;
        c
    }
}

/// Runtime ISA tier selected once per comparison.
#[derive(Clone, Copy)]
enum IsaTier {
    #[cfg(target_arch = "x86_64")]
    Avx2Bmi2Popcnt,
    #[cfg(target_arch = "x86_64")]
    Avx2Popcnt,
    Portable,
}

/// Execute the streaming canonical-index comparison over `query` against a pre-built [`RefIndex`],
/// populating the `.MAT` [`CoverageMatrix`].
pub fn compare_streaming(
    query: &Genome,
    reference: &Genome,
    ref_index: &RefIndex,
    rev_ref: Option<&Genome>,
    config: &CompareConfig,
) -> Result<CompareResult> {
    compare_streaming_with_coverage(query, reference, ref_index, rev_ref, config, true)
}

/// Execute the streaming comparison over `query` against a pre-built [`RefIndex`] with optional
/// `.MAT` [`CoverageMatrix`] collection (`collect_coverage = false` skips per-worker 800 KiB
/// matrix allocations and rejected-hit FP64 similarity divisions while keeping fragments and
/// stats 100% bit-identical).
pub fn compare_streaming_with_coverage(
    query: &Genome,
    reference: &Genome,
    ref_index: &RefIndex,
    rev_ref: Option<&Genome>,
    config: &CompareConfig,
    collect_coverage: bool,
) -> Result<CompareResult> {
    config.validate()?;

    let compare_reverse = config.compare_reverse && rev_ref.is_some();
    if ref_index.forward_only && compare_reverse {
        return Err(GeckoError::InvalidParameter(
            "forward-only RefIndex cannot be queried with compare_reverse=true".into(),
        ));
    }

    let wl = config.word_length;
    if query.seq.len < wl || ref_index.kept_entries == 0 {
        return Ok(CompareResult {
            fragments: Vec::new(),
            coverage: CoverageMatrix::new(),
            stats: CompareStats::default(),
        });
    }

    let use_compact_diag =
        query.seq.len <= COMPACT_DIAG_MAX_LEN && reference.seq.len <= COMPACT_DIAG_MAX_LEN;

    match (ref_index.pos_strand_vec(), use_compact_diag) {
        (PosStrandVec::U32(v), true) => compare_streaming_typed::<&[u32], DiagSlot32>(
            query,
            reference,
            ref_index,
            v.as_slice(),
            rev_ref,
            config,
            compare_reverse,
            collect_coverage,
        ),
        (PosStrandVec::U32(v), false) => compare_streaming_typed::<&[u32], DiagSlot64>(
            query,
            reference,
            ref_index,
            v.as_slice(),
            rev_ref,
            config,
            compare_reverse,
            collect_coverage,
        ),
        (PosStrandVec::U64(v), _) => compare_streaming_typed::<&[u64], DiagSlot64>(
            query,
            reference,
            ref_index,
            v.as_slice(),
            rev_ref,
            config,
            compare_reverse,
            collect_coverage,
        ),
    }
}

#[allow(clippy::too_many_arguments)]
fn compare_streaming_typed<PS: PsSlice, S: DiagSlotOps>(
    query: &Genome,
    reference: &Genome,
    ref_index: &RefIndex,
    ps_slice: PS,
    rev_ref: Option<&Genome>,
    config: &CompareConfig,
    compare_reverse: bool,
    collect_coverage: bool,
) -> Result<CompareResult> {
    let wl = config.word_length;
    let isa = {
        #[cfg(target_arch = "x86_64")]
        {
            if has_avx2() && has_popcnt() {
                if has_bmi2() {
                    IsaTier::Avx2Bmi2Popcnt
                } else {
                    IsaTier::Avx2Popcnt
                }
            } else {
                IsaTier::Portable
            }
        }
        #[cfg(not(target_arch = "x86_64"))]
        {
            let _ = (has_avx2(), has_bmi2(), has_popcnt());
            IsaTier::Portable
        }
    };

    let mask = if wl == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * wl)) - 1
    };
    let rc_shift = (2 * (wl - 1)) as u32;
    let ny = reference.seq.len as u64;

    let fwd_opts = ExtendOptions {
        min_length: config.min_length,
        min_similarity: config.min_similarity,
        word_length: wl,
        fixed_length: config.fixed_length,
        strand: 'f',
        scoring: config.scoring,
        spaced_seed: ref_index.spaced_mask.is_some(),
    };
    let rev_opts = ExtendOptions {
        strand: 'r',
        ..fwd_opts
    };

    let max_starts = query.seq.len - wl + 1;
    let mut chunks = Vec::with_capacity(max_starts.div_ceil(QUERY_CHUNK_BYTES));
    let mut s = 0usize;
    while s < max_starts {
        let e = (s + QUERY_CHUNK_BYTES).min(max_starts);
        chunks.push((s, e));
        s = e;
    }

    let has_unique_bm = !ref_index.unique_bm.is_empty();
    let two_hit_window = config.two_hit_window;
    let t_query = Instant::now();
    let worker_results: Vec<WorkerState<S>> = chunks
        .into_par_iter()
        .fold(
            || {
                WorkerState::<S>::new(
                    compare_reverse,
                    collect_coverage,
                    has_unique_bm,
                    two_hit_window,
                )
            },
            |mut state, (chunk_start, chunk_end)| {
                if chunk_start < state.last_chunk_end {
                    state.fwd_cache.clear();
                    if let Some(rc) = state.rev_cache.as_mut() {
                        rc.clear();
                    }
                    if let Some(ft) = state.fwd_two_hit.as_mut() {
                        ft.clear();
                    }
                    if let Some(rt) = state.rev_two_hit.as_mut() {
                        rt.clear();
                    }
                    state.active_ring.clear();
                    state.sp_x_cursor = if query.seq.special_runs.is_empty() {
                        0
                    } else {
                        query
                            .seq
                            .special_runs
                            .partition_point(|&(_, e)| e <= chunk_start as u64)
                    };
                } else if state.last_chunk_end == 0 && chunk_start > 0 {
                    state.sp_x_cursor = if query.seq.special_runs.is_empty() {
                        0
                    } else {
                        query
                            .seq
                            .special_runs
                            .partition_point(|&(_, e)| e <= chunk_start as u64)
                    };
                }

                scan_query_chunk(
                    query,
                    reference,
                    rev_ref,
                    ref_index,
                    ps_slice,
                    chunk_start,
                    chunk_end,
                    wl,
                    ny,
                    mask,
                    rc_shift,
                    config.soft_mask,
                    config.max_rep.max(1),
                    fwd_opts,
                    rev_opts,
                    isa,
                    &mut state,
                );

                state.last_chunk_end = chunk_end;
                state
            },
        )
        .collect();
    let query_secs = t_query.elapsed().as_secs_f64();
    let t_dedup = Instant::now();

    let total_frags: usize = worker_results.iter().map(|w| w.fragments.len()).sum();
    let mut fragments = Vec::with_capacity(total_frags);
    let mut coverage = CoverageMatrix::new();
    let mut stats = CompareStats::default();
    let mut bypass_hits = 0u64;
    let mut bloom_rejected = 0u64;
    let mut fast_rejects = 0u64;
    let mut two_hit_dropped = 0u64;

    for mut w in worker_results {
        stats.word_matches += w.word_matches;
        stats.raw_hits += w.raw_hits;
        stats.hits_used += w.hits_used;
        bypass_hits += w.bypass_hits;
        bloom_rejected += w.bloom_rejected;
        fast_rejects += w.fast_rejects;
        two_hit_dropped += w.two_hit_dropped;
        if let Some(cov) = w.coverage.as_ref() {
            coverage.merge_from(cov);
        }
        fragments.append(&mut w.fragments);
    }
    stats.filtered_hits = stats.hits_used;

    dedup_fragments(&mut fragments);

    for f in &fragments {
        if f.strand == 'f' {
            stats.forward_fragments += 1;
        } else {
            stats.reverse_fragments += 1;
        }
    }
    let dedup_secs = t_dedup.elapsed().as_secs_f64();

    if profile_enabled() {
        let kind_label = match &ref_index.kind {
            RefIndexKind::DirectCsr { .. } => "DirectCsr",
            RefIndexKind::HashCsr { .. } => "HashCsr",
        };
        eprintln!(
            "[GECKO_PROFILE] stream_query kind={} k={} parts=[{}..{}) kept_entries={} word_matches={} raw_hits={} hits_used={} frags={} bypass={} bloom_reject={} fast_reject={} two_hit_drop={} query_ms={:.2} dedup_ms={:.2}",
            kind_label,
            wl,
            ref_index.part_start,
            ref_index.part_end,
            ref_index.kept_entries,
            stats.word_matches,
            stats.raw_hits,
            stats.hits_used,
            fragments.len(),
            bypass_hits,
            bloom_rejected,
            fast_rejects,
            two_hit_dropped,
            query_secs * 1e3,
            dedup_secs * 1e3,
        );
    }

    Ok(CompareResult {
        fragments,
        coverage,
        stats,
    })
}

/// Inline AVX2 8-wide sorted `u32` tag probe expanded directly inside `flush_hash`'s Stage 3a
/// loop. For `len <= 8` (~98% of all buckets), executes a single branchless 256-bit load +
/// `vpcmpeqd` + `vmovmskps` + length mask with ZERO `vzeroupper` transitions.
#[cfg(target_arch = "x86_64")]
macro_rules! probe_bucket_avx2_inline {
    ($tags:expr, $tags_ptr:expr, $start:expr, $end:expr, $tag:expr) => {{
        use core::arch::x86_64::*;

        let start = $start;
        let end = $end;
        let tag = $tag;
        let len = end - start;
        let ptr = $tags_ptr;
        'probe: {
            if len <= 8 {
                let vtag = _mm256_set1_epi32(tag as i32);
                let vcand = unsafe { _mm256_loadu_si256(ptr.add(start) as *const __m256i) };
                let cmp = _mm256_cmpeq_epi32(vcand, vtag);
                let valid_mask = (1u32 << len) - 1;
                let mask = (_mm256_movemask_ps(_mm256_castsi256_ps(cmp)) as u32) & valid_mask;
                if mask == 0 {
                    break 'probe (start, start);
                }
                let lo = start + mask.trailing_zeros() as usize;
                let hi = start + (32 - mask.leading_zeros()) as usize;
                break 'probe (lo, hi);
            }
            if len <= 32 {
                let vtag = _mm256_set1_epi32(tag as i32);
                let mut j = start;
                while j + 8 <= end {
                    let vcand = unsafe { _mm256_loadu_si256(ptr.add(j) as *const __m256i) };
                    let cmp = _mm256_cmpeq_epi32(vcand, vtag);
                    let mask = _mm256_movemask_ps(_mm256_castsi256_ps(cmp)) as u32;
                    if mask != 0 {
                        let lo = j + mask.trailing_zeros() as usize;
                        if (mask & 0x80) == 0 {
                            let hi = j + (32 - mask.leading_zeros()) as usize;
                            break 'probe (lo, hi);
                        }
                        let mut hi = j + 8;
                        while hi < end && unsafe { *ptr.add(hi) } == tag {
                            hi += 1;
                        }
                        break 'probe (lo, hi);
                    }
                    if unsafe { *ptr.add(j + 7) } > tag {
                        break 'probe (start, start);
                    }
                    j += 8;
                }
                let rem = end - j;
                if rem > 0 {
                    let vcand = unsafe { _mm256_loadu_si256(ptr.add(j) as *const __m256i) };
                    let cmp = _mm256_cmpeq_epi32(vcand, vtag);
                    let valid_mask = (1u32 << rem) - 1;
                    let mask = (_mm256_movemask_ps(_mm256_castsi256_ps(cmp)) as u32) & valid_mask;
                    if mask != 0 {
                        let lo = j + mask.trailing_zeros() as usize;
                        let hi = j + (32 - mask.leading_zeros()) as usize;
                        break 'probe (lo, hi);
                    }
                }
                break 'probe (start, start);
            }
            probe_bucket_binary_u32($tags, start, end, tag)
        }
    }};
}

macro_rules! probe_bucket_scalar_inline {
    ($tags:expr, $_tags_ptr:expr, $start:expr, $end:expr, $tag:expr) => {
        probe_sorted_tags_scalar($tags, $start, $end, $tag)
    };
}

/// Generate a complete chunk scanner + 5-stage Bloom/CSR prefetched batch flusher + zero-`vzeroupper`
/// pure-scalar Stage-3c miss drainer + cold AVX2 slow extension handler.
macro_rules! impl_stream_tier {
    (
        $scan_hash:ident,
        $scan_direct:ident,
        $flush_hash:ident,
        $flush_direct:ident,
        $drain_misses:ident,
        $miss_slow:ident,
        [$($hash_attr:meta),*],
        [$($scalar_attr:meta),*],
        $probe_macro:ident,
        $pack_mm32:ident,
        $extend_call:expr
    ) => {
        $(#[$hash_attr])*
        #[allow(clippy::too_many_arguments)]
        unsafe fn $scan_hash<
            const WL32: bool,
            const PART_FILTER: bool,
            const SPACED: bool,
            const FWD_ONLY: bool,
            P: PsRead,
            S: DiagSlotOps,
        >(
            query: &Genome,
            reference: &Genome,
            rev_ref: Option<&Genome>,
            offsets: &[u32],
            tags: &[u32],
            ps_ptr: P,
            bloom: Option<&BlockedBloom>,
            unique_bm: &[u64],
            sub_bits: u32,
            hash_shift: u32,
            part_start: u32,
            part_end: u32,
            chunk_start: usize,
            chunk_end: usize,
            wl: usize,
            ny: u64,
            mask: u64,
            seed_mask: u64,
            rc_shift: u32,
            soft_mask: bool,
            max_rep: usize,
            fwd_opts: ExtendOptions,
            rev_opts: ExtendOptions,
            state: &mut WorkerState<S>,
        ) {
            let seq_x = &query.seq;
            let seq_y = &reference.seq;
            let rev_seq_y = rev_ref.map(|r| &r.seq);
            let max_y_start = ny - wl as u64;
            let scan_end = (chunk_end + wl - 1).min(seq_x.len);
            let part_bucket_base = part_start << sub_bits;
            let bloom_shift = bloom.map_or(56u32, |b| b.block_shift);
            let bloom_base = bloom.map_or(0u32, |b| part_start << b.sub_bits);
            let mut batch = [QueryItem::default(); BATCH_SIZE];
            let mut batch_len = 0usize;

            seq_x.for_each_valid_run(chunk_start, scan_end, soft_mask, |rs, re| {
                if re - rs < wl {
                    return;
                }
                let mut kf = 0u64;
                let mut kr = 0u64;
                let mut idx = rs;
                let prime_end = rs + wl - 1;
                while idx < prime_end {
                    let w_idx = idx >> 5;
                    let lane = idx & 31;
                    let avail = (32 - lane).min(prime_end - idx);
                    let mut word = unsafe { *seq_x.bits.get_unchecked(w_idx) } >> (lane << 1);
                    for _ in 0..avail {
                        let b = word & 3;
                        word >>= 2;
                        kf = if WL32 { (kf << 2) | b } else { ((kf << 2) | b) & mask };
                        if !FWD_ONLY {
                            kr = (kr >> 2) | ((b ^ 3) << if WL32 { 62 } else { rc_shift });
                        }
                    }
                    idx += avail;
                }

                while idx < re {
                    let w_idx = idx >> 5;
                    let lane = idx & 31;
                    let avail = (32 - lane).min(re - idx);
                    let mut word = unsafe { *seq_x.bits.get_unchecked(w_idx) } >> (lane << 1);
                    let start_pos = (idx + 1 - wl) as u64;
                    for pos_x in start_pos..start_pos + avail as u64 {
                        let b = word & 3;
                        word >>= 2;
                        kf = if WL32 { (kf << 2) | b } else { ((kf << 2) | b) & mask };
                        if !FWD_ONLY {
                            kr = (kr >> 2) | ((b ^ 3) << if WL32 { 62 } else { rc_shift });
                        }
                        if unsafe {
                            try_active_diag_bypass::<true, S>(
                                pos_x,
                                seq_x,
                                seq_y,
                                rev_seq_y,
                                unique_bm,
                                max_y_start,
                                seed_mask,
                                state,
                            )
                        } {
                            continue;
                        }
                        let (canon, q_kind) = if FWD_ONLY {
                            let sf = if SPACED { kf & seed_mask } else { kf };
                            (sf, 0u8)
                        } else {
                            let (sf, sr) = if SPACED {
                                (kf & seed_mask, kr & seed_mask)
                            } else {
                                (kf, kr)
                            };
                            if sf < sr {
                                (sf, 0u8)
                            } else if sr < sf {
                                (sr, 1u8)
                            } else {
                                (sf, 2u8)
                            }
                        };
                        let h = fmix64(canon);
                        if PART_FILTER {
                            let p = (h >> 56) as u32;
                            if p < part_start || p >= part_end {
                                continue;
                            }
                        }
                        let b_global = (h >> hash_shift) as u32;
                        let b_bucket = if PART_FILTER {
                            b_global - part_bucket_base
                        } else {
                            b_global
                        };
                        let b_idx = b_bucket + (b_bucket >> sub_bits);
                        let bloom_blk = if PART_FILTER {
                            ((h >> bloom_shift) as u32) - bloom_base
                        } else {
                            (h >> bloom_shift) as u32
                        };
                        unsafe {
                            *batch.get_unchecked_mut(batch_len) = QueryItem {
                                pos_x,
                                b_idx,
                                tag: h as u32,
                                start: 0,
                                end: 0,
                                bloom_blk,
                                q_kind,
                            };
                        }
                        batch_len += 1;
                        if batch_len == BATCH_SIZE {
                            unsafe {
                                $flush_hash::<FWD_ONLY, P, S>(
                                    &mut batch[..BATCH_SIZE],
                                    query,
                                    reference,
                                    rev_ref,
                                    offsets,
                                    tags,
                                    ps_ptr,
                                    bloom,
                                    wl as u64,
                                    seed_mask,
                                    ny,
                                    max_rep,
                                    fwd_opts,
                                    rev_opts,
                                    state,
                                );
                            }
                            batch_len = 0;
                        }
                    }
                    idx += avail;
                }
            });

            if batch_len > 0 {
                unsafe {
                    $flush_hash::<FWD_ONLY, P, S>(
                        &mut batch[..batch_len],
                        query,
                        reference,
                        rev_ref,
                        offsets,
                        tags,
                        ps_ptr,
                        bloom,
                        wl as u64,
                        seed_mask,
                        ny,
                        max_rep,
                        fwd_opts,
                        rev_opts,
                        state,
                    );
                }
            }
        }

        $(#[$scalar_attr])*
        #[allow(clippy::too_many_arguments)]
        unsafe fn $scan_direct<
            const PART_FILTER: bool,
            const FWD_ONLY: bool,
            P: PsRead,
            S: DiagSlotOps,
        >(
            query: &Genome,
            reference: &Genome,
            rev_ref: Option<&Genome>,
            offsets: &[u32],
            ps_ptr: P,
            unique_bm: &[u64],
            sub_bits: u32,
            part_start: u32,
            part_end: u32,
            chunk_start: usize,
            chunk_end: usize,
            wl: usize,
            ny: u64,
            mask: u64,
            rc_shift: u32,
            soft_mask: bool,
            max_rep: usize,
            fwd_opts: ExtendOptions,
            rev_opts: ExtendOptions,
            state: &mut WorkerState<S>,
        ) {
            let seq_x = &query.seq;
            let seq_y = &reference.seq;
            let rev_seq_y = rev_ref.map(|r| &r.seq);
            let max_y_start = ny - wl as u64;
            let scan_end = (chunk_end + wl - 1).min(seq_x.len);
            let part_bucket_base = part_start << sub_bits;
            let mut batch = [QueryItem::default(); BATCH_SIZE];
            let mut batch_len = 0usize;

            seq_x.for_each_valid_run(chunk_start, scan_end, soft_mask, |rs, re| {
                if re - rs < wl {
                    return;
                }
                let mut kf = 0u64;
                let mut kr = 0u64;
                let mut idx = rs;
                let prime_end = rs + wl - 1;
                while idx < prime_end {
                    let w_idx = idx >> 5;
                    let lane = idx & 31;
                    let avail = (32 - lane).min(prime_end - idx);
                    let mut word = unsafe { *seq_x.bits.get_unchecked(w_idx) } >> (lane << 1);
                    for _ in 0..avail {
                        let b = word & 3;
                        word >>= 2;
                        kf = ((kf << 2) | b) & mask;
                        if !FWD_ONLY {
                            kr = (kr >> 2) | ((b ^ 3) << rc_shift);
                        }
                    }
                    idx += avail;
                }

                while idx < re {
                    let w_idx = idx >> 5;
                    let lane = idx & 31;
                    let avail = (32 - lane).min(re - idx);
                    let mut word = unsafe { *seq_x.bits.get_unchecked(w_idx) } >> (lane << 1);
                    let start_pos = (idx + 1 - wl) as u64;
                    for pos_x in start_pos..start_pos + avail as u64 {
                        let b = word & 3;
                        word >>= 2;
                        kf = ((kf << 2) | b) & mask;
                        if !FWD_ONLY {
                            kr = (kr >> 2) | ((b ^ 3) << rc_shift);
                        }
                        if unsafe {
                            try_active_diag_bypass::<false, S>(
                                pos_x,
                                seq_x,
                                seq_y,
                                rev_seq_y,
                                unique_bm,
                                max_y_start,
                                mask,
                                state,
                            )
                        } {
                            continue;
                        }
                        let (canon, q_kind) = if FWD_ONLY {
                            (kf, 0u8)
                        } else if kf < kr {
                            (kf, 0u8)
                        } else if kr < kf {
                            (kr, 1u8)
                        } else {
                            (kf, 2u8)
                        };
                        let b_global = canon as u32;
                        if PART_FILTER {
                            let p = b_global >> sub_bits;
                            if p < part_start || p >= part_end {
                                continue;
                            }
                        }
                        let b_bucket = if PART_FILTER {
                            b_global - part_bucket_base
                        } else {
                            b_global
                        };
                        let b_idx = b_bucket + (b_bucket >> sub_bits);
                        unsafe {
                            *batch.get_unchecked_mut(batch_len) = QueryItem {
                                pos_x,
                                b_idx,
                                tag: 0,
                                start: 0,
                                end: 0,
                                bloom_blk: 0,
                                q_kind,
                            };
                        }
                        batch_len += 1;
                        if batch_len == BATCH_SIZE {
                            unsafe {
                                $flush_direct::<FWD_ONLY, P, S>(
                                    &mut batch[..BATCH_SIZE],
                                    query,
                                    reference,
                                    rev_ref,
                                    offsets,
                                    ps_ptr,
                                    wl as u64,
                                    mask,
                                    ny,
                                    max_rep,
                                    fwd_opts,
                                    rev_opts,
                                    state,
                                );
                            }
                            batch_len = 0;
                        }
                    }
                    idx += avail;
                }
            });

            if batch_len > 0 {
                unsafe {
                    $flush_direct::<FWD_ONLY, P, S>(
                        &mut batch[..batch_len],
                        query,
                        reference,
                        rev_ref,
                        offsets,
                        ps_ptr,
                        wl as u64,
                        mask,
                        ny,
                        max_rep,
                        fwd_opts,
                        rev_opts,
                        state,
                    );
                }
            }
        }

        $(#[$hash_attr])*
        #[inline(never)]
        #[allow(clippy::too_many_arguments)]
        unsafe fn $flush_hash<const FWD_ONLY: bool, P: PsRead, S: DiagSlotOps>(
            batch: &mut [QueryItem],
            query: &Genome,
            reference: &Genome,
            rev_ref: Option<&Genome>,
            offsets: &[u32],
            tags: &[u32],
            ps_ptr: P,
            bloom: Option<&BlockedBloom>,
            wl: u64,
            seed_mask: u64,
            ny: u64,
            max_rep: usize,
            fwd_opts: ExtendOptions,
            rev_opts: ExtendOptions,
            state: &mut WorkerState<S>,
        ) {
            let off_ptr = offsets.as_ptr();
            let tags_ptr = tags.as_ptr();
            let mut survivors = [0u8; BATCH_SIZE];
            let n_surv: usize;

            if let Some(bloom) = bloom {
                let bloom_ptr = bloom.blocks.as_ptr();
                // Stage 0: Prefetch 64 L3-resident Bloom cache lines into L2 in parallel.
                for item in batch.iter() {
                    prefetch_t1(unsafe { bloom_ptr.add(item.bloom_blk as usize) });
                }
                // Stage 1: Branchless Bloom test; only survivors prefetch `offsets[b_idx]`.
                let mut ns = 0usize;
                for (i, item) in batch.iter().enumerate() {
                    if unsafe { bloom.contains_unchecked(item.bloom_blk as usize, item.tag) } {
                        prefetch_t2(unsafe { off_ptr.add(item.b_idx as usize) });
                        unsafe { *survivors.get_unchecked_mut(ns) = i as u8 };
                        ns += 1;
                    }
                }
                state.bloom_rejected += (batch.len() - ns) as u64;
                if ns == 0 {
                    return;
                }
                n_surv = ns;
            } else {
                for (i, item) in batch.iter().enumerate() {
                    prefetch_t2(unsafe { off_ptr.add(item.b_idx as usize) });
                    unsafe { *survivors.get_unchecked_mut(i) = i as u8 };
                }
                n_surv = batch.len();
            }

            // Stage 2: Load `[start, end)` and prefetch `tags[start]` ONLY for Bloom survivors.
            for &idx in &survivors[..n_surv] {
                let item = unsafe { batch.get_unchecked_mut(idx as usize) };
                let b_idx = item.b_idx as usize;
                let s = unsafe { *off_ptr.add(b_idx) };
                let e = unsafe { *off_ptr.add(b_idx + 1) };
                item.start = s;
                item.end = e;
                if s < e {
                    prefetch_t2(unsafe { tags_ptr.add(s as usize) });
                }
            }

            // Stage 3a: Probe 32-bit tag buckets ONLY for survivors, compacting hitting items.
            let mut hits = [0u8; BATCH_SIZE];
            let mut n_hits = 0usize;
            for &idx in &survivors[..n_surv] {
                let item = unsafe { batch.get_unchecked_mut(idx as usize) };
                let s = item.start as usize;
                let e = item.end as usize;
                if s == e {
                    continue;
                }
                let (lo, hi) = $probe_macro!(tags, tags_ptr, s, e, item.tag);
                item.start = lo as u32;
                item.end = hi as u32;
                if lo < hi {
                    unsafe {
                        ps_ptr.prefetch_at(lo);
                        *hits.get_unchecked_mut(n_hits) = idx;
                    }
                    n_hits += 1;
                }
            }
            if n_hits == 0 {
                return;
            }

            // Stage 3b + 3c: Unpack prefetched `pos_strand`, filter via diagonal cache, prefetch
            // `seq_y.bits` for uncovered misses into `miss_buf`, and drain in FIFO order.
            let y_fwd_ptr = reference.seq.bits.as_ptr();
            let y_rev_ptr = rev_ref.map_or(core::ptr::null(), |r| r.seq.bits.as_ptr());
            let mut miss_buf = [MissItem::default(); MISS_BUF_CAP];
            let mut miss_len = 0usize;
            let collect_cov = state.coverage.is_some();

            macro_rules! flush_misses_hash {
                ($len:expr) => {
                    if collect_cov {
                        unsafe {
                            $drain_misses::<true, true, S>(
                                &miss_buf[..$len],
                                query,
                                reference,
                                rev_ref,
                                wl,
                                seed_mask,
                                fwd_opts,
                                rev_opts,
                                state,
                            );
                        }
                    } else {
                        unsafe {
                            $drain_misses::<true, false, S>(
                                &miss_buf[..$len],
                                query,
                                reference,
                                rev_ref,
                                wl,
                                seed_mask,
                                fwd_opts,
                                rev_opts,
                                state,
                            );
                        }
                    }
                };
            }

            for &hit_idx in &hits[..n_hits] {
                let item = unsafe { *batch.get_unchecked(hit_idx as usize) };
                let lo = item.start as usize;
                let hi = item.end as usize;
                let total = hi - lo;
                let step = if total > max_rep {
                    (total / max_rep).max(1)
                } else {
                    1
                };
                let pos_x = item.pos_x;
                let mut idx = lo;

                if FWD_ONLY {
                    state.word_matches += 1;
                    state.raw_hits += ((total - 1) / step + 1) as u64;
                    while idx < hi {
                        let (pos_y, _) = unsafe { ps_ptr.unpack_at(idx) };
                        idx += step;
                        let diag_f = pos_x as i64 - pos_y as i64;
                        if !state.fwd_cache.is_covered(diag_f, pos_x) {
                            prefetch_t2(unsafe { y_fwd_ptr.add((pos_y as usize) >> 5) });
                            unsafe {
                                *miss_buf.get_unchecked_mut(miss_len) = MissItem {
                                    pos_x,
                                    pos_y_rev: pos_y,
                                };
                            }
                            miss_len += 1;
                            if miss_len == MISS_BUF_CAP {
                                flush_misses_hash!(MISS_BUF_CAP);
                                miss_len = 0;
                            }
                        }
                    }
                } else {
                    let q_kind = item.q_kind;
                    let mut matched_any = false;
                    while idx < hi {
                        let (pos_y, r_strand) = unsafe { ps_ptr.unpack_at(idx) };
                        idx += step;

                        if q_kind == 2 {
                            matched_any = true;
                            state.raw_hits += 1;
                            let diag_f = pos_x as i64 - pos_y as i64;
                            if !state.fwd_cache.is_covered(diag_f, pos_x) {
                                prefetch_t2(unsafe { y_fwd_ptr.add((pos_y as usize) >> 5) });
                                unsafe {
                                    *miss_buf.get_unchecked_mut(miss_len) = MissItem {
                                        pos_x,
                                        pos_y_rev: pos_y,
                                    };
                                }
                                miss_len += 1;
                                if miss_len == MISS_BUF_CAP {
                                    flush_misses_hash!(MISS_BUF_CAP);
                                    miss_len = 0;
                                }
                            }
                            if rev_ref.is_some() {
                                state.raw_hits += 1;
                                let rev_pos_y = ny - pos_y - wl;
                                let diag_r = pos_x as i64 - rev_pos_y as i64;
                                if !state.rev_cache.as_ref().unwrap().is_covered(diag_r, pos_x) {
                                    prefetch_t2(unsafe {
                                        y_rev_ptr.add((rev_pos_y as usize) >> 5)
                                    });
                                    unsafe {
                                        *miss_buf.get_unchecked_mut(miss_len) = MissItem {
                                            pos_x,
                                            pos_y_rev: rev_pos_y | STRAND_BIT_MASK,
                                        };
                                    }
                                    miss_len += 1;
                                    if miss_len == MISS_BUF_CAP {
                                        flush_misses_hash!(MISS_BUF_CAP);
                                        miss_len = 0;
                                    }
                                }
                            }
                        } else if r_strand == q_kind {
                            matched_any = true;
                            state.raw_hits += 1;
                            let diag_f = pos_x as i64 - pos_y as i64;
                            if !state.fwd_cache.is_covered(diag_f, pos_x) {
                                prefetch_t2(unsafe { y_fwd_ptr.add((pos_y as usize) >> 5) });
                                unsafe {
                                    *miss_buf.get_unchecked_mut(miss_len) = MissItem {
                                        pos_x,
                                        pos_y_rev: pos_y,
                                    };
                                }
                                miss_len += 1;
                                if miss_len == MISS_BUF_CAP {
                                    flush_misses_hash!(MISS_BUF_CAP);
                                    miss_len = 0;
                                }
                            }
                        } else if rev_ref.is_some() {
                            matched_any = true;
                            state.raw_hits += 1;
                            let rev_pos_y = ny - pos_y - wl;
                            let diag_r = pos_x as i64 - rev_pos_y as i64;
                            if !state.rev_cache.as_ref().unwrap().is_covered(diag_r, pos_x) {
                                prefetch_t2(unsafe { y_rev_ptr.add((rev_pos_y as usize) >> 5) });
                                unsafe {
                                    *miss_buf.get_unchecked_mut(miss_len) = MissItem {
                                        pos_x,
                                        pos_y_rev: rev_pos_y | STRAND_BIT_MASK,
                                    };
                                }
                                miss_len += 1;
                                if miss_len == MISS_BUF_CAP {
                                    flush_misses_hash!(MISS_BUF_CAP);
                                    miss_len = 0;
                                }
                            }
                        }
                    }

                    if matched_any {
                        state.word_matches += 1;
                    }
                }
            }

            if miss_len > 0 {
                flush_misses_hash!(miss_len);
            }
        }

        $(#[$scalar_attr])*
        #[inline(never)]
        #[allow(clippy::too_many_arguments)]
        unsafe fn $flush_direct<const FWD_ONLY: bool, P: PsRead, S: DiagSlotOps>(
            batch: &mut [QueryItem],
            query: &Genome,
            reference: &Genome,
            rev_ref: Option<&Genome>,
            offsets: &[u32],
            ps_ptr: P,
            wl: u64,
            seed_mask: u64,
            ny: u64,
            max_rep: usize,
            fwd_opts: ExtendOptions,
            rev_opts: ExtendOptions,
            state: &mut WorkerState<S>,
        ) {
            let off_ptr = offsets.as_ptr();

            for item in batch.iter() {
                prefetch_t2(unsafe { off_ptr.add(item.b_idx as usize) });
            }

            for item in batch.iter_mut() {
                let b_idx = item.b_idx as usize;
                let s = unsafe { *off_ptr.add(b_idx) };
                let e = unsafe { *off_ptr.add(b_idx + 1) };
                item.start = s;
                item.end = e;
                if s < e {
                    unsafe { ps_ptr.prefetch_at(s as usize) };
                }
            }

            let y_fwd_ptr = reference.seq.bits.as_ptr();
            let y_rev_ptr = rev_ref.map_or(core::ptr::null(), |r| r.seq.bits.as_ptr());
            let mut miss_buf = [MissItem::default(); MISS_BUF_CAP];
            let mut miss_len = 0usize;
            let collect_cov = state.coverage.is_some();

            macro_rules! flush_misses_direct {
                ($len:expr) => {
                    if collect_cov {
                        unsafe {
                            $drain_misses::<false, true, S>(
                                &miss_buf[..$len],
                                query,
                                reference,
                                rev_ref,
                                wl,
                                seed_mask,
                                fwd_opts,
                                rev_opts,
                                state,
                            );
                        }
                    } else {
                        unsafe {
                            $drain_misses::<false, false, S>(
                                &miss_buf[..$len],
                                query,
                                reference,
                                rev_ref,
                                wl,
                                seed_mask,
                                fwd_opts,
                                rev_opts,
                                state,
                            );
                        }
                    }
                };
            }

            for item in batch.iter() {
                let lo = item.start as usize;
                let hi = item.end as usize;
                if lo == hi {
                    continue;
                }
                let total = hi - lo;
                let step = if total > max_rep {
                    (total / max_rep).max(1)
                } else {
                    1
                };
                let pos_x = item.pos_x;
                let mut idx = lo;

                if FWD_ONLY {
                    state.word_matches += 1;
                    state.raw_hits += ((total - 1) / step + 1) as u64;
                    while idx < hi {
                        let (pos_y, _) = unsafe { ps_ptr.unpack_at(idx) };
                        idx += step;
                        let diag_f = pos_x as i64 - pos_y as i64;
                        if !state.fwd_cache.is_covered(diag_f, pos_x) {
                            prefetch_t2(unsafe { y_fwd_ptr.add((pos_y as usize) >> 5) });
                            unsafe {
                                *miss_buf.get_unchecked_mut(miss_len) = MissItem {
                                    pos_x,
                                    pos_y_rev: pos_y,
                                };
                            }
                            miss_len += 1;
                            if miss_len == MISS_BUF_CAP {
                                flush_misses_direct!(MISS_BUF_CAP);
                                miss_len = 0;
                            }
                        }
                    }
                } else {
                    let q_kind = item.q_kind;
                    let mut matched_any = false;
                    while idx < hi {
                        let (pos_y, r_strand) = unsafe { ps_ptr.unpack_at(idx) };
                        idx += step;

                        if q_kind == 2 {
                            matched_any = true;
                            state.raw_hits += 1;
                            let diag_f = pos_x as i64 - pos_y as i64;
                            if !state.fwd_cache.is_covered(diag_f, pos_x) {
                                prefetch_t2(unsafe { y_fwd_ptr.add((pos_y as usize) >> 5) });
                                unsafe {
                                    *miss_buf.get_unchecked_mut(miss_len) = MissItem {
                                        pos_x,
                                        pos_y_rev: pos_y,
                                    };
                                }
                                miss_len += 1;
                                if miss_len == MISS_BUF_CAP {
                                    flush_misses_direct!(MISS_BUF_CAP);
                                    miss_len = 0;
                                }
                            }
                            if rev_ref.is_some() {
                                state.raw_hits += 1;
                                let rev_pos_y = ny - pos_y - wl;
                                let diag_r = pos_x as i64 - rev_pos_y as i64;
                                if !state.rev_cache.as_ref().unwrap().is_covered(diag_r, pos_x) {
                                    prefetch_t2(unsafe {
                                        y_rev_ptr.add((rev_pos_y as usize) >> 5)
                                    });
                                    unsafe {
                                        *miss_buf.get_unchecked_mut(miss_len) = MissItem {
                                            pos_x,
                                            pos_y_rev: rev_pos_y | STRAND_BIT_MASK,
                                        };
                                    }
                                    miss_len += 1;
                                    if miss_len == MISS_BUF_CAP {
                                        flush_misses_direct!(MISS_BUF_CAP);
                                        miss_len = 0;
                                    }
                                }
                            }
                        } else if r_strand == q_kind {
                            matched_any = true;
                            state.raw_hits += 1;
                            let diag_f = pos_x as i64 - pos_y as i64;
                            if !state.fwd_cache.is_covered(diag_f, pos_x) {
                                prefetch_t2(unsafe { y_fwd_ptr.add((pos_y as usize) >> 5) });
                                unsafe {
                                    *miss_buf.get_unchecked_mut(miss_len) = MissItem {
                                        pos_x,
                                        pos_y_rev: pos_y,
                                    };
                                }
                                miss_len += 1;
                                if miss_len == MISS_BUF_CAP {
                                    flush_misses_direct!(MISS_BUF_CAP);
                                    miss_len = 0;
                                }
                            }
                        } else if rev_ref.is_some() {
                            matched_any = true;
                            state.raw_hits += 1;
                            let rev_pos_y = ny - pos_y - wl;
                            let diag_r = pos_x as i64 - rev_pos_y as i64;
                            if !state.rev_cache.as_ref().unwrap().is_covered(diag_r, pos_x) {
                                prefetch_t2(unsafe { y_rev_ptr.add((rev_pos_y as usize) >> 5) });
                                unsafe {
                                    *miss_buf.get_unchecked_mut(miss_len) = MissItem {
                                        pos_x,
                                        pos_y_rev: rev_pos_y | STRAND_BIT_MASK,
                                    };
                                }
                                miss_len += 1;
                                if miss_len == MISS_BUF_CAP {
                                    flush_misses_direct!(MISS_BUF_CAP);
                                    miss_len = 0;
                                }
                            }
                        }
                    }

                    if matched_any {
                        state.word_matches += 1;
                    }
                }
            }

            if miss_len > 0 {
                flush_misses_direct!(miss_len);
            }
        }

        /// Pure-scalar (NO AVX2, 0 `%ymm`, 0 `vzeroupper`) Stage-3c miss drainer.
        ///
        /// Re-checks `is_covered` for same-diagonal duplicates in `miss_buf`, verifies exact
        /// 2-bit seeds only on `HashCsr` (`VERIFY_SEED == true`), reuses the seed XOR for spaced
        /// seeds, rejects 97%-99.8% of noise hits inline via `try_fast_short_reject_64`, and
        /// isolates full AVX2 extension + `Fragment` push inside cold `$miss_slow`.
        $(#[$scalar_attr])*
        #[inline(never)]
        #[allow(clippy::too_many_arguments)]
        unsafe fn $drain_misses<
            const VERIFY_SEED: bool,
            const COLLECT_COV: bool,
            S: DiagSlotOps,
        >(
            miss_buf: &[MissItem],
            query: &Genome,
            reference: &Genome,
            rev_ref: Option<&Genome>,
            wl: u64,
            seed_mask: u64,
            fwd_opts: ExtendOptions,
            rev_opts: ExtendOptions,
            state: &mut WorkerState<S>,
        ) {
            let seq_x = &query.seq;
            let nx = seq_x.len as u64;
            let default_fixed = fwd_opts.fixed_length && fwd_opts.scoring.is_default_gecko();
            let noise_add = if fwd_opts.spaced_seed {
                wl - 1
            } else {
                2 * wl - 1
            };
            let min_len = fwd_opts.min_length;
            let contig_seed_score = wl as i64 * fwd_opts.scoring.match_score;

            for &m in miss_buf {
                let pos_x = m.pos_x;
                let is_rev = (m.pos_y_rev & STRAND_BIT_MASK) != 0;
                let pos_y = m.pos_y_rev & !STRAND_BIT_MASK;
                let diag = pos_x as i64 - pos_y as i64;

                // Mandatory re-check: an earlier miss in `miss_buf` on the same diagonal may have
                // just updated `fwd_cache` / `rev_cache` to cover `pos_x`.
                if !is_rev {
                    if state.fwd_cache.is_covered(diag, pos_x) {
                        continue;
                    }
                } else if state.rev_cache.as_ref().unwrap().is_covered(diag, pos_x) {
                    continue;
                }

                let (target, opts) = if !is_rev {
                    (reference, fwd_opts)
                } else {
                    (rev_ref.unwrap(), rev_opts)
                };
                let seq_y = &target.seq;

                let (seed_ident, seed_score) = if VERIFY_SEED {
                    let seed_diff = seq_x.load_32(pos_x as usize) ^ seq_y.load_32(pos_y as usize);
                    if (seed_diff & seed_mask) != 0 {
                        continue;
                    }
                    seed_init_stats_from_diff(seed_diff, wl, opts.scoring, opts.spaced_seed)
                } else {
                    (wl, contig_seed_score)
                };

                if state.two_hit_window > 0 {
                    let th_ok = if !is_rev {
                        state
                            .fwd_two_hit
                            .as_mut()
                            .unwrap()
                            .accept(diag, pos_x, wl, state.two_hit_window as u64)
                    } else {
                        state
                            .rev_two_hit
                            .as_mut()
                            .unwrap()
                            .accept(diag, pos_x, wl, state.two_hit_window as u64)
                    };
                    if !th_ok {
                        state.two_hit_dropped += 1;
                        continue;
                    }
                }

                state.hits_used += 1;

                let ny = seq_y.len as u64;
                if default_fixed
                    && (0..128).contains(&seed_score)
                    && pos_x >= 64
                    && pos_y >= 64
                    && pos_x + wl + 64 <= nx
                    && pos_y + wl + 64 <= ny
                    && !seq_x
                        .has_special_in_short_range((pos_x - 64) as usize, (pos_x + wl + 64) as usize)
                    && !seq_y
                        .has_special_in_short_range((pos_y - 64) as usize, (pos_y + wl + 64) as usize)
                    && let Some((length, u_max)) = unsafe {
                        try_fast_short_reject_64(
                            seq_x,
                            seq_y,
                            pos_x,
                            pos_y,
                            wl,
                            seed_score,
                            min_len,
                            |mm_even| $pack_mm32(mm_even),
                        )
                    }
                {
                    if COLLECT_COV && let Some(cov) = state.coverage.as_mut() {
                        cov.record_unit_exact(length, u_max);
                    }
                    state.fast_rejects += 1;
                    let skip = pos_x + noise_add;
                    if !is_rev {
                        state.fwd_cache.update_noise(diag, skip);
                    } else {
                        state.rev_cache.as_mut().unwrap().update_noise(diag, skip);
                    }
                    continue;
                }

                unsafe {
                    $miss_slow::<COLLECT_COV, S>(
                        is_rev,
                        diag,
                        pos_x,
                        pos_y,
                        seed_ident,
                        seed_score,
                        query,
                        target,
                        opts,
                        state,
                    );
                }
            }
        }

        $(#[$hash_attr])*
        #[cold]
        #[inline(never)]
        #[allow(clippy::too_many_arguments)]
        unsafe fn $miss_slow<const COLLECT_COV: bool, S: DiagSlotOps>(
            is_rev: bool,
            diag: i64,
            pos_x: u64,
            pos_y: u64,
            seed_ident: u64,
            seed_score: i64,
            query: &Genome,
            target: &Genome,
            opts: ExtendOptions,
            state: &mut WorkerState<S>,
        ) {
            let cx_x = state.sp_x_cursor_for(&query.seq.special_runs, pos_x);
            let cov_opt = if COLLECT_COV {
                state.coverage.as_mut()
            } else {
                None
            };
            let outcome = ($extend_call)(
                Hit { diag, pos_x },
                query,
                target,
                opts,
                seed_ident,
                seed_score,
                cx_x,
                cov_opt,
            );
            let strand = is_rev as u8;
            if let Some(frag) = outcome.fragment {
                if !is_rev {
                    state.fwd_cache.update_hsp(diag, outcome.xfil_max);
                } else {
                    state
                        .rev_cache
                        .as_mut()
                        .unwrap()
                        .update_hsp(diag, outcome.xfil_max);
                }
                if state.has_unique_bm {
                    state.active_ring.record(diag, outcome.xfil_max, strand);
                }
                state.fragments.push(frag);
            } else if outcome.xfil_max > pos_x + 256 {
                let skip_end = exact_rejected_skip(
                    &query.seq,
                    &target.seq,
                    pos_x,
                    pos_y,
                    outcome.xfil_max,
                    opts.word_length,
                    opts.scoring,
                    &mut state.skip_scratch,
                );
                if !is_rev {
                    state.fwd_cache.update_hsp(diag, skip_end);
                } else {
                    state
                        .rev_cache
                        .as_mut()
                        .unwrap()
                        .update_hsp(diag, skip_end);
                }
                if state.has_unique_bm {
                    state.active_ring.record(diag, skip_end, strand);
                }
            } else {
                let skip = if opts.spaced_seed {
                    pos_x + opts.word_length as u64 - 1
                } else {
                    pos_x + 2 * opts.word_length as u64 - 1
                };
                if !is_rev {
                    state.fwd_cache.update_noise(diag, skip);
                } else {
                    state.rev_cache.as_mut().unwrap().update_noise(diag, skip);
                }
            }
        }
    };
}

#[cfg(target_arch = "x86_64")]
impl_stream_tier!(
    scan_hash_avx2_bmi2,
    scan_direct_bmi2_popcnt,
    flush_hash_avx2_bmi2,
    flush_direct_bmi2_popcnt,
    drain_misses_bmi2_popcnt,
    miss_slow_avx2_bmi2,
    [target_feature(enable = "avx2,bmi2,popcnt")],
    [target_feature(enable = "bmi2,popcnt")],
    probe_bucket_avx2_inline,
    pack_mm32_bmi2,
    |h, gx, gy, opts, si, ss, cx_x, cov| unsafe {
        extend_single_hit_with_seed_avx2_bmi2(h, gx, gy, opts, si, ss, cx_x, cov)
    }
);

#[cfg(target_arch = "x86_64")]
impl_stream_tier!(
    scan_hash_avx2_swar,
    scan_direct_popcnt,
    flush_hash_avx2_swar,
    flush_direct_popcnt,
    drain_misses_popcnt,
    miss_slow_avx2_swar,
    [target_feature(enable = "avx2,popcnt")],
    [target_feature(enable = "popcnt")],
    probe_bucket_avx2_inline,
    pack_mm32_swar,
    |h, gx, gy, opts, si, ss, cx_x, cov| unsafe {
        extend_single_hit_with_seed_avx2_swar(h, gx, gy, opts, si, ss, cx_x, cov)
    }
);

impl_stream_tier!(
    scan_hash_portable,
    scan_direct_portable,
    flush_hash_portable,
    flush_direct_portable,
    drain_misses_portable,
    miss_slow_portable,
    [],
    [],
    probe_bucket_scalar_inline,
    pack_mm32_swar,
    |h, gx, gy, opts, si, ss, _cx_x, cov| {
        extend_single_hit_with_seed_portable(h, gx, gy, opts, si, ss, cov)
    }
);

macro_rules! dispatch_hash_tier_inner {
    (
        $fn_name:ident, $wl32:expr, $part_filter:expr, $spaced_const:expr, $fwd_const:expr,
        $ps_ptr:expr, $bloom:expr, $unique_bm:expr,
        $query:expr, $reference:expr, $rev_ref:expr, $offsets:expr, $tags:expr,
        $sub_bits:expr, $hash_shift:expr, $p_start:expr, $p_end:expr,
        $chunk_start:expr, $chunk_end:expr, $wl:expr, $ny:expr, $mask:expr, $seed_mask:expr,
        $rc_shift:expr, $soft_mask:expr, $max_rep:expr, $fwd_opts:expr, $rev_opts:expr, $state:expr
    ) => {
        match ($wl32, $part_filter) {
            (true, false) => unsafe {
                $fn_name::<true, false, $spaced_const, $fwd_const, _, _>(
                    $query,
                    $reference,
                    $rev_ref,
                    $offsets,
                    $tags,
                    $ps_ptr,
                    $bloom,
                    $unique_bm,
                    $sub_bits,
                    $hash_shift,
                    $p_start,
                    $p_end,
                    $chunk_start,
                    $chunk_end,
                    32,
                    $ny,
                    $mask,
                    $seed_mask,
                    $rc_shift,
                    $soft_mask,
                    $max_rep,
                    $fwd_opts,
                    $rev_opts,
                    $state,
                );
            },
            (true, true) => unsafe {
                $fn_name::<true, true, $spaced_const, $fwd_const, _, _>(
                    $query,
                    $reference,
                    $rev_ref,
                    $offsets,
                    $tags,
                    $ps_ptr,
                    $bloom,
                    $unique_bm,
                    $sub_bits,
                    $hash_shift,
                    $p_start,
                    $p_end,
                    $chunk_start,
                    $chunk_end,
                    32,
                    $ny,
                    $mask,
                    $seed_mask,
                    $rc_shift,
                    $soft_mask,
                    $max_rep,
                    $fwd_opts,
                    $rev_opts,
                    $state,
                );
            },
            (false, false) => unsafe {
                $fn_name::<false, false, $spaced_const, $fwd_const, _, _>(
                    $query,
                    $reference,
                    $rev_ref,
                    $offsets,
                    $tags,
                    $ps_ptr,
                    $bloom,
                    $unique_bm,
                    $sub_bits,
                    $hash_shift,
                    $p_start,
                    $p_end,
                    $chunk_start,
                    $chunk_end,
                    $wl,
                    $ny,
                    $mask,
                    $seed_mask,
                    $rc_shift,
                    $soft_mask,
                    $max_rep,
                    $fwd_opts,
                    $rev_opts,
                    $state,
                );
            },
            (false, true) => unsafe {
                $fn_name::<false, true, $spaced_const, $fwd_const, _, _>(
                    $query,
                    $reference,
                    $rev_ref,
                    $offsets,
                    $tags,
                    $ps_ptr,
                    $bloom,
                    $unique_bm,
                    $sub_bits,
                    $hash_shift,
                    $p_start,
                    $p_end,
                    $chunk_start,
                    $chunk_end,
                    $wl,
                    $ny,
                    $mask,
                    $seed_mask,
                    $rc_shift,
                    $soft_mask,
                    $max_rep,
                    $fwd_opts,
                    $rev_opts,
                    $state,
                );
            },
        }
    };
}

macro_rules! dispatch_hash_tier {
    (
        $fn_name:ident, $wl32:expr, $part_filter:expr, $is_spaced:expr, $fwd_only:expr,
        $ps_ptr:expr, $bloom:expr, $unique_bm:expr,
        $query:expr, $reference:expr, $rev_ref:expr, $offsets:expr, $tags:expr,
        $sub_bits:expr, $hash_shift:expr, $p_start:expr, $p_end:expr,
        $chunk_start:expr, $chunk_end:expr, $wl:expr, $ny:expr, $mask:expr, $seed_mask:expr,
        $rc_shift:expr, $soft_mask:expr, $max_rep:expr, $fwd_opts:expr, $rev_opts:expr, $state:expr
    ) => {
        match ($is_spaced, $fwd_only) {
            (false, false) => dispatch_hash_tier_inner!(
                $fn_name,
                $wl32,
                $part_filter,
                false,
                false,
                $ps_ptr,
                $bloom,
                $unique_bm,
                $query,
                $reference,
                $rev_ref,
                $offsets,
                $tags,
                $sub_bits,
                $hash_shift,
                $p_start,
                $p_end,
                $chunk_start,
                $chunk_end,
                $wl,
                $ny,
                $mask,
                $seed_mask,
                $rc_shift,
                $soft_mask,
                $max_rep,
                $fwd_opts,
                $rev_opts,
                $state
            ),
            (true, false) => dispatch_hash_tier_inner!(
                $fn_name,
                $wl32,
                $part_filter,
                true,
                false,
                $ps_ptr,
                $bloom,
                $unique_bm,
                $query,
                $reference,
                $rev_ref,
                $offsets,
                $tags,
                $sub_bits,
                $hash_shift,
                $p_start,
                $p_end,
                $chunk_start,
                $chunk_end,
                $wl,
                $ny,
                $mask,
                $seed_mask,
                $rc_shift,
                $soft_mask,
                $max_rep,
                $fwd_opts,
                $rev_opts,
                $state
            ),
            (false, true) => dispatch_hash_tier_inner!(
                $fn_name,
                $wl32,
                $part_filter,
                false,
                true,
                $ps_ptr,
                $bloom,
                $unique_bm,
                $query,
                $reference,
                $rev_ref,
                $offsets,
                $tags,
                $sub_bits,
                $hash_shift,
                $p_start,
                $p_end,
                $chunk_start,
                $chunk_end,
                $wl,
                $ny,
                $mask,
                $seed_mask,
                $rc_shift,
                $soft_mask,
                $max_rep,
                $fwd_opts,
                $rev_opts,
                $state
            ),
            (true, true) => dispatch_hash_tier_inner!(
                $fn_name,
                $wl32,
                $part_filter,
                true,
                true,
                $ps_ptr,
                $bloom,
                $unique_bm,
                $query,
                $reference,
                $rev_ref,
                $offsets,
                $tags,
                $sub_bits,
                $hash_shift,
                $p_start,
                $p_end,
                $chunk_start,
                $chunk_end,
                $wl,
                $ny,
                $mask,
                $seed_mask,
                $rc_shift,
                $soft_mask,
                $max_rep,
                $fwd_opts,
                $rev_opts,
                $state
            ),
        }
    };
}

macro_rules! dispatch_direct_tier {
    (
        $fn_name:ident, $part_filter:expr, $fwd_only:expr, $ps_ptr:expr, $unique_bm:expr,
        $query:expr, $reference:expr, $rev_ref:expr, $offsets:expr,
        $sub_bits:expr, $p_start:expr, $p_end:expr,
        $chunk_start:expr, $chunk_end:expr, $wl:expr, $ny:expr, $mask:expr,
        $rc_shift:expr, $soft_mask:expr, $max_rep:expr, $fwd_opts:expr, $rev_opts:expr, $state:expr
    ) => {
        match ($part_filter, $fwd_only) {
            (false, false) => unsafe {
                $fn_name::<false, false, _, _>(
                    $query,
                    $reference,
                    $rev_ref,
                    $offsets,
                    $ps_ptr,
                    $unique_bm,
                    $sub_bits,
                    $p_start,
                    $p_end,
                    $chunk_start,
                    $chunk_end,
                    $wl,
                    $ny,
                    $mask,
                    $rc_shift,
                    $soft_mask,
                    $max_rep,
                    $fwd_opts,
                    $rev_opts,
                    $state,
                );
            },
            (true, false) => unsafe {
                $fn_name::<true, false, _, _>(
                    $query,
                    $reference,
                    $rev_ref,
                    $offsets,
                    $ps_ptr,
                    $unique_bm,
                    $sub_bits,
                    $p_start,
                    $p_end,
                    $chunk_start,
                    $chunk_end,
                    $wl,
                    $ny,
                    $mask,
                    $rc_shift,
                    $soft_mask,
                    $max_rep,
                    $fwd_opts,
                    $rev_opts,
                    $state,
                );
            },
            (false, true) => unsafe {
                $fn_name::<false, true, _, _>(
                    $query,
                    $reference,
                    $rev_ref,
                    $offsets,
                    $ps_ptr,
                    $unique_bm,
                    $sub_bits,
                    $p_start,
                    $p_end,
                    $chunk_start,
                    $chunk_end,
                    $wl,
                    $ny,
                    $mask,
                    $rc_shift,
                    $soft_mask,
                    $max_rep,
                    $fwd_opts,
                    $rev_opts,
                    $state,
                );
            },
            (true, true) => unsafe {
                $fn_name::<true, true, _, _>(
                    $query,
                    $reference,
                    $rev_ref,
                    $offsets,
                    $ps_ptr,
                    $unique_bm,
                    $sub_bits,
                    $p_start,
                    $p_end,
                    $chunk_start,
                    $chunk_end,
                    $wl,
                    $ny,
                    $mask,
                    $rc_shift,
                    $soft_mask,
                    $max_rep,
                    $fwd_opts,
                    $rev_opts,
                    $state,
                );
            },
        }
    };
}

#[allow(clippy::too_many_arguments)]
fn scan_query_chunk<PS: PsSlice, S: DiagSlotOps>(
    query: &Genome,
    reference: &Genome,
    rev_ref: Option<&Genome>,
    ref_index: &RefIndex,
    ps_slice: PS,
    chunk_start: usize,
    chunk_end: usize,
    wl: usize,
    ny: u64,
    mask: u64,
    rc_shift: u32,
    soft_mask: bool,
    max_rep: usize,
    fwd_opts: ExtendOptions,
    rev_opts: ExtendOptions,
    isa: IsaTier,
    state: &mut WorkerState<S>,
) {
    let p_start = ref_index.part_start as u32;
    let p_end = ref_index.part_end as u32;
    let part_filter = p_start != 0 || p_end != 256;
    let is_spaced = ref_index.spaced_mask.is_some();
    let seed_mask = ref_index.spaced_mask.unwrap_or(mask);
    let fwd_only = ref_index.forward_only;
    let unique_bm = &ref_index.unique_bm;
    let ps_ptr = ps_slice.as_ps_ptr();
    debug_assert_eq!(
        ref_index.unique_bm_xor,
        if ref_index.is_direct() { 0 } else { u64::MAX }
    );

    match (&ref_index.kind, isa) {
        #[cfg(target_arch = "x86_64")]
        (
            RefIndexKind::HashCsr {
                sub_bits,
                hash_shift,
                offsets,
                tags,
                bloom,
                ..
            },
            IsaTier::Avx2Bmi2Popcnt,
        ) => dispatch_hash_tier!(
            scan_hash_avx2_bmi2,
            wl == 32,
            part_filter,
            is_spaced,
            fwd_only,
            ps_ptr,
            bloom.as_ref(),
            unique_bm,
            query,
            reference,
            rev_ref,
            offsets,
            tags,
            *sub_bits,
            *hash_shift,
            p_start,
            p_end,
            chunk_start,
            chunk_end,
            wl,
            ny,
            mask,
            seed_mask,
            rc_shift,
            soft_mask,
            max_rep,
            fwd_opts,
            rev_opts,
            state
        ),
        #[cfg(target_arch = "x86_64")]
        (
            RefIndexKind::HashCsr {
                sub_bits,
                hash_shift,
                offsets,
                tags,
                bloom,
                ..
            },
            IsaTier::Avx2Popcnt,
        ) => dispatch_hash_tier!(
            scan_hash_avx2_swar,
            wl == 32,
            part_filter,
            is_spaced,
            fwd_only,
            ps_ptr,
            bloom.as_ref(),
            unique_bm,
            query,
            reference,
            rev_ref,
            offsets,
            tags,
            *sub_bits,
            *hash_shift,
            p_start,
            p_end,
            chunk_start,
            chunk_end,
            wl,
            ny,
            mask,
            seed_mask,
            rc_shift,
            soft_mask,
            max_rep,
            fwd_opts,
            rev_opts,
            state
        ),
        (
            RefIndexKind::HashCsr {
                sub_bits,
                hash_shift,
                offsets,
                tags,
                bloom,
                ..
            },
            IsaTier::Portable,
        ) => dispatch_hash_tier!(
            scan_hash_portable,
            wl == 32,
            part_filter,
            is_spaced,
            fwd_only,
            ps_ptr,
            bloom.as_ref(),
            unique_bm,
            query,
            reference,
            rev_ref,
            offsets,
            tags,
            *sub_bits,
            *hash_shift,
            p_start,
            p_end,
            chunk_start,
            chunk_end,
            wl,
            ny,
            mask,
            seed_mask,
            rc_shift,
            soft_mask,
            max_rep,
            fwd_opts,
            rev_opts,
            state
        ),
        #[cfg(target_arch = "x86_64")]
        (
            RefIndexKind::DirectCsr {
                sub_bits, offsets, ..
            },
            IsaTier::Avx2Bmi2Popcnt,
        ) => dispatch_direct_tier!(
            scan_direct_bmi2_popcnt,
            part_filter,
            fwd_only,
            ps_ptr,
            unique_bm,
            query,
            reference,
            rev_ref,
            offsets,
            *sub_bits,
            p_start,
            p_end,
            chunk_start,
            chunk_end,
            wl,
            ny,
            mask,
            rc_shift,
            soft_mask,
            max_rep,
            fwd_opts,
            rev_opts,
            state
        ),
        #[cfg(target_arch = "x86_64")]
        (
            RefIndexKind::DirectCsr {
                sub_bits, offsets, ..
            },
            IsaTier::Avx2Popcnt,
        ) => dispatch_direct_tier!(
            scan_direct_popcnt,
            part_filter,
            fwd_only,
            ps_ptr,
            unique_bm,
            query,
            reference,
            rev_ref,
            offsets,
            *sub_bits,
            p_start,
            p_end,
            chunk_start,
            chunk_end,
            wl,
            ny,
            mask,
            rc_shift,
            soft_mask,
            max_rep,
            fwd_opts,
            rev_opts,
            state
        ),
        (
            RefIndexKind::DirectCsr {
                sub_bits, offsets, ..
            },
            IsaTier::Portable,
        ) => dispatch_direct_tier!(
            scan_direct_portable,
            part_filter,
            fwd_only,
            ps_ptr,
            unique_bm,
            query,
            reference,
            rev_ref,
            offsets,
            *sub_bits,
            p_start,
            p_end,
            chunk_start,
            chunk_end,
            wl,
            ny,
            mask,
            rc_shift,
            soft_mask,
            max_rep,
            fwd_opts,
            rev_opts,
            state
        ),
    }
}

/// Sort merged fragments and remove any fragment fully contained within a preceding fragment
/// on the same `(strand, diag, seq_x, seq_y)`.
pub fn dedup_fragments(fragments: &mut Vec<Fragment>) {
    if fragments.len() <= 1 {
        return;
    }
    fragments.sort_unstable_by_key(|f| {
        (
            f.strand,
            f.diag,
            f.seq_x,
            f.seq_y,
            f.x_start,
            Reverse(f.x_end),
            Reverse(f.score),
        )
    });

    let mut write = 1usize;
    for read in 1..fragments.len() {
        let prev = &fragments[write - 1];
        let cur = &fragments[read];
        let contained = cur.strand == prev.strand
            && cur.diag == prev.diag
            && cur.seq_x == prev.seq_x
            && cur.seq_y == prev.seq_y
            && cur.x_end <= prev.x_end;
        if !contained {
            if write != read {
                fragments[write] = fragments[read].clone();
            }
            write += 1;
        }
    }
    fragments.truncate(write);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::extend::{Scoring, extend_single_hit_outcome};
    use crate::spaced::SpacedSeed;

    #[test]
    fn dedup_removes_identical_and_contained_fragments_only() {
        let mk = |strand, diag, seq_x, seq_y, x_start, x_end, score| Fragment {
            diag,
            x_start,
            y_start: x_start,
            x_end,
            y_end: x_end,
            length: x_end - x_start + 1,
            ident: x_end - x_start + 1,
            score,
            similarity: 100.0,
            seq_x,
            seq_y,
            block: 0,
            strand,
        };

        let mut frags = vec![
            mk('f', 10, 0, 0, 100, 200, 400),
            mk('f', 10, 0, 0, 100, 200, 400), // exact duplicate
            mk('f', 10, 0, 0, 120, 180, 240), // strictly contained
            mk('f', 10, 0, 0, 150, 250, 400), // overlapping but extends further -> keep
            mk('f', 10, 1, 0, 120, 180, 240), // different seq_x -> keep
            mk('r', 10, 0, 0, 100, 200, 400), // different strand -> keep
        ];

        dedup_fragments(&mut frags);
        assert_eq!(frags.len(), 4);
    }

    #[test]
    fn fast_short_reject_64_fuzz_10k_bit_exact_against_scalar_oracle() {
        let mut rng = 0x5EED_CAFE_1234_5678u64;
        let mut next_u64 = || {
            rng = rng
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            rng
        };

        let alphabet = *b"ACGT";
        let wls = [4usize, 8, 12, 16, 20, 24, 28];
        let mut fast_hits = 0usize;
        let mut cov_fast = CoverageMatrix::new();
        let mut cov_oracle = CoverageMatrix::new();

        for trial in 0..10_000usize {
            let wl = wls[trial % wls.len()];
            let spaced = (trial % 3) == 0 && wl >= 8;
            let len = 260 + (next_u64() as usize % 80);
            let mut sx = vec![0u8; len];
            let mut sy = vec![0u8; len];
            let mismatch_pct = [15, 35, 55, 70, 75, 85][trial % 6];

            for i in 0..len {
                let r = next_u64();
                let bx = alphabet[((r >> 16) & 3) as usize];
                sx[i] = bx;
                let r2 = next_u64();
                sy[i] = if (r2 % 100) < mismatch_pct {
                    alphabet[(((r >> 16) + 1 + (r2 % 3)) & 3) as usize]
                } else {
                    bx
                };
            }

            let pos_x = 64 + (next_u64() as usize % 32);
            let pos_y = 64 + (next_u64() as usize % 32);
            for k in 0..wl {
                if !spaced || (k % 3) != 1 {
                    sy[pos_y + k] = sx[pos_x + k];
                }
            }

            // Inject special bases at exact +/- 64 boundaries and inside/outside window on subset.
            if trial % 11 == 0 {
                let target_pos = match trial % 4 {
                    0 => pos_x - 64,
                    1 => pos_x + wl + 63,
                    2 => pos_x - 65,
                    _ => pos_x + wl + 64,
                };
                if target_pos < len {
                    sx[target_pos] = b'N';
                }
            }

            let gx = Genome::from_records("x", [("x0", sx.as_slice())]);
            let gy = Genome::from_records("y", [("y0", sy.as_slice())]);
            let min_len = (wl as u64) + (next_u64() % 45);
            let opts = ExtendOptions {
                min_length: min_len,
                min_similarity: 60.0,
                word_length: wl,
                fixed_length: true,
                strand: 'f',
                scoring: Scoring::default(),
                spaced_seed: spaced,
            };

            let seed_diff = gx.seq.load_32(pos_x) ^ gy.seq.load_32(pos_y);
            let (_si, ss) = seed_init_stats_from_diff(seed_diff, wl as u64, opts.scoring, spaced);
            let diag = pos_x as i64 - pos_y as i64;
            let hit = Hit {
                diag,
                pos_x: pos_x as u64,
            };

            let oracle =
                extend_single_hit_outcome(hit, &gx, &gy, opts, false, false, &mut cov_oracle);

            let can_fast = (0..128).contains(&ss)
                && !gx
                    .seq
                    .has_special_in_short_range(pos_x - 64, pos_x + wl + 64)
                && !gy
                    .seq
                    .has_special_in_short_range(pos_y - 64, pos_y + wl + 64);

            let fast = if can_fast {
                unsafe {
                    try_fast_short_reject_64(
                        &gx.seq,
                        &gy.seq,
                        pos_x as u64,
                        pos_y as u64,
                        wl as u64,
                        ss,
                        min_len,
                        |mm| pack_mm32_swar(mm),
                    )
                }
            } else {
                None
            };

            if let Some((flen, u_max)) = fast {
                fast_hits += 1;
                assert!(
                    oracle.fragment.is_none(),
                    "fast path rejected an accepted fragment on trial {trial}"
                );
                assert!(
                    oracle.xfil_max <= pos_x as u64 + 256,
                    "fast path rejected a >256 bp span on trial {trial}"
                );
                cov_fast.record_unit_exact(flen, u_max);
            } else {
                // Slow path falls through to exact extension.
                let _ = extend_single_hit_outcome(hit, &gx, &gy, opts, false, false, &mut cov_fast);
            }
        }

        assert!(
            fast_hits > 3_000,
            "expected substantial fast-reject rate, got {fast_hits}"
        );
        assert_eq!(
            cov_fast.bins, cov_oracle.bins,
            "CoverageMatrix diverged between fast and oracle"
        );
    }

    #[test]
    fn diag_slot_32_and_64_are_bit_exact_contiguous_and_spaced() {
        let alphabet = *b"ACGT";
        let mut rng = 0x1122_3344_5566_7788u64;
        let mut make_dna = |n: usize| -> Vec<u8> {
            (0..n)
                .map(|_| {
                    rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
                    alphabet[((rng >> 32) & 3) as usize]
                })
                .collect()
        };
        let mut sx = make_dna(8_000);
        let mut sy = make_dna(8_000);
        // Embed homologous forward and reverse-complement blocks.
        sy[1_200..1_500].copy_from_slice(&sx[500..800]);
        let rc_blk: Vec<u8> = sx[2_000..2_300]
            .iter()
            .rev()
            .map(|&b| match b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                _ => b'A',
            })
            .collect();
        sy[4_000..4_300].copy_from_slice(&rc_blk);
        sx[600] = b'C'; // single SNP inside forward block

        let gx = Genome::from_records("x", [("x0", sx.as_slice())]);
        let gy = Genome::from_records("y", [("y0", sy.as_slice())]);
        let rev_y = gy.reverse_complement();

        for spaced in [None, Some(SpacedSeed::parse("w12s18").unwrap())] {
            let wl = spaced.as_ref().map_or(12, |s| s.span);
            let cfg = CompareConfig {
                word_length: wl,
                spaced_seed: spaced,
                min_length: 40,
                min_similarity: 60.0,
                freq_limit: 50,
                ..Default::default()
            };
            let idx = RefIndex::build_with_mask(
                &gy,
                wl,
                cfg.spaced_mask_2bit(),
                cfg.freq_limit,
                cfg.soft_mask,
            )
            .unwrap();
            let PosStrandVec::U32(ref v) = *idx.pos_strand_vec() else {
                panic!("expected U32 pos_strand");
            };
            let r32 = compare_streaming_typed::<&[u32], DiagSlot32>(
                &gx,
                &gy,
                &idx,
                v.as_slice(),
                Some(&rev_y),
                &cfg,
                true,
                true,
            )
            .unwrap();
            let r64 = compare_streaming_typed::<&[u32], DiagSlot64>(
                &gx,
                &gy,
                &idx,
                v.as_slice(),
                Some(&rev_y),
                &cfg,
                true,
                true,
            )
            .unwrap();

            assert_eq!(r32.fragments, r64.fragments);
            assert_eq!(r32.coverage.bins, r64.coverage.bins);
            assert_eq!(r32.stats.word_matches, r64.stats.word_matches);
            assert_eq!(r32.stats.raw_hits, r64.stats.raw_hits);
            assert_eq!(r32.stats.hits_used, r64.stats.hits_used);
            assert!(r32.fragments.len() >= 2);
        }
    }
}
