use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use flate2::read::MultiGzDecoder;

use crate::error::{GeckoError, Result};

/// Separator byte placed between concatenated multi-FASTA records in coordinate space.
/// Stored as `0x00` so uppercasing (`b & 0xDF`) preserves separators identically.
pub const SEQ_SEP: u8 = 0u8;

/// Block shift for `PackedSeq::special_block_bm` (1 bit per 512 bp = 256 B/Mb).
pub(crate) const SPECIAL_BLOCK_SHIFT: usize = 9;
/// Block size in bases for `PackedSeq::special_block_bm` (512 bp).
pub(crate) const SPECIAL_BLOCK_SIZE: usize = 1 << SPECIAL_BLOCK_SHIFT;

const BASE_DECODE_UPPER: [u8; 4] = *b"ACGT";
const BASE_DECODE_LOWER: [u8; 4] = *b"acgt";

/// Reverse the order of the 32 2-bit lanes inside a `u64` (lane 0 ↔ lane 31, ..., lane 15 ↔ lane 16).
#[inline(always)]
pub(crate) fn reverse_32_2bit(mut x: u64) -> u64 {
    x = x.swap_bytes();
    x = ((x & 0x0F0F_0F0F_0F0F_0F0F) << 4) | ((x >> 4) & 0x0F0F_0F0F_0F0F_0F0F);
    x = ((x & 0x3333_3333_3333_3333) << 2) | ((x >> 2) & 0x3333_3333_3333_3333);
    x
}

/// 2-bit/bp packed DNA sequence storage (`0.25 bytes/bp`) with compact sparse annotations for
/// `'N'` runs, soft-masked lowercase runs, multi-FASTA separators, and rare IUPAC ambiguity codes.
///
/// Coordinate positions match GECKO's concatenated multi-FASTA layout (including 1 `SEQ_SEP`
/// position between adjacent records). Base `i` is stored in `bits[i >> 5]` at bit offset
/// `2 * (i & 31)`: `00 = A`, `01 = C`, `10 = G`, `11 = T/U`. Two trailing `0u64` sentinel words
/// ensure branchless 32-base window loads never read out of bounds.
#[derive(Debug, Clone)]
pub struct PackedSeq {
    /// Total length in coordinate positions (bases + `SEQ_SEP` slots).
    pub len: usize,
    /// Packed 2-bit bases (32 bases/word, LSB-first), padded with 2 trailing `0u64` sentinels.
    pub bits: Vec<u64>,
    /// Sorted, non-overlapping half-open `[start, end)` runs of `'N'`/`'n'` bases.
    pub n_runs: Vec<(u64, u64)>,
    /// Sorted, non-overlapping half-open `[start, end)` runs of soft-masked lowercase `a/c/g/t/u`.
    pub soft_runs: Vec<(u64, u64)>,
    /// Sorted `SEQ_SEP` positions (between multi-FASTA records or literal `'*'`).
    pub sep_pos: Vec<u64>,
    /// Rare non-ACGTN IUPAC ambiguity codes `(pos, upper_ascii_byte)` (empty on standard RefSeq).
    pub iupac: Vec<(u64, u8)>,
    /// Sorted, coalesced `[start, end)` union of `n_runs`, `sep_pos`, and `iupac` — allows the
    /// extension loop to compute exact `N`/separator-free spans in a single lookup.
    pub special_runs: Vec<(u64, u64)>,
    /// 1-bit-per-512-bp summary bitmap of `special_runs` (padded with 1 trailing `0u64` word,
    /// or empty when `special_runs` is empty) for O(1) branchless special-free window checks.
    pub(crate) special_block_bm: Vec<u64>,
}

impl PackedSeq {
    /// Total number of coordinate positions (bases + separators).
    #[inline(always)]
    pub fn len(&self) -> usize {
        self.len
    }

    #[inline(always)]
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Approximate heap bytes currently held by this packed sequence.
    pub fn heap_bytes(&self) -> usize {
        self.bits.capacity() * std::mem::size_of::<u64>()
            + self.n_runs.capacity() * std::mem::size_of::<(u64, u64)>()
            + self.soft_runs.capacity() * std::mem::size_of::<(u64, u64)>()
            + self.sep_pos.capacity() * std::mem::size_of::<u64>()
            + self.iupac.capacity() * std::mem::size_of::<(u64, u8)>()
            + self.special_runs.capacity() * std::mem::size_of::<(u64, u64)>()
            + self.special_block_bm.capacity() * std::mem::size_of::<u64>()
    }

    /// Branchless O(1) check whether any 512-bp block overlapping `[lo, hi)` contains a special
    /// base (`'N'`, `SEQ_SEP`, or IUPAC). Zero false negatives.
    ///
    /// Caller must ensure `lo < hi && hi <= self.len && hi - lo <= SPECIAL_BLOCK_SIZE`.
    #[inline(always)]
    pub(crate) fn has_special_in_short_range(&self, lo: usize, hi: usize) -> bool {
        if self.special_block_bm.is_empty() {
            return false;
        }
        debug_assert!(lo < hi && hi <= self.len && hi - lo <= SPECIAL_BLOCK_SIZE);
        let b0 = lo >> SPECIAL_BLOCK_SHIFT;
        let b1 = (hi - 1) >> SPECIAL_BLOCK_SHIFT;
        // Safety: `b0 <= b1 <= (self.len - 1) >> 9`, so `b0 >> 6` and `b1 >> 6` are < `bm.len()`.
        let w0 = unsafe { *self.special_block_bm.get_unchecked(b0 >> 6) };
        let w1 = unsafe { *self.special_block_bm.get_unchecked(b1 >> 6) };
        (((w0 >> (b0 & 63)) | (w1 >> (b1 & 63))) & 1) != 0
    }

    /// Check whether any 512-bp block overlapping `[lo, hi)` contains a special base (`'N'`,
    /// `SEQ_SEP`, or IUPAC). Zero false negatives.
    #[inline]
    pub fn has_special_in_range(&self, lo: usize, hi: usize) -> bool {
        if self.special_block_bm.is_empty() {
            return false;
        }
        let lo = lo.min(self.len);
        let hi = hi.min(self.len);
        if lo >= hi {
            return false;
        }
        let b0 = lo >> SPECIAL_BLOCK_SHIFT;
        let b1 = (hi - 1) >> SPECIAL_BLOCK_SHIFT;
        let w0 = b0 >> 6;
        let w1 = b1 >> 6;
        let bit0 = b0 & 63;
        let bit1 = b1 & 63;
        if w0 == w1 {
            let mask = (u64::MAX << bit0) & (u64::MAX >> (63 - bit1));
            (self.special_block_bm[w0] & mask) != 0
        } else {
            if (self.special_block_bm[w0] & (u64::MAX << bit0)) != 0 {
                return true;
            }
            if self.special_block_bm[w0 + 1..w1].iter().any(|&w| w != 0) {
                return true;
            }
            (self.special_block_bm[w1] & (u64::MAX >> (63 - bit1))) != 0
        }
    }

    /// Return the 2-bit base value (`0..=3`) at coordinate `pos`.
    #[inline(always)]
    pub fn base_2bit(&self, pos: usize) -> u8 {
        let w = unsafe { *self.bits.get_unchecked(pos >> 5) };
        ((w >> ((pos & 31) << 1)) & 3) as u8
    }

    /// Branchless load of 32 consecutive 2-bit bases `[pos .. pos + 32)` packed LSB-first
    /// (lane 0 = `pos`, lane 31 = `pos + 31`).
    ///
    /// Caller must ensure `pos <= self.len`.
    #[inline(always)]
    pub fn load_32(&self, pos: usize) -> u64 {
        let w_idx = pos >> 5;
        let shift = (pos & 31) << 1;
        // Safety: `bits` is padded with 2 trailing sentinel u64 words so `w_idx + 1 < bits.len()`.
        let lo = unsafe { *self.bits.get_unchecked(w_idx) };
        let hi = unsafe { *self.bits.get_unchecked(w_idx + 1) };
        let pair = (lo as u128) | ((hi as u128) << 64);
        (pair >> shift) as u64
    }

    /// Branchless load of 32 consecutive 2-bit bases `[pos - 31 ..= pos]` in **reverse** order
    /// (lane 0 = `pos`, lane 1 = `pos - 1`, ..., lane 31 = `pos - 31`).
    ///
    /// Caller must ensure `pos >= 31 && pos < self.len`.
    #[inline(always)]
    pub fn load_32_rev(&self, pos: usize) -> u64 {
        reverse_32_2bit(self.load_32(pos - 31))
    }

    /// Return the uppercase ASCII byte at `pos` (`b'A'`, `b'C'`, `b'G'`, `b'T'`, `b'N'`,
    /// `SEQ_SEP`, or rare uppercase IUPAC code).
    #[inline]
    pub fn byte_upper_at(&self, pos: usize) -> u8 {
        if !self.special_runs.is_empty() {
            let p = pos as u64;
            if contains_run(&self.n_runs, p) {
                return b'N';
            }
            if !self.sep_pos.is_empty() && self.sep_pos.binary_search(&p).is_ok() {
                return SEQ_SEP;
            }
            if !self.iupac.is_empty()
                && let Ok(idx) = self.iupac.binary_search_by_key(&p, |&(q, _)| q)
            {
                return self.iupac[idx].1;
            }
        }
        BASE_DECODE_UPPER[self.base_2bit(pos) as usize]
    }

    /// Reconstruct the raw ASCII byte (preserving lowercase soft-masking) at `pos`.
    pub fn byte_raw_at(&self, pos: usize) -> u8 {
        let up = self.byte_upper_at(pos);
        if up != SEQ_SEP
            && up != b'N'
            && !self.soft_runs.is_empty()
            && contains_run(&self.soft_runs, pos as u64)
        {
            return BASE_DECODE_LOWER[self.base_2bit(pos) as usize];
        }
        up
    }

    /// Bulk-decode the entire sequence into raw ASCII bytes (preserving case, `'N'`, `SEQ_SEP`,
    /// and rare IUPAC codes) in linear time (~20 ms for 250 Mb; used for GPU buffer uploads).
    pub fn to_ascii_vec(&self) -> Vec<u8> {
        self.to_ascii_vec_padded(self.len)
    }

    /// Bulk-decode into raw ASCII bytes zero-padded to a 4-byte boundary plus 2 trailing zero
    /// `u32` sentinel words (for direct `u32` GPU storage buffer upload).
    pub fn to_ascii_vec_padded_u32(&self) -> Vec<u8> {
        let padded = (self.len.div_ceil(4) + 2).max(4) * 4;
        self.to_ascii_vec_padded(padded)
    }

    fn to_ascii_vec_padded(&self, alloc_len: usize) -> Vec<u8> {
        let mut out = vec![0u8; alloc_len.max(self.len)];
        let full_words = self.len / 32;
        for (w_idx, &word) in self.bits.iter().enumerate().take(full_words) {
            let mut w = word;
            let base_off = w_idx * 32;
            let dst = &mut out[base_off..base_off + 32];
            for slot in dst.iter_mut() {
                *slot = BASE_DECODE_UPPER[(w & 3) as usize];
                w >>= 2;
            }
        }
        let rem_start = full_words * 32;
        if rem_start < self.len {
            let mut w = self.bits[full_words];
            for slot in &mut out[rem_start..self.len] {
                *slot = BASE_DECODE_UPPER[(w & 3) as usize];
                w >>= 2;
            }
        }
        for &(s, e) in &self.soft_runs {
            for b in &mut out[s as usize..e as usize] {
                *b |= 0x20;
            }
        }
        for &(s, e) in &self.n_runs {
            out[s as usize..e as usize].fill(b'N');
        }
        for &p in &self.sep_pos {
            out[p as usize] = SEQ_SEP;
        }
        for &(p, b) in &self.iupac {
            out[p as usize] = b;
        }
        out
    }

    /// Invoke `f(run_start, run_end)` for every maximal contiguous run of valid A/C/G/T bases
    /// within `[start, end)`, jumping over `N`s, `SEQ_SEP`, IUPAC codes, and (when `soft_mask`
    /// is true) soft-masked lowercase runs in O(1) per gap.
    #[inline]
    pub fn for_each_valid_run<F>(&self, start: usize, end: usize, soft_mask: bool, mut f: F)
    where
        F: FnMut(usize, usize),
    {
        let start = start.min(self.len);
        let end = end.min(self.len);
        if start >= end {
            return;
        }

        let sp = &self.special_runs;
        let sm: &[(u64, u64)] = if soft_mask { &self.soft_runs } else { &[] };

        let mut si = sp.partition_point(|&(_, e)| (e as usize) <= start);
        let mut mi = sm.partition_point(|&(_, e)| (e as usize) <= start);
        let mut cur = start;

        while cur < end {
            while si < sp.len() && (sp[si].1 as usize) <= cur {
                si += 1;
            }
            while mi < sm.len() && (sm[mi].1 as usize) <= cur {
                mi += 1;
            }
            let in_sp = si < sp.len() && (sp[si].0 as usize) <= cur;
            let in_sm = mi < sm.len() && (sm[mi].0 as usize) <= cur;
            if in_sp || in_sm {
                let jump_sp = if in_sp { sp[si].1 as usize } else { cur };
                let jump_sm = if in_sm { sm[mi].1 as usize } else { cur };
                cur = jump_sp.max(jump_sm);
                continue;
            }
            let next_sp = if si < sp.len() {
                sp[si].0 as usize
            } else {
                end
            };
            let next_sm = if mi < sm.len() {
                sm[mi].0 as usize
            } else {
                end
            };
            let run_end = end.min(next_sp).min(next_sm);
            if run_end > cur {
                f(cur, run_end);
            }
            cur = run_end;
        }
    }

    /// Fast 2-bit reverse complement of the entire sequence (exact universal mirror
    /// `rev_pos = len - 1 - pos` across multi-FASTA records and separators).
    pub fn reverse_complement(&self) -> Self {
        let len = self.len;
        let n_words = len.div_ceil(32);
        let mut bits = Vec::with_capacity(n_words + 2);
        let full_words = len / 32;

        {
            let spare = &mut bits.spare_capacity_mut()[..n_words + 2];
            if full_words >= 65_536 && rayon::current_num_threads() > 1 {
                use rayon::prelude::*;
                const CHUNK_WORDS: usize = 16_384;
                spare[..full_words]
                    .par_chunks_mut(CHUNK_WORDS)
                    .enumerate()
                    .for_each(|(chunk_idx, dst)| {
                        let base_w = chunk_idx * CHUNK_WORDS;
                        for (off, slot) in dst.iter_mut().enumerate() {
                            let w = base_w + off;
                            let src_start = len - (w + 1) * 32;
                            // Complementing 2-bit bases (A=00<->T=11, C=01<->G=10) is bitwise NOT!
                            slot.write(!reverse_32_2bit(self.load_32(src_start)));
                        }
                    });
            } else {
                for (w, slot) in spare[..full_words].iter_mut().enumerate() {
                    let src_start = len - (w + 1) * 32;
                    slot.write(!reverse_32_2bit(self.load_32(src_start)));
                }
            }
            let rem = len % 32;
            if rem > 0 {
                let mut last = 0u64;
                for k in 0..rem {
                    let src_pos = rem - 1 - k;
                    let comp = (self.base_2bit(src_pos) ^ 3) as u64;
                    last |= comp << (2 * k);
                }
                spare[full_words].write(last);
            }
            spare[n_words].write(0);
            spare[n_words + 1].write(0);
        }
        // Safety: all `n_words + 2` slots in `spare` were initialized above.
        unsafe { bits.set_len(n_words + 2) };

        let len_u64 = len as u64;
        let mirror_runs = |runs: &[(u64, u64)]| -> Vec<(u64, u64)> {
            runs.iter()
                .rev()
                .map(|&(s, e)| (len_u64 - e, len_u64 - s))
                .collect()
        };

        let n_runs = mirror_runs(&self.n_runs);
        let soft_runs = mirror_runs(&self.soft_runs);
        let special_runs = mirror_runs(&self.special_runs);
        let special_block_bm = build_special_block_bm(len, &special_runs);
        let sep_pos: Vec<u64> = self
            .sep_pos
            .iter()
            .rev()
            .map(|&p| len_u64 - 1 - p)
            .collect();
        let iupac: Vec<(u64, u8)> = self
            .iupac
            .iter()
            .rev()
            .map(|&(p, b)| (len_u64 - 1 - p, complement_base(b)))
            .collect();

        Self {
            len,
            bits,
            n_runs,
            soft_runs,
            sep_pos,
            iupac,
            special_runs,
            special_block_bm,
        }
    }

    /// Build a 1-bit/bp bitmap (`u32` words, padded with 2 trailing `0u32` sentinels) marking
    /// every special position (`'N'`, `SEQ_SEP`, or rare IUPAC code). Returns `None` if empty.
    pub fn special_bitmap(&self) -> Option<Vec<u32>> {
        if self.special_runs.is_empty() || self.len == 0 {
            return None;
        }
        Some(build_runs_bitmap(self.len, &[&self.special_runs]))
    }

    /// Build a 1-bit/bp bitmap (`u32` words, padded with 2 trailing `0u32` sentinels) marking
    /// every position that cannot participate in a seed k-mer (special bases, plus soft-masked
    /// lowercase runs when `soft_mask` is true).
    ///
    /// Returns `None` if all positions are valid seed bases. Callers can also check
    /// `self.invalid_bitmap_aliases_special(soft_mask)` to reuse `special_bitmap` directly
    /// without allocating or uploading a second bitmap.
    pub fn invalid_seed_bitmap(&self, soft_mask: bool) -> Option<Vec<u32>> {
        if self.len == 0 {
            return None;
        }
        let has_soft = soft_mask && !self.soft_runs.is_empty();
        if self.special_runs.is_empty() && !has_soft {
            return None;
        }
        if has_soft {
            Some(build_runs_bitmap(
                self.len,
                &[&self.special_runs, &self.soft_runs],
            ))
        } else {
            Some(build_runs_bitmap(self.len, &[&self.special_runs]))
        }
    }

    /// Returns `true` when the invalid-seed bitmap is identical to `special_bitmap`
    /// (i.e. `!soft_mask || self.soft_runs.is_empty()`), allowing GPU backends to alias both.
    #[inline]
    pub fn invalid_bitmap_aliases_special(&self, soft_mask: bool) -> bool {
        !soft_mask || self.soft_runs.is_empty()
    }
}

/// Fill a 1-bit/bp `u32` bitmap (with 2 trailing zero sentinel words for unaligned 32-bit loads)
/// from one or more sorted half-open `[start, end)` run lists.
fn build_runs_bitmap(len: usize, run_slices: &[&[(u64, u64)]]) -> Vec<u32> {
    let n_words = len.div_ceil(32);
    let mut words = vec![0u32; n_words + 2];
    for runs in run_slices {
        for &(s, e) in *runs {
            let start = (s as usize).min(len);
            let end = (e as usize).min(len);
            if start >= end {
                continue;
            }
            let w_s = start >> 5;
            let b_s = start & 31;
            let last = end - 1;
            let w_e = last >> 5;
            let b_e = last & 31;
            if w_s == w_e {
                let mask = (u32::MAX << b_s) & (u32::MAX >> (31 - b_e));
                words[w_s] |= mask;
            } else {
                words[w_s] |= u32::MAX << b_s;
                if w_s + 1 < w_e {
                    words[w_s + 1..w_e].fill(u32::MAX);
                }
                words[w_e] |= u32::MAX >> (31 - b_e);
            }
        }
    }
    words
}

#[inline(always)]
fn contains_run(runs: &[(u64, u64)], pos: u64) -> bool {
    let idx = runs.partition_point(|&(s, _)| s <= pos);
    idx > 0 && pos < runs[idx - 1].1
}

/// Compact LUT for fast scalar FASTA byte packing:
/// `0..=3` = uppercase A/C/G/T(U), `4..=7` = lowercase a/c/g/t(u), `8` = N/n, `0xFF` = other.
const PACK_LUT: [u8; 256] = {
    let mut lut = [0xFFu8; 256];
    lut[b'A' as usize] = 0;
    lut[b'C' as usize] = 1;
    lut[b'G' as usize] = 2;
    lut[b'T' as usize] = 3;
    lut[b'U' as usize] = 3;
    lut[b'a' as usize] = 4;
    lut[b'c' as usize] = 5;
    lut[b'g' as usize] = 6;
    lut[b't' as usize] = 7;
    lut[b'u' as usize] = 7;
    lut[b'N' as usize] = 8;
    lut[b'n' as usize] = 8;
    lut
};

/// Spread 32 bits of `x` into even bits (`0, 2, 4, ..., 62`) of a `u64` (portable Morton).
#[inline(always)]
fn spread_32_to_even_64(x: u32) -> u64 {
    let mut v = x as u64;
    v = (v | (v << 16)) & 0x0000_FFFF_0000_FFFF;
    v = (v | (v << 8)) & 0x00FF_00FF_00FF_00FF;
    v = (v | (v << 4)) & 0x0F0F_0F0F_0F0F_0F0F;
    v = (v | (v << 2)) & 0x3333_3333_3333_3333;
    v = (v | (v << 1)) & 0x5555_5555_5555_5555;
    v
}

/// Spread 16 bits of `x` into even bits (`0, 2, 4, ..., 30`) of a `u32` (portable Morton).
#[inline(always)]
fn spread_16_to_even_32(x: u16) -> u32 {
    let mut v = x as u32;
    v = (v | (v << 8)) & 0x00FF_00FF;
    v = (v | (v << 4)) & 0x0F0F_0F0F;
    v = (v | (v << 2)) & 0x3333_3333;
    v = (v | (v << 1)) & 0x5555_5555;
    v
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "bmi2")]
#[inline]
unsafe fn interleave_2bit_32_bmi2(lo_mask: u32, hi_mask: u32) -> u64 {
    core::arch::x86_64::_pdep_u64(lo_mask as u64, 0x5555_5555_5555_5555)
        | core::arch::x86_64::_pdep_u64(hi_mask as u64, 0xAAAA_AAAA_AAAA_AAAA)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "bmi2")]
#[inline]
unsafe fn interleave_2bit_16_bmi2(lo_mask: u16, hi_mask: u16) -> u32 {
    core::arch::x86_64::_pdep_u32(lo_mask as u32, 0x5555_5555)
        | core::arch::x86_64::_pdep_u32(hi_mask as u32, 0xAAAA_AAAA)
}

#[inline(always)]
fn interleave_2bit_32(lo_mask: u32, hi_mask: u32, use_bmi2: bool) -> u64 {
    #[cfg(target_arch = "x86_64")]
    if use_bmi2 {
        // Safety: caller verified BMI2 feature detection.
        return unsafe { interleave_2bit_32_bmi2(lo_mask, hi_mask) };
    }
    let _ = use_bmi2;
    spread_32_to_even_64(lo_mask) | (spread_32_to_even_64(hi_mask) << 1)
}

#[inline(always)]
fn interleave_2bit_16(lo_mask: u16, hi_mask: u16, use_bmi2: bool) -> u32 {
    #[cfg(target_arch = "x86_64")]
    if use_bmi2 {
        // Safety: caller verified BMI2 feature detection.
        return unsafe { interleave_2bit_16_bmi2(lo_mask, hi_mask) };
    }
    let _ = use_bmi2;
    spread_16_to_even_32(lo_mask) | (spread_16_to_even_32(hi_mask) << 1)
}

/// Streaming builder that packs ASCII FASTA bytes directly into 2-bit [`PackedSeq`] words
/// without allocating an intermediate 1-byte/bp buffer.
pub(crate) struct PackedSeqBuilder {
    len: usize,
    bits: Vec<u64>,
    cur_word: u64,
    bit_pos: u32,
    n_runs: Vec<(u64, u64)>,
    soft_runs: Vec<(u64, u64)>,
    sep_pos: Vec<u64>,
    iupac: Vec<(u64, u8)>,
    use_avx2: bool,
    use_bmi2: bool,
}

impl PackedSeqBuilder {
    pub fn with_capacity_bases(cap_bases: usize) -> Self {
        #[cfg(target_arch = "x86_64")]
        let (use_avx2, use_bmi2) = (
            std::is_x86_feature_detected!("avx2"),
            std::is_x86_feature_detected!("bmi2"),
        );
        #[cfg(not(target_arch = "x86_64"))]
        let (use_avx2, use_bmi2) = (false, false);

        Self {
            len: 0,
            bits: Vec::with_capacity(cap_bases.div_ceil(32) + 2),
            cur_word: 0,
            bit_pos: 0,
            n_runs: Vec::new(),
            soft_runs: Vec::new(),
            sep_pos: Vec::new(),
            iupac: Vec::new(),
            use_avx2,
            use_bmi2,
        }
    }

    #[inline(always)]
    fn push_2bit(&mut self, v: u64) {
        self.cur_word |= v << self.bit_pos;
        self.bit_pos += 2;
        if self.bit_pos == 64 {
            self.bits.push(self.cur_word);
            self.cur_word = 0;
            self.bit_pos = 0;
        }
        self.len += 1;
    }

    #[inline(always)]
    fn push_32_2bit(&mut self, word: u64) {
        if self.bit_pos == 0 {
            self.bits.push(word);
        } else {
            self.cur_word |= word << self.bit_pos;
            self.bits.push(self.cur_word);
            self.cur_word = word >> (64 - self.bit_pos);
        }
        self.len += 32;
    }

    #[inline(always)]
    fn push_16_2bit(&mut self, word: u32) {
        let w = word as u64;
        self.cur_word |= w << self.bit_pos;
        let next_pos = self.bit_pos + 32;
        if next_pos >= 64 {
            self.bits.push(self.cur_word);
            let rem = next_pos - 64;
            self.cur_word = if rem == 0 { 0 } else { w >> (32 - rem) };
            self.bit_pos = rem;
        } else {
            self.bit_pos = next_pos;
        }
        self.len += 16;
    }

    #[inline(always)]
    fn record_n_range(&mut self, start: u64, end: u64) {
        if let Some(last) = self.n_runs.last_mut()
            && last.1 == start
        {
            last.1 = end;
            return;
        }
        self.n_runs.push((start, end));
    }

    #[inline(always)]
    fn record_soft_range(&mut self, start: u64, end: u64) {
        if let Some(last) = self.soft_runs.last_mut()
            && last.1 == start
        {
            last.1 = end;
            return;
        }
        self.soft_runs.push((start, end));
    }

    #[inline(always)]
    fn record_soft_bitmask(&mut self, base_pos: u64, mut mask: u32) {
        while mask != 0 {
            let s = mask.trailing_zeros();
            let shifted = mask >> s;
            let run_len = shifted.trailing_ones();
            self.record_soft_range(base_pos + s as u64, base_pos + (s + run_len) as u64);
            if s + run_len >= 32 {
                break;
            }
            mask &= !(((1u32 << run_len) - 1) << s);
        }
    }

    #[inline]
    pub fn push_sep(&mut self) {
        let pos = self.len as u64;
        self.sep_pos.push(pos);
        self.push_2bit(0);
    }

    /// Push a raw FASTA sequence line/slice, returning the number of nucleotide bases added
    /// (excluding separators). Uses AVX2 (32 bp/step) + SSE2 (16 bp/step) on x86_64.
    #[inline]
    pub fn push_fasta_bytes(&mut self, bytes: &[u8]) -> u64 {
        #[cfg(target_arch = "x86_64")]
        if self.use_avx2 && bytes.len() >= 16 {
            // Safety: `self.use_avx2` verified by runtime `is_x86_feature_detected!("avx2")`.
            return unsafe { self.push_fasta_bytes_avx2(bytes) };
        }
        self.push_fasta_bytes_scalar(bytes)
    }

    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "avx2")]
    unsafe fn push_fasta_bytes_avx2(&mut self, bytes: &[u8]) -> u64 {
        use core::arch::x86_64::*;

        let mut added = 0u64;
        let mut i = 0usize;
        let n = bytes.len();
        let use_bmi2 = self.use_bmi2;

        let va = _mm256_set1_epi8(b'A' as i8);
        let vc = _mm256_set1_epi8(b'C' as i8);
        let vg = _mm256_set1_epi8(b'G' as i8);
        let vt = _mm256_set1_epi8(b'T' as i8);
        let vn = _mm256_set1_epi8(b'N' as i8);
        let v_upper_mask = _mm256_set1_epi8(0xDFu8 as i8);
        let v_lower_bit = _mm256_set1_epi8(0x20);

        while i + 32 <= n {
            let v = unsafe { _mm256_loadu_si256(bytes.as_ptr().add(i) as *const __m256i) };
            let eq_a = _mm256_cmpeq_epi8(v, va);
            let eq_c = _mm256_cmpeq_epi8(v, vc);
            let eq_g = _mm256_cmpeq_epi8(v, vg);
            let eq_t = _mm256_cmpeq_epi8(v, vt);
            let valid_up =
                _mm256_or_si256(_mm256_or_si256(eq_a, eq_c), _mm256_or_si256(eq_g, eq_t));

            if _mm256_movemask_epi8(valid_up) == -1 {
                let lo = _mm256_movemask_epi8(_mm256_or_si256(eq_c, eq_t)) as u32;
                let hi = _mm256_movemask_epi8(_mm256_or_si256(eq_g, eq_t)) as u32;
                self.push_32_2bit(interleave_2bit_32(lo, hi, use_bmi2));
                added += 32;
                i += 32;
                continue;
            }

            // Case-insensitive A/C/G/T check (handles soft-masked repeats & mixed boundaries).
            let vu = _mm256_and_si256(v, v_upper_mask);
            let ua = _mm256_cmpeq_epi8(vu, va);
            let uc = _mm256_cmpeq_epi8(vu, vc);
            let ug = _mm256_cmpeq_epi8(vu, vg);
            let ut = _mm256_cmpeq_epi8(vu, vt);
            let valid_ci = _mm256_or_si256(_mm256_or_si256(ua, uc), _mm256_or_si256(ug, ut));

            if _mm256_movemask_epi8(valid_ci) == -1 {
                let base_pos = self.len as u64;
                let lo = _mm256_movemask_epi8(_mm256_or_si256(uc, ut)) as u32;
                let hi = _mm256_movemask_epi8(_mm256_or_si256(ug, ut)) as u32;
                let lower = _mm256_movemask_epi8(_mm256_cmpeq_epi8(
                    _mm256_and_si256(v, v_lower_bit),
                    v_lower_bit,
                )) as u32;
                if lower == u32::MAX {
                    self.record_soft_range(base_pos, base_pos + 32);
                } else if lower != 0 {
                    self.record_soft_bitmask(base_pos, lower);
                }
                self.push_32_2bit(interleave_2bit_32(lo, hi, use_bmi2));
                added += 32;
                i += 32;
                continue;
            }

            // All-N chunk fast path (centromeres / assembly gaps).
            if _mm256_movemask_epi8(_mm256_cmpeq_epi8(vu, vn)) == -1 {
                let base_pos = self.len as u64;
                self.record_n_range(base_pos, base_pos + 32);
                self.push_32_2bit(0);
                added += 32;
                i += 32;
                continue;
            }

            added += self.push_fasta_bytes_scalar(&bytes[i..i + 32]);
            i += 32;
        }

        // 16-byte SSE2 path covers the tail of 80-bp NCBI lines (32+32+16) with zero scalar work
        // and half the tail of 60-bp Ensembl lines (32+16+12).
        if i + 16 <= n {
            let va16 = _mm_set1_epi8(b'A' as i8);
            let vc16 = _mm_set1_epi8(b'C' as i8);
            let vg16 = _mm_set1_epi8(b'G' as i8);
            let vt16 = _mm_set1_epi8(b'T' as i8);
            let vn16 = _mm_set1_epi8(b'N' as i8);
            let vum16 = _mm_set1_epi8(0xDFu8 as i8);
            let vlb16 = _mm_set1_epi8(0x20);

            let v = unsafe { _mm_loadu_si128(bytes.as_ptr().add(i) as *const __m128i) };
            let eq_a = _mm_cmpeq_epi8(v, va16);
            let eq_c = _mm_cmpeq_epi8(v, vc16);
            let eq_g = _mm_cmpeq_epi8(v, vg16);
            let eq_t = _mm_cmpeq_epi8(v, vt16);
            let valid_up = _mm_or_si128(_mm_or_si128(eq_a, eq_c), _mm_or_si128(eq_g, eq_t));

            if _mm_movemask_epi8(valid_up) == 0xFFFF {
                let lo = _mm_movemask_epi8(_mm_or_si128(eq_c, eq_t)) as u16;
                let hi = _mm_movemask_epi8(_mm_or_si128(eq_g, eq_t)) as u16;
                self.push_16_2bit(interleave_2bit_16(lo, hi, use_bmi2));
                added += 16;
                i += 16;
            } else {
                let vu = _mm_and_si128(v, vum16);
                let ua = _mm_cmpeq_epi8(vu, va16);
                let uc = _mm_cmpeq_epi8(vu, vc16);
                let ug = _mm_cmpeq_epi8(vu, vg16);
                let ut = _mm_cmpeq_epi8(vu, vt16);
                let valid_ci = _mm_or_si128(_mm_or_si128(ua, uc), _mm_or_si128(ug, ut));

                if _mm_movemask_epi8(valid_ci) == 0xFFFF {
                    let base_pos = self.len as u64;
                    let lo = _mm_movemask_epi8(_mm_or_si128(uc, ut)) as u16;
                    let hi = _mm_movemask_epi8(_mm_or_si128(ug, ut)) as u16;
                    let lower =
                        _mm_movemask_epi8(_mm_cmpeq_epi8(_mm_and_si128(v, vlb16), vlb16)) as u32;
                    if lower == 0xFFFF {
                        self.record_soft_range(base_pos, base_pos + 16);
                    } else if lower != 0 {
                        self.record_soft_bitmask(base_pos, lower);
                    }
                    self.push_16_2bit(interleave_2bit_16(lo, hi, use_bmi2));
                    added += 16;
                    i += 16;
                } else if _mm_movemask_epi8(_mm_cmpeq_epi8(vu, vn16)) == 0xFFFF {
                    let base_pos = self.len as u64;
                    self.record_n_range(base_pos, base_pos + 16);
                    self.push_16_2bit(0);
                    added += 16;
                    i += 16;
                } else {
                    added += self.push_fasta_bytes_scalar(&bytes[i..i + 16]);
                    i += 16;
                }
            }
        }

        if i < n {
            added += self.push_fasta_bytes_scalar(&bytes[i..]);
        }
        added
    }

    #[inline]
    fn push_fasta_bytes_scalar(&mut self, bytes: &[u8]) -> u64 {
        let mut added = 0u64;
        for &b in bytes {
            let code = PACK_LUT[b as usize];
            if code < 4 {
                self.push_2bit(code as u64);
                added += 1;
            } else if code < 8 {
                let pos = self.len as u64;
                self.record_soft_range(pos, pos + 1);
                self.push_2bit((code & 3) as u64);
                added += 1;
            } else if code == 8 {
                let pos = self.len as u64;
                self.record_n_range(pos, pos + 1);
                self.push_2bit(0);
                added += 1;
            } else if b == b'*' || b == SEQ_SEP {
                self.push_sep();
            } else if b.is_ascii_alphabetic() {
                let pos = self.len as u64;
                self.iupac.push((pos, b.to_ascii_uppercase()));
                self.push_2bit(0);
                added += 1;
            }
        }
        added
    }

    /// Stitch another [`PackedSeqBuilder`] onto the end of `self` with exact 2-bit word alignment
    /// and shifted/merged sparse metadata runs.
    fn append_builder(&mut self, other: &PackedSeqBuilder) {
        if other.len == 0 {
            return;
        }
        let base_off = self.len as u64;
        if self.bit_pos == 0 {
            self.bits.extend_from_slice(&other.bits);
            self.cur_word = other.cur_word;
            self.bit_pos = other.bit_pos;
            self.len += other.len;
        } else {
            let shift = self.bit_pos;
            let inv = 64 - shift;
            self.bits.reserve(other.bits.len() + 1);
            for &w in &other.bits {
                self.cur_word |= w << shift;
                self.bits.push(self.cur_word);
                self.cur_word = w >> inv;
            }
            if other.bit_pos > 0 {
                self.cur_word |= other.cur_word << shift;
                let next_pos = shift + other.bit_pos;
                if next_pos >= 64 {
                    self.bits.push(self.cur_word);
                    let rem = next_pos - 64;
                    self.cur_word = if rem == 0 { 0 } else { other.cur_word >> inv };
                    self.bit_pos = rem;
                } else {
                    self.bit_pos = next_pos;
                }
            }
            self.len += other.len;
        }

        for &(s, e) in &other.n_runs {
            self.record_n_range(s + base_off, e + base_off);
        }
        for &(s, e) in &other.soft_runs {
            self.record_soft_range(s + base_off, e + base_off);
        }
        for &p in &other.sep_pos {
            self.sep_pos.push(p + base_off);
        }
        for &(p, b) in &other.iupac {
            self.iupac.push((p + base_off, b));
        }
    }

    pub fn finish(mut self) -> PackedSeq {
        if self.bit_pos > 0 {
            self.bits.push(self.cur_word);
        }
        // Two trailing zero sentinel words for branchless `load_32`.
        self.bits.push(0);
        self.bits.push(0);

        let special_runs = build_special_runs(&self.n_runs, &self.sep_pos, &self.iupac);
        let special_block_bm = build_special_block_bm(self.len, &special_runs);

        PackedSeq {
            len: self.len,
            bits: self.bits,
            n_runs: self.n_runs,
            soft_runs: self.soft_runs,
            sep_pos: self.sep_pos,
            iupac: self.iupac,
            special_runs,
            special_block_bm,
        }
    }
}

/// Build a 1-bit-per-512-bp summary bitmap over `special_runs` (with 1 trailing `0u64` word),
/// returning an empty `Vec` when `special_runs` is empty or `len == 0`.
fn build_special_block_bm(len: usize, special_runs: &[(u64, u64)]) -> Vec<u64> {
    if special_runs.is_empty() || len == 0 {
        return Vec::new();
    }
    let n_blocks = len.div_ceil(SPECIAL_BLOCK_SIZE);
    let mut words = vec![0u64; n_blocks.div_ceil(64) + 1];
    for &(s, e) in special_runs {
        let start = (s as usize).min(len);
        let end = (e as usize).min(len);
        if start >= end {
            continue;
        }
        let b_s = start >> SPECIAL_BLOCK_SHIFT;
        let b_e = (end - 1) >> SPECIAL_BLOCK_SHIFT;
        let w_s = b_s >> 6;
        let w_e = b_e >> 6;
        let bit_s = b_s & 63;
        let bit_e = b_e & 63;
        if w_s == w_e {
            words[w_s] |= (u64::MAX << bit_s) & (u64::MAX >> (63 - bit_e));
        } else {
            words[w_s] |= u64::MAX << bit_s;
            words[w_s + 1..w_e].fill(u64::MAX);
            words[w_e] |= u64::MAX >> (63 - bit_e);
        }
    }
    words
}

fn build_special_runs(
    n_runs: &[(u64, u64)],
    sep_pos: &[u64],
    iupac: &[(u64, u8)],
) -> Vec<(u64, u64)> {
    let total = n_runs.len() + sep_pos.len() + iupac.len();
    if total == 0 {
        return Vec::new();
    }
    let mut raw = Vec::with_capacity(total);
    raw.extend_from_slice(n_runs);
    for &p in sep_pos {
        raw.push((p, p + 1));
    }
    for &(p, _) in iupac {
        raw.push((p, p + 1));
    }
    raw.sort_unstable_by_key(|&(s, _)| s);

    let mut merged: Vec<(u64, u64)> = Vec::with_capacity(raw.len());
    for (s, e) in raw {
        if let Some(last) = merged.last_mut()
            && s <= last.1
        {
            if e > last.1 {
                last.1 = e;
            }
        } else {
            merged.push((s, e));
        }
    }
    merged
}

/// Metadata for a single sequence record inside a [`Genome`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SeqRecord {
    /// Short sequence identifier (first whitespace-delimited token after `>`).
    pub id: String,
    /// Full header line without leading `>` or trailing newline.
    pub header: String,
    /// Start offset in `Genome::seq` (includes `SEQ_SEP` slots from preceding records).
    pub start_with_sep: u64,
    /// Start offset in concatenated nucleotide coordinates (excluding `SEQ_SEP`).
    pub start_no_sep: u64,
    /// Sequence length in bases.
    pub len: u64,
}

/// An in-memory FASTA collection (single or multi-FASTA) stored at 2 bits/bp in [`PackedSeq`].
#[derive(Debug, Clone)]
pub struct Genome {
    /// Display name (typically the file stem or first record ID).
    pub name: String,
    /// 2-bit/bp packed concatenated sequence with `SEQ_SEP` slots between records.
    pub seq: PackedSeq,
    /// Per-record metadata in file order.
    pub records: Vec<SeqRecord>,
    /// Total number of bases across all records (excluding separators).
    pub total_bases: u64,
}

impl Genome {
    /// Load a plain or gzip-compressed FASTA file directly into 2-bit/bp storage.
    ///
    /// For uncompressed FASTAs >= 8 MiB on Unix when multiple Rayon threads are available,
    /// partitions the file along newline boundaries and packs 2-bit sequences in parallel via
    /// positional reads (`pread`), then stitches chunks bit-identically.
    pub fn from_fasta_path(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let file = File::open(path)?;
        let file_len = file.metadata().map(|m| m.len() as usize).unwrap_or(0);
        let name = file_stem_name(path);

        // Auto-detect gzip via 2-byte magic header (0x1f 0x8b) without a 4 MiB upfront read.
        #[cfg(unix)]
        let is_gz = {
            use std::os::unix::fs::FileExt;
            let mut magic = [0u8; 2];
            file.read_at(&mut magic, 0).unwrap_or(0) >= 2 && magic == [0x1f, 0x8b]
        };
        #[cfg(not(unix))]
        let is_gz = {
            use std::io::{Read, Seek};
            let mut f = &file;
            let mut magic = [0u8; 2];
            let n = f.read(&mut magic).unwrap_or(0);
            f.seek(std::io::SeekFrom::Start(0))?;
            n >= 2 && magic == [0x1f, 0x8b]
        };

        if is_gz {
            let reader = BufReader::with_capacity(4 << 20, file);
            let gz = MultiGzDecoder::new(reader);
            Self::from_fasta_reader_with_hint(
                BufReader::with_capacity(4 << 20, gz),
                name,
                path,
                file_len.saturating_mul(3),
            )
        } else {
            #[cfg(unix)]
            {
                const PAR_THRESHOLD: usize = 8 * 1024 * 1024;
                if file_len >= PAR_THRESHOLD && rayon::current_num_threads() > 1 {
                    return Self::from_plain_fasta_parallel(&file, file_len, name, path);
                }
            }
            let reader = BufReader::with_capacity(4 << 20, file);
            Self::from_fasta_reader_with_hint(reader, name, path, file_len)
        }
    }

    #[cfg(unix)]
    fn from_plain_fasta_parallel(
        file: &File,
        file_len: usize,
        name: String,
        path: &Path,
    ) -> Result<Self> {
        use rayon::prelude::*;
        use std::os::unix::fs::FileExt;

        let n_chunks = (rayon::current_num_threads() * 2).clamp(2, 32);
        let mut boundaries = Vec::with_capacity(n_chunks + 1);
        boundaries.push(0usize);

        let mut probe = [0u8; 4096];
        for c in 1..n_chunks {
            let approx = file_len * c / n_chunks;
            let mut pos = approx;
            let mut found = None;
            while pos < file_len {
                let to_read = probe.len().min(file_len - pos);
                file.read_exact_at(&mut probe[..to_read], pos as u64)?;
                if let Some(rel) = probe[..to_read].iter().position(|&b| b == b'\n') {
                    found = Some(pos + rel + 1);
                    break;
                }
                pos += to_read;
            }
            let next_b = found.unwrap_or(file_len);
            if next_b > *boundaries.last().unwrap() && next_b < file_len {
                boundaries.push(next_b);
            }
        }
        boundaries.push(file_len);

        if boundaries.len() <= 2 {
            let mut f = file.try_clone()?;
            use std::io::Seek;
            f.seek(std::io::SeekFrom::Start(0))?;
            return Self::from_fasta_reader_with_hint(
                BufReader::with_capacity(4 << 20, f),
                name,
                path,
                file_len,
            );
        }

        struct ChunkParsed {
            leading: PackedSeqBuilder,
            records: Vec<(String, String, PackedSeqBuilder)>,
        }

        let chunk_ranges: Vec<(usize, usize)> =
            boundaries.windows(2).map(|w| (w[0], w[1])).collect();

        const READ_BUF_CAP: usize = 256 * 1024;

        let parsed: Vec<ChunkParsed> = chunk_ranges
            .par_iter()
            .map(|&(start, end)| -> Result<ChunkParsed> {
                let total_len = end - start;
                let mut leading = PackedSeqBuilder::with_capacity_bases(0);
                let mut records: Vec<(String, String, PackedSeqBuilder)> = Vec::new();

                let mut buf = vec![0u8; READ_BUF_CAP];
                let mut file_off = start;
                let mut carry = 0usize;
                let mut pred_w = 0usize;

                while file_off < end || carry > 0 {
                    let rem_file = end - file_off;
                    let to_read = if rem_file > 0 {
                        let space = buf.len() - carry;
                        if space == 0 {
                            buf.resize(buf.len() * 2, 0);
                        }
                        let n = (buf.len() - carry).min(rem_file);
                        file.read_exact_at(&mut buf[carry..carry + n], file_off as u64)?;
                        file_off += n;
                        n
                    } else {
                        0
                    };
                    let valid = carry + to_read;
                    let at_end = file_off == end;

                    // Process all complete lines in buf[..valid] (or final partial line if at_end).
                    let limit = if at_end {
                        valid
                    } else if let Some(last_nl) = buf[..valid].iter().rposition(|&b| b == b'\n') {
                        last_nl + 1
                    } else {
                        carry = valid;
                        continue;
                    };

                    let mut pos = 0usize;
                    while pos < limit {
                        if pred_w > 0
                            && pos + pred_w < limit
                            && buf[pos + pred_w] == b'\n'
                            && buf[pos] != b'>'
                        {
                            let target = if let Some((_, _, b)) = records.last_mut() {
                                b
                            } else {
                                if leading.len == 0 {
                                    leading.bits.reserve(total_len / 32 + 2);
                                }
                                &mut leading
                            };
                            target.push_fasta_bytes(&buf[pos..pos + pred_w]);
                            pos += pred_w + 1;
                            continue;
                        }

                        let nl = buf[pos..limit]
                            .iter()
                            .position(|&b| b == b'\n')
                            .map_or(limit, |rel| pos + rel);
                        let has_cr = nl > pos && buf[nl - 1] == b'\r';
                        let line_end = if has_cr { nl - 1 } else { nl };
                        let line = &buf[pos..line_end];

                        if !line.is_empty() {
                            if line[0] == b'>' {
                                let header =
                                    String::from_utf8_lossy(&line[1..]).trim_start().to_string();
                                let id = header.split_whitespace().next().unwrap_or("").to_string();
                                let rem_est = (end - file_off) + (limit - nl);
                                records.push((
                                    id,
                                    header,
                                    PackedSeqBuilder::with_capacity_bases(rem_est),
                                ));
                            } else {
                                if !has_cr && line.len() >= 16 {
                                    pred_w = line.len();
                                }
                                let target = if let Some((_, _, b)) = records.last_mut() {
                                    b
                                } else {
                                    if leading.len == 0 {
                                        leading.bits.reserve(total_len / 32 + 2);
                                    }
                                    &mut leading
                                };
                                target.push_fasta_bytes(line);
                            }
                        }
                        pos = (nl + 1).min(limit);
                    }

                    carry = valid - limit;
                    if carry > 0 {
                        buf.copy_within(limit..valid, 0);
                    }
                    if at_end {
                        break;
                    }
                }

                Ok(ChunkParsed { leading, records })
            })
            .collect::<Result<Vec<_>>>()?;

        let mut merged: Option<PackedSeqBuilder> = None;
        let mut records: Vec<SeqRecord> = Vec::new();
        let mut total_bases: u64 = 0;

        for (ci, chunk) in parsed.into_iter().enumerate() {
            if chunk.leading.len > 0 {
                if ci == 0 || records.is_empty() {
                    return Err(GeckoError::InvalidFasta {
                        path: path.to_path_buf(),
                        message: "sequence data encountered before first '>' header".into(),
                    });
                }
                let added = chunk.leading.len as u64;
                merged.as_mut().unwrap().append_builder(&chunk.leading);
                records.last_mut().unwrap().len += added;
                total_bases += added;
            }
            for (id, header, mut rec_builder) in chunk.records {
                if let Some(m) = merged.as_mut() {
                    m.push_sep();
                    let start_with_sep = m.len as u64;
                    let added = rec_builder.len as u64;
                    m.append_builder(&rec_builder);
                    records.push(SeqRecord {
                        id,
                        header,
                        start_with_sep,
                        start_no_sep: total_bases,
                        len: added,
                    });
                    total_bases += added;
                } else {
                    // Reuse the very first record's builder directly without copying its bits.
                    let extra_words =
                        (file_len / 32 + 4).saturating_sub(rec_builder.bits.capacity());
                    if extra_words > 0 {
                        rec_builder.bits.reserve(extra_words);
                    }
                    let added = rec_builder.len as u64;
                    merged = Some(rec_builder);
                    records.push(SeqRecord {
                        id,
                        header,
                        start_with_sep: 0,
                        start_no_sep: 0,
                        len: added,
                    });
                    total_bases += added;
                }
            }
        }

        let Some(merged) = merged else {
            return Err(GeckoError::InvalidFasta {
                path: path.to_path_buf(),
                message: "no FASTA records found".into(),
            });
        };

        Ok(Self {
            name,
            seq: merged.finish(),
            records,
            total_bases,
        })
    }

    /// Load query and reference FASTAs concurrently across Rayon workers.
    pub fn load_pair_parallel(
        query_path: impl AsRef<Path>,
        ref_path: impl AsRef<Path>,
    ) -> Result<(Self, Self)> {
        let qp = query_path.as_ref();
        let rp = ref_path.as_ref();
        let (q, r) = rayon::join(|| Self::from_fasta_path(qp), || Self::from_fasta_path(rp));
        Ok((q?, r?))
    }

    /// Parse FASTA records from any [`BufRead`] source.
    pub fn from_fasta_reader<R: BufRead>(reader: R, name: String, path: &Path) -> Result<Self> {
        Self::from_fasta_reader_with_hint(reader, name, path, 0)
    }

    fn from_fasta_reader_with_hint<R: BufRead>(
        mut reader: R,
        name: String,
        path: &Path,
        capacity_hint: usize,
    ) -> Result<Self> {
        const CHUNK_CAP: usize = 4 * 1024 * 1024;
        let mut builder = PackedSeqBuilder::with_capacity_bases(capacity_hint);
        let mut records: Vec<SeqRecord> = Vec::new();
        let mut total_bases: u64 = 0;
        let mut in_record = false;

        let mut buf = vec![0u8; CHUNK_CAP];
        let mut carry = 0usize;
        let mut pred_w = 0usize;

        loop {
            let nr = reader.read(&mut buf[carry..])?;
            if nr == 0 {
                break;
            }
            let tot = carry + nr;
            let mut pos = 0usize;

            while pos < tot {
                // Fast path: constant FASTA line width (no memchr scan needed).
                if in_record
                    && pred_w > 0
                    && pos + pred_w < tot
                    && buf[pos + pred_w] == b'\n'
                    && buf[pos] != b'>'
                {
                    let added = builder.push_fasta_bytes(&buf[pos..pos + pred_w]);
                    records.last_mut().unwrap().len += added;
                    total_bases += added;
                    pos += pred_w + 1;
                    continue;
                }

                let Some(rel_nl) = buf[pos..tot].iter().position(|&b| b == b'\n') else {
                    break;
                };
                let nl = pos + rel_nl;
                let has_cr = nl > pos && buf[nl - 1] == b'\r';
                let line_end = if has_cr { nl - 1 } else { nl };
                let line = &buf[pos..line_end];

                if !line.is_empty() {
                    if line[0] == b'>' {
                        if in_record {
                            builder.push_sep();
                        }
                        let header = String::from_utf8_lossy(&line[1..]).trim_start().to_string();
                        let id = header.split_whitespace().next().unwrap_or("").to_string();
                        let start_with_sep = builder.len as u64;
                        records.push(SeqRecord {
                            id,
                            header,
                            start_with_sep,
                            start_no_sep: total_bases,
                            len: 0,
                        });
                        in_record = true;
                    } else {
                        if !in_record {
                            return Err(GeckoError::InvalidFasta {
                                path: path.to_path_buf(),
                                message: "sequence data encountered before first '>' header".into(),
                            });
                        }
                        if !has_cr && line.len() >= 16 {
                            pred_w = line.len();
                        }
                        let added = builder.push_fasta_bytes(line);
                        records.last_mut().unwrap().len += added;
                        total_bases += added;
                    }
                }
                pos = nl + 1;
            }

            carry = tot - pos;
            if carry > 0 {
                if carry == CHUNK_CAP {
                    if !in_record {
                        return Err(GeckoError::InvalidFasta {
                            path: path.to_path_buf(),
                            message: "sequence data encountered before first '>' header".into(),
                        });
                    }
                    if buf[0] != b'>' {
                        let added = builder.push_fasta_bytes(&buf[..carry]);
                        records.last_mut().unwrap().len += added;
                        total_bases += added;
                    }
                    carry = 0;
                } else {
                    buf.copy_within(pos..tot, 0);
                }
            }
        }

        if carry > 0 {
            let line_end = if buf[carry - 1] == b'\r' {
                carry - 1
            } else {
                carry
            };
            let line = &buf[..line_end];
            if !line.is_empty() {
                if line[0] == b'>' {
                    if in_record {
                        builder.push_sep();
                    }
                    let header = String::from_utf8_lossy(&line[1..]).trim_start().to_string();
                    let id = header.split_whitespace().next().unwrap_or("").to_string();
                    let start_with_sep = builder.len as u64;
                    records.push(SeqRecord {
                        id,
                        header,
                        start_with_sep,
                        start_no_sep: total_bases,
                        len: 0,
                    });
                } else {
                    if !in_record {
                        return Err(GeckoError::InvalidFasta {
                            path: path.to_path_buf(),
                            message: "sequence data encountered before first '>' header".into(),
                        });
                    }
                    let added = builder.push_fasta_bytes(line);
                    records.last_mut().unwrap().len += added;
                    total_bases += added;
                }
            }
        }

        if records.is_empty() {
            return Err(GeckoError::InvalidFasta {
                path: path.to_path_buf(),
                message: "no FASTA records found".into(),
            });
        }

        Ok(Self {
            name,
            seq: builder.finish(),
            records,
            total_bases,
        })
    }

    /// Build a `Genome` directly from `(id, sequence)` pairs (useful for tests).
    pub fn from_records<S1: Into<String>, S2: AsRef<[u8]>>(
        name: impl Into<String>,
        seqs: impl IntoIterator<Item = (S1, S2)>,
    ) -> Self {
        let mut builder = PackedSeqBuilder::with_capacity_bases(0);
        let mut records = Vec::new();
        let mut total_bases = 0u64;

        for (i, (id, seq)) in seqs.into_iter().enumerate() {
            if i > 0 {
                builder.push_sep();
            }
            let id_str = id.into();
            let seq_bytes = seq.as_ref();
            let start_with_sep = builder.len as u64;
            let len = builder.push_fasta_bytes(seq_bytes);
            records.push(SeqRecord {
                id: id_str.clone(),
                header: id_str,
                start_with_sep,
                start_no_sep: total_bases,
                len,
            });
            total_bases += len;
        }

        Self {
            name: name.into(),
            seq: builder.finish(),
            records,
            total_bases,
        }
    }

    /// Total coordinate length (bases + `SEQ_SEP` slots).
    #[inline(always)]
    pub fn data_len(&self) -> usize {
        self.seq.len
    }

    /// Return the reverse-complement `Genome` following GECKO's multi-FASTA convention:
    /// record order is reversed (last input record becomes record 0), and each sequence
    /// is reverse-complemented while preserving soft-masking.
    pub fn reverse_complement(&self) -> Self {
        let seq = self.seq.reverse_complement();
        let mut records = Vec::with_capacity(self.records.len());
        let mut total_bases = 0u64;
        let mut cursor = 0u64;

        for (out_idx, rec) in self.records.iter().rev().enumerate() {
            if out_idx > 0 {
                cursor += 1; // SEQ_SEP
            }
            records.push(SeqRecord {
                id: rec.id.clone(),
                header: rec.header.clone(),
                start_with_sep: cursor,
                start_no_sep: total_bases,
                len: rec.len,
            });
            cursor += rec.len;
            total_bases += rec.len;
        }

        Self {
            name: format!("{}-revercomp", self.name),
            seq,
            records,
            total_bases,
        }
    }

    /// Number of FASTA records in this genome.
    #[inline]
    pub fn num_records(&self) -> usize {
        self.records.len()
    }

    /// Primary sequence identifier (first record's ID, matching GECKO `.INF` `SeqX name`).
    pub fn primary_id(&self) -> &str {
        self.records
            .first()
            .map(|r| r.id.as_str())
            .unwrap_or(&self.name)
    }

    /// Map a position in `self.seq` (with separators) to its 0-based record index.
    #[inline]
    pub fn seq_index_for_data_pos(&self, pos: u64) -> usize {
        if self.records.len() == 1 {
            return 0;
        }
        let idx = self.records.partition_point(|r| r.start_with_sep <= pos);
        idx.saturating_sub(1)
    }

    /// Length of the sequence containing `pos` (in `seq` coordinates).
    #[inline]
    pub fn seq_len_at(&self, pos: u64) -> u64 {
        self.records[self.seq_index_for_data_pos(pos)].len
    }

    /// Extract a nucleotide slice `[start_no_sep .. start_no_sep + len]` belonging to `seq_idx`
    /// in uppercase ASCII.
    pub fn slice_bases(&self, seq_idx: usize, start_no_sep: u64, len: u64) -> Vec<u8> {
        let rec = &self.records[seq_idx];
        let local = start_no_sep.saturating_sub(rec.start_no_sep) as usize;
        let data_start = rec.start_with_sep as usize + local;
        let data_end =
            (data_start + len as usize).min(rec.start_with_sep as usize + rec.len as usize);
        (data_start..data_end)
            .map(|i| self.seq.byte_upper_at(i))
            .collect()
    }
}

/// Complement a single DNA/RNA base while preserving case.
#[inline]
pub fn complement_base(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' | b'U' => b'A',
        b'a' => b't',
        b'c' => b'g',
        b'g' => b'c',
        b't' | b'u' => b'a',
        other => other,
    }
}

/// Derive a clean stem from a FASTA path, stripping `.gz` and the sequence extension.
pub fn file_stem_name(path: &Path) -> String {
    let fname = path
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("sequence");
    let stripped = fname.strip_suffix(".gz").unwrap_or(fname);
    match stripped.rsplit_once('.') {
        Some((stem, _)) if !stem.is_empty() => stem.to_string(),
        _ => stripped.to_string(),
    }
}

/// Helper to read a file that might be gzip-compressed.
pub fn open_maybe_gz(path: &Path) -> Result<Box<dyn Read>> {
    let file = File::open(path)?;
    let mut reader = BufReader::with_capacity(1 << 16, file);
    let is_gz = {
        let buf = reader.fill_buf()?;
        buf.len() >= 2 && buf[0] == 0x1f && buf[1] == 0x8b
    };
    if is_gz {
        Ok(Box::new(MultiGzDecoder::new(reader)))
    } else {
        Ok(Box::new(reader))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn parses_multi_fasta_and_reverse_complements() {
        let fasta = b">seqA desc1\nACGTac\nGT\n>seqB desc2\nTTAA\n";
        let g = Genome::from_fasta_reader(Cursor::new(fasta), "test".into(), Path::new("test.fa"))
            .unwrap();
        assert_eq!(g.num_records(), 2);
        assert_eq!(g.total_bases, 12);
        assert_eq!(g.records[0].id, "seqA");
        assert_eq!(g.records[0].len, 8);
        assert_eq!(g.records[0].start_with_sep, 0);
        assert_eq!(g.records[0].start_no_sep, 0);
        assert_eq!(g.records[1].id, "seqB");
        assert_eq!(g.records[1].len, 4);
        assert_eq!(g.records[1].start_with_sep, 9);
        assert_eq!(g.records[1].start_no_sep, 8);
        assert_eq!(g.seq.to_ascii_vec(), b"ACGTacGT\0TTAA");

        let rc = g.reverse_complement();
        assert_eq!(rc.num_records(), 2);
        assert_eq!(rc.records[0].id, "seqB");
        assert_eq!(rc.records[1].id, "seqA");
        assert_eq!(rc.seq.to_ascii_vec(), b"TTAA\0ACgtACGT");
    }

    #[test]
    fn maps_data_pos_to_seq_index() {
        let g = Genome::from_records("g", [("s0", "AAAA"), ("s1", "CCCC"), ("s2", "GG")]);
        // coords: AAAA\0CCCC\0GG (len 12)
        assert_eq!(g.seq_index_for_data_pos(0), 0);
        assert_eq!(g.seq_index_for_data_pos(3), 0);
        assert_eq!(g.seq_index_for_data_pos(5), 1);
        assert_eq!(g.seq_index_for_data_pos(8), 1);
        assert_eq!(g.seq_index_for_data_pos(10), 2);
        assert_eq!(g.seq_index_for_data_pos(11), 2);
    }

    #[test]
    fn packed_seq_load32_rev_and_valid_runs_are_exact() {
        let s0 = "ACGTACGTNNNacgtTTGGCCAAAACCCGGGTTTACGTACGTACGTACGTACGTACGTACGTACGT";
        let s1 = "GGCCttaaNNRYSWACGTACGTACGTACGTACGTACGTACGTACGT";
        let g = Genome::from_records("g", [("s0", s0), ("s1", s1)]);
        let raw = g.seq.to_ascii_vec();
        let expected_raw: Vec<u8> = s0
            .bytes()
            .chain(std::iter::once(SEQ_SEP))
            .chain(s1.bytes())
            .collect();
        assert_eq!(raw, expected_raw);

        // Verify reverse_complement byte-for-byte
        let rc = g.reverse_complement();
        let expected_rc: Vec<u8> = expected_raw
            .iter()
            .rev()
            .map(|&b| complement_base(b))
            .collect();
        assert_eq!(rc.seq.to_ascii_vec(), expected_rc);

        // Verify load_32 and load_32_rev at every offset
        for pos in 0..g.seq.len.saturating_sub(32) {
            let w = g.seq.load_32(pos);
            for k in 0..32 {
                assert_eq!(((w >> (2 * k)) & 3) as u8, g.seq.base_2bit(pos + k));
            }
            let wr = g.seq.load_32_rev(pos + 31);
            for k in 0..32 {
                assert_eq!(((wr >> (2 * k)) & 3) as u8, g.seq.base_2bit(pos + 31 - k));
            }
        }

        // Verify valid runs match manual scan for both soft_mask = true and false
        for soft_mask in [true, false] {
            let mut got = vec![false; g.seq.len];
            g.seq.for_each_valid_run(0, g.seq.len, soft_mask, |s, e| {
                for slot in &mut got[s..e] {
                    *slot = true;
                }
            });
            for (i, &b) in expected_raw.iter().enumerate() {
                let want = if soft_mask {
                    matches!(b, b'A' | b'C' | b'G' | b'T' | b'U')
                } else {
                    matches!(
                        b,
                        b'A' | b'C' | b'G' | b'T' | b'U' | b'a' | b'c' | b'g' | b't' | b'u'
                    )
                };
                assert_eq!(
                    got[i], want,
                    "valid_run mismatch at i={i}, soft_mask={soft_mask}"
                );
            }
            let inv_bm = g.seq.invalid_seed_bitmap(soft_mask).unwrap();
            for (i, &is_valid) in got.iter().enumerate() {
                let bit_set = ((inv_bm[i >> 5] >> (i & 31)) & 1) != 0;
                assert_eq!(
                    bit_set, !is_valid,
                    "invalid_seed_bitmap mismatch at i={i}, soft_mask={soft_mask}"
                );
            }
        }
        let sp_bm = g.seq.special_bitmap().unwrap();
        for (i, &b) in expected_raw.iter().enumerate() {
            let bit_set = ((sp_bm[i >> 5] >> (i & 31)) & 1) != 0;
            let want_sp = !matches!(
                b,
                b'A' | b'C' | b'G' | b'T' | b'U' | b'a' | b'c' | b'g' | b't' | b'u'
            );
            assert_eq!(bit_set, want_sp, "special_bitmap mismatch at i={i}");
        }
        // Multi-word run (> 100 bp) bitmap verification
        let long_n = format!("AAAA{}TTTT{}CCCC", "N".repeat(100), "a".repeat(75));
        let g2 = Genome::from_records("g2", [("s0", long_n.as_bytes())]);
        let inv2 = g2.seq.invalid_seed_bitmap(true).unwrap();
        for (i, &b) in long_n.as_bytes().iter().enumerate() {
            let bit_set = ((inv2[i >> 5] >> (i & 31)) & 1) != 0;
            let want_inv = !matches!(b, b'A' | b'C' | b'G' | b'T');
            assert_eq!(bit_set, want_inv, "multi-word bitmap mismatch at i={i}");
        }
    }

    #[test]
    fn simd_fasta_packing_matches_scalar_exactly() {
        // Exercise 32-bp AVX2, 16-bp SSE2, and scalar paths across every starting bit_pos
        // alignment (0..32) and all line types (pure upper, pure lower, mixed, all-N, U/u, IUPAC).
        let lines: [&[u8]; 9] = [
            b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT", // 80 bp pure upper
            b"acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgt", // 80 bp pure lower
            b"ACGTacgtACGTacgtTTTTggggCCCCaaaaACGTacgtACGTacgtTTTTggggCCCCaaaa", // 64 bp mixed
            b"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN",     // 60 bp all-N
            b"ACGTCCGGTTAAACGTTTGGCCAAAACCCGGGTTTACGTACGGTTAACCGGAATTCCGGAA",    // 61 bp upper
            b"acgtccggttaaacgtttggccaaAACCCGGGtttacgtacggttaaccggaattccgg",      // 59 bp mixed
            b"ACGTNNNNacgtNNNNACGTUUUUacgtuuuuRYSWKM*ACGT", // special + U + IUPAC + sep
            b"GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGgggggggggggggggg", // 32 upper + 16 lower
            b"NNNNNNNNNNNNNNNNACGTACGTACGTACGT",            // 16 N + 16 upper
        ];

        for prefix_len in 0..32usize {
            let prefix = vec![b'C'; prefix_len];
            let mut b_simd = PackedSeqBuilder::with_capacity_bases(1024);
            let mut b_scalar = PackedSeqBuilder::with_capacity_bases(1024);
            b_scalar.use_avx2 = false;

            let a1 = b_simd.push_fasta_bytes(&prefix);
            let a2 = b_scalar.push_fasta_bytes(&prefix);
            assert_eq!(a1, a2);
            for line in &lines {
                assert_eq!(
                    b_simd.push_fasta_bytes(line),
                    b_scalar.push_fasta_bytes(line),
                    "added mismatch for prefix_len={prefix_len}"
                );
            }

            let s_simd = b_simd.finish();
            let s_scalar = b_scalar.finish();
            assert_eq!(s_simd.len, s_scalar.len);
            assert_eq!(
                s_simd.bits, s_scalar.bits,
                "bits mismatch for prefix_len={prefix_len}"
            );
            assert_eq!(s_simd.n_runs, s_scalar.n_runs);
            assert_eq!(s_simd.soft_runs, s_scalar.soft_runs);
            assert_eq!(s_simd.sep_pos, s_scalar.sep_pos);
            assert_eq!(s_simd.iupac, s_scalar.iupac);
            assert_eq!(s_simd.special_runs, s_scalar.special_runs);
        }
    }

    #[cfg(unix)]
    #[test]
    fn parallel_fasta_loader_matches_serial_bit_exact() {
        use std::io::Write;
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        let line_a = "ACGTacgtNNNNACGTCCGGTTAAACGTTTGGCCAAAACCCGGGTTTACGTACGGTTAAC\n";
        let line_b = "NNNNNNNNNNNNNNNNggggccccaaaaTTTTACGTACGTACGTACGTACGTACGTACGT\r\n";
        let line_c = "RYSWKMACGTacgtACGTacgtACGTacgtACGTacgtACGTacgtACGTacgtACGTac\n";
        for rec in 0..5 {
            writeln!(tmp, ">chr{rec} synthetic multi-fasta record {rec}").unwrap();
            for rep in 0..400 {
                let s = match (rec + rep) % 3 {
                    0 => line_a,
                    1 => line_b,
                    _ => line_c,
                };
                tmp.write_all(s.as_bytes()).unwrap();
            }
        }
        tmp.flush().unwrap();

        let path = tmp.path();
        let file = File::open(path).unwrap();
        let file_len = file.metadata().unwrap().len() as usize;
        let par = Genome::from_plain_fasta_parallel(&file, file_len, "t".into(), path).unwrap();

        let file2 = File::open(path).unwrap();
        let ser = Genome::from_fasta_reader_with_hint(
            BufReader::with_capacity(64 * 1024, file2),
            "t".into(),
            path,
            file_len,
        )
        .unwrap();

        assert_eq!(par.total_bases, ser.total_bases);
        assert_eq!(par.records.len(), ser.records.len());
        for (a, b) in par.records.iter().zip(ser.records.iter()) {
            assert_eq!(a.id, b.id);
            assert_eq!(a.header, b.header);
            assert_eq!(a.start_with_sep, b.start_with_sep);
            assert_eq!(a.start_no_sep, b.start_no_sep);
            assert_eq!(a.len, b.len);
        }
        assert_eq!(par.seq.len, ser.seq.len);
        assert_eq!(par.seq.bits, ser.seq.bits);
        assert_eq!(par.seq.n_runs, ser.seq.n_runs);
        assert_eq!(par.seq.soft_runs, ser.seq.soft_runs);
        assert_eq!(par.seq.sep_pos, ser.seq.sep_pos);
        assert_eq!(par.seq.iupac, ser.seq.iupac);
        assert_eq!(par.seq.special_runs, ser.seq.special_runs);
        assert_eq!(par.seq.special_block_bm, ser.seq.special_block_bm);
    }

    #[test]
    fn special_block_bm_exact_forward_and_reverse_complement() {
        // Empty and clean sequences have empty bitmap.
        let empty = Genome::from_records("e", [("e0", b"".as_slice())]);
        assert!(empty.seq.special_block_bm.is_empty());
        assert!(!empty.seq.has_special_in_range(0, 0));

        let clean = Genome::from_records("c", [("c0", vec![b'A'; 2000])]);
        assert!(clean.seq.special_block_bm.is_empty());
        assert!(!clean.seq.has_special_in_short_range(0, 512));
        assert!(!clean.seq.has_special_in_range(0, 2000));

        // Non-512-multiple length with specials at exact block boundaries 511, 512, 513,
        // plus a multi-word run spanning > 64 blocks (32,768+ bp).
        let len = 40_013usize;
        let mut raw = vec![b'A'; len];
        raw[511] = b'N';
        raw[512] = b'R';
        raw[513] = b'N';
        raw[1024] = b'*';
        raw[3_000..37_000].fill(b'N'); // spans multiple u64 bitmap words (> 64 * 512)
        raw[len - 1] = b'N';

        let g = Genome::from_records("s", [("s0", raw.as_slice())]);
        let grc = g.reverse_complement();

        for seq in [&g.seq, &grc.seq] {
            let mut exact_blocks = vec![false; seq.len.div_ceil(SPECIAL_BLOCK_SIZE)];
            for &(s, e) in &seq.special_runs {
                for p in s as usize..e as usize {
                    exact_blocks[p >> SPECIAL_BLOCK_SHIFT] = true;
                }
            }
            for (b, &expected) in exact_blocks.iter().enumerate() {
                let lo = b * SPECIAL_BLOCK_SIZE;
                let hi = ((b + 1) * SPECIAL_BLOCK_SIZE).min(seq.len);
                assert_eq!(
                    seq.has_special_in_short_range(lo, hi),
                    expected,
                    "short-range mismatch on block {b}"
                );
                assert_eq!(
                    seq.has_special_in_range(lo, hi),
                    expected,
                    "range mismatch on block {b}"
                );
            }
            // Cross-boundary checks around 511/512/513 and mirrored positions.
            for lo in [0usize, 500, 511, 512, 513, 1020, 2990, 36_990, len - 64] {
                for span in [1usize, 16, 160, 512] {
                    let hi = (lo + span).min(seq.len);
                    if lo >= hi {
                        continue;
                    }
                    let b0 = lo >> SPECIAL_BLOCK_SHIFT;
                    let b1 = (hi - 1) >> SPECIAL_BLOCK_SHIFT;
                    let expected = exact_blocks[b0..=b1].iter().any(|&x| x);
                    assert_eq!(seq.has_special_in_short_range(lo, hi), expected);
                    assert_eq!(seq.has_special_in_range(lo, hi), expected);
                }
            }
            // Wide multi-word range checks.
            assert!(seq.has_special_in_range(0, seq.len));
            let expected_mid = exact_blocks[(1500 >> 9)..=((2500 - 1) >> 9)]
                .iter()
                .any(|&x| x);
            assert_eq!(seq.has_special_in_range(1500, 2500), expected_mid);
        }
    }

    #[test]
    fn reverse_complement_exact_across_lengths_including_parallel_threshold() {
        let alphabet = *b"ACGT";
        let comp = |b: u8| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => unreachable!(),
        };
        for &len in &[0usize, 1, 31, 32, 33, 64, 65, 2_100_005] {
            let mut raw = Vec::with_capacity(len);
            let mut x = 0x1234_5678_9ABC_DEF0u64 ^ (len as u64);
            for _ in 0..len {
                x = x.wrapping_mul(6364136223846793005).wrapping_add(1);
                raw.push(alphabet[((x >> 32) & 3) as usize]);
            }
            let g = Genome::from_records("rc", [("r0", raw.as_slice())]);
            let rc = g.seq.reverse_complement();
            assert_eq!(rc.len, len);
            assert_eq!(rc.bits.len(), len.div_ceil(32) + 2);
            assert_eq!(rc.bits[rc.bits.len() - 1], 0);
            assert_eq!(rc.bits[rc.bits.len() - 2], 0);
            if len <= 1024 {
                for i in 0..len {
                    assert_eq!(rc.byte_upper_at(i), comp(raw[len - 1 - i]));
                }
            } else {
                // Verify start, boundaries around parallel 16384-word chunks, and tail.
                for &i in &[
                    0,
                    1,
                    31,
                    32,
                    33,
                    524_287,
                    524_288,
                    524_289,
                    1_048_576,
                    len - 33,
                    len - 32,
                    len - 1,
                ] {
                    assert_eq!(rc.byte_upper_at(i), comp(raw[len - 1 - i]), "pos {i}");
                }
                // Double reverse complement must be bit-exact with original.
                let rc2 = rc.reverse_complement();
                assert_eq!(rc2.bits, g.seq.bits);
            }
        }
    }
}
