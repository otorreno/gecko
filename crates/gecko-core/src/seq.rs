use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use flate2::read::MultiGzDecoder;

use crate::error::{GeckoError, Result};

/// Separator byte placed between concatenated multi-FASTA records in coordinate space.
/// Stored as `0x00` so uppercasing (`b & 0xDF`) preserves separators identically.
pub const SEQ_SEP: u8 = 0u8;

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
        for w in 0..full_words {
            let src_start = len - (w + 1) * 32;
            let fwd = self.load_32(src_start);
            // Complementing 2-bit bases (A=00<->T=11, C=01<->G=10) is bitwise NOT!
            bits.push(!reverse_32_2bit(fwd));
        }
        let rem = len % 32;
        if rem > 0 {
            let mut last = 0u64;
            for k in 0..rem {
                let src_pos = rem - 1 - k;
                let comp = (self.base_2bit(src_pos) ^ 3) as u64;
                last |= comp << (2 * k);
            }
            bits.push(last);
        }
        bits.push(0);
        bits.push(0);

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
}

impl PackedSeqBuilder {
    pub fn with_capacity_bases(cap_bases: usize) -> Self {
        Self {
            len: 0,
            bits: Vec::with_capacity(cap_bases.div_ceil(32) + 2),
            cur_word: 0,
            bit_pos: 0,
            n_runs: Vec::new(),
            soft_runs: Vec::new(),
            sep_pos: Vec::new(),
            iupac: Vec::new(),
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
    fn record_n(&mut self, pos: u64) {
        if let Some(last) = self.n_runs.last_mut()
            && last.1 == pos
        {
            last.1 = pos + 1;
            return;
        }
        self.n_runs.push((pos, pos + 1));
    }

    #[inline(always)]
    fn record_soft(&mut self, pos: u64) {
        if let Some(last) = self.soft_runs.last_mut()
            && last.1 == pos
        {
            last.1 = pos + 1;
            return;
        }
        self.soft_runs.push((pos, pos + 1));
    }

    #[inline]
    pub fn push_sep(&mut self) {
        let pos = self.len as u64;
        self.sep_pos.push(pos);
        self.push_2bit(0);
    }

    /// Push a raw FASTA sequence line/slice, returning the number of nucleotide bases added
    /// (excluding separators).
    #[inline]
    pub fn push_fasta_bytes(&mut self, bytes: &[u8]) -> u64 {
        let mut added = 0u64;
        for &b in bytes {
            match b {
                b'A' => {
                    self.push_2bit(0);
                    added += 1;
                }
                b'C' => {
                    self.push_2bit(1);
                    added += 1;
                }
                b'G' => {
                    self.push_2bit(2);
                    added += 1;
                }
                b'T' | b'U' => {
                    self.push_2bit(3);
                    added += 1;
                }
                b'a' => {
                    let pos = self.len as u64;
                    self.record_soft(pos);
                    self.push_2bit(0);
                    added += 1;
                }
                b'c' => {
                    let pos = self.len as u64;
                    self.record_soft(pos);
                    self.push_2bit(1);
                    added += 1;
                }
                b'g' => {
                    let pos = self.len as u64;
                    self.record_soft(pos);
                    self.push_2bit(2);
                    added += 1;
                }
                b't' | b'u' => {
                    let pos = self.len as u64;
                    self.record_soft(pos);
                    self.push_2bit(3);
                    added += 1;
                }
                b'N' | b'n' => {
                    let pos = self.len as u64;
                    self.record_n(pos);
                    self.push_2bit(0);
                    added += 1;
                }
                b'*' | SEQ_SEP => {
                    self.push_sep();
                }
                other if other.is_ascii_alphabetic() => {
                    let pos = self.len as u64;
                    self.iupac.push((pos, other.to_ascii_uppercase()));
                    self.push_2bit(0);
                    added += 1;
                }
                _ => {}
            }
        }
        added
    }

    pub fn finish(mut self) -> PackedSeq {
        if self.bit_pos > 0 {
            self.bits.push(self.cur_word);
        }
        // Two trailing zero sentinel words for branchless `load_32`.
        self.bits.push(0);
        self.bits.push(0);

        let special_runs = build_special_runs(&self.n_runs, &self.sep_pos, &self.iupac);

        PackedSeq {
            len: self.len,
            bits: self.bits,
            n_runs: self.n_runs,
            soft_runs: self.soft_runs,
            sep_pos: self.sep_pos,
            iupac: self.iupac,
            special_runs,
        }
    }
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
    /// Load a plain or gzip-compressed FASTA file directly into 2-bit/bp storage using 4 MB I/O.
    pub fn from_fasta_path(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let file = File::open(path)?;
        let file_len = file.metadata().map(|m| m.len() as usize).unwrap_or(0);
        let mut reader = BufReader::with_capacity(4 << 20, file);

        // Auto-detect gzip via magic bytes 0x1f 0x8b
        let is_gz = {
            let buf = reader.fill_buf()?;
            buf.len() >= 2 && buf[0] == 0x1f && buf[1] == 0x8b
        };

        let name = file_stem_name(path);
        if is_gz {
            let gz = MultiGzDecoder::new(reader);
            Self::from_fasta_reader_with_hint(
                BufReader::with_capacity(4 << 20, gz),
                name,
                path,
                file_len.saturating_mul(3),
            )
        } else {
            Self::from_fasta_reader_with_hint(reader, name, path, file_len)
        }
    }

    /// Load query and reference FASTAs concurrently across two threads.
    pub fn load_pair_parallel(
        query_path: impl AsRef<Path>,
        ref_path: impl AsRef<Path>,
    ) -> Result<(Self, Self)> {
        let qp = query_path.as_ref();
        let rp = ref_path.as_ref();
        std::thread::scope(|s| {
            let hq = s.spawn(|| Self::from_fasta_path(qp));
            let r = Self::from_fasta_path(rp);
            let q = hq.join().expect("FASTA loader thread panicked");
            Ok((q?, r?))
        })
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
}
