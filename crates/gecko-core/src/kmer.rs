use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;

use crate::error::{GeckoError, Result};
use crate::radix::parallel_radix_sort_u128_key;
use crate::seq::Genome;

/// Canonical k-mer word size stored in GECKO dictionaries (32 bp packed into 64 bits).
pub const INDEX_WORD_SIZE: usize = 32;

const GIDX_MAGIC: &[u8; 8] = b"GECKOIDX";
const GIDX_VERSION: u32 = 1;

/// Uncompacted k-mer occurrence entry (16 bytes vs 24 bytes in C `wentry`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(C)]
pub struct WordEntry {
    pub kmer: u64,
    pub pos: u64,
}

/// Compact sorted k-mer dictionary supporting prefix queries for any word length `k <= index_k`.
#[derive(Debug, Clone)]
pub struct KmerDict {
    /// Indexed word length (32 for standard sequences, or shorter if genome < 32 bp).
    pub index_k: usize,
    /// Sorted unique 2-bit encoded k-mers (MSB-aligned in `u64`).
    pub words: Vec<u64>,
    /// Start index into `positions` for each word (`offsets.len() == words.len() + 1`).
    pub offsets: Vec<u64>,
    /// Positions in `Genome::data` (with `SEQ_SEP` separators) grouped per unique word.
    pub positions: Vec<u64>,
}

/// Options controlling dictionary construction.
#[derive(Debug, Clone, Copy)]
pub struct DictOptions {
    /// If true (GECKO default), lowercase bases break k-mer windows (soft-masking).
    pub soft_mask: bool,
    /// Maximum memory budget in bytes for sorting before spilling sorted runs to disk.
    pub memory_budget: usize,
    /// Word size to index (default [`INDEX_WORD_SIZE`] = 32).
    pub index_k: usize,
}

impl Default for DictOptions {
    fn default() -> Self {
        Self {
            soft_mask: true,
            memory_budget: 2 * 1024 * 1024 * 1024, // 2 GiB
            index_k: INDEX_WORD_SIZE,
        }
    }
}

/// Summary statistics for a [`KmerDict`].
#[derive(Debug, Clone, serde::Serialize)]
pub struct DictStats {
    pub index_k: usize,
    pub unique_words: usize,
    pub total_occurrences: usize,
    pub max_frequency: usize,
}

impl KmerDict {
    /// Build a sorted, compacted k-mer dictionary from a [`Genome`].
    ///
    /// Automatically sorts in memory when entries fit in `opts.memory_budget`, or falls
    /// back to chunked external sort + k-way merge when exceeding the budget.
    pub fn build(genome: &Genome, opts: DictOptions) -> Result<Self> {
        let index_k = opts.index_k;
        if !(4..=32).contains(&index_k) {
            return Err(GeckoError::InvalidParameter(format!(
                "index word size must be in 4..=32, got {index_k}"
            )));
        }

        // Reserve half of memory_budget for the entries Vec and half for radix aux buffer.
        let max_entries_in_mem =
            (opts.memory_budget / (2 * std::mem::size_of::<WordEntry>())).max(1024);
        let shift_to_msb = 64 - 2 * index_k;
        let kmer_mask: u64 = if index_k == 32 {
            u64::MAX
        } else {
            (1u64 << (2 * index_k)) - 1
        };

        let mut entries: Vec<WordEntry> =
            Vec::with_capacity(genome.seq.len.min(max_entries_in_mem));
        let mut spill_files: Vec<File> = Vec::new();
        let mut io_err: Option<GeckoError> = None;

        genome
            .seq
            .for_each_valid_run(0, genome.seq.len, opts.soft_mask, |rs, re| {
                if io_err.is_some() || re - rs < index_k {
                    return;
                }
                let mut cur_word: u64 = 0;
                for idx in rs..re {
                    let bits = genome.seq.base_2bit(idx) as u64;
                    cur_word = ((cur_word << 2) | bits) & kmer_mask;
                    if idx + 1 >= rs + index_k {
                        let start_pos = (idx + 1 - index_k) as u64;
                        entries.push(WordEntry {
                            kmer: cur_word << shift_to_msb,
                            pos: start_pos,
                        });
                        if entries.len() >= max_entries_in_mem {
                            match sort_and_spill_entries(&mut entries) {
                                Ok(f) => spill_files.push(f),
                                Err(e) => {
                                    io_err = Some(e);
                                    return;
                                }
                            }
                        }
                    }
                }
            });
        if let Some(e) = io_err {
            return Err(e);
        }

        if spill_files.is_empty() {
            parallel_radix_sort_u128_key(&mut entries, |e| (e.kmer, e.pos));
            Ok(Self::compact_sorted(index_k, &entries))
        } else {
            if !entries.is_empty() {
                spill_files.push(sort_and_spill_entries(&mut entries)?);
            }
            drop(entries);
            Self::merge_spilled(index_k, spill_files)
        }
    }

    fn compact_sorted(index_k: usize, sorted: &[WordEntry]) -> Self {
        if sorted.is_empty() {
            return Self {
                index_k,
                words: Vec::new(),
                offsets: vec![0],
                positions: Vec::new(),
            };
        }

        let mut unique_count = 1usize;
        for pair in sorted.windows(2) {
            if pair[0].kmer != pair[1].kmer {
                unique_count += 1;
            }
        }

        let mut words = Vec::with_capacity(unique_count);
        let mut offsets = Vec::with_capacity(unique_count + 1);
        let mut positions = Vec::with_capacity(sorted.len());

        let mut prev_kmer = sorted[0].kmer;
        words.push(prev_kmer);
        offsets.push(0);

        for entry in sorted {
            if entry.kmer != prev_kmer {
                offsets.push(positions.len() as u64);
                words.push(entry.kmer);
                prev_kmer = entry.kmer;
            }
            positions.push(entry.pos);
        }
        offsets.push(positions.len() as u64);

        Self {
            index_k,
            words,
            offsets,
            positions,
        }
    }

    fn merge_spilled(index_k: usize, mut files: Vec<File>) -> Result<Self> {
        use std::io::Seek;

        let mut readers: Vec<BufReader<File>> = Vec::with_capacity(files.len());
        for mut f in files.drain(..) {
            f.rewind()?;
            readers.push(BufReader::with_capacity(256 * 1024, f));
        }

        let mut heap: BinaryHeap<Reverse<(u64, u64, usize)>> = BinaryHeap::new();
        for (i, r) in readers.iter_mut().enumerate() {
            if let Some(entry) = read_word_entry(r)? {
                heap.push(Reverse((entry.kmer, entry.pos, i)));
            }
        }

        let mut words = Vec::new();
        let mut offsets = Vec::new();
        let mut positions = Vec::new();
        let mut current_kmer: Option<u64> = None;

        while let Some(Reverse((kmer, pos, run_idx))) = heap.pop() {
            match current_kmer {
                None => {
                    words.push(kmer);
                    offsets.push(0);
                    current_kmer = Some(kmer);
                }
                Some(prev) if prev != kmer => {
                    offsets.push(positions.len() as u64);
                    words.push(kmer);
                    current_kmer = Some(kmer);
                }
                Some(_) => {}
            }
            positions.push(pos);

            if let Some(next) = read_word_entry(&mut readers[run_idx])? {
                heap.push(Reverse((next.kmer, next.pos, run_idx)));
            }
        }

        if offsets.is_empty() {
            offsets.push(0);
        } else {
            offsets.push(positions.len() as u64);
        }

        Ok(Self {
            index_k,
            words,
            offsets,
            positions,
        })
    }

    /// Number of unique k-mers in the dictionary.
    #[inline]
    pub fn len(&self) -> usize {
        self.words.len()
    }

    /// Whether the dictionary contains zero k-mers.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.words.is_empty()
    }

    /// Occurrence slice for the `i`-th unique word.
    #[inline]
    pub fn positions_for(&self, i: usize) -> &[u64] {
        let start = self.offsets[i] as usize;
        let end = self.offsets[i + 1] as usize;
        &self.positions[start..end]
    }

    /// Compute summary statistics for this dictionary.
    pub fn stats(&self) -> DictStats {
        let max_frequency = (0..self.words.len())
            .map(|i| (self.offsets[i + 1] - self.offsets[i]) as usize)
            .max()
            .unwrap_or(0);
        DictStats {
            index_k: self.index_k,
            unique_words: self.words.len(),
            total_occurrences: self.positions.len(),
            max_frequency,
        }
    }

    /// Persist the dictionary in compact binary `.gidx` format for reuse across runs.
    pub fn save_to_path(&self, path: impl AsRef<Path>) -> Result<()> {
        let file = File::create(path)?;
        let mut w = BufWriter::with_capacity(1 << 20, file);
        w.write_all(GIDX_MAGIC)?;
        w.write_all(&GIDX_VERSION.to_le_bytes())?;
        w.write_all(&(self.index_k as u32).to_le_bytes())?;
        w.write_all(&(self.words.len() as u64).to_le_bytes())?;
        w.write_all(&(self.positions.len() as u64).to_le_bytes())?;

        for &word in &self.words {
            w.write_all(&word.to_le_bytes())?;
        }
        for &off in &self.offsets {
            w.write_all(&off.to_le_bytes())?;
        }
        for &pos in &self.positions {
            w.write_all(&pos.to_le_bytes())?;
        }
        w.flush()?;
        Ok(())
    }

    /// Load a previously saved `.gidx` dictionary from disk.
    pub fn load_from_path(path: impl AsRef<Path>) -> Result<Self> {
        let file = File::open(path)?;
        let mut r = BufReader::with_capacity(1 << 20, file);

        let mut magic = [0u8; 8];
        r.read_exact(&mut magic)?;
        if &magic != GIDX_MAGIC {
            return Err(GeckoError::InvalidFormat(
                "not a valid GECKO .gidx dictionary file".into(),
            ));
        }

        let version = read_u32_le(&mut r)?;
        if version != GIDX_VERSION {
            return Err(GeckoError::InvalidFormat(format!(
                "unsupported .gidx version {version}"
            )));
        }
        let index_k = read_u32_le(&mut r)? as usize;
        let num_words = read_u64_le(&mut r)? as usize;
        let num_positions = read_u64_le(&mut r)? as usize;

        let mut words = vec![0u64; num_words];
        for slot in &mut words {
            *slot = read_u64_le(&mut r)?;
        }

        let mut offsets = vec![0u64; num_words + 1];
        for slot in &mut offsets {
            *slot = read_u64_le(&mut r)?;
        }

        let mut positions = vec![0u64; num_positions];
        for slot in &mut positions {
            *slot = read_u64_le(&mut r)?;
        }

        Ok(Self {
            index_k,
            words,
            offsets,
            positions,
        })
    }
}

/// Extract the top `word_length` bases (as a right-aligned `u64` prefix) from an MSB-aligned k-mer.
#[inline]
pub fn kmer_prefix(msb_kmer: u64, word_length: usize) -> u64 {
    if word_length == 32 {
        msb_kmer
    } else {
        msb_kmer >> (64 - 2 * word_length)
    }
}

/// Decode an MSB-aligned 2-bit k-mer into its ASCII string of length `k`.
pub fn decode_kmer(msb_kmer: u64, k: usize) -> String {
    const ALF: [u8; 4] = *b"ACGT";
    let mut out = Vec::with_capacity(k);
    for i in 0..k {
        let shift = 62 - 2 * i;
        let bits = ((msb_kmer >> shift) & 0b11) as usize;
        out.push(ALF[bits]);
    }
    String::from_utf8(out).unwrap()
}

fn sort_and_spill_entries(entries: &mut Vec<WordEntry>) -> Result<File> {
    parallel_radix_sort_u128_key(entries, |e| (e.kmer, e.pos));
    let mut tmp = tempfile::tempfile()?;
    {
        let mut w = BufWriter::with_capacity(1 << 20, &mut tmp);
        for e in entries.iter() {
            w.write_all(&e.kmer.to_le_bytes())?;
            w.write_all(&e.pos.to_le_bytes())?;
        }
        w.flush()?;
    }
    entries.clear();
    Ok(tmp)
}

fn read_word_entry<R: Read>(r: &mut R) -> Result<Option<WordEntry>> {
    let mut buf = [0u8; 16];
    match r.read_exact(&mut buf) {
        Ok(()) => Ok(Some(WordEntry {
            kmer: u64::from_le_bytes(buf[0..8].try_into().unwrap()),
            pos: u64::from_le_bytes(buf[8..16].try_into().unwrap()),
        })),
        Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => Ok(None),
        Err(e) => Err(GeckoError::Io(e)),
    }
}

fn read_u32_le<R: Read>(r: &mut R) -> Result<u32> {
    let mut b = [0u8; 4];
    r.read_exact(&mut b)?;
    Ok(u32::from_le_bytes(b))
}

fn read_u64_le<R: Read>(r: &mut R) -> Result<u64> {
    let mut b = [0u8; 8];
    r.read_exact(&mut b)?;
    Ok(u64::from_le_bytes(b))
}

/// On-disk sorted k-mer stream storing 16-byte `(kmer: u64, pos: u64)` entries in a
/// temporary file (auto-removed when dropped) so 200+ Mb chromosomes compare in bounded RAM.
#[derive(Debug)]
pub struct DiskKmerIndex {
    pub index_k: usize,
    file: File,
}

/// A k-mer dictionary stored either fully in RAM (fast path) or streamed from a sorted
/// temporary file on disk (out-of-core path for huge chromosomes).
#[derive(Debug)]
pub enum KmerIndex {
    InMemory(KmerDict),
    OnDisk(DiskKmerIndex),
}

impl KmerIndex {
    /// Build a [`KmerIndex`] respecting `opts.memory_budget`: stays in memory if all k-mers
    /// fit in RAM, or spills sorted chunks and k-way merges into a streaming temporary file.
    pub fn build(genome: &Genome, opts: DictOptions) -> Result<Self> {
        let index_k = opts.index_k;
        if !(4..=32).contains(&index_k) {
            return Err(GeckoError::InvalidParameter(format!(
                "index word size must be in 4..=32, got {index_k}"
            )));
        }

        let max_entries_in_mem =
            (opts.memory_budget / (2 * std::mem::size_of::<WordEntry>())).max(1024);
        let shift_to_msb = 64 - 2 * index_k;
        let kmer_mask: u64 = if index_k == 32 {
            u64::MAX
        } else {
            (1u64 << (2 * index_k)) - 1
        };

        let mut entries: Vec<WordEntry> =
            Vec::with_capacity(genome.seq.len.min(max_entries_in_mem));
        let mut spill_files: Vec<File> = Vec::new();
        let mut io_err: Option<GeckoError> = None;

        genome
            .seq
            .for_each_valid_run(0, genome.seq.len, opts.soft_mask, |rs, re| {
                if io_err.is_some() || re - rs < index_k {
                    return;
                }
                let mut cur_word: u64 = 0;
                for idx in rs..re {
                    let bits = genome.seq.base_2bit(idx) as u64;
                    cur_word = ((cur_word << 2) | bits) & kmer_mask;
                    if idx + 1 >= rs + index_k {
                        let start_pos = (idx + 1 - index_k) as u64;
                        entries.push(WordEntry {
                            kmer: cur_word << shift_to_msb,
                            pos: start_pos,
                        });
                        if entries.len() >= max_entries_in_mem {
                            match sort_and_spill_entries(&mut entries) {
                                Ok(f) => spill_files.push(f),
                                Err(e) => {
                                    io_err = Some(e);
                                    return;
                                }
                            }
                        }
                    }
                }
            });
        if let Some(e) = io_err {
            return Err(e);
        }

        if spill_files.is_empty() {
            parallel_radix_sort_u128_key(&mut entries, |e| (e.kmer, e.pos));
            Ok(KmerIndex::InMemory(KmerDict::compact_sorted(
                index_k, &entries,
            )))
        } else {
            if !entries.is_empty() {
                spill_files.push(sort_and_spill_entries(&mut entries)?);
            }
            drop(entries);
            let merged = merge_spilled_to_disk(spill_files)?;
            Ok(KmerIndex::OnDisk(DiskKmerIndex {
                index_k,
                file: merged,
            }))
        }
    }

    #[inline]
    pub fn index_k(&self) -> usize {
        match self {
            KmerIndex::InMemory(d) => d.index_k,
            KmerIndex::OnDisk(d) => d.index_k,
        }
    }

    /// Create a sequential cursor yielding prefix-grouped occurrence lists (skipping words
    /// whose frequency exceeds `freq_limit`).
    pub fn cursor(&self, word_length: usize, freq_limit: usize) -> Result<KmerCursor<'_>> {
        use std::io::Seek;
        match self {
            KmerIndex::InMemory(dict) => Ok(KmerCursor::InMemory {
                dict,
                pos: 0,
                word_length,
                freq_limit,
            }),
            KmerIndex::OnDisk(disk) => {
                let mut cloned = disk.file.try_clone()?;
                cloned.rewind()?;
                let mut reader = BufReader::with_capacity(1 << 20, cloned);
                let peeked = read_word_entry(&mut reader)?;
                Ok(KmerCursor::OnDisk {
                    reader,
                    peeked,
                    word_length,
                    freq_limit,
                })
            }
        }
    }
}

/// Reusable flat buffer holding all distinct 32-mer occurrence lists (passing `freq_limit`)
/// sharing a `word_length` prefix—eliminating per-k-mer heap allocations during streaming merge.
#[derive(Debug, Default)]
pub struct KmerPrefixGroup {
    pub prefix: u64,
    positions: Vec<u64>,
    offsets: Vec<usize>,
}

impl KmerPrefixGroup {
    #[inline]
    pub fn clear(&mut self) {
        self.positions.clear();
        self.offsets.clear();
        self.offsets.push(0);
    }

    #[inline]
    pub fn num_words(&self) -> usize {
        self.offsets.len().saturating_sub(1)
    }

    #[inline]
    pub fn word_positions(&self, idx: usize) -> &[u64] {
        &self.positions[self.offsets[idx]..self.offsets[idx + 1]]
    }
}

/// Sequential prefix-group cursor over either an in-memory [`KmerDict`] or an on-disk
/// [`DiskKmerIndex`].
pub enum KmerCursor<'a> {
    InMemory {
        dict: &'a KmerDict,
        pos: usize,
        word_length: usize,
        freq_limit: usize,
    },
    OnDisk {
        reader: BufReader<File>,
        peeked: Option<WordEntry>,
        word_length: usize,
        freq_limit: usize,
    },
}

impl<'a> KmerCursor<'a> {
    /// Fill `out` with the next prefix group containing at least one word with `freq <= freq_limit`.
    /// Returns `Ok(true)` if a group was produced, or `Ok(false)` at EOF.
    pub fn fill_next_prefix_group(&mut self, out: &mut KmerPrefixGroup) -> Result<bool> {
        match self {
            KmerCursor::InMemory {
                dict,
                pos,
                word_length,
                freq_limit,
            } => {
                let n = dict.words.len();
                let wl = *word_length;
                let fl = *freq_limit;
                while *pos < n {
                    let prefix = kmer_prefix(dict.words[*pos], wl);
                    out.clear();
                    out.prefix = prefix;
                    while *pos < n && kmer_prefix(dict.words[*pos], wl) == prefix {
                        let slice = dict.positions_for(*pos);
                        if slice.len() <= fl {
                            out.positions.extend_from_slice(slice);
                            out.offsets.push(out.positions.len());
                        }
                        *pos += 1;
                    }
                    if out.num_words() > 0 {
                        return Ok(true);
                    }
                }
                Ok(false)
            }
            KmerCursor::OnDisk {
                reader,
                peeked,
                word_length,
                freq_limit,
            } => {
                let wl = *word_length;
                let fl = *freq_limit;
                while let Some(first) = *peeked {
                    let prefix = kmer_prefix(first.kmer, wl);
                    out.clear();
                    out.prefix = prefix;

                    while let Some(cur) = *peeked {
                        if kmer_prefix(cur.kmer, wl) != prefix {
                            break;
                        }
                        let cur_kmer = cur.kmer;
                        let start_len = out.positions.len();
                        let mut count = 0usize;
                        while let Some(e) = *peeked {
                            if e.kmer != cur_kmer {
                                break;
                            }
                            count += 1;
                            if count <= fl {
                                out.positions.push(e.pos);
                            } else if out.positions.len() > start_len {
                                out.positions.truncate(start_len);
                            }
                            *peeked = read_word_entry(reader)?;
                        }
                        if count <= fl {
                            out.offsets.push(out.positions.len());
                        }
                    }

                    if out.num_words() > 0 {
                        return Ok(true);
                    }
                }
                Ok(false)
            }
        }
    }
}

fn merge_spilled_to_disk(mut files: Vec<File>) -> Result<File> {
    use std::io::Seek;

    let mut readers: Vec<BufReader<File>> = Vec::with_capacity(files.len());
    for mut f in files.drain(..) {
        f.rewind()?;
        readers.push(BufReader::with_capacity(256 * 1024, f));
    }

    let mut heap: BinaryHeap<Reverse<(u64, u64, usize)>> = BinaryHeap::new();
    for (i, r) in readers.iter_mut().enumerate() {
        if let Some(entry) = read_word_entry(r)? {
            heap.push(Reverse((entry.kmer, entry.pos, i)));
        }
    }

    let mut out = tempfile::tempfile()?;
    {
        let mut w = BufWriter::with_capacity(1 << 20, &mut out);
        while let Some(Reverse((kmer, pos, run_idx))) = heap.pop() {
            w.write_all(&kmer.to_le_bytes())?;
            w.write_all(&pos.to_le_bytes())?;
            if let Some(next) = read_word_entry(&mut readers[run_idx])? {
                heap.push(Reverse((next.kmer, next.pos, run_idx)));
            }
        }
        w.flush()?;
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn builds_and_roundtrips_dictionary_in_mem_and_spilled() {
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let genome = Genome::from_records("g", [("s0", seq)]);

        let d_mem = KmerDict::build(&genome, DictOptions::default()).unwrap();
        assert!(!d_mem.is_empty());
        assert_eq!(d_mem.stats().total_occurrences, seq.len() - 32 + 1);

        // Force out-of-core spilling with a tiny memory budget
        let d_spill = KmerDict::build(
            &genome,
            DictOptions {
                memory_budget: 16 * 1024,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(d_mem.words, d_spill.words);
        assert_eq!(d_mem.offsets, d_spill.offsets);
        assert_eq!(d_mem.positions, d_spill.positions);

        // Roundtrip .gidx serialization
        let tmp = tempfile::NamedTempFile::new().unwrap();
        d_mem.save_to_path(tmp.path()).unwrap();
        let d_loaded = KmerDict::load_from_path(tmp.path()).unwrap();
        assert_eq!(d_mem.words, d_loaded.words);
        assert_eq!(d_mem.offsets, d_loaded.offsets);
        assert_eq!(d_mem.positions, d_loaded.positions);
    }

    #[test]
    fn soft_masking_skips_lowercase_seeds() {
        let seq = "ACGTACGTACGTACGTacgtACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let g = Genome::from_records("g", [("s0", seq)]);
        let masked = KmerDict::build(
            &g,
            DictOptions {
                soft_mask: true,
                ..Default::default()
            },
        )
        .unwrap();
        let unmasked = KmerDict::build(
            &g,
            DictOptions {
                soft_mask: false,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(masked.stats().total_occurrences < unmasked.stats().total_occurrences);
    }
}
