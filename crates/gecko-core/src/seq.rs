use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use flate2::read::MultiGzDecoder;

use crate::error::{GeckoError, Result};

/// Separator byte placed between concatenated multi-FASTA records in `Genome::data`.
/// Stored as `0x00` so SIMD uppercasing (`b & 0xDF`) preserves separators identically.
pub const SEQ_SEP: u8 = 0u8;

/// Metadata for a single sequence record inside a [`Genome`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SeqRecord {
    /// Short sequence identifier (first whitespace-delimited token after `>`).
    pub id: String,
    /// Full header line without leading `>` or trailing newline.
    pub header: String,
    /// Start offset in `Genome::data` (includes `SEQ_SEP` bytes from preceding records).
    pub start_with_sep: u64,
    /// Start offset in concatenated nucleotide coordinates (excluding `SEQ_SEP`).
    pub start_no_sep: u64,
    /// Sequence length in bases.
    pub len: u64,
}

/// An in-memory FASTA collection (single or multi-FASTA).
///
/// Sequences are stored contiguously in `data`, separated by [`SEQ_SEP`] (`0x00`),
/// preserving original ASCII case so soft-masked lowercase bases can be skipped
/// during k-mer seeding while still participating in ungapped extension.
#[derive(Debug, Clone)]
pub struct Genome {
    /// Display name (typically the file stem or first record ID).
    pub name: String,
    /// Concatenated sequence bytes with `SEQ_SEP` between records.
    pub data: Vec<u8>,
    /// Per-record metadata in file order.
    pub records: Vec<SeqRecord>,
    /// Total number of bases across all records (excluding separators).
    pub total_bases: u64,
}

impl Genome {
    /// Load a plain or gzip-compressed FASTA file.
    pub fn from_fasta_path(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let file = File::open(path)?;
        let file_len = file.metadata().map(|m| m.len() as usize).unwrap_or(0);
        let mut reader = BufReader::with_capacity(1 << 20, file);

        // Auto-detect gzip via magic bytes 0x1f 0x8b
        let is_gz = {
            let buf = reader.fill_buf()?;
            buf.len() >= 2 && buf[0] == 0x1f && buf[1] == 0x8b
        };

        let name = file_stem_name(path);
        if is_gz {
            let gz = MultiGzDecoder::new(reader);
            Self::from_fasta_reader_with_hint(
                BufReader::with_capacity(1 << 20, gz),
                name,
                path,
                file_len.saturating_mul(3),
            )
        } else {
            Self::from_fasta_reader_with_hint(reader, name, path, file_len)
        }
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
        let mut data = Vec::with_capacity(capacity_hint);
        let mut records: Vec<SeqRecord> = Vec::new();
        let mut total_bases: u64 = 0;
        let mut line = Vec::with_capacity(256);
        let mut in_record = false;

        loop {
            line.clear();
            let n = reader.read_until(b'\n', &mut line)?;
            if n == 0 {
                break;
            }
            let mut end = line.len();
            while end > 0 && (line[end - 1] == b'\n' || line[end - 1] == b'\r') {
                end -= 1;
            }
            if end == 0 {
                continue;
            }
            let trimmed = &line[..end];
            if trimmed[0] == b'>' {
                if in_record {
                    data.push(SEQ_SEP);
                }
                let header_str = String::from_utf8_lossy(&trimmed[1..]);
                let header = header_str.trim_start().to_string();
                let id = header.split_whitespace().next().unwrap_or("").to_string();
                let start_with_sep = data.len() as u64;
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
                let rec = records.last_mut().unwrap();
                // Fast path: most FASTA lines consist purely of ASCII letters.
                if trimmed.iter().all(|b| b.is_ascii_alphabetic()) {
                    data.extend_from_slice(trimmed);
                    let len = trimmed.len() as u64;
                    rec.len += len;
                    total_bases += len;
                } else {
                    for &b in trimmed {
                        if b.is_ascii_alphabetic() {
                            data.push(b);
                            rec.len += 1;
                            total_bases += 1;
                        } else if b == b'*' || b == SEQ_SEP {
                            data.push(SEQ_SEP);
                        }
                    }
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
            data,
            records,
            total_bases,
        })
    }

    /// Build a `Genome` directly from `(id, sequence)` pairs (useful for tests).
    pub fn from_records<S1: Into<String>, S2: AsRef<[u8]>>(
        name: impl Into<String>,
        seqs: impl IntoIterator<Item = (S1, S2)>,
    ) -> Self {
        let mut data = Vec::new();
        let mut records = Vec::new();
        let mut total_bases = 0u64;

        for (i, (id, seq)) in seqs.into_iter().enumerate() {
            if i > 0 {
                data.push(SEQ_SEP);
            }
            let id_str = id.into();
            let seq_bytes = seq.as_ref();
            let start_with_sep = data.len() as u64;
            data.extend_from_slice(seq_bytes);
            let len = seq_bytes.len() as u64;
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
            data,
            records,
            total_bases,
        }
    }

    /// Return the reverse-complement `Genome` following GECKO's multi-FASTA convention:
    /// record order is reversed (last input record becomes record 0), and each sequence
    /// is reverse-complemented while preserving ASCII case.
    pub fn reverse_complement(&self) -> Self {
        let mut data = Vec::with_capacity(self.data.len());
        let mut records = Vec::with_capacity(self.records.len());
        let mut total_bases = 0u64;

        for (out_idx, rec) in self.records.iter().rev().enumerate() {
            if out_idx > 0 {
                data.push(SEQ_SEP);
            }
            let start_with_sep = data.len() as u64;
            let s = rec.start_with_sep as usize;
            let e = s + rec.len as usize;
            for &b in self.data[s..e].iter().rev() {
                data.push(complement_base(b));
            }
            records.push(SeqRecord {
                id: rec.id.clone(),
                header: rec.header.clone(),
                start_with_sep,
                start_no_sep: total_bases,
                len: rec.len,
            });
            total_bases += rec.len;
        }

        Self {
            name: format!("{}-revercomp", self.name),
            data,
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

    /// Map a position in `self.data` (with separators) to its 0-based record index.
    #[inline]
    pub fn seq_index_for_data_pos(&self, pos: u64) -> usize {
        if self.records.len() == 1 {
            return 0;
        }
        let idx = self.records.partition_point(|r| r.start_with_sep <= pos);
        idx.saturating_sub(1)
    }

    /// Length of the sequence containing `pos` (in `data` coordinates).
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
        self.data[data_start..data_end]
            .iter()
            .map(|b| b.to_ascii_uppercase())
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
        assert_eq!(&g.data, b"ACGTacGT\0TTAA");

        let rc = g.reverse_complement();
        assert_eq!(rc.num_records(), 2);
        assert_eq!(rc.records[0].id, "seqB");
        assert_eq!(rc.records[1].id, "seqA");
        assert_eq!(&rc.data, b"TTAA\0ACgtACGT");
    }

    #[test]
    fn maps_data_pos_to_seq_index() {
        let g = Genome::from_records("g", [("s0", "AAAA"), ("s1", "CCCC"), ("s2", "GG")]);
        // data: AAAA*CCCC*GG (len 12)
        assert_eq!(g.seq_index_for_data_pos(0), 0);
        assert_eq!(g.seq_index_for_data_pos(3), 0);
        assert_eq!(g.seq_index_for_data_pos(5), 1);
        assert_eq!(g.seq_index_for_data_pos(8), 1);
        assert_eq!(g.seq_index_for_data_pos(10), 2);
        assert_eq!(g.seq_index_for_data_pos(11), 2);
    }
}
