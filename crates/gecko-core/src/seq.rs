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

const COMP_LUT: [u8; 256] = {
    let mut t = [0u8; 256];
    let mut i = 0usize;
    while i < 256 {
        t[i] = i as u8;
        i += 1;
    }
    t[b'A' as usize] = b'T';
    t[b'C' as usize] = b'G';
    t[b'G' as usize] = b'C';
    t[b'T' as usize] = b'A';
    t[b'U' as usize] = b'A';
    t[b'a' as usize] = b't';
    t[b'c' as usize] = b'g';
    t[b'g' as usize] = b'c';
    t[b't' as usize] = b'a';
    t[b'u' as usize] = b'a';
    t
};

#[inline]
fn line_is_pure_alpha(line: &[u8]) -> bool {
    let mut bad = 0u8;
    for &b in line {
        bad |= (((b & 0xDF).wrapping_sub(b'A')) >= 26) as u8;
    }
    bad == 0
}

#[inline]
fn append_seq_line(line: &[u8], data: &mut Vec<u8>, rec_len: &mut u64, total_bases: &mut u64) {
    if line_is_pure_alpha(line) {
        data.extend_from_slice(line);
        let n = line.len() as u64;
        *rec_len += n;
        *total_bases += n;
    } else {
        for &b in line {
            if (b & 0xDF).wrapping_sub(b'A') < 26 {
                data.push(b);
                *rec_len += 1;
                *total_bases += 1;
            } else if b == b'*' || b == SEQ_SEP {
                data.push(SEQ_SEP);
            }
        }
    }
}

impl Genome {
    /// Load a plain or gzip-compressed FASTA file using a 4 MB chunked stride-predicted parser.
    pub fn from_fasta_path(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let file = File::open(path)?;
        let file_len = file.metadata().map(|m| m.len() as usize).unwrap_or(0);
        let mut reader = BufReader::with_capacity(64 * 1024, file);

        // Auto-detect gzip via magic bytes 0x1f 0x8b
        let is_gz = {
            let buf = reader.fill_buf()?;
            buf.len() >= 2 && buf[0] == 0x1f && buf[1] == 0x8b
        };

        let name = file_stem_name(path);
        if is_gz {
            let gz = MultiGzDecoder::new(reader);
            Self::from_raw_reader(gz, name, path, file_len.saturating_mul(3))
        } else {
            Self::from_raw_reader(reader, name, path, file_len)
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
        Self::from_raw_reader(reader, name, path, 0)
    }

    /// High-speed 4 MB chunked + constant-line-width stride-predicted FASTA parser.
    fn from_raw_reader<R: Read>(
        mut reader: R,
        name: String,
        path: &Path,
        capacity_hint: usize,
    ) -> Result<Self> {
        const CHUNK_CAP: usize = 4 * 1024 * 1024;
        let mut buf = vec![0u8; CHUNK_CAP];
        let mut data = Vec::with_capacity(capacity_hint.saturating_add(64));
        let mut records: Vec<SeqRecord> = Vec::new();
        let mut total_bases: u64 = 0;
        let mut in_record = false;
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
                // Ultra-fast path: constant FASTA line width with pure alphabetic bytes.
                if in_record
                    && pred_w > 0
                    && pos + pred_w < tot
                    && buf[pos + pred_w] == b'\n'
                    && line_is_pure_alpha(&buf[pos..pos + pred_w])
                {
                    data.extend_from_slice(&buf[pos..pos + pred_w]);
                    let n = pred_w as u64;
                    records.last_mut().unwrap().len += n;
                    total_bases += n;
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
                            data.push(SEQ_SEP);
                        }
                        let rest = &line[1..];
                        let header = String::from_utf8_lossy(rest).trim_start().to_string();
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
                        if !has_cr && line.len() >= 16 {
                            pred_w = line.len();
                        }
                        let rec = records.last_mut().unwrap();
                        append_seq_line(line, &mut data, &mut rec.len, &mut total_bases);
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
                        let rec = records.last_mut().unwrap();
                        append_seq_line(&buf[..carry], &mut data, &mut rec.len, &mut total_bases);
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
                        data.push(SEQ_SEP);
                    }
                    let rest = &line[1..];
                    let header = String::from_utf8_lossy(rest).trim_start().to_string();
                    let id = header.split_whitespace().next().unwrap_or("").to_string();
                    let start_with_sep = data.len() as u64;
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
                    let rec = records.last_mut().unwrap();
                    append_seq_line(line, &mut data, &mut rec.len, &mut total_bases);
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
            data.extend(self.data[s..e].iter().rev().map(|&b| COMP_LUT[b as usize]));
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
