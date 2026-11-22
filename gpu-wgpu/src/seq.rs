//! Multi-FASTA sequence loading, soft-masking preservation, reverse-complement
//! generation, and 4-byte-packed u32 view for WGSL storage buffers.

use std::fs::File;
use std::io::{BufReader, Read};
use std::path::Path;

pub struct FastSeq {
    pub ident: String,
    pub data: Vec<u8>,       // raw ASCII bytes, 0x00 between multi-FASTA records
    pub buf_len: usize,      // unpadded length including (n_seqs - 1) separator bytes
    pub n_bases: u64,        // nucleotide count excluding separators
    pub n_seqs: u64,
    pub seq_starts: Vec<u64>,
}

#[inline]
pub fn base_to_2bit(c: u8) -> u8 {
    match c {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 4,
    }
}

#[inline]
pub fn comp_base(c: u8) -> u8 {
    match c {
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

impl FastSeq {
    pub fn load_fasta(path: &Path) -> Result<Self, String> {
        let file = File::open(path).map_err(|e| format!("Cannot open {}: {e}", path.display()))?;
        let fsize = file
            .metadata()
            .map_err(|e| format!("Cannot stat {}: {e}", path.display()))?
            .len() as usize;
        if fsize == 0 {
            return Err(format!("Empty FASTA file: {}", path.display()));
        }

        let mut reader = BufReader::with_capacity(4 * 1024 * 1024, file);
        let mut raw = Vec::with_capacity(fsize);
        reader
            .read_to_end(&mut raw)
            .map_err(|e| format!("Read error {}: {e}", path.display()))?;

        let mut data = Vec::with_capacity(fsize + 16);
        let mut seq_starts = Vec::new();
        let mut ident = String::new();
        let mut in_header = false;
        let mut capture_first_id = false;
        let mut seen_first_gt = false;

        for &c in &raw {
            if in_header {
                if c == b'\n' {
                    in_header = false;
                    capture_first_id = false;
                } else if capture_first_id {
                    if ident.is_empty() && c == b' ' {
                        continue;
                    } else if c == b' ' || c == b'\r' || c == b'\t' {
                        capture_first_id = false;
                    } else if ident.len() < 200 {
                        ident.push(c as char);
                    }
                }
                continue;
            }
            if c == b'>' {
                in_header = true;
                if !seen_first_gt {
                    seen_first_gt = true;
                    capture_first_id = true;
                    seq_starts.push(0);
                } else {
                    data.push(0x00);
                    seq_starts.push(data.len() as u64);
                }
                continue;
            }
            if !seen_first_gt {
                continue;
            }
            if c.is_ascii_alphabetic() {
                data.push(c);
            } else if c == b'*' {
                data.push(0x00);
            }
        }

        let buf_len = data.len();
        let n_seqs = seq_starts.len() as u64;
        let n_bases = if n_seqs > 0 {
            buf_len as u64 - (n_seqs - 1)
        } else {
            0
        };
        seq_starts.push(buf_len as u64 + 1);

        Ok(Self {
            ident,
            data,
            buf_len,
            n_bases,
            n_seqs,
            seq_starts,
        })
    }

    pub fn reverse_complement(&self) -> Self {
        let mut rc_data = vec![0u8; self.buf_len];
        let mut seq_starts = Vec::with_capacity(self.seq_starts.len());
        seq_starts.push(0);
        for i in 0..self.buf_len {
            let b = comp_base(self.data[self.buf_len - 1 - i]);
            rc_data[i] = b;
            if b == 0x00 {
                seq_starts.push((i + 1) as u64);
            }
        }
        seq_starts.push(self.buf_len as u64 + 1);

        Self {
            ident: self.ident.clone(),
            data: rc_data,
            buf_len: self.buf_len,
            n_bases: self.n_bases,
            n_seqs: self.n_seqs,
            seq_starts,
        }
    }

    /// Pack ASCII bytes into little-endian u32 words (4 bases per u32) for WGSL storage buffers.
    pub fn pack_u32_words(&self) -> Vec<u32> {
        let n_words = (self.buf_len + 3) / 4 + 2;
        let mut words = vec![0u32; n_words];
        for (i, &b) in self.data.iter().enumerate() {
            words[i >> 2] |= (b as u32) << ((i & 3) * 8);
        }
        words
    }

    #[inline]
    pub fn find_seq_id(&self, pos: u64) -> u64 {
        if self.n_seqs <= 1 {
            return 0;
        }
        let mut lo = 0usize;
        let mut hi = (self.n_seqs - 1) as usize;
        while lo < hi {
            let mid = lo + ((hi - lo + 1) >> 1);
            if self.seq_starts[mid] <= pos {
                lo = mid;
            } else {
                hi = mid - 1;
            }
        }
        lo as u64
    }
}
