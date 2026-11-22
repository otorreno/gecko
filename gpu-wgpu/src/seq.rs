//! Multi-FASTA sequence loading, soft-masking preservation, reverse-complement
//! generation, and 4-byte-packed u32 view for WGSL storage buffers.

use std::fs::File;
use std::io::Read;
use std::path::Path;

pub struct FastSeq {
    pub ident: String,
    pub data: Vec<u8>,       // raw ASCII bytes, 0x00 between multi-FASTA records
    pub buf_len: usize,      // unpadded length including (n_seqs - 1) separator bytes
    pub n_bases: u64,        // nucleotide count excluding separators
    pub n_seqs: u64,
    pub seq_starts: Vec<u64>,
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
fn line_is_pure_alpha(line: &[u8]) -> bool {
    let mut bad = 0u8;
    for &b in line {
        bad |= (((b & 0xDF).wrapping_sub(b'A')) >= 26) as u8;
    }
    bad == 0
}

#[inline]
fn append_seq_line(line: &[u8], data: &mut Vec<u8>) {
    if line_is_pure_alpha(line) {
        data.extend_from_slice(line);
    } else {
        for &c in line {
            if (c & 0xDF).wrapping_sub(b'A') < 26 {
                data.push(c);
            } else if c == b'*' {
                data.push(0x00);
            }
        }
    }
}

impl FastSeq {
    pub fn load_fasta(path: &Path) -> Result<Self, String> {
        let mut file =
            File::open(path).map_err(|e| format!("Cannot open {}: {e}", path.display()))?;
        let fsize = file
            .metadata()
            .map_err(|e| format!("Cannot stat {}: {e}", path.display()))?
            .len() as usize;
        if fsize == 0 {
            return Err(format!("Empty FASTA file: {}", path.display()));
        }

        const CHUNK_CAP: usize = 4 * 1024 * 1024;
        let mut buf = vec![0u8; CHUNK_CAP];
        let mut data = Vec::with_capacity(fsize + 16);
        let mut seq_starts = Vec::new();
        let mut ident = String::new();
        let mut seen_first_gt = false;
        let mut carry = 0usize;
        let mut pred_w = 0usize;

        loop {
            let nr = file
                .read(&mut buf[carry..])
                .map_err(|e| format!("Read error {}: {e}", path.display()))?;
            if nr == 0 {
                break;
            }
            let tot = carry + nr;
            let mut pos = 0usize;

            while pos < tot {
                // Fast path: constant FASTA line width skips newline scan on 99.99% of lines
                if seen_first_gt
                    && pred_w > 0
                    && pos + pred_w < tot
                    && buf[pos + pred_w] == b'\n'
                    && line_is_pure_alpha(&buf[pos..pos + pred_w])
                {
                    data.extend_from_slice(&buf[pos..pos + pred_w]);
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

                if !line.is_empty() && line[0] == b'>' {
                    if !seen_first_gt {
                        seen_first_gt = true;
                        seq_starts.push(0);
                        let mut k = 1usize;
                        while k < line.len() && line[k] == b' ' {
                            k += 1;
                        }
                        let start = k;
                        while k < line.len()
                            && line[k] != b' '
                            && line[k] != b'\t'
                            && (k - start) < 200
                        {
                            k += 1;
                        }
                        ident.push_str(&String::from_utf8_lossy(&line[start..k]));
                    } else {
                        data.push(0x00);
                        seq_starts.push(data.len() as u64);
                    }
                } else if seen_first_gt && !line.is_empty() {
                    if !has_cr && line.len() >= 16 {
                        pred_w = line.len();
                    }
                    append_seq_line(line, &mut data);
                }
                pos = nl + 1;
            }

            carry = tot - pos;
            if carry > 0 {
                if carry == CHUNK_CAP {
                    if seen_first_gt && buf[0] != b'>' {
                        append_seq_line(&buf[..carry], &mut data);
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
            if seen_first_gt && line_end > 0 && buf[0] != b'>' {
                append_seq_line(&buf[..line_end], &mut data);
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
        let mut rc_data = Vec::with_capacity(self.buf_len);
        rc_data.extend(self.data.iter().rev().map(|&b| COMP_LUT[b as usize]));

        let mut seq_starts = Vec::with_capacity(self.seq_starts.len());
        seq_starts.push(0);
        for k in (2..=self.n_seqs as usize).rev() {
            let orig_start = self.seq_starts[k - 1];
            seq_starts.push(self.buf_len as u64 - orig_start + 1);
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
        let mut words = Vec::with_capacity(n_words);
        let chunks = self.data.chunks_exact(4);
        let rem = chunks.remainder();
        for ch in chunks {
            words.push(u32::from_le_bytes([ch[0], ch[1], ch[2], ch[3]]));
        }
        if !rem.is_empty() {
            let mut last = [0u8; 4];
            last[..rem.len()].copy_from_slice(rem);
            words.push(u32::from_le_bytes(last));
        }
        words.resize(n_words, 0);
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
