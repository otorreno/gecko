use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::error::{GeckoError, Result};

/// Exact wire size in bytes of a single fragment record in GECKO `.frags` files.
///
/// Layout (big-endian):
/// `diag(i64) + x_start(u64) + y_start(u64) + x_end(u64) + y_end(u64) + length(u64)`
/// `+ ident(u64) + score(u64) + similarity(f32) + seq_x(u64) + seq_y(u64) + block(i64) + strand(u8)`
/// = 11 * 8 + 4 + 1 = 93 bytes.
pub const FRAG_WIRE_SIZE: usize = 93;

/// A single detected High-scoring Segment Pair (HSP), binary-compatible with C `struct FragFile`.
///
/// Note on coordinate conventions (matching GECKO):
/// - `x_start`, `x_end`: 0-based inclusive coordinates in concatenated `query` bases (no separators).
/// - `y_start`, `y_end`: 0-based inclusive coordinates in concatenated `reference` bases (if `strand == 'f'`)
///   or in concatenated reverse-complemented `reference` bases (if `strand == 'r'`).
///   CSV export transforms reverse-strand `y` via `y_total - y - 1`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Fragment {
    pub diag: i64,
    pub x_start: u64,
    pub y_start: u64,
    pub x_end: u64,
    pub y_end: u64,
    pub length: u64,
    pub ident: u64,
    pub score: u64,
    pub similarity: f32,
    pub seq_x: u64,
    pub seq_y: u64,
    pub block: i64,
    pub strand: char,
}

impl Fragment {
    /// Percentage identity (`0.0 ..= 100.0`).
    #[inline]
    pub fn pct_identity(&self) -> f64 {
        if self.length == 0 {
            0.0
        } else {
            100.0 * (self.ident as f64) / (self.length as f64)
        }
    }

    /// Encode into the 93-byte big-endian GECKO `.frags` wire format.
    pub fn to_be_bytes(&self) -> [u8; FRAG_WIRE_SIZE] {
        let mut buf = [0u8; FRAG_WIRE_SIZE];
        buf[0..8].copy_from_slice(&self.diag.to_be_bytes());
        buf[8..16].copy_from_slice(&self.x_start.to_be_bytes());
        buf[16..24].copy_from_slice(&self.y_start.to_be_bytes());
        buf[24..32].copy_from_slice(&self.x_end.to_be_bytes());
        buf[32..40].copy_from_slice(&self.y_end.to_be_bytes());
        buf[40..48].copy_from_slice(&self.length.to_be_bytes());
        buf[48..56].copy_from_slice(&self.ident.to_be_bytes());
        buf[56..64].copy_from_slice(&self.score.to_be_bytes());
        buf[64..68].copy_from_slice(&self.similarity.to_be_bytes());
        buf[68..76].copy_from_slice(&self.seq_x.to_be_bytes());
        buf[76..84].copy_from_slice(&self.seq_y.to_be_bytes());
        buf[84..92].copy_from_slice(&self.block.to_be_bytes());
        buf[92] = self.strand as u8;
        buf
    }

    /// Decode from the 93-byte big-endian GECKO `.frags` wire format.
    pub fn from_be_bytes(buf: &[u8; FRAG_WIRE_SIZE]) -> Result<Self> {
        let strand = buf[92] as char;
        if strand != 'f' && strand != 'r' {
            return Err(GeckoError::InvalidFormat(format!(
                "invalid fragment strand byte 0x{:02x}",
                buf[92]
            )));
        }
        Ok(Self {
            diag: i64::from_be_bytes(buf[0..8].try_into().unwrap()),
            x_start: u64::from_be_bytes(buf[8..16].try_into().unwrap()),
            y_start: u64::from_be_bytes(buf[16..24].try_into().unwrap()),
            x_end: u64::from_be_bytes(buf[24..32].try_into().unwrap()),
            y_end: u64::from_be_bytes(buf[32..40].try_into().unwrap()),
            length: u64::from_be_bytes(buf[40..48].try_into().unwrap()),
            ident: u64::from_be_bytes(buf[48..56].try_into().unwrap()),
            score: u64::from_be_bytes(buf[56..64].try_into().unwrap()),
            similarity: f32::from_be_bytes(buf[64..68].try_into().unwrap()),
            seq_x: u64::from_be_bytes(buf[68..76].try_into().unwrap()),
            seq_y: u64::from_be_bytes(buf[76..84].try_into().unwrap()),
            block: i64::from_be_bytes(buf[84..92].try_into().unwrap()),
            strand,
        })
    }
}

/// Decoded contents of a GECKO `.frags` binary file.
#[derive(Debug, Clone)]
pub struct FragsFile {
    pub x_total: u64,
    pub y_total: u64,
    pub fragments: Vec<Fragment>,
}

/// Write a `.frags` file in exact GECKO big-endian binary format.
pub fn write_frags_path(
    path: impl AsRef<Path>,
    x_total: u64,
    y_total: u64,
    fragments: &[Fragment],
) -> Result<()> {
    let file = File::create(path)?;
    let mut w = BufWriter::with_capacity(1 << 20, file);
    write_frags(&mut w, x_total, y_total, fragments)?;
    w.flush()?;
    Ok(())
}

/// Write a `.frags` stream in exact GECKO big-endian binary format.
pub fn write_frags<W: Write>(
    w: &mut W,
    x_total: u64,
    y_total: u64,
    fragments: &[Fragment],
) -> Result<()> {
    w.write_all(&x_total.to_be_bytes())?;
    w.write_all(&y_total.to_be_bytes())?;
    for frag in fragments {
        w.write_all(&frag.to_be_bytes())?;
    }
    Ok(())
}

/// Read all fragments from a GECKO `.frags` binary file (fixing C `fragmentv2.c`'s
/// 104-byte vs 93-byte truncation bug).
pub fn read_frags_path(path: impl AsRef<Path>) -> Result<FragsFile> {
    let file = File::open(path)?;
    let mut r = BufReader::with_capacity(1 << 20, file);
    read_frags(&mut r)
}

/// Read all fragments from any `.frags` binary stream.
pub fn read_frags<R: Read>(r: &mut R) -> Result<FragsFile> {
    let mut header = [0u8; 16];
    r.read_exact(&mut header)?;
    let x_total = u64::from_be_bytes(header[0..8].try_into().unwrap());
    let y_total = u64::from_be_bytes(header[8..16].try_into().unwrap());

    let mut fragments = Vec::new();
    let mut rec = [0u8; FRAG_WIRE_SIZE];
    while read_exact_or_eof(r, &mut rec)? {
        fragments.push(Fragment::from_be_bytes(&rec)?);
    }

    Ok(FragsFile {
        x_total,
        y_total,
        fragments,
    })
}

fn read_exact_or_eof<R: Read>(r: &mut R, buf: &mut [u8]) -> Result<bool> {
    let mut read = 0;
    while read < buf.len() {
        match r.read(&mut buf[read..])? {
            0 if read == 0 => return Ok(false),
            0 => {
                return Err(GeckoError::InvalidFormat(format!(
                    "truncated .frags record ({read} of {} bytes)",
                    buf.len()
                )));
            }
            n => read += n,
        }
    }
    Ok(true)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn frags_binary_roundtrip_is_exact_93_bytes() {
        let frags = vec![
            Fragment {
                diag: -123,
                x_start: 10,
                y_start: 133,
                x_end: 59,
                y_end: 182,
                length: 50,
                ident: 48,
                score: 184,
                similarity: 92.0,
                seq_x: 0,
                seq_y: 1,
                block: 0,
                strand: 'f',
            },
            Fragment {
                diag: 42,
                x_start: 100,
                y_start: 58,
                x_end: 199,
                y_end: 157,
                length: 100,
                ident: 95,
                score: 360,
                similarity: 90.0,
                seq_x: 2,
                seq_y: 0,
                block: 0,
                strand: 'r',
            },
        ];

        let mut bytes = Vec::new();
        write_frags(&mut bytes, 1000, 2000, &frags).unwrap();
        assert_eq!(bytes.len(), 16 + 2 * FRAG_WIRE_SIZE);

        let decoded = read_frags(&mut Cursor::new(&bytes)).unwrap();
        assert_eq!(decoded.x_total, 1000);
        assert_eq!(decoded.y_total, 2000);
        assert_eq!(decoded.fragments, frags);
    }
}
