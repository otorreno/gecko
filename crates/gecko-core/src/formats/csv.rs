use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use crate::error::{GeckoError, Result};
use crate::extend::CoverageMatrix;
use crate::formats::frags::{Fragment, FragsFile};

/// Metadata written into GECKO `.INF` and `.csv` headers (compatible with GECKO-MGV).
#[derive(Debug, Clone)]
pub struct ComparisonMetadata {
    pub seq_x_filename: String,
    pub seq_y_filename: String,
    pub seq_x_name: String,
    pub seq_y_name: String,
    pub seq_x_len: u64,
    pub seq_y_len: u64,
    pub min_length: u64,
    pub min_similarity: u64,
    pub total_hits: u64,
    pub total_hits_used: u64,
    pub total_fragments: u64,
}

/// Write the legacy `.INF` metadata block.
pub fn write_inf<W: Write>(w: &mut W, meta: &ComparisonMetadata) -> Result<()> {
    writeln!(
        w,
        "All by-Identity Ungapped Fragments (Hits based approach)"
    )?;
    writeln!(w, "[GECKO v2 (Rust) -- compatible with GECKO-MGV]")?;
    writeln!(w, "SeqX filename        : {}", meta.seq_x_filename)?;
    writeln!(w, "SeqY filename        : {}", meta.seq_y_filename)?;
    writeln!(w, "SeqX name            : {}", meta.seq_x_name)?;
    writeln!(w, "SeqY name            : {}", meta.seq_y_name)?;
    writeln!(w, "SeqX length          : {}", meta.seq_x_len)?;
    writeln!(w, "SeqY length          : {}", meta.seq_y_len)?;
    writeln!(w, "Min.fragment.length  : {}", meta.min_length)?;
    writeln!(w, "Min.Identity         : {}", meta.min_similarity)?;
    writeln!(w, "Tot Hits (seeds)     : {}", meta.total_hits)?;
    writeln!(w, "Tot Hits (seeds) used: {}", meta.total_hits_used)?;
    writeln!(w, "Total fragments      : {}", meta.total_fragments)?;
    writeln!(
        w,
        "========================================================"
    )?;
    Ok(())
}

/// Write `.INF` to a file path.
pub fn write_inf_path(path: impl AsRef<Path>, meta: &ComparisonMetadata) -> Result<()> {
    let file = File::create(path)?;
    let mut w = BufWriter::new(file);
    write_inf(&mut w, meta)?;
    w.flush()?;
    Ok(())
}

#[inline]
fn push_u64(buf: &mut Vec<u8>, mut v: u64) {
    let mut tmp = [0u8; 20];
    let mut n = 0usize;
    loop {
        tmp[n] = b'0' + (v % 10) as u8;
        v /= 10;
        n += 1;
        if v == 0 {
            break;
        }
    }
    while n > 0 {
        n -= 1;
        buf.push(tmp[n]);
    }
}

#[inline]
fn push_i64(buf: &mut Vec<u8>, v: i64) {
    if v < 0 {
        buf.push(b'-');
        push_u64(buf, v.unsigned_abs());
    } else {
        push_u64(buf, v as u64);
    }
}

#[inline]
fn push_fixed2(buf: &mut Vec<u8>, mut val: f64) {
    if val < 0.0 {
        buf.push(b'-');
        val = -val;
    }
    let scaled = (val * 100.0 + 0.5) as u64;
    push_u64(buf, scaled / 100);
    buf.push(b'.');
    let frac = scaled % 100;
    buf.push(b'0' + (frac / 10) as u8);
    buf.push(b'0' + (frac % 10) as u8);
}

/// Write the 1000 x 100 coverage matrix `.MAT` file.
pub fn write_mat_path(path: impl AsRef<Path>, mat: &CoverageMatrix) -> Result<()> {
    let file = File::create(path)?;
    let mut w = BufWriter::with_capacity(4 * 1024 * 1024, file);
    let mut buf = Vec::with_capacity(256 * 1024);
    for row in &mat.bins {
        if buf.len() + 2048 > buf.capacity() {
            w.write_all(&buf)?;
            buf.clear();
        }
        for &val in row {
            push_u64(&mut buf, val);
            buf.push(b'\t');
        }
        buf.push(b'\n');
    }
    if !buf.is_empty() {
        w.write_all(&buf)?;
    }
    w.flush()?;
    Ok(())
}

/// Write a drop-in GECKO-MGV compatible `.csv` file (metadata header + fragment rows).
pub fn write_csv_path(
    path: impl AsRef<Path>,
    meta: &ComparisonMetadata,
    fragments: &[Fragment],
) -> Result<()> {
    let file = File::create(path)?;
    let mut w = BufWriter::with_capacity(4 * 1024 * 1024, file);
    write_csv(&mut w, meta, fragments)?;
    w.flush()?;
    Ok(())
}

/// Write a drop-in GECKO-MGV compatible `.csv` stream.
pub fn write_csv<W: Write>(
    w: &mut W,
    meta: &ComparisonMetadata,
    fragments: &[Fragment],
) -> Result<()> {
    write_inf(w, meta)?;
    w.write_all(
        b"Total CSB: 0\n\
          ========================================================\n\
          Type,xStart,yStart,xEnd,yEnd,strand(f/r),block,length,score,ident,similarity,%ident,SeqX,SeqY\n",
    )?;

    const FLUSH_CAP: usize = 4 * 1024 * 1024;
    let mut buf = Vec::with_capacity(FLUSH_CAP.min(fragments.len().saturating_mul(96).max(4096)));
    let y_total = meta.seq_y_len;

    for f in fragments {
        if buf.len() + 256 > buf.capacity() && buf.len() >= 1024 * 1024 {
            w.write_all(&buf)?;
            buf.clear();
        }
        let (y_start, y_end) = if f.strand == 'r' {
            (
                y_total.saturating_sub(f.y_start + 1),
                y_total.saturating_sub(f.y_end + 1),
            )
        } else {
            (f.y_start, f.y_end)
        };
        buf.extend_from_slice(b"Frag,");
        push_u64(&mut buf, f.x_start);
        buf.push(b',');
        push_u64(&mut buf, y_start);
        buf.push(b',');
        push_u64(&mut buf, f.x_end);
        buf.push(b',');
        push_u64(&mut buf, y_end);
        buf.push(b',');
        buf.push(f.strand as u8);
        buf.push(b',');
        push_i64(&mut buf, f.block);
        buf.push(b',');
        push_u64(&mut buf, f.length);
        buf.push(b',');
        push_u64(&mut buf, f.score);
        buf.push(b',');
        push_u64(&mut buf, f.ident);
        buf.push(b',');
        push_fixed2(&mut buf, f.similarity as f64);
        buf.push(b',');
        push_fixed2(&mut buf, f.pct_identity());
        buf.push(b',');
        push_u64(&mut buf, f.seq_x);
        buf.push(b',');
        push_u64(&mut buf, f.seq_y);
        buf.push(b'\n');
    }
    if !buf.is_empty() {
        w.write_all(&buf)?;
    }
    Ok(())
}

/// Parse a GECKO `.csv` file back into [`FragsFile`] (converting reverse-strand `y` coordinates
/// back into internal `.frags` reverse-complement coordinates so `align` and `borders` work
/// identically on both `.frags` and `.csv` inputs).
pub fn read_csv_path(path: impl AsRef<Path>) -> Result<FragsFile> {
    let file = File::open(path)?;
    read_csv(BufReader::with_capacity(1 << 20, file))
}

pub fn read_csv<R: BufRead>(reader: R) -> Result<FragsFile> {
    let mut x_total: Option<u64> = None;
    let mut y_total: Option<u64> = None;
    let mut fragments = Vec::new();

    for line_res in reader.lines() {
        let line = line_res?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if let Some(rest) = trimmed.strip_prefix("SeqX length") {
            if let Some((_, val)) = rest.split_once(':') {
                x_total = val.trim().parse::<u64>().ok();
            }
            continue;
        }
        if let Some(rest) = trimmed.strip_prefix("SeqY length") {
            if let Some((_, val)) = rest.split_once(':') {
                y_total = val.trim().parse::<u64>().ok();
            }
            continue;
        }
        if !trimmed.starts_with("Frag,") && !trimmed.starts_with("Frag\t") {
            continue;
        }

        let cols: Vec<&str> = trimmed
            .split([',', '\t', ' '])
            .filter(|s| !s.is_empty())
            .collect();
        if cols.len() < 14 {
            return Err(GeckoError::InvalidFormat(format!(
                "expected 14 columns in CSV fragment row, got {}: {trimmed}",
                cols.len()
            )));
        }

        let x_start: u64 = parse_field(cols[1], "xStart")?;
        let mut y_start: u64 = parse_field(cols[2], "yStart")?;
        let x_end: u64 = parse_field(cols[3], "xEnd")?;
        let mut y_end: u64 = parse_field(cols[4], "yEnd")?;
        let strand = cols[5].chars().next().unwrap_or('f');
        let block: i64 = parse_field(cols[6], "block")?;
        let length: u64 = parse_field(cols[7], "length")?;
        let score: u64 = parse_field(cols[8], "score")?;
        let ident: u64 = parse_field(cols[9], "ident")?;
        let similarity: f32 = parse_field(cols[10], "similarity")?;
        let seq_x: u64 = parse_field(cols[12], "SeqX")?;
        let seq_y: u64 = parse_field(cols[13], "SeqY")?;

        // Convert CSV reverse-strand coordinates (where y_start >= y_end in forward coords)
        // back into internal reverse-complement coordinates.
        if strand == 'r'
            && let Some(yt) = y_total
        {
            y_start = yt.saturating_sub(y_start + 1);
            y_end = yt.saturating_sub(y_end + 1);
        }

        let diag = x_start as i64 - y_start as i64;
        fragments.push(Fragment {
            diag,
            x_start,
            y_start,
            x_end,
            y_end,
            length,
            ident,
            score,
            similarity,
            seq_x,
            seq_y,
            block,
            strand,
        });
    }

    Ok(FragsFile {
        x_total: x_total.unwrap_or(0),
        y_total: y_total.unwrap_or(0),
        fragments,
    })
}

fn parse_field<T: std::str::FromStr>(s: &str, name: &str) -> Result<T> {
    s.parse::<T>().map_err(|_| {
        GeckoError::InvalidFormat(format!("failed to parse CSV field '{name}': '{s}'"))
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn csv_roundtrips_forward_and_reverse_coordinates() {
        let meta = ComparisonMetadata {
            seq_x_filename: "q.fa".into(),
            seq_y_filename: "r.fa".into(),
            seq_x_name: "q0".into(),
            seq_y_name: "r0".into(),
            seq_x_len: 500,
            seq_y_len: 1000,
            min_length: 40,
            min_similarity: 60,
            total_hits: 100,
            total_hits_used: 20,
            total_fragments: 2,
        };
        let frags = vec![
            Fragment {
                diag: 0,
                x_start: 10,
                y_start: 10,
                x_end: 59,
                y_end: 59,
                length: 50,
                ident: 50,
                score: 200,
                similarity: 100.0,
                seq_x: 0,
                seq_y: 0,
                block: 0,
                strand: 'f',
            },
            Fragment {
                diag: -10,
                x_start: 20,
                y_start: 30,
                x_end: 69,
                y_end: 79,
                length: 50,
                ident: 45,
                score: 160,
                similarity: 80.0,
                seq_x: 0,
                seq_y: 0,
                block: 0,
                strand: 'r',
            },
        ];

        let mut buf = Vec::new();
        write_csv(&mut buf, &meta, &frags).unwrap();
        let text = String::from_utf8(buf).unwrap();
        // Reverse-strand yStart=30, yEnd=79 on y_total=1000 => 969, 920 in CSV
        assert!(text.contains("Frag,20,969,69,920,r,0,50,160,45,80.00,90.00,0,0"));

        let parsed = read_csv(Cursor::new(text.as_bytes())).unwrap();
        assert_eq!(parsed.x_total, 500);
        assert_eq!(parsed.y_total, 1000);
        assert_eq!(parsed.fragments, frags);
    }
}
