use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::error::Result;
use crate::formats::frags::Fragment;
use crate::seq::Genome;

const WRAP_WIDTH: usize = 70;

/// Extract formatted pairwise text alignments from fragments (replaces `frags2align.sh` /
/// `frags2text` / `csvFrags2text`).
pub fn write_alignments_path(
    path: impl AsRef<Path>,
    query: &Genome,
    reference: &Genome,
    fragments: &[Fragment],
) -> Result<()> {
    let file = File::create(path)?;
    let mut w = BufWriter::with_capacity(1 << 20, file);
    write_alignments(
        &mut w,
        query,
        reference,
        fragments,
        0,
        None::<&mut std::io::Sink>,
    )?;
    w.flush()?;
    Ok(())
}

/// Extract flanking border alignments (extended by `border_size` bp on both ends, clamped
/// safely to sequence boundaries) and optionally emit extracted query segments as FASTA
/// (replaces `frags2borders.sh` / `csvExtractBorders`).
pub fn write_borders_path(
    path: impl AsRef<Path>,
    query: &Genome,
    reference: &Genome,
    fragments: &[Fragment],
    border_size: u64,
    query_fasta_path: Option<&Path>,
) -> Result<()> {
    let file = File::create(path)?;
    let mut w = BufWriter::with_capacity(1 << 20, file);

    let mut q_writer = match query_fasta_path {
        Some(p) => Some(BufWriter::with_capacity(1 << 20, File::create(p)?)),
        None => None,
    };

    write_alignments(
        &mut w,
        query,
        reference,
        fragments,
        border_size,
        q_writer.as_mut(),
    )?;

    w.flush()?;
    if let Some(ref mut qw) = q_writer {
        qw.flush()?;
    }
    Ok(())
}

/// Core alignment & border writer.
pub fn write_alignments<W: Write, QW: Write>(
    w: &mut W,
    query: &Genome,
    reference: &Genome,
    fragments: &[Fragment],
    border_size: u64,
    mut query_fasta: Option<&mut QW>,
) -> Result<()> {
    let needs_rev = fragments.iter().any(|f| f.strand == 'r');
    let rev_ref = if needs_rev {
        Some(reference.reverse_complement())
    } else {
        None
    };

    for (idx, frag) in fragments.iter().enumerate() {
        let q_idx = frag.seq_x as usize;
        let q_rec = &query.records[q_idx];

        let (y_genome, y_idx) = if frag.strand == 'f' {
            (reference, frag.seq_y as usize)
        } else {
            let rr = rev_ref.as_ref().unwrap();
            let rev_idx = rr.num_records() - 1 - frag.seq_y as usize;
            (rr, rev_idx)
        };
        let y_rec = &y_genome.records[y_idx];

        // Compute safe flank extension clamped to both X and Y record boundaries.
        let q_local_start = frag.x_start.saturating_sub(q_rec.start_no_sep);
        let q_local_end = frag.x_end.saturating_sub(q_rec.start_no_sep);
        let y_local_start = frag.y_start.saturating_sub(y_rec.start_no_sep);
        let y_local_end = frag.y_end.saturating_sub(y_rec.start_no_sep);

        let left_x = q_local_start.min(border_size);
        let left_y = y_local_start.min(border_size);
        let left_ext = left_x.min(left_y);

        let right_x = q_rec.len.saturating_sub(q_local_end + 1).min(border_size);
        let right_y = y_rec.len.saturating_sub(y_local_end + 1).min(border_size);
        let right_ext = right_x.min(right_y);

        let ext_len = frag.length + left_ext + right_ext;
        let seq_x = query.slice_bases(q_idx, frag.x_start - left_ext, ext_len);
        let seq_y = y_genome.slice_bases(y_idx, frag.y_start - left_ext, ext_len);
        let n = seq_x.len().min(seq_y.len());

        let header_line = format!(
            ">{} ALIGNED WITH >{} LENGTH: {} IDENT: {} STRAND: {} @({}, {})",
            q_rec.header,
            y_rec.header,
            frag.length,
            frag.ident,
            frag.strand,
            frag.x_start,
            frag.y_start
        );
        writeln!(w, "{header_line}")?;

        if let Some(ref mut qw) = query_fasta {
            let sanitized = header_line.replace('>', "#");
            writeln!(qw, ">S{idx}-{sanitized}")?;
            for chunk in seq_x[..n].chunks(WRAP_WIDTH) {
                qw.write_all(chunk)?;
                writeln!(qw)?;
            }
        }

        let mut offset = 0usize;
        while offset < n {
            let end = (offset + WRAP_WIDTH).min(n);
            let cx = &seq_x[offset..end];
            let cy = &seq_y[offset..end];

            write!(w, "X:\t")?;
            w.write_all(cx)?;
            write!(w, "\n\t")?;
            for (&a, &b) in cx.iter().zip(cy.iter()) {
                w.write_all(if a == b { b"|" } else { b" " })?;
            }
            write!(w, "\nY:\t")?;
            w.write_all(cy)?;
            write!(w, "\n\n")?;
            offset = end;
        }
        writeln!(w)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn extracts_alignment_and_borders() {
        let q = Genome::from_records("q", [("q0", "TTTTACGTACGTACGTGGGG")]);
        let r = Genome::from_records("r", [("r0", "AAAAACGTACGTACGTCCCC")]);
        let frag = Fragment {
            diag: 0,
            x_start: 4,
            y_start: 4,
            x_end: 15,
            y_end: 15,
            length: 12,
            ident: 12,
            score: 48,
            similarity: 100.0,
            seq_x: 0,
            seq_y: 0,
            block: 0,
            strand: 'f',
        };

        let mut out = Vec::new();
        let mut qfa = Vec::new();
        write_alignments(&mut out, &q, &r, &[frag], 2, Some(&mut qfa)).unwrap();

        let text = String::from_utf8(out).unwrap();
        assert!(text.contains("X:\tTTACGTACGTACGTGG"));
        assert!(text.contains("Y:\tAAACGTACGTACGTCC"));
        assert!(text.contains("  ||||||||||||  "));

        let qfa_str = String::from_utf8(qfa).unwrap();
        assert!(qfa_str.starts_with(">S0-#q0"));
        assert!(qfa_str.contains("TTACGTACGTACGTGG"));
    }
}
