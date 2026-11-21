use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use serde::Serialize;

use crate::error::Result;
use crate::formats::frags::Fragment;
use crate::seq::Genome;

/// Per-sequence local 0-based half-open `[start, end)` coordinates for a fragment.
#[derive(Debug, Clone, Serialize)]
pub struct LocalCoords<'a> {
    pub query_id: &'a str,
    pub query_len: u64,
    pub query_start: u64,
    pub query_end: u64,
    pub target_id: &'a str,
    pub target_len: u64,
    pub target_start: u64,
    pub target_end: u64,
    pub strand: char, // '+' or '-'
}

/// Resolve a GECKO [`Fragment`] into per-sequence 0-based half-open coordinates on both
/// `query` and `reference` (converting reverse-strand coordinates onto the forward reference
/// sequence so PAF, BLAST, BEDPE, and JSONL consumers get standard genomic intervals).
pub fn resolve_local_coords<'a>(
    frag: &Fragment,
    query: &'a Genome,
    reference: &'a Genome,
) -> LocalCoords<'a> {
    let q_rec = &query.records[frag.seq_x as usize];
    let t_rec = &reference.records[frag.seq_y as usize];

    let query_start = frag.x_start.saturating_sub(q_rec.start_no_sep);
    let query_end = (frag.x_end.saturating_sub(q_rec.start_no_sep) + 1).min(q_rec.len);

    let (target_start, target_end, strand) = if frag.strand == 'r' {
        let y_total = reference.total_bases;
        let concat_start = y_total.saturating_sub(frag.y_end + 1);
        let concat_end_incl = y_total.saturating_sub(frag.y_start + 1);
        let local_start = concat_start.saturating_sub(t_rec.start_no_sep);
        let local_end = (concat_end_incl.saturating_sub(t_rec.start_no_sep) + 1).min(t_rec.len);
        (local_start, local_end, '-')
    } else {
        let local_start = frag.y_start.saturating_sub(t_rec.start_no_sep);
        let local_end = (frag.y_end.saturating_sub(t_rec.start_no_sep) + 1).min(t_rec.len);
        (local_start, local_end, '+')
    };

    LocalCoords {
        query_id: &q_rec.id,
        query_len: q_rec.len,
        query_start,
        query_end,
        target_id: &t_rec.id,
        target_len: t_rec.len,
        target_start,
        target_end,
        strand,
    }
}

/// Write alignments in **PAF** (Pairwise mApping Format, minimap2 standard) with GECKO tags:
/// `AS:i:<raw_score>`, `si:f:<similarity>`, `cg:Z:<len>M`.
pub fn write_paf_path(
    path: impl AsRef<Path>,
    query: &Genome,
    reference: &Genome,
    fragments: &[Fragment],
) -> Result<()> {
    let file = File::create(path)?;
    let mut w = BufWriter::with_capacity(1 << 20, file);
    write_paf(&mut w, query, reference, fragments)?;
    w.flush()?;
    Ok(())
}

pub fn write_paf<W: Write>(
    w: &mut W,
    query: &Genome,
    reference: &Genome,
    fragments: &[Fragment],
) -> Result<()> {
    for f in fragments {
        let c = resolve_local_coords(f, query, reference);
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t255\tAS:i:{}\tsi:f:{:.2}\tcg:Z:{}M",
            c.query_id,
            c.query_len,
            c.query_start,
            c.query_end,
            c.strand,
            c.target_id,
            c.target_len,
            c.target_start,
            c.target_end,
            f.ident,
            f.length,
            f.score,
            f.similarity,
            f.length
        )?;
    }
    Ok(())
}

/// Write alignments in **BLAST tabular (`outfmt 6`)** format:
/// `qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore`
/// (1-based inclusive coordinates; `sstart > send` on reverse strand per BLAST convention).
pub fn write_blast_tabular_path(
    path: impl AsRef<Path>,
    query: &Genome,
    reference: &Genome,
    fragments: &[Fragment],
) -> Result<()> {
    let file = File::create(path)?;
    let mut w = BufWriter::with_capacity(1 << 20, file);
    write_blast_tabular(&mut w, query, reference, fragments)?;
    w.flush()?;
    Ok(())
}

pub fn write_blast_tabular<W: Write>(
    w: &mut W,
    query: &Genome,
    reference: &Genome,
    fragments: &[Fragment],
) -> Result<()> {
    // Ungapped Karlin-Altschul constants for +4/-4 nucleotide scoring on uniform ACGT background:
    // sum(p_i p_j exp(lambda * s_ij)) = 0.25*e^(4*lambda) + 0.75*e^(-4*lambda) = 1 => lambda = ln(3)/4.
    let lambda = 3.0_f64.ln() / 4.0;
    let k_param = 0.1_f64;
    let search_space = (query.total_bases as f64) * (reference.total_bases as f64);

    for f in fragments {
        let c = resolve_local_coords(f, query, reference);
        let qstart = c.query_start + 1;
        let qend = c.query_end;
        let (sstart, send) = if c.strand == '-' {
            (c.target_end, c.target_start + 1)
        } else {
            (c.target_start + 1, c.target_end)
        };
        let mismatch = f.length.saturating_sub(f.ident);
        let raw = f.score as f64;
        let bitscore = (lambda * raw - k_param.ln()) / std::f64::consts::LN_2;
        let evalue = (k_param * search_space * (-lambda * raw).exp()).max(0.0);

        writeln!(
            w,
            "{}\t{}\t{:.2}\t{}\t{}\t0\t{}\t{}\t{}\t{}\t{:.2e}\t{:.1}",
            c.query_id,
            c.target_id,
            f.pct_identity(),
            f.length,
            mismatch,
            qstart,
            qend,
            sstart,
            send,
            evalue,
            bitscore.max(0.0)
        )?;
    }
    Ok(())
}

/// Write pairwise synteny blocks in **BEDPE** format:
/// `chrom1 start1 end1 chrom2 start2 end2 name score strand1 strand2`
pub fn write_bedpe_path(
    path: impl AsRef<Path>,
    query: &Genome,
    reference: &Genome,
    fragments: &[Fragment],
) -> Result<()> {
    let file = File::create(path)?;
    let mut w = BufWriter::with_capacity(1 << 20, file);
    write_bedpe(&mut w, query, reference, fragments)?;
    w.flush()?;
    Ok(())
}

pub fn write_bedpe<W: Write>(
    w: &mut W,
    query: &Genome,
    reference: &Genome,
    fragments: &[Fragment],
) -> Result<()> {
    for (idx, f) in fragments.iter().enumerate() {
        let c = resolve_local_coords(f, query, reference);
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\tfrag_{}\t{}\t+\t{}",
            c.query_id,
            c.query_start,
            c.query_end,
            c.target_id,
            c.target_start,
            c.target_end,
            idx,
            f.score,
            c.strand
        )?;
    }
    Ok(())
}

/// Write fragments as newline-delimited JSON (**JSONL**) with resolved sequence IDs & coordinates.
pub fn write_jsonl_path(
    path: impl AsRef<Path>,
    query: &Genome,
    reference: &Genome,
    fragments: &[Fragment],
) -> Result<()> {
    let file = File::create(path)?;
    let mut w = BufWriter::with_capacity(1 << 20, file);
    write_jsonl(&mut w, query, reference, fragments)?;
    w.flush()?;
    Ok(())
}

#[derive(Serialize)]
struct JsonlRecord<'a> {
    query_id: &'a str,
    query_start: u64,
    query_end: u64,
    query_len: u64,
    target_id: &'a str,
    target_start: u64,
    target_end: u64,
    target_len: u64,
    strand: char,
    length: u64,
    ident: u64,
    pct_identity: f64,
    score: u64,
    similarity: f32,
    seq_x: u64,
    seq_y: u64,
}

pub fn write_jsonl<W: Write>(
    w: &mut W,
    query: &Genome,
    reference: &Genome,
    fragments: &[Fragment],
) -> Result<()> {
    for f in fragments {
        let c = resolve_local_coords(f, query, reference);
        let rec = JsonlRecord {
            query_id: c.query_id,
            query_start: c.query_start,
            query_end: c.query_end,
            query_len: c.query_len,
            target_id: c.target_id,
            target_start: c.target_start,
            target_end: c.target_end,
            target_len: c.target_len,
            strand: c.strand,
            length: f.length,
            ident: f.ident,
            pct_identity: f.pct_identity(),
            score: f.score,
            similarity: f.similarity,
            seq_x: f.seq_x,
            seq_y: f.seq_y,
        };
        serde_json::to_writer(&mut *w, &rec)?;
        writeln!(w)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn resolves_multi_fasta_reverse_strand_local_coords() {
        let q = Genome::from_records("q", [("q0", "ACGTACGTAC")]);
        // r0 len=10 (0..10), r1 len=20 (10..30), total=30
        let r = Genome::from_records("r", [("r0", "AAAAAAAAAA"), ("r1", "CCCCCCCCCCCCCCCCCCCC")]);
        // In rev_r: rev_r0 = rc(r1) at [0..20), rev_r1 = rc(r0) at [20..30)
        // A match at rev_r0 [2..=5] (len 4) => orig seq_y=1 (r1), local [14..18)
        let f = Fragment {
            diag: -2,
            x_start: 0,
            y_start: 2,
            x_end: 3,
            y_end: 5,
            length: 4,
            ident: 4,
            score: 16,
            similarity: 100.0,
            seq_x: 0,
            seq_y: 1,
            block: 0,
            strand: 'r',
        };
        let c = resolve_local_coords(&f, &q, &r);
        assert_eq!(c.target_id, "r1");
        assert_eq!(c.target_start, 14);
        assert_eq!(c.target_end, 18);
        assert_eq!(c.strand, '-');

        let mut paf = Vec::new();
        write_paf(&mut paf, &q, &r, &[f]).unwrap();
        let line = String::from_utf8(paf).unwrap();
        assert!(line.starts_with("q0\t10\t0\t4\t-\tr1\t20\t14\t18\t4\t4\t255\tAS:i:16"));
    }
}
