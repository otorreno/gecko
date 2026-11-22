//! Fragment deduplication and GECKO-compatible output writers
//! (.frags exact 93-byte big-endian format, .frags.INF, .frags.MAT, .csv).

use std::fs::File;
use std::io::{BufWriter, Write};

use crate::gpu::FragOut;
use crate::seq::FastSeq;

#[derive(Clone, Debug)]
pub struct GeckoFrag {
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

pub fn finalize_and_dedup(
    raw: &[FragOut],
    s_x: &FastSeq,
    s_y: &FastSeq,
    s_y_rc: &FastSeq,
) -> (Vec<GeckoFrag>, Vec<GeckoFrag>, Vec<[i64; 100]>) {
    let mut frags_f = Vec::new();
    let mut frags_r = Vec::new();
    let mut coverage = vec![[0i64; 100]; 1000];

    for gf in raw {
        let is_rev = gf.strand != 0;
        let seq_x = s_x.find_seq_id(gf.x_start as u64);
        let seq_y_buf = if is_rev {
            s_y_rc.find_seq_id(gf.y_start as u64)
        } else {
            s_y.find_seq_id(gf.y_start as u64)
        };
        let length = (gf.x_end - gf.x_start + 1) as u64;
        let sim = ((gf.score as f64) * 100.0 / ((length * 4) as f64)) as f32;

        let frag = GeckoFrag {
            diag: (gf.x_start as i64) - (gf.y_start as i64),
            x_start: (gf.x_start as u64) - seq_x,
            y_start: (gf.y_start as u64) - seq_y_buf,
            x_end: (gf.x_end as u64) - seq_x,
            y_end: (gf.y_end as u64) - seq_y_buf,
            length,
            ident: gf.ident as u64,
            score: gf.score as u64,
            similarity: sim,
            seq_x,
            seq_y: if is_rev {
                s_y.n_seqs - seq_y_buf - 1
            } else {
                seq_y_buf
            },
            block: 0,
            strand: if is_rev { 'r' } else { 'f' },
        };

        let lb = (length.min(999)) as usize;
        let sb = (sim as i32).clamp(0, 99) as usize;
        coverage[lb][sb] += 1;

        if is_rev {
            frags_r.push(frag);
        } else {
            frags_f.push(frag);
        }
    }

    dedup_strand(&mut frags_f);
    dedup_strand(&mut frags_r);
    (frags_f, frags_r, coverage)
}

fn dedup_strand(v: &mut Vec<GeckoFrag>) {
    if v.len() <= 1 {
        return;
    }
    v.sort_unstable_by(|a, b| {
        a.diag
            .cmp(&b.diag)
            .then_with(|| a.x_start.cmp(&b.x_start))
            .then_with(|| b.x_end.cmp(&a.x_end))
            .then_with(|| b.score.cmp(&a.score))
    });
    let mut out = 0usize;
    for i in 0..v.len() {
        if out > 0
            && v[i].diag == v[out - 1].diag
            && v[i].seq_x == v[out - 1].seq_x
            && v[i].seq_y == v[out - 1].seq_y
            && v[i].x_end <= v[out - 1].x_end
        {
            continue;
        }
        v.swap(out, i);
        out += 1;
    }
    v.truncate(out);
}

fn write_frag_be<W: Write>(w: &mut W, f: &GeckoFrag) -> std::io::Result<()> {
    w.write_all(&f.diag.to_be_bytes())?;
    w.write_all(&f.x_start.to_be_bytes())?;
    w.write_all(&f.y_start.to_be_bytes())?;
    w.write_all(&f.x_end.to_be_bytes())?;
    w.write_all(&f.y_end.to_be_bytes())?;
    w.write_all(&f.length.to_be_bytes())?;
    w.write_all(&f.ident.to_be_bytes())?;
    w.write_all(&f.score.to_be_bytes())?;
    w.write_all(&f.similarity.to_be_bytes())?;
    w.write_all(&f.seq_x.to_be_bytes())?;
    w.write_all(&f.seq_y.to_be_bytes())?;
    w.write_all(&f.block.to_be_bytes())?;
    w.write_all(&[f.strand as u8])?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub fn write_outputs(
    out_prefix: &str,
    seq_x_path: &str,
    seq_y_path: &str,
    s_x: &FastSeq,
    s_y: &FastSeq,
    l_min: u64,
    sim_th: u64,
    tot_hits: u64,
    tot_hits_used: u64,
    frags_f: &[GeckoFrag],
    frags_r: &[GeckoFrag],
    coverage: &[[i64; 100]],
) -> Result<(), String> {
    let n0 = s_x.n_bases;
    let n1 = s_y.n_bases;
    let total_frags = frags_f.len() + frags_r.len();

    // 1. Binary .frags (exact 93-byte big-endian records)
    {
        let f = File::create(format!("{out_prefix}.frags"))
            .map_err(|e| format!("creating .frags: {e}"))?;
        let mut w = BufWriter::new(f);
        w.write_all(&n0.to_be_bytes()).map_err(|e| e.to_string())?;
        w.write_all(&n1.to_be_bytes()).map_err(|e| e.to_string())?;
        for fr in frags_f {
            write_frag_be(&mut w, fr).map_err(|e| e.to_string())?;
        }
        for fr in frags_r {
            write_frag_be(&mut w, fr).map_err(|e| e.to_string())?;
        }
    }

    // 2. .frags.MAT
    {
        let f = File::create(format!("{out_prefix}.frags.MAT"))
            .map_err(|e| format!("creating .frags.MAT: {e}"))?;
        let mut w = BufWriter::new(f);
        for row in coverage {
            for &val in row {
                write!(w, "{val}\t").map_err(|e| e.to_string())?;
            }
            writeln!(w).map_err(|e| e.to_string())?;
        }
    }

    // 3. .frags.INF
    let header = format!(
        "All by-Identity Ungapped Fragments (gecko_wgpu GPU engine)\n\
         [GECKO-compatible output]\n\
         SeqX filename        : {seq_x_path}\n\
         SeqY filename        : {seq_y_path}\n\
         SeqX name            : {}\n\
         SeqY name            : {}\n\
         SeqX length          : {n0}\n\
         SeqY length          : {n1}\n\
         Min.fragment.length  : {l_min}\n\
         Min.Identity         : {sim_th}\n\
         Tot Hits (seeds)     : {tot_hits}\n\
         Tot Hits (seeds) used: {tot_hits_used}\n\
         Total fragments      : {total_frags}\n\
         ========================================================\n",
        s_x.ident, s_y.ident
    );
    std::fs::write(format!("{out_prefix}.frags.INF"), &header)
        .map_err(|e| format!("writing .frags.INF: {e}"))?;

    // 4. .csv
    {
        let f = File::create(format!("{out_prefix}.csv"))
            .map_err(|e| format!("creating .csv: {e}"))?;
        let mut w = BufWriter::new(f);
        w.write_all(header.as_bytes()).map_err(|e| e.to_string())?;
        writeln!(w, "Total CSB: 0").map_err(|e| e.to_string())?;
        writeln!(w, "========================================================")
            .map_err(|e| e.to_string())?;
        writeln!(
            w,
            "Type,xStart,yStart,xEnd,yEnd,strand(f/r),block,length,score,ident,similarity,%ident,SeqX,SeqY"
        )
        .map_err(|e| e.to_string())?;

        for fr in frags_f {
            let pct_ident = 100.0 * (fr.ident as f64) / (fr.length as f64);
            writeln!(
                w,
                "Frag,{},{},{},{},{},{},{},{},{},{:.2},{:.2},{},{}",
                fr.x_start,
                fr.y_start,
                fr.x_end,
                fr.y_end,
                fr.strand,
                fr.block,
                fr.length,
                fr.score,
                fr.ident,
                fr.similarity,
                pct_ident,
                fr.seq_x,
                fr.seq_y
            )
            .map_err(|e| e.to_string())?;
        }
        for fr in frags_r {
            let pct_ident = 100.0 * (fr.ident as f64) / (fr.length as f64);
            let ys = n1 - fr.y_start - 1;
            let ye = n1 - fr.y_end - 1;
            writeln!(
                w,
                "Frag,{},{},{},{},{},{},{},{},{},{:.2},{:.2},{},{}",
                fr.x_start,
                ys,
                fr.x_end,
                ye,
                fr.strand,
                fr.block,
                fr.length,
                fr.score,
                fr.ident,
                fr.similarity,
                pct_ident,
                fr.seq_x,
                fr.seq_y
            )
            .map_err(|e| e.to_string())?;
        }
    }

    Ok(())
}
