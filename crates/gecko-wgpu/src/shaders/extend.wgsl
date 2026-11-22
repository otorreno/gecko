// extend.wgsl
// WGSL compute shader for exact GECKO bidirectional ungapped extension
// supporting configurable scoring and '*' (42) / 0x00 multi-FASTA separators.

struct ExtendParams {
    len_x: u32,
    len_y: u32,
    wl: u32,
    min_len: u32,
    sim_th_x100: u32,
    n_hits: u32,
    max_frags: u32,
    match_score: i32,
    mismatch_penalty: i32,
    n_penalty: i32,
    _pad0: u32,
    _pad1: u32,
};

struct HitRec {
    pos_x: u32,
    pos_y: u32,
    strand: u32,
    _pad: u32,
};

struct FragOut {
    x_start: u32,
    y_start: u32,
    x_end: u32,
    y_end: u32,
    ident: u32,
    score: u32,
    strand: u32,
    _pad: u32,
};

@group(0) @binding(0) var<uniform> params: ExtendParams;
@group(0) @binding(1) var<storage, read> seq_x: array<u32>;
@group(0) @binding(2) var<storage, read> seq_y: array<u32>;
@group(0) @binding(3) var<storage, read> seq_y_rc: array<u32>;
@group(0) @binding(4) var<storage, read> hits: array<HitRec>;
@group(0) @binding(5) var<storage, read_write> out_frags: array<FragOut>;
@group(0) @binding(6) var<storage, read_write> frag_count: atomic<u32>;

fn read_raw_x(idx: u32) -> u32 {
    return (seq_x[idx >> 2u] >> ((idx & 3u) * 8u)) & 0xFFu;
}

fn read_raw_y(idx: u32, is_rev: u32) -> u32 {
    var word: u32 = 0u;
    if (is_rev == 0u) {
        word = seq_y[idx >> 2u];
    } else {
        word = seq_y_rc[idx >> 2u];
    }
    return (word >> ((idx & 3u) * 8u)) & 0xFFu;
}

@compute @workgroup_size(128)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let idx = gid.x;
    if (idx >= params.n_hits) {
        return;
    }

    let h = hits[idx];
    let pos_x = h.pos_x;
    let pos_y = h.pos_y;
    let is_rev = h.strand;
    let wl = params.wl;
    let match_sc = params.match_score;
    let mm_pen = params.mismatch_penalty;
    let n_pen = params.n_penalty;

    // Forward pass
    var xfil = pos_x + wl;
    var ycol = pos_y + wl;
    let rem_fwd = min(params.len_x - xfil, params.len_y - ycol);

    var xfilmax = pos_x + wl - 1u;
    var ycolmax = pos_y + wl - 1u;
    var fscore: i32 = match_sc * i32(wl);
    var fscore_max: i32 = fscore;
    var n_ident: u32 = wl;
    var max_ident: u32 = wl;

    for (var step: u32 = 0u; step < rem_fwd; step = step + 1u) {
        let rx = read_raw_x(xfil + step);
        let ry = read_raw_y(ycol + step, is_rev);
        if (rx == 0u || ry == 0u || rx == 42u || ry == 42u) {
            break;
        }
        let cx = rx & 0xDFu;
        let cy = ry & 0xDFu;
        if (cx == 78u || cy == 78u) { // 'N'
            fscore = fscore - n_pen;
        } else if (cx == cy) {
            fscore = fscore + match_sc;
            n_ident = n_ident + 1u;
            if (fscore_max < fscore) {
                fscore_max = fscore;
                xfilmax = xfil + step;
                ycolmax = ycol + step;
                max_ident = n_ident;
            }
        } else {
            fscore = fscore - mm_pen;
        }
        if (fscore < 0) {
            break;
        }
    }

    // Backward pass: starts with fscore = fscore_max (exact GECKO FragHits.c semantics)
    let rem_bwd = min(pos_x, pos_y);
    var xfilmax2 = pos_x;
    var ycolmax2 = pos_y;
    fscore = fscore_max;
    n_ident = max_ident;

    for (var step: u32 = 1u; step <= rem_bwd; step = step + 1u) {
        let rx = read_raw_x(pos_x - step);
        let ry = read_raw_y(pos_y - step, is_rev);
        if (rx == 0u || ry == 0u || rx == 42u || ry == 42u) {
            break;
        }
        let cx = rx & 0xDFu;
        let cy = ry & 0xDFu;
        if (cx == 78u || cy == 78u) {
            fscore = fscore - n_pen;
        } else if (cx == cy) {
            fscore = fscore + match_sc;
            n_ident = n_ident + 1u;
            if (fscore_max < fscore) {
                fscore_max = fscore;
                xfilmax2 = pos_x - step;
                ycolmax2 = pos_y - step;
                max_ident = n_ident;
            }
        } else {
            fscore = fscore - mm_pen;
        }
        if (fscore < 0) {
            break;
        }
    }

    let length = xfilmax - xfilmax2 + 1u;
    if (length > params.min_len) {
        let sim = (f32(fscore_max) * 100.0) / f32(length * u32(match_sc));
        let sim_th = f32(params.sim_th_x100) / 100.0;
        if (sim > sim_th) {
            let slot = atomicAdd(&frag_count, 1u);
            if (slot < params.max_frags) {
                out_frags[slot] = FragOut(
                    xfilmax2, ycolmax2, xfilmax, ycolmax,
                    max_ident, u32(fscore_max), is_rev, 0u
                );
            }
        }
    }
}
