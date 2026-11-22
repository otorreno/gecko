// extend.wgsl
// WGSL compute shader for exact GECKO bidirectional ungapped extension
// (+4 match, -4 mismatch, -1 for N, stop on 0x00 separator or fscore < 0).

struct ExtendParams {
    len_x: u32,
    len_y: u32,
    wl: u32,
    min_len: u32,
    sim_th_x100: u32,   // SimTh * 100 (e.g. 60 -> 6000) for integer comparison
    n_hits: u32,
    max_frags: u32,
    _pad: u32,
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

fn read_upper_x(idx: u32) -> u32 {
    return ((seq_x[idx >> 2u] >> ((idx & 3u) * 8u)) & 0xFFu) & 0xDFu;
}

fn read_upper_y(idx: u32, is_rev: u32) -> u32 {
    var word: u32 = 0u;
    if (is_rev == 0u) {
        word = seq_y[idx >> 2u];
    } else {
        word = seq_y_rc[idx >> 2u];
    }
    return ((word >> ((idx & 3u) * 8u)) & 0xFFu) & 0xDFu;
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

    // Forward pass
    var xfil = pos_x + wl;
    var ycol = pos_y + wl;
    let rem_fwd = min(params.len_x - xfil, params.len_y - ycol);

    var xfilmax = pos_x + wl - 1u;
    var ycolmax = pos_y + wl - 1u;
    var fscore: i32 = i32(4u * wl);
    var fscore_max: i32 = fscore;
    var n_ident: u32 = wl;
    var max_ident: u32 = wl;

    for (var step: u32 = 0u; step < rem_fwd; step = step + 1u) {
        let cx = read_upper_x(xfil + step);
        let cy = read_upper_y(ycol + step, is_rev);
        if (cx == 0u || cy == 0u) {
            break;
        }
        if (cx == 78u || cy == 78u) { // 'N'
            fscore = fscore - 1;
        } else if (cx == cy) {
            fscore = fscore + 4;
            n_ident = n_ident + 1u;
            if (fscore_max < fscore) {
                fscore_max = fscore;
                xfilmax = xfil + step;
                ycolmax = ycol + step;
                max_ident = n_ident;
            }
        } else {
            fscore = fscore - 4;
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
        let cx = read_upper_x(pos_x - step);
        let cy = read_upper_y(pos_y - step, is_rev);
        if (cx == 0u || cy == 0u) {
            break;
        }
        if (cx == 78u || cy == 78u) {
            fscore = fscore - 1;
        } else if (cx == cy) {
            fscore = fscore + 4;
            n_ident = n_ident + 1u;
            if (fscore_max < fscore) {
                fscore_max = fscore;
                xfilmax2 = pos_x - step;
                ycolmax2 = pos_y - step;
                max_ident = n_ident;
            }
        } else {
            fscore = fscore - 4;
        }
        if (fscore < 0) {
            break;
        }
    }

    let length = xfilmax - xfilmax2 + 1u;
    if (length > params.min_len) {
        let sim = (f32(fscore_max) * 100.0) / f32(length * 4u);
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
