// extend.wgsl
// 2-bit PackedSeq WGSL compute shader for exact GECKO bidirectional ungapped extension
// supporting configurable scoring, 'N' penalties, multi-FASTA SEQ_SEP, rare IUPAC codes,
// palindromic spaced seeds, and 16-bp word-at-a-time exact-match fast paths.

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
    sp_x_off: u32,       // word offset into special_bm, or 0xFFFFFFFF if none
    sp_y_off: u32,
    sp_yrc_off: u32,
    sep_x_off: u32,      // u32 offset into meta_buf
    sep_x_len: u32,
    sep_y_off: u32,
    sep_y_len: u32,
    sep_yrc_off: u32,
    sep_yrc_len: u32,
    iupac_x_off: u32,    // u32 offset into meta_buf (pairs of [pos, upper_ascii])
    iupac_x_len: u32,
    iupac_y_off: u32,
    iupac_y_len: u32,
    iupac_yrc_off: u32,
    iupac_yrc_len: u32,
    is_spaced: u32,
    _pad1: u32,
    _pad2: u32,
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
@group(0) @binding(1) var<storage, read> seq_x_bits: array<u32>;
@group(0) @binding(2) var<storage, read> seq_y_bits: array<u32>;
@group(0) @binding(3) var<storage, read> seq_y_rc_bits: array<u32>;
@group(0) @binding(4) var<storage, read> special_bm: array<u32>;
@group(0) @binding(5) var<storage, read> meta_buf: array<u32>;
@group(0) @binding(6) var<storage, read> hits: array<HitRec>;
@group(0) @binding(7) var<storage, read_write> out_frags: array<FragOut>;
@group(0) @binding(8) var<storage, read_write> frag_count: atomic<u32>;

const NONE_OFF: u32 = 0xFFFFFFFFu;

fn base_2bit_x(pos: u32) -> u32 {
    return (seq_x_bits[pos >> 4u] >> ((pos & 15u) * 2u)) & 3u;
}

fn base_2bit_y(pos: u32, is_rev: u32) -> u32 {
    var word: u32 = 0u;
    if (is_rev == 0u) {
        word = seq_y_bits[pos >> 4u];
    } else {
        word = seq_y_rc_bits[pos >> 4u];
    }
    return (word >> ((pos & 15u) * 2u)) & 3u;
}

fn load_16_2bit_x(pos: u32) -> u32 {
    let w = pos >> 4u;
    let sh = (pos & 15u) * 2u;
    let lo = seq_x_bits[w];
    if (sh == 0u) {
        return lo;
    }
    return (lo >> sh) | (seq_x_bits[w + 1u] << (32u - sh));
}

fn load_16_2bit_y(pos: u32, is_rev: u32) -> u32 {
    let w = pos >> 4u;
    let sh = (pos & 15u) * 2u;
    if (is_rev == 0u) {
        let lo = seq_y_bits[w];
        if (sh == 0u) {
            return lo;
        }
        return (lo >> sh) | (seq_y_bits[w + 1u] << (32u - sh));
    } else {
        let lo = seq_y_rc_bits[w];
        if (sh == 0u) {
            return lo;
        }
        return (lo >> sh) | (seq_y_rc_bits[w + 1u] << (32u - sh));
    }
}

fn is_special(bm_off: u32, pos: u32) -> bool {
    if (bm_off == NONE_OFF) {
        return false;
    }
    return ((special_bm[bm_off + (pos >> 5u)] >> (pos & 31u)) & 1u) != 0u;
}

fn load_16_special(bm_off: u32, pos: u32) -> u32 {
    if (bm_off == NONE_OFF) {
        return 0u;
    }
    let w = bm_off + (pos >> 5u);
    let sh = pos & 31u;
    let lo = special_bm[w];
    if (sh <= 16u) {
        return (lo >> sh) & 0xFFFFu;
    }
    let hi = special_bm[w + 1u];
    return ((lo >> sh) | (hi << (32u - sh))) & 0xFFFFu;
}

fn base_2bit_to_ascii(b: u32) -> u32 {
    if (b == 0u) { return 65u; } // 'A'
    if (b == 1u) { return 67u; } // 'C'
    if (b == 2u) { return 71u; } // 'G'
    return 84u;                  // 'T'
}

fn decode_special(pos: u32, sep_off: u32, sep_len: u32, iupac_off: u32, iupac_len: u32) -> u32 {
    if (sep_len > 0u) {
        var lo: u32 = 0u;
        var hi: u32 = sep_len;
        loop {
            if (lo >= hi) { break; }
            let mid = lo + ((hi - lo) >> 1u);
            if (meta_buf[sep_off + mid] < pos) {
                lo = mid + 1u;
            } else {
                hi = mid;
            }
        }
        if (lo < sep_len && meta_buf[sep_off + lo] == pos) {
            return 0u; // SEQ_SEP
        }
    }
    if (iupac_len > 0u) {
        var lo: u32 = 0u;
        var hi: u32 = iupac_len;
        loop {
            if (lo >= hi) { break; }
            let mid = lo + ((hi - lo) >> 1u);
            if (meta_buf[iupac_off + mid * 2u] < pos) {
                lo = mid + 1u;
            } else {
                hi = mid;
            }
        }
        if (lo < iupac_len && meta_buf[iupac_off + lo * 2u] == pos) {
            return meta_buf[iupac_off + lo * 2u + 1u];
        }
    }
    return 78u; // 'N'
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

    let sp_x_off = params.sp_x_off;
    let sp_y_off = select(params.sp_y_off, params.sp_yrc_off, is_rev != 0u);
    let sep_y_off = select(params.sep_y_off, params.sep_yrc_off, is_rev != 0u);
    let sep_y_len = select(params.sep_y_len, params.sep_yrc_len, is_rev != 0u);
    let iupac_y_off = select(params.iupac_y_off, params.iupac_yrc_off, is_rev != 0u);
    let iupac_y_len = select(params.iupac_y_len, params.iupac_yrc_len, is_rev != 0u);

    // Initial seed identity and score (exact for both contiguous and spaced seeds)
    var n_ident: u32 = wl;
    var fscore: i32 = match_sc * i32(wl);
    if (params.is_spaced != 0u) {
        let bits = 2u * wl;
        var m_lo: u32 = 0xFFFFFFFFu;
        var m_hi: u32 = 0xFFFFFFFFu;
        if (bits < 32u) {
            m_lo = (1u << bits) - 1u;
            m_hi = 0u;
        } else if (bits == 32u) {
            m_lo = 0xFFFFFFFFu;
            m_hi = 0u;
        } else if (bits < 64u) {
            m_lo = 0xFFFFFFFFu;
            m_hi = (1u << (bits - 32u)) - 1u;
        }
        let d_lo = (load_16_2bit_x(pos_x) ^ load_16_2bit_y(pos_y, is_rev)) & m_lo;
        let d_hi = (load_16_2bit_x(pos_x + 16u) ^ load_16_2bit_y(pos_y + 16u, is_rev)) & m_hi;
        let mm_lo = (d_lo | (d_lo >> 1u)) & 0x55555555u;
        let mm_hi = (d_hi | (d_hi >> 1u)) & 0x55555555u;
        let mm = countOneBits(mm_lo) + countOneBits(mm_hi);
        n_ident = wl - mm;
        fscore = i32(n_ident) * match_sc - i32(mm) * mm_pen;
    }

    // Forward pass with 16-bp exact-match word fast path
    let xfil = pos_x + wl;
    let ycol = pos_y + wl;
    let rem_fwd = min(params.len_x - xfil, params.len_y - ycol);

    var xfilmax = pos_x + wl - 1u;
    var ycolmax = pos_y + wl - 1u;
    var fscore_max: i32 = fscore;
    var max_ident: u32 = n_ident;

    var step: u32 = 0u;
    loop {
        if (step >= rem_fwd) { break; }
        let px = xfil + step;
        let py = ycol + step;

        if (step + 16u <= rem_fwd) {
            let sp16 = load_16_special(sp_x_off, px) | load_16_special(sp_y_off, py);
            if (sp16 == 0u) {
                let wx = load_16_2bit_x(px);
                let wy = load_16_2bit_y(py, is_rev);
                if (wx == wy) {
                    fscore = fscore + 16 * match_sc;
                    n_ident = n_ident + 16u;
                    if (fscore_max < fscore) {
                        fscore_max = fscore;
                        xfilmax = px + 15u;
                        ycolmax = py + 15u;
                        max_ident = n_ident;
                    }
                    step = step + 16u;
                    continue;
                }
                var diff = wx ^ wy;
                var terminated = false;
                for (var k: u32 = 0u; k < 16u; k = k + 1u) {
                    if ((diff & 3u) == 0u) {
                        fscore = fscore + match_sc;
                        n_ident = n_ident + 1u;
                        if (fscore_max < fscore) {
                            fscore_max = fscore;
                            xfilmax = px + k;
                            ycolmax = py + k;
                            max_ident = n_ident;
                        }
                    } else {
                        fscore = fscore - mm_pen;
                        if (fscore < 0) {
                            terminated = true;
                            break;
                        }
                    }
                    diff = diff >> 2u;
                }
                if (terminated) { break; }
                step = step + 16u;
                continue;
            }
        }

        let sp_x = is_special(sp_x_off, px);
        let sp_y = is_special(sp_y_off, py);
        if (!sp_x && !sp_y) {
            if (base_2bit_x(px) == base_2bit_y(py, is_rev)) {
                fscore = fscore + match_sc;
                n_ident = n_ident + 1u;
                if (fscore_max < fscore) {
                    fscore_max = fscore;
                    xfilmax = px;
                    ycolmax = py;
                    max_ident = n_ident;
                }
            } else {
                fscore = fscore - mm_pen;
            }
        } else {
            let cx = select(
                base_2bit_to_ascii(base_2bit_x(px)),
                decode_special(px, params.sep_x_off, params.sep_x_len, params.iupac_x_off, params.iupac_x_len),
                sp_x
            );
            let cy = select(
                base_2bit_to_ascii(base_2bit_y(py, is_rev)),
                decode_special(py, sep_y_off, sep_y_len, iupac_y_off, iupac_y_len),
                sp_y
            );
            if (cx == 0u || cy == 0u) {
                break;
            }
            if (cx == 78u || cy == 78u) {
                fscore = fscore - n_pen;
            } else if (cx == cy) {
                fscore = fscore + match_sc;
                n_ident = n_ident + 1u;
                if (fscore_max < fscore) {
                    fscore_max = fscore;
                    xfilmax = px;
                    ycolmax = py;
                    max_ident = n_ident;
                }
            } else {
                fscore = fscore - mm_pen;
            }
        }
        if (fscore < 0) {
            break;
        }
        step = step + 1u;
    }

    // Backward pass: starts with fscore = fscore_max (exact GECKO FragHits.c semantics)
    let rem_bwd = min(pos_x, pos_y);
    var xfilmax2 = pos_x;
    var ycolmax2 = pos_y;
    fscore = fscore_max;
    n_ident = max_ident;

    var bstep: u32 = 1u;
    loop {
        if (bstep > rem_bwd) { break; }

        if (bstep + 15u <= rem_bwd) {
            let px_lo = pos_x - bstep - 15u;
            let py_lo = pos_y - bstep - 15u;
            let sp16 = load_16_special(sp_x_off, px_lo) | load_16_special(sp_y_off, py_lo);
            if (sp16 == 0u) {
                let wx = load_16_2bit_x(px_lo);
                let wy = load_16_2bit_y(py_lo, is_rev);
                if (wx == wy) {
                    fscore = fscore + 16 * match_sc;
                    n_ident = n_ident + 16u;
                    if (fscore_max < fscore) {
                        fscore_max = fscore;
                        xfilmax2 = px_lo;
                        ycolmax2 = py_lo;
                        max_ident = n_ident;
                    }
                    bstep = bstep + 16u;
                    continue;
                }
                var diff = wx ^ wy;
                var terminated = false;
                for (var k: u32 = 0u; k < 16u; k = k + 1u) {
                    let lane = 15u - k;
                    if (((diff >> (lane * 2u)) & 3u) == 0u) {
                        fscore = fscore + match_sc;
                        n_ident = n_ident + 1u;
                        if (fscore_max < fscore) {
                            fscore_max = fscore;
                            xfilmax2 = px_lo + lane;
                            ycolmax2 = py_lo + lane;
                            max_ident = n_ident;
                        }
                    } else {
                        fscore = fscore - mm_pen;
                        if (fscore < 0) {
                            terminated = true;
                            break;
                        }
                    }
                }
                if (terminated) { break; }
                bstep = bstep + 16u;
                continue;
            }
        }

        let px = pos_x - bstep;
        let py = pos_y - bstep;
        let sp_x = is_special(sp_x_off, px);
        let sp_y = is_special(sp_y_off, py);
        if (!sp_x && !sp_y) {
            if (base_2bit_x(px) == base_2bit_y(py, is_rev)) {
                fscore = fscore + match_sc;
                n_ident = n_ident + 1u;
                if (fscore_max < fscore) {
                    fscore_max = fscore;
                    xfilmax2 = px;
                    ycolmax2 = py;
                    max_ident = n_ident;
                }
            } else {
                fscore = fscore - mm_pen;
            }
        } else {
            let cx = select(
                base_2bit_to_ascii(base_2bit_x(px)),
                decode_special(px, params.sep_x_off, params.sep_x_len, params.iupac_x_off, params.iupac_x_len),
                sp_x
            );
            let cy = select(
                base_2bit_to_ascii(base_2bit_y(py, is_rev)),
                decode_special(py, sep_y_off, sep_y_len, iupac_y_off, iupac_y_len),
                sp_y
            );
            if (cx == 0u || cy == 0u) {
                break;
            }
            if (cx == 78u || cy == 78u) {
                fscore = fscore - n_pen;
            } else if (cx == cy) {
                fscore = fscore + match_sc;
                n_ident = n_ident + 1u;
                if (fscore_max < fscore) {
                    fscore_max = fscore;
                    xfilmax2 = px;
                    ycolmax2 = py;
                    max_ident = n_ident;
                }
            } else {
                fscore = fscore - mm_pen;
            }
        }
        if (fscore < 0) {
            break;
        }
        bstep = bstep + 1u;
    }

    let length = xfilmax - xfilmax2 + 1u;
    // Strict `>` matches GECKO CPU/CUDA filter semantics.
    let lhs = u32(max(fscore_max, 0)) * 10000u;
    let rhs = length * u32(match_sc) * params.sim_th_x100;
    if (length > params.min_len && lhs > rhs) {
        let slot = atomicAdd(&frag_count, 1u);
        if (slot < params.max_frags) {
            out_frags[slot] = FragOut(
                xfilmax2,
                ycolmax2,
                xfilmax,
                ycolmax,
                max_ident,
                u32(fscore_max),
                is_rev,
                0u
            );
        }
    }
}
