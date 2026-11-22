// query.wgsl
// Tiled WGSL compute shader for canonical k-mer lookup against a Hash-CSR index
// with bit-exact 64-bit MurmurHash3 fmix64 emulation (vec2<u32>) and lock-free
// Left-Maximal Seed Filtering.

struct QueryParams {
    len_x: u32,
    len_y: u32,
    wl: u32,
    shift_bits: u32,   // 64 - bucket_bits (in [40..52], so bucket = hi >> (shift_bits - 32))
    tile_start: u32,
    tile_end: u32,
    max_hits: u32,
    stride: u32,
};

struct KeyPair {
    lo: u32,
    hi: u32,
};

struct HitRec {
    pos_x: u32,
    pos_y: u32,
    strand: u32,
    _pad: u32,
};

@group(0) @binding(0) var<uniform> params: QueryParams;
@group(0) @binding(1) var<storage, read> seq_x: array<u32>;
@group(0) @binding(2) var<storage, read> seq_y: array<u32>;
@group(0) @binding(3) var<storage, read> seq_y_rc: array<u32>;
@group(0) @binding(4) var<storage, read> bucket_off: array<u32>;
@group(0) @binding(5) var<storage, read> idx_keys: array<KeyPair>;
@group(0) @binding(6) var<storage, read> idx_vals: array<u32>;
@group(0) @binding(7) var<storage, read_write> out_hits: array<HitRec>;
@group(0) @binding(8) var<storage, read_write> hit_count: atomic<u32>;

fn read_byte_x(idx: u32) -> u32 {
    return (seq_x[idx >> 2u] >> ((idx & 3u) * 8u)) & 0xFFu;
}

fn read_byte_y(idx: u32) -> u32 {
    return (seq_y[idx >> 2u] >> ((idx & 3u) * 8u)) & 0xFFu;
}

fn read_byte_y_rc(idx: u32) -> u32 {
    return (seq_y_rc[idx >> 2u] >> ((idx & 3u) * 8u)) & 0xFFu;
}

fn base_to_2bit(c: u32) -> u32 {
    if (c == 65u) { return 0u; } // 'A'
    if (c == 67u) { return 1u; } // 'C'
    if (c == 71u) { return 2u; } // 'G'
    if (c == 84u) { return 3u; } // 'T'
    return 4u;
}

// Exact 64-bit wrapping multiplication returning (lo, hi)
fn mul64(a: vec2<u32>, b: vec2<u32>) -> vec2<u32> {
    let a0 = a.x & 0xFFFFu;
    let a1 = a.x >> 16u;
    let b0 = b.x & 0xFFFFu;
    let b1 = b.x >> 16u;

    let m00 = a0 * b0;
    let m01 = a0 * b1;
    let m10 = a1 * b0;
    let m11 = a1 * b1;

    let mid = (m00 >> 16u) + (m01 & 0xFFFFu) + (m10 & 0xFFFFu);
    let p_lo = (m00 & 0xFFFFu) | ((mid & 0xFFFFu) << 16u);
    let p_hi = m11 + (m01 >> 16u) + (m10 >> 16u) + (mid >> 16u);

    let res_hi = p_hi + a.x * b.y + a.y * b.x;
    return vec2<u32>(p_lo, res_hi);
}

// Exact 64-bit MurmurHash3 fmix64 bijection
fn fmix64(x_in: vec2<u32>) -> vec2<u32> {
    var x = vec2<u32>(x_in.x ^ (x_in.y >> 1u), x_in.y);
    x = mul64(x, vec2<u32>(0xed558ccdu, 0xff51afd7u));
    x = vec2<u32>(x.x ^ (x.y >> 1u), x.y);
    x = mul64(x, vec2<u32>(0x1a85ec53u, 0xc4ceb9feu));
    x = vec2<u32>(x.x ^ (x.y >> 1u), x.y);
    return x;
}

fn lt64(a: vec2<u32>, b: vec2<u32>) -> bool {
    return (a.y < b.y) || (a.y == b.y && a.x < b.x);
}

fn le64(a: vec2<u32>, b: vec2<u32>) -> bool {
    return (a.y < b.y) || (a.y == b.y && a.x <= b.x);
}

fn eq64(a: vec2<u32>, b: vec2<u32>) -> bool {
    return a.x == b.x && a.y == b.y;
}

fn shl2_mask(v: vec2<u32>, b: u32, mask_lo: u32, mask_hi: u32) -> vec2<u32> {
    let hi = ((v.y << 2u) | (v.x >> 30u)) & mask_hi;
    let lo = ((v.x << 2u) | b) & mask_lo;
    return vec2<u32>(lo, hi);
}

fn shr2_insert(v: vec2<u32>, cb: u32, rc_shift: u32) -> vec2<u32> {
    var lo = (v.x >> 2u) | ((v.y & 3u) << 30u);
    var hi = v.y >> 2u;
    if (rc_shift < 32u) {
        lo = lo | (cb << rc_shift);
    } else {
        hi = hi | (cb << (rc_shift - 32u));
    }
    return vec2<u32>(lo, hi);
}

@compute @workgroup_size(128)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let t_start = params.tile_start + gid.x * params.stride;
    if (t_start >= params.tile_end) {
        return;
    }
    let t_end = min(t_start + params.stride, params.tile_end);

    let wl = params.wl;
    let bits = 2u * wl;
    var mask_lo: u32 = 0xFFFFFFFFu;
    var mask_hi: u32 = 0xFFFFFFFFu;
    if (bits < 32u) {
        mask_lo = (1u << bits) - 1u;
        mask_hi = 0u;
    } else if (bits == 32u) {
        mask_lo = 0xFFFFFFFFu;
        mask_hi = 0u;
    } else if (bits < 64u) {
        mask_lo = 0xFFFFFFFFu;
        mask_hi = (1u << (bits - 32u)) - 1u;
    }
    let rc_shift = 2u * (wl - 1u);
    let hi_shift = params.shift_bits - 32u;

    var scan_from: u32 = 0u;
    if (t_start >= wl - 1u) {
        scan_from = t_start - (wl - 1u);
    }

    var kf = vec2<u32>(0u, 0u);
    var kr = vec2<u32>(0u, 0u);
    var valid: u32 = 0u;

    for (var i = scan_from; i < t_end; i = i + 1u) {
        let b = base_to_2bit(read_byte_x(i));
        if (b < 4u) {
            kf = shl2_mask(kf, b, mask_lo, mask_hi);
            kr = shr2_insert(kr, b ^ 3u, rc_shift);
            valid = valid + 1u;

            if (valid >= wl && i >= t_start) {
                let pos_x = i + 1u - wl;
                let x_strand = select(1u, 0u, le64(kf, kr));
                let canon = select(kr, kf, x_strand == 0u);
                let is_pal = eq64(kf, kr);

                let h = fmix64(canon);
                let bkt = h.y >> hi_shift;
                let st = bucket_off[bkt];
                let en = bucket_off[bkt + 1u];
                if (st == en) {
                    continue;
                }

                var j = st;
                if (en - st <= 16u) {
                    loop {
                        if (j >= en) { break; }
                        let kp = idx_keys[j];
                        if (!lt64(vec2<u32>(kp.lo, kp.hi), h)) { break; }
                        j = j + 1u;
                    }
                } else {
                    var lo = st;
                    var hi = en;
                    loop {
                        if (lo >= hi) { break; }
                        let mid = lo + ((hi - lo) >> 1u);
                        let kp = idx_keys[mid];
                        if (lt64(vec2<u32>(kp.lo, kp.hi), h)) {
                            lo = mid + 1u;
                        } else {
                            hi = mid;
                        }
                    }
                    j = lo;
                }

                var prev_x: u32 = 0u;
                var prev_x_valid = false;
                if (pos_x > 0u) {
                    prev_x = read_byte_x(pos_x - 1u);
                    prev_x_valid = (base_to_2bit(prev_x) < 4u);
                }
                var next_x: u32 = 0u;
                var next_x_valid = false;
                if (pos_x + wl < params.len_x) {
                    next_x = read_byte_x(pos_x + wl);
                    next_x_valid = (base_to_2bit(next_x) < 4u);
                }

                loop {
                    if (j >= en) { break; }
                    let kp = idx_keys[j];
                    if (!eq64(vec2<u32>(kp.lo, kp.hi), h)) { break; }

                    let ps = idx_vals[j];
                    let pos_y_orig = ps >> 1u;
                    let y_strand = ps & 1u;

                    if (is_pal || x_strand == y_strand) {
                        var left_ext = false;
                        if (prev_x_valid && pos_y_orig > 0u) {
                            left_ext = (read_byte_y(pos_y_orig - 1u) == prev_x);
                        }
                        var right_ext = false;
                        if (next_x_valid && pos_y_orig + wl < params.len_y) {
                            right_ext = (read_byte_y(pos_y_orig + wl) == next_x);
                        }
                        if (!left_ext || !right_ext) {
                            let slot = atomicAdd(&hit_count, 1u);
                            if (slot < params.max_hits) {
                                out_hits[slot] = HitRec(pos_x, pos_y_orig, 0u, 0u);
                            }
                        }
                    }
                    if (is_pal || x_strand != y_strand) {
                        let pos_y_rc = params.len_y - pos_y_orig - wl;
                        var left_ext = false;
                        if (prev_x_valid && pos_y_rc > 0u) {
                            left_ext = (read_byte_y_rc(pos_y_rc - 1u) == prev_x);
                        }
                        var right_ext = false;
                        if (next_x_valid && pos_y_rc + wl < params.len_y) {
                            right_ext = (read_byte_y_rc(pos_y_rc + wl) == next_x);
                        }
                        if (!left_ext || !right_ext) {
                            let slot = atomicAdd(&hit_count, 1u);
                            if (slot < params.max_hits) {
                                out_hits[slot] = HitRec(pos_x, pos_y_rc, 1u, 0u);
                            }
                        }
                    }
                    j = j + 1u;
                }
            }
        } else {
            valid = 0u;
            kf = vec2<u32>(0u, 0u);
            kr = vec2<u32>(0u, 0u);
        }
    }
}
