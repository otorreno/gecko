//! Parallel canonical k-mer Hash-CSR reference index builder for WGPU.
//! Stores bijective MurmurHash3 fmix64(canonical_kmer) keys split into `[lo, hi]` u32 pairs
//! so WGSL compute shaders can probe exact 64-bit keys with zero hash collisions.

use bytemuck::{Pod, Zeroable};
use gecko_core::ref_index::{base_lut, fmix64};
use gecko_core::seq::Genome;
use rayon::prelude::*;

#[repr(C)]
#[derive(Copy, Clone, Debug, Pod, Zeroable)]
pub struct KeyU64Pair {
    pub lo: u32,
    pub hi: u32,
}

pub struct WgpuRefIndex {
    pub shift_bits: u32,
    pub bucket_off: Vec<u32>,
    pub keys: Vec<KeyU64Pair>,
    pub vals: Vec<u32>,
}

impl WgpuRefIndex {
    pub fn build(reference: &Genome, wl: usize, freq_limit: u32, soft_mask: bool) -> Self {
        let lut = base_lut(soft_mask);
        let data = &reference.data;
        let len = data.len();
        let mask: u64 = if wl == 32 {
            u64::MAX
        } else {
            (1u64 << (2 * wl)) - 1
        };
        let rc_shift = 2 * (wl - 1);

        let n_threads = rayon::current_num_threads().max(1);
        let chunk_sz = len.div_ceil(n_threads).max(65_536);
        let n_chunks = if len == 0 { 0 } else { len.div_ceil(chunk_sz) };

        let per_chunk_parts: Vec<Vec<Vec<(u64, u32)>>> = (0..n_chunks)
            .into_par_iter()
            .map(|cid| {
                let c_start = cid * chunk_sz;
                let c_end = (c_start + chunk_sz).min(len);
                let scan_from = c_start.saturating_sub(wl - 1);
                let mut parts: Vec<Vec<(u64, u32)>> = (0..256).map(|_| Vec::new()).collect();

                let mut kf: u64 = 0;
                let mut kr: u64 = 0;
                let mut valid: usize = 0;

                for (off, &byte) in data[scan_from..c_end].iter().enumerate() {
                    let i = scan_from + off;
                    let b = lut[byte as usize];
                    if b < 4 {
                        kf = ((kf << 2) | b as u64) & mask;
                        kr = (kr >> 2) | (((b ^ 3) as u64) << rc_shift);
                        valid += 1;
                        if valid >= wl && i >= c_start {
                            let strand: u32 = if kf <= kr { 0 } else { 1 };
                            let key = if strand == 1 { kr } else { kf };
                            let h = fmix64(key);
                            let p = (h >> 56) as usize;
                            let start_pos = (i + 1 - wl) as u32;
                            parts[p].push((h, (start_pos << 1) | strand));
                        }
                    } else {
                        valid = 0;
                        kf = 0;
                        kr = 0;
                    }
                }
                parts
            })
            .collect();

        let mut part_vecs: Vec<Vec<(u64, u32)>> = (0..256)
            .map(|p| {
                let cap: usize = per_chunk_parts.iter().map(|c| c[p].len()).sum();
                let mut v = Vec::with_capacity(cap);
                for c in &per_chunk_parts {
                    v.extend_from_slice(&c[p]);
                }
                v
            })
            .collect();
        drop(per_chunk_parts);

        part_vecs.par_iter_mut().for_each(|v| {
            v.sort_unstable_by_key(|&(h, ps)| (h, ps));
            if freq_limit > 0 && !v.is_empty() {
                let fl = freq_limit as usize;
                let mut write = 0usize;
                let mut i = 0usize;
                let n = v.len();
                while i < n {
                    let h = v[i].0;
                    let mut j = i + 1;
                    while j < n && v[j].0 == h {
                        j += 1;
                    }
                    let run_len = j - i;
                    if run_len <= fl {
                        v.copy_within(i..j, write);
                        write += run_len;
                    } else {
                        let mut c0 = 0usize;
                        let mut c1 = 0usize;
                        for &(_, ps) in &v[i..j] {
                            if (ps & 1) == 0 {
                                c0 += 1;
                            } else {
                                c1 += 1;
                            }
                        }
                        let k0 = c0 <= fl;
                        let k1 = c1 <= fl;
                        for m in i..j {
                            let keep = if (v[m].1 & 1) == 0 { k0 } else { k1 };
                            if keep {
                                v[write] = v[m];
                                write += 1;
                            }
                        }
                    }
                    i = j;
                }
                v.truncate(write);
            }
        });

        let n_kept: u32 = part_vecs.iter().map(|v| v.len() as u32).sum();

        let mut bucket_bits: u32 = 12;
        while bucket_bits < 24 && (1u64 << bucket_bits) * 4 < (n_kept as u64) {
            bucket_bits += 1;
        }
        let num_buckets = 1usize << bucket_bits;
        let shift_bits = 64 - bucket_bits;

        let mut keys = Vec::with_capacity(n_kept.max(1) as usize);
        let mut vals = Vec::with_capacity(n_kept.max(1) as usize);
        let mut bucket_off = vec![0u32; num_buckets + 1];

        for v in part_vecs {
            for (h, ps) in v {
                let b = (h >> shift_bits) as usize;
                bucket_off[b + 1] += 1;
                keys.push(KeyU64Pair {
                    lo: h as u32,
                    hi: (h >> 32) as u32,
                });
                vals.push(ps);
            }
        }
        for b in 1..=num_buckets {
            bucket_off[b] += bucket_off[b - 1];
        }
        if keys.is_empty() {
            keys.push(KeyU64Pair { lo: 0, hi: 0 });
            vals.push(0);
        }

        Self {
            shift_bits,
            bucket_off,
            keys,
            vals,
        }
    }
}
