use rayon::prelude::*;

/// Threshold below which standard unstable sort beats radix setup overhead.
const SMALL_SORT_THRESHOLD: usize = 2048;

/// Parallel radix sort for items with a `(primary: u64, secondary: u64)` key pair.
///
/// Uses a top-level 256-way MSD scatter on the highest non-constant byte of `primary`
/// (parallelized across buckets with `rayon`), followed by an in-bucket LSD radix sort
/// that skips constant bytes. Total complexity is `O(n)` with linear cache access.
pub fn parallel_radix_sort_u128_key<T, F>(data: &mut [T], key_fn: F)
where
    T: Copy + Send + Sync,
    F: Fn(&T) -> (u64, u64) + Sync,
{
    let n = data.len();
    if n <= 1 {
        return;
    }
    if n <= SMALL_SORT_THRESHOLD {
        data.sort_unstable_by_key(|x| key_fn(x));
        return;
    }

    // Find which bytes actually vary across the dataset so we skip constant passes.
    let (p_first, s_first) = key_fn(&data[0]);
    let (p_diff, s_diff) = data
        .par_chunks(65536)
        .map(|chunk| {
            let mut pd = 0u64;
            let mut sd = 0u64;
            for item in chunk {
                let (p, s) = key_fn(item);
                pd |= p ^ p_first;
                sd |= s ^ s_first;
            }
            (pd, sd)
        })
        .reduce(|| (0, 0), |(a1, b1), (a2, b2)| (a1 | a2, b1 | b2));

    if p_diff == 0 && s_diff == 0 {
        return;
    }

    // Determine the highest varying byte (15..=8 => primary, 7..=0 => secondary)
    let top_byte = highest_varying_byte(p_diff, s_diff);

    // Parallel count per 256-bucket on top_byte
    let chunk_size = (n / (rayon::current_num_threads() * 4).max(1)).clamp(4096, 131_072);
    let chunk_counts: Vec<[usize; 256]> = data
        .par_chunks(chunk_size)
        .map(|chunk| {
            let mut counts = [0usize; 256];
            for item in chunk {
                let b = extract_byte(key_fn(item), top_byte);
                counts[b] += 1;
            }
            counts
        })
        .collect();

    let mut total_counts = [0usize; 256];
    for c in &chunk_counts {
        for i in 0..256 {
            total_counts[i] += c[i];
        }
    }

    let mut bucket_starts = [0usize; 256];
    let mut sum = 0usize;
    for i in 0..256 {
        bucket_starts[i] = sum;
        sum += total_counts[i];
    }

    // Scatter into auxiliary buffer
    let mut aux = vec![data[0]; n];
    {
        let mut cursor = bucket_starts;
        for item in data.iter() {
            let b = extract_byte(key_fn(item), top_byte);
            aux[cursor[b]] = *item;
            cursor[b] += 1;
        }
    }

    // Copy partitioned buckets back into `data` and sort each bucket in parallel.
    data.copy_from_slice(&aux);

    // Build non-empty bucket slices and sort each bucket in parallel via LSD radix sort.
    let mut slices: Vec<&mut [T]> = Vec::with_capacity(256);
    let mut rest = data;
    for &count in &total_counts {
        let (head, tail) = rest.split_at_mut(count);
        if count > 1 {
            slices.push(head);
        }
        rest = tail;
    }

    slices.into_par_iter().for_each(|bucket| {
        lsd_radix_sort_bucket(bucket, top_byte, &key_fn);
    });
}

fn lsd_radix_sort_bucket<T, F>(bucket: &mut [T], max_byte_exclusive: usize, key_fn: &F)
where
    T: Copy,
    F: Fn(&T) -> (u64, u64),
{
    let n = bucket.len();
    if n <= 1 || max_byte_exclusive == 0 {
        return;
    }
    if n <= SMALL_SORT_THRESHOLD {
        bucket.sort_unstable_by_key(|x| key_fn(x));
        return;
    }

    let (p0, s0) = key_fn(&bucket[0]);
    let mut p_diff = 0u64;
    let mut s_diff = 0u64;
    for item in bucket.iter() {
        let (p, s) = key_fn(item);
        p_diff |= p ^ p0;
        s_diff |= s ^ s0;
    }

    let mut active_bytes = Vec::with_capacity(max_byte_exclusive);
    for byte_idx in 0..max_byte_exclusive {
        let mask = 0xffu64 << ((byte_idx % 8) * 8);
        let diff = if byte_idx >= 8 { p_diff } else { s_diff };
        if (diff & mask) != 0 {
            active_bytes.push(byte_idx);
        }
    }
    if active_bytes.is_empty() {
        return;
    }

    let mut buf = vec![bucket[0]; n];
    let mut src_is_bucket = true;

    for byte_idx in active_bytes {
        let (src, dst): (&[T], &mut [T]) = if src_is_bucket {
            (&*bucket, &mut buf)
        } else {
            (&buf, &mut *bucket)
        };

        let mut counts = [0usize; 256];
        for item in src {
            counts[extract_byte(key_fn(item), byte_idx)] += 1;
        }
        let mut offsets = [0usize; 256];
        let mut sum = 0usize;
        for i in 0..256 {
            offsets[i] = sum;
            sum += counts[i];
        }
        for item in src {
            let b = extract_byte(key_fn(item), byte_idx);
            dst[offsets[b]] = *item;
            offsets[b] += 1;
        }

        src_is_bucket = !src_is_bucket;
    }

    if !src_is_bucket {
        bucket.copy_from_slice(&buf);
    }
}

#[inline]
fn highest_varying_byte(p_diff: u64, s_diff: u64) -> usize {
    if p_diff != 0 {
        8 + (63 - p_diff.leading_zeros() as usize) / 8
    } else {
        (63 - s_diff.leading_zeros() as usize) / 8
    }
}

#[inline]
fn extract_byte((primary, secondary): (u64, u64), byte_idx: usize) -> usize {
    if byte_idx >= 8 {
        ((primary >> ((byte_idx - 8) * 8)) & 0xff) as usize
    } else {
        ((secondary >> (byte_idx * 8)) & 0xff) as usize
    }
}

/// Convert a signed `i64` into an order-preserving `u64` key for radix sort.
#[inline]
pub fn i64_to_sortable_u64(v: i64) -> u64 {
    (v as u64) ^ (1u64 << 63)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sorts_pairs_matching_std() {
        let mut items: Vec<(u64, u64)> = (0..25_000u64)
            .map(|i| {
                let p = i
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                let s = i.wrapping_mul(2862933555777941757).wrapping_add(3037000493);
                (p, s)
            })
            .collect();
        let mut expected = items.clone();
        expected.sort_unstable();

        parallel_radix_sort_u128_key(&mut items, |&(p, s)| (p, s));
        assert_eq!(items, expected);
    }

    #[test]
    fn sorts_signed_diagonals() {
        let mut vals = vec![i64::MAX, -1, 0, i64::MIN, 42, -999_999, 12345];
        let mut expected = vals.clone();
        expected.sort_unstable();
        parallel_radix_sort_u128_key(&mut vals, |&v| (i64_to_sortable_u64(v), 0));
        assert_eq!(vals, expected);
    }
}
