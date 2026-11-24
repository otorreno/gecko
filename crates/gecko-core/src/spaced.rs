use std::fmt;
use std::str::FromStr;

use crate::error::{GeckoError, Result};

/// Symmetric (palindromic) spaced-seed pattern for high-sensitivity divergent genome comparison.
///
/// Care positions are marked with `'1'` (or `'#'`) and don't-care positions with `'0'` (or `'_'`).
/// Requiring palindrome symmetry (`pattern[i] == pattern[span - 1 - i]`) guarantees that the
/// forward-strand MSB-first mask, reverse-complement MSB-first mask, and LSB-first `load_32`
/// verification mask are all identical (`mask_2bit`), enabling single-instruction `kf & mask_2bit`
/// canonical spaced-seed hashing on both CPU and GPU with zero extra memory.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SpacedSeed {
    /// Canonical `'1'`/`'0'` pattern string (length = `span`).
    pub pattern: String,
    /// Total window width in bp (`4..=32`).
    pub span: usize,
    /// Number of care (`'1'`) bases (`4..=span`).
    pub weight: usize,
    /// 2-bit per base care mask (`0b11` at every `'1'` position, `0b00` at every `'0'` position).
    /// Because `pattern` is symmetric, this mask applies identically to MSB-first rolling k-mers
    /// (`kf`, `kr`) and LSB-first `PackedSeq::load_32` words.
    pub mask_2bit: u64,
}

impl SpacedSeed {
    /// Named presets tuned for divergent genomic comparisons (all strictly palindromic).
    pub const PRESETS: &'static [(&'static str, &'static str)] = &[
        ("w12s18", "111010110011010111"),
        ("w14s20", "11110101100110101111"),
        ("w16s22", "1111011010110101101111"),
        ("w18s26", "11110110101100110101101111"),
        ("w20s28", "1111011011011001101101101111"),
    ];

    /// Parse a named preset (e.g. `"w14s20"`, `"w18s26"`) or an explicit `'1'`/`'0'` (`'#'`/`'_'`)
    /// palindromic spaced-seed mask.
    pub fn parse(spec: &str) -> Result<Self> {
        let trimmed = spec.trim();
        let raw = Self::PRESETS
            .iter()
            .find(|(name, _)| name.eq_ignore_ascii_case(trimmed))
            .map(|&(_, pat)| pat)
            .unwrap_or(trimmed);

        let mut canonical = String::with_capacity(raw.len());
        for ch in raw.chars() {
            match ch {
                '1' | '#' => canonical.push('1'),
                '0' | '_' | '-' => canonical.push('0'),
                _ => {
                    return Err(GeckoError::InvalidParameter(format!(
                        "invalid character '{ch}' in spaced seed '{spec}' (use '1'/'0', '#'/'_', or preset: w12s18, w14s20, w16s22, w18s26, w20s28)"
                    )));
                }
            }
        }

        let span = canonical.len();
        if !(4..=32).contains(&span) {
            return Err(GeckoError::InvalidParameter(format!(
                "spaced seed span must be between 4 and 32 bp, got {span} ('{spec}')"
            )));
        }
        let bytes = canonical.as_bytes();
        if bytes[0] != b'1' || bytes[span - 1] != b'1' {
            return Err(GeckoError::InvalidParameter(format!(
                "spaced seed must start and end with a care ('1') position: '{canonical}'"
            )));
        }
        for i in 0..span / 2 {
            if bytes[i] != bytes[span - 1 - i] {
                return Err(GeckoError::InvalidParameter(format!(
                    "spaced seed must be palindromic (symmetric) for canonical dual-strand indexing: '{canonical}' differs at positions {i} and {}",
                    span - 1 - i
                )));
            }
        }

        let mut weight = 0usize;
        let mut mask_2bit = 0u64;
        for (i, &b) in bytes.iter().enumerate() {
            if b == b'1' {
                weight += 1;
                mask_2bit |= 3u64 << (2 * i);
            }
        }
        if weight < 4 {
            return Err(GeckoError::InvalidParameter(format!(
                "spaced seed weight (number of '1' care positions) must be at least 4, got {weight}"
            )));
        }

        Ok(Self {
            pattern: canonical,
            span,
            weight,
            mask_2bit,
        })
    }

    /// Returns `true` if the pattern contains at least one `'0'` don't-care position.
    #[inline]
    pub fn is_gapped(&self) -> bool {
        self.weight < self.span
    }

    /// Returns `Some(mask_2bit)` when gapped, or `None` if all positions are `'1'` (contiguous).
    #[inline]
    pub fn gapped_mask_2bit(&self) -> Option<u64> {
        if self.is_gapped() {
            Some(self.mask_2bit)
        } else {
            None
        }
    }
}

impl FromStr for SpacedSeed {
    type Err = GeckoError;

    fn from_str(s: &str) -> Result<Self> {
        Self::parse(s)
    }
}

impl fmt::Display for SpacedSeed {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} (w={}/s={})", self.pattern, self.weight, self.span)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn all_presets_parse_and_are_symmetric() {
        for &(name, pat) in SpacedSeed::PRESETS {
            let s = SpacedSeed::parse(name).expect(name);
            assert_eq!(s.pattern, pat);
            assert!(s.is_gapped());
            assert!(s.weight >= 12 && s.span <= 32);
            // Verify MSB-first and LSB-first 2-bit masks are identical.
            let mut msb_mask = 0u64;
            for (i, b) in s.pattern.bytes().enumerate() {
                if b == b'1' {
                    msb_mask |= 3u64 << (2 * (s.span - 1 - i));
                }
            }
            assert_eq!(s.mask_2bit, msb_mask, "preset {name} symmetry mismatch");
        }
    }

    #[test]
    fn rejects_asymmetric_or_invalid_masks() {
        assert!(SpacedSeed::parse("111011").is_err());
        assert!(SpacedSeed::parse("011110").is_err());
        assert!(SpacedSeed::parse("111").is_err());
        assert!(SpacedSeed::parse("100001").is_err());
        let contig = SpacedSeed::parse("11111111").unwrap();
        assert!(!contig.is_gapped());
        assert_eq!(contig.gapped_mask_2bit(), None);
        let alt = SpacedSeed::parse("##_##").unwrap();
        assert_eq!(alt.pattern, "11011");
        assert_eq!(alt.weight, 4);
        assert_eq!(alt.span, 5);
    }
}
