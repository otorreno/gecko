//! Opaque-box 4-tier E2E test suite for GECKO v2 (Wave 3: CPU & CUDA Profiling,
//! Lossless Parity & Next-Level Optimization).
//!
//! Design rules:
//! - Strictly opaque-box: invokes the compiled `gecko` binary (and `cargo`/`ar`/`cuobjdump`
//!   for CUDA static verification) with zero imports from `gecko_core` or `gecko_cuda`.
//! - Deterministic on CPU-only VMs: all comparison tests pass `--backend cpu` and strip
//!   inherited `GECKO_PROFILE` / `GECKO_BACKEND` env vars.
//! - Isolated hygiene: every test writes solely inside its own RAII `tempfile::TempDir`
//!   and never modifies repository files, `.agents/teamwork/`, or untracked benchmark FASTAs.
//! - 4-tier naming (33 tests total):
//!   - `tier1_*` (18 tests: 12 CPU feature coverage + 6 CUDA static SASS/PTX/build contracts)
//!   - `tier2_*` (7 tests: boundary, corner & negative cases)
//!   - `tier3_*` (3 tests: cross-feature combinations)
//!   - `tier4_*` (5 tests: real-world primate & mammalian chromosome workloads)

use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::OnceLock;

// ============================================================================
// Common Opaque-Box Process, DNA Synthesis & FASTA Helpers
// ============================================================================

fn gecko_bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_gecko"))
}

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("crates dir")
        .parent()
        .expect("repo root")
        .to_path_buf()
}

fn run_gecko(args: &[&str]) -> Output {
    Command::new(gecko_bin())
        .args(args)
        .env_remove("GECKO_PROFILE")
        .env_remove("GECKO_BACKEND")
        .output()
        .expect("failed to spawn gecko binary")
}

fn run_gecko_with_env(args: &[&str], envs: &[(&str, Option<&str>)]) -> Output {
    let mut cmd = Command::new(gecko_bin());
    cmd.args(args)
        .env_remove("GECKO_PROFILE")
        .env_remove("GECKO_BACKEND");
    for &(k, v) in envs {
        match v {
            Some(val) => {
                cmd.env(k, val);
            }
            None => {
                cmd.env_remove(k);
            }
        }
    }
    cmd.output().expect("failed to spawn gecko binary")
}

fn assert_success(out: &Output, ctx: &str) {
    assert!(
        out.status.success(),
        "{ctx}: expected exit 0, got {:?}\nstdout:\n{}\nstderr:\n{}",
        out.status.code(),
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr),
    );
}

fn assert_clean_failure(out: &Output, ctx: &str) {
    let code = out.status.code();
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        !out.status.success(),
        "{ctx}: expected non-zero exit, got 0\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&out.stdout),
        stderr,
    );
    assert_ne!(
        code,
        Some(101),
        "{ctx}: binary panicked (exit 101)\nstderr:\n{stderr}"
    );
    assert!(
        !stderr.contains("panicked at") && !stderr.contains("RUST_BACKTRACE"),
        "{ctx}: stderr contains panic markers:\n{stderr}"
    );
}

fn pseudo_dna(seed: u64, len: usize) -> String {
    const BASES: [u8; 4] = *b"ACGT";
    let mut x = seed;
    let mut out = Vec::with_capacity(len);
    for _ in 0..len {
        x = x
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        out.push(BASES[((x >> 32) & 3) as usize]);
    }
    String::from_utf8(out).unwrap()
}

fn complement_base(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' | b'U' => b'A',
        b'a' => b't',
        b'c' => b'g',
        b'g' => b'c',
        b't' | b'u' => b'a',
        other => other,
    }
}

fn revcomp(s: &str) -> String {
    s.bytes()
        .rev()
        .map(complement_base)
        .map(|b| b as char)
        .collect()
}

fn alt_base(b: u8) -> u8 {
    match b {
        b'A' => b'C',
        b'C' => b'G',
        b'G' => b'T',
        b'T' => b'A',
        b'a' => b'c',
        b'c' => b'g',
        b'g' => b't',
        b't' => b'a',
        other => other,
    }
}

fn all_mismatch(orig: &str) -> String {
    orig.bytes().map(alt_base).map(|b| b as char).collect()
}

fn mutate_every(s: &str, every: usize, offset: usize) -> String {
    s.bytes()
        .enumerate()
        .map(|(i, b)| {
            if i >= offset && (i - offset).is_multiple_of(every) {
                alt_base(b)
            } else {
                b
            }
        })
        .map(|b| b as char)
        .collect()
}

fn mutate_positions(s: &str, positions: &[usize]) -> String {
    let set: HashSet<usize> = positions.iter().copied().collect();
    s.bytes()
        .enumerate()
        .map(|(i, b)| if set.contains(&i) { alt_base(b) } else { b })
        .map(|b| b as char)
        .collect()
}

/// Repeat `(match_run matches, mismatch_run mismatches)` with period `match_run + mismatch_run`.
fn pattern_mutate(orig: &str, match_run: usize, mismatch_run: usize) -> String {
    let period = match_run + mismatch_run;
    assert!(period > 0);
    orig.bytes()
        .enumerate()
        .map(|(i, b)| {
            if i % period < match_run {
                b
            } else {
                alt_base(b)
            }
        })
        .map(|b| b as char)
        .collect()
}

/// Mirror of `pattern_mutate` for backward-entering flanks: repeats `(mismatch_run, match_run)`.
fn pattern_mutate_rev(orig: &str, mismatch_run: usize, match_run: usize) -> String {
    let period = match_run + mismatch_run;
    assert!(period > 0);
    orig.bytes()
        .enumerate()
        .map(|(i, b)| {
            if i % period < mismatch_run {
                alt_base(b)
            } else {
                b
            }
        })
        .map(|b| b as char)
        .collect()
}

/// Encode `v` in `0..1024` as a unique 5-base ACGT string.
fn encode_base4_5mer(v: usize) -> String {
    const BASES: [u8; 4] = *b"ACGT";
    assert!(v < 1024, "encode_base4_5mer expects v < 1024, got {v}");
    let mut out = [0u8; 5];
    for (k, slot) in out.iter_mut().enumerate() {
        *slot = BASES[(v >> (2 * k)) & 3];
    }
    String::from_utf8(out.to_vec()).unwrap()
}

fn write_fasta(path: &Path, records: &[(&str, &str)]) {
    let mut f = File::create(path).expect("create fasta");
    for (header, seq) in records {
        writeln!(f, ">{header}").unwrap();
        for chunk in seq.as_bytes().chunks(80) {
            f.write_all(chunk).unwrap();
            f.write_all(b"\n").unwrap();
        }
        if seq.is_empty() {
            f.write_all(b"\n").unwrap();
        }
    }
}

/// Streamingly slice `[start_base..end_base)` (0-based) from a single-record FASTA file
/// into a temporary FASTA file without ever modifying the source file.
fn slice_real_fasta_to_temp(
    src_path: &Path,
    dst_path: &Path,
    header: &str,
    start_base: usize,
    end_base: usize,
) -> usize {
    let file = File::open(src_path).expect("open real fasta");
    let mut reader = BufReader::with_capacity(256 * 1024, file);
    let mut out = File::create(dst_path).expect("create sliced fasta");
    writeln!(out, ">{header}").unwrap();

    let mut pos = 0usize;
    let mut written = 0usize;
    let mut line = String::new();
    let mut buf = Vec::with_capacity(60);

    loop {
        line.clear();
        let n = reader.read_line(&mut line).expect("read fasta line");
        if n == 0 {
            break;
        }
        if line.starts_with('>') {
            continue;
        }
        let trimmed = line.trim_end_matches(['\r', '\n']);
        let line_len = trimmed.len();
        let line_end = pos + line_len;
        if line_end <= start_base {
            pos = line_end;
            continue;
        }
        if pos >= end_base {
            break;
        }
        let s = start_base.saturating_sub(pos);
        let e = (end_base - pos).min(line_len);
        for &b in &trimmed.as_bytes()[s..e] {
            buf.push(b);
            written += 1;
            if buf.len() == 60 {
                out.write_all(&buf).unwrap();
                out.write_all(b"\n").unwrap();
                buf.clear();
            }
        }
        pos = line_end;
        if pos >= end_base {
            break;
        }
    }
    if !buf.is_empty() {
        out.write_all(&buf).unwrap();
        out.write_all(b"\n").unwrap();
    }
    written
}

// ============================================================================
// Pure-Rust FIPS 180-4 SHA-256 Helper (Zero External Crates)
// ============================================================================

fn sha256_bytes(data: &[u8]) -> [u8; 32] {
    const K: [u32; 64] = [
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4,
        0xab1c5ed5, 0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe,
        0x9bdc06a7, 0xc19bf174, 0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f,
        0x4a7484aa, 0x5cb0a9dc, 0x76f988da, 0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
        0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967, 0x27b70a85, 0x2e1b2138, 0x4d2c6dfc,
        0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85, 0xa2bfe8a1, 0xa81a664b,
        0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070, 0x19a4c116,
        0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
        0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7,
        0xc67178f2,
    ];

    let mut h: [u32; 8] = [
        0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab,
        0x5be0cd19,
    ];

    let bit_len = (data.len() as u64).wrapping_mul(8);
    let mut padded = Vec::with_capacity(data.len() + 72);
    padded.extend_from_slice(data);
    padded.push(0x80);
    while (padded.len() % 64) != 56 {
        padded.push(0);
    }
    padded.extend_from_slice(&bit_len.to_be_bytes());

    for chunk in padded.as_chunks::<64>().0 {
        let mut w = [0u32; 64];
        for (i, word) in chunk.as_chunks::<4>().0.iter().enumerate() {
            w[i] = u32::from_be_bytes(*word);
        }
        for i in 16..64 {
            let s0 = w[i - 15].rotate_right(7) ^ w[i - 15].rotate_right(18) ^ (w[i - 15] >> 3);
            let s1 = w[i - 2].rotate_right(17) ^ w[i - 2].rotate_right(19) ^ (w[i - 2] >> 10);
            w[i] = w[i - 16]
                .wrapping_add(s0)
                .wrapping_add(w[i - 7])
                .wrapping_add(s1);
        }

        let [mut a, mut b, mut c, mut d, mut e, mut f, mut g, mut hh] = h;
        for i in 0..64 {
            let s1 = e.rotate_right(6) ^ e.rotate_right(11) ^ e.rotate_right(25);
            let ch = (e & f) ^ ((!e) & g);
            let temp1 = hh
                .wrapping_add(s1)
                .wrapping_add(ch)
                .wrapping_add(K[i])
                .wrapping_add(w[i]);
            let s0 = a.rotate_right(2) ^ a.rotate_right(13) ^ a.rotate_right(22);
            let maj = (a & b) ^ (a & c) ^ (b & c);
            let temp2 = s0.wrapping_add(maj);

            hh = g;
            g = f;
            f = e;
            e = d.wrapping_add(temp1);
            d = c;
            c = b;
            b = a;
            a = temp1.wrapping_add(temp2);
        }

        h[0] = h[0].wrapping_add(a);
        h[1] = h[1].wrapping_add(b);
        h[2] = h[2].wrapping_add(c);
        h[3] = h[3].wrapping_add(d);
        h[4] = h[4].wrapping_add(e);
        h[5] = h[5].wrapping_add(f);
        h[6] = h[6].wrapping_add(g);
        h[7] = h[7].wrapping_add(hh);
    }

    let mut digest = [0u8; 32];
    for (i, word) in h.iter().enumerate() {
        digest[i * 4..(i + 1) * 4].copy_from_slice(&word.to_be_bytes());
    }
    digest
}

fn sha256_hex(data: &[u8]) -> String {
    const HEX: &[u8; 16] = b"0123456789abcdef";
    let digest = sha256_bytes(data);
    let mut out = String::with_capacity(64);
    for b in digest {
        out.push(HEX[(b >> 4) as usize] as char);
        out.push(HEX[(b & 0xf) as usize] as char);
    }
    out
}

fn verify_sha256_self_test() {
    static VERIFIED: OnceLock<()> = OnceLock::new();
    VERIFIED.get_or_init(|| {
        assert_eq!(
            sha256_hex(b"abc"),
            "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad",
            "FIPS 180-4 SHA-256 self-test failed for b\"abc\""
        );
        assert_eq!(
            sha256_hex(b""),
            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
            "FIPS 180-4 SHA-256 self-test failed for empty input"
        );
    });
}

fn assert_frags_sha256_identical(a: &Path, b: &Path, ctx: &str) {
    verify_sha256_self_test();
    let ba = fs::read(a).unwrap_or_else(|e| panic!("{ctx}: read {}: {e}", a.display()));
    let bb = fs::read(b).unwrap_or_else(|e| panic!("{ctx}: read {}: {e}", b.display()));
    let ha = sha256_hex(&ba);
    let hb = sha256_hex(&bb);
    if ha != hb || ba != bb {
        let fa = parse_frags_file(a);
        let fb = parse_frags_file(b);
        panic!(
            "{ctx}: .frags SHA-256 mismatch:\n  A ({}): sha256={ha} ({} bytes, {} frags)\n  B ({}): sha256={hb} ({} bytes, {} frags)",
            a.display(),
            ba.len(),
            fa.fragments.len(),
            b.display(),
            bb.len(),
            fb.fragments.len(),
        );
    }
}

// ============================================================================
// Standalone Opaque-Box Parsers (`.frags`, `.MAT`, All 9 Formats)
// ============================================================================

const FRAG_WIRE_SIZE: usize = 93;

#[derive(Debug, Clone)]
struct ParsedFrag {
    diag: i64,
    x_start: u64,
    y_start: u64,
    x_end: u64,
    y_end: u64,
    length: u64,
    ident: u64,
    score: u64,
    similarity: f32,
    seq_x: u64,
    seq_y: u64,
    block: i64,
    strand: char,
}

impl PartialEq for ParsedFrag {
    fn eq(&self, other: &Self) -> bool {
        self.diag == other.diag
            && self.x_start == other.x_start
            && self.y_start == other.y_start
            && self.x_end == other.x_end
            && self.y_end == other.y_end
            && self.length == other.length
            && self.ident == other.ident
            && self.score == other.score
            && self.similarity.to_bits() == other.similarity.to_bits()
            && self.seq_x == other.seq_x
            && self.seq_y == other.seq_y
            && self.block == other.block
            && self.strand == other.strand
    }
}

impl Eq for ParsedFrag {}

impl ParsedFrag {
    fn pct_identity(&self) -> f64 {
        if self.length == 0 {
            0.0
        } else {
            (self.ident as f64) * 100.0 / (self.length as f64)
        }
    }
}

#[derive(Debug, Clone)]
struct ParsedFragsFile {
    x_total: u64,
    y_total: u64,
    fragments: Vec<ParsedFrag>,
}

fn parse_frags_file(path: &Path) -> ParsedFragsFile {
    let bytes = fs::read(path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
    assert!(
        bytes.len() >= 16,
        ".frags file {} shorter than 16-byte header ({} bytes)",
        path.display(),
        bytes.len()
    );
    assert_eq!(
        (bytes.len() - 16) % FRAG_WIRE_SIZE,
        0,
        ".frags file {} body size {} is not a multiple of {FRAG_WIRE_SIZE}",
        path.display(),
        bytes.len() - 16
    );
    let x_total = u64::from_be_bytes(bytes[0..8].try_into().unwrap());
    let y_total = u64::from_be_bytes(bytes[8..16].try_into().unwrap());
    let n = (bytes.len() - 16) / FRAG_WIRE_SIZE;
    let mut fragments = Vec::with_capacity(n);
    for i in 0..n {
        let b = &bytes[16 + i * FRAG_WIRE_SIZE..16 + (i + 1) * FRAG_WIRE_SIZE];
        let strand = match b[92] {
            b'f' => 'f',
            b'r' => 'r',
            other => panic!("invalid strand byte 0x{other:02x} in {}", path.display()),
        };
        fragments.push(ParsedFrag {
            diag: i64::from_be_bytes(b[0..8].try_into().unwrap()),
            x_start: u64::from_be_bytes(b[8..16].try_into().unwrap()),
            y_start: u64::from_be_bytes(b[16..24].try_into().unwrap()),
            x_end: u64::from_be_bytes(b[24..32].try_into().unwrap()),
            y_end: u64::from_be_bytes(b[32..40].try_into().unwrap()),
            length: u64::from_be_bytes(b[40..48].try_into().unwrap()),
            ident: u64::from_be_bytes(b[48..56].try_into().unwrap()),
            score: u64::from_be_bytes(b[56..64].try_into().unwrap()),
            similarity: f32::from_be_bytes(b[64..68].try_into().unwrap()),
            seq_x: u64::from_be_bytes(b[68..76].try_into().unwrap()),
            seq_y: u64::from_be_bytes(b[76..84].try_into().unwrap()),
            block: i64::from_be_bytes(b[84..92].try_into().unwrap()),
            strand,
        });
    }
    ParsedFragsFile {
        x_total,
        y_total,
        fragments,
    }
}

fn parse_mat_file(path: &Path) -> Vec<Vec<u64>> {
    let text = fs::read_to_string(path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(
        lines.len(),
        1000,
        ".MAT {} must have 1000 rows, got {}",
        path.display(),
        lines.len()
    );
    let mut rows = Vec::with_capacity(1000);
    for (r, line) in lines.iter().enumerate() {
        assert!(
            line.ends_with('\t'),
            ".MAT {} row {r} missing trailing tab",
            path.display()
        );
        let vals: Vec<u64> = line
            .split_whitespace()
            .map(|s| s.parse::<u64>().expect("valid u64 in .MAT"))
            .collect();
        assert_eq!(
            vals.len(),
            100,
            ".MAT {} row {r} must have 100 columns, got {}",
            path.display(),
            vals.len()
        );
        rows.push(vals);
    }
    rows
}

fn assert_mat_covers_frags(mat: &[Vec<u64>], frags: &[ParsedFrag]) {
    let total: u64 = mat.iter().flat_map(|r| r.iter()).copied().sum();
    assert!(
        total >= frags.len() as u64,
        ".MAT total ({total}) < fragment count ({})",
        frags.len()
    );
    for f in frags {
        let row = (f.length as usize).min(999);
        let col = (f.similarity.max(0.0) as usize).min(99);
        assert!(
            mat[row][col] >= 1,
            ".MAT[{row}][{col}] is 0 for accepted frag (len={}, sim={})",
            f.length,
            f.similarity
        );
    }
}

fn validate_all_9_formats(
    out_dir: &Path,
    stem: &str,
    expect_forward_only: bool,
) -> ParsedFragsFile {
    for suffix in [
        ".frags",
        ".csv",
        ".frags.INF",
        ".frags.MAT",
        ".paf",
        ".blast.tsv",
        ".bedpe",
        ".jsonl",
        ".png",
    ] {
        let p = out_dir.join(format!("{stem}{suffix}"));
        assert!(
            p.exists(),
            "missing expected --format all output {}",
            p.display()
        );
    }

    let ff = parse_frags_file(&out_dir.join(format!("{stem}.frags")));
    let n = ff.fragments.len();

    // .frags.INF (14 lines, 56 '=' separator)
    let inf_text = fs::read_to_string(out_dir.join(format!("{stem}.frags.INF"))).unwrap();
    let inf_lines: Vec<&str> = inf_text.lines().collect();
    assert_eq!(inf_lines.len(), 14, ".frags.INF must have 14 lines");
    assert_eq!(inf_lines[13], "=".repeat(56));
    assert!(inf_lines[6].ends_with(&ff.x_total.to_string()));
    assert!(inf_lines[7].ends_with(&ff.y_total.to_string()));
    assert!(inf_lines[12].ends_with(&n.to_string()));

    // .csv (17 header lines + N rows, 14 columns)
    let csv_text = fs::read_to_string(out_dir.join(format!("{stem}.csv"))).unwrap();
    let csv_lines: Vec<&str> = csv_text.lines().collect();
    assert_eq!(csv_lines.len(), 17 + n);
    assert_eq!(&csv_lines[..14], &inf_lines[..14]);
    assert_eq!(csv_lines[14], "Total CSB: 0");
    assert_eq!(csv_lines[15], "=".repeat(56));
    assert_eq!(
        csv_lines[16],
        "Type,xStart,yStart,xEnd,yEnd,strand(f/r),block,length,score,ident,similarity,%ident,SeqX,SeqY"
    );
    for (f, row) in ff.fragments.iter().zip(&csv_lines[17..]) {
        let cols: Vec<&str> = row.split(',').collect();
        assert_eq!(cols.len(), 14, "CSV row must have 14 cols: {row}");
        assert_eq!(cols[0], "Frag");
        assert_eq!(cols[5], f.strand.to_string());
        assert_eq!(cols[7].parse::<u64>().unwrap(), f.length);
        assert_eq!(cols[8].parse::<u64>().unwrap(), f.score);
        assert_eq!(cols[9].parse::<u64>().unwrap(), f.ident);
        assert_eq!(cols[12].parse::<u64>().unwrap(), f.seq_x);
        assert_eq!(cols[13].parse::<u64>().unwrap(), f.seq_y);
    }

    // .frags.MAT (1000 x 100 covering all fragments)
    let mat = parse_mat_file(&out_dir.join(format!("{stem}.frags.MAT")));
    assert_mat_covers_frags(&mat, &ff.fragments);

    // .paf (15 columns)
    let paf_lines: Vec<String> = fs::read_to_string(out_dir.join(format!("{stem}.paf")))
        .unwrap()
        .lines()
        .map(str::to_string)
        .collect();
    assert_eq!(paf_lines.len(), n);
    for (f, line) in ff.fragments.iter().zip(&paf_lines) {
        let cols: Vec<&str> = line.split('\t').collect();
        assert_eq!(cols.len(), 15, "PAF must have 15 cols: {line}");
        let expected_strand = if f.strand == 'f' { "+" } else { "-" };
        assert_eq!(cols[4], expected_strand);
        if expect_forward_only {
            assert_eq!(cols[4], "+", "PAF must be '+' under --forward-only");
        }
        assert_eq!(cols[9].parse::<u64>().unwrap(), f.ident);
        assert_eq!(cols[10].parse::<u64>().unwrap(), f.length);
        assert_eq!(cols[11], "255");
        assert_eq!(cols[12], format!("AS:i:{}", f.score));
        assert!(cols[13].starts_with("si:f:"));
        assert_eq!(cols[14], format!("cg:Z:{}M", f.length));
    }

    // .blast.tsv (12 columns)
    let blast_lines: Vec<String> = fs::read_to_string(out_dir.join(format!("{stem}.blast.tsv")))
        .unwrap()
        .lines()
        .map(str::to_string)
        .collect();
    assert_eq!(blast_lines.len(), n);
    for (f, line) in ff.fragments.iter().zip(&blast_lines) {
        let cols: Vec<&str> = line.split('\t').collect();
        assert_eq!(cols.len(), 12, "BLAST must have 12 cols: {line}");
        assert_eq!(cols[2], format!("{:.2}", f.pct_identity()));
        assert_eq!(cols[3].parse::<u64>().unwrap(), f.length);
        assert_eq!(cols[4].parse::<u64>().unwrap(), f.length - f.ident);
        assert_eq!(cols[5], "0");
        let sstart: u64 = cols[8].parse().unwrap();
        let send: u64 = cols[9].parse().unwrap();
        if f.strand == 'f' {
            assert!(sstart <= send);
        } else {
            assert!(sstart >= send);
        }
    }

    // .bedpe (10 columns)
    let bedpe_lines: Vec<String> = fs::read_to_string(out_dir.join(format!("{stem}.bedpe")))
        .unwrap()
        .lines()
        .map(str::to_string)
        .collect();
    assert_eq!(bedpe_lines.len(), n);
    for (i, (f, line)) in ff.fragments.iter().zip(&bedpe_lines).enumerate() {
        let cols: Vec<&str> = line.split('\t').collect();
        assert_eq!(cols.len(), 10, "BEDPE must have 10 cols: {line}");
        assert_eq!(cols[6], format!("frag_{i}"));
        assert_eq!(cols[7].parse::<u64>().unwrap(), f.score);
        assert_eq!(cols[8], "+");
        assert_eq!(cols[9], if f.strand == 'f' { "+" } else { "-" });
    }

    // .jsonl (exact 16 keys)
    let expected_keys: HashSet<&str> = [
        "query_id",
        "query_start",
        "query_end",
        "query_len",
        "target_id",
        "target_start",
        "target_end",
        "target_len",
        "strand",
        "length",
        "ident",
        "pct_identity",
        "score",
        "similarity",
        "seq_x",
        "seq_y",
    ]
    .into_iter()
    .collect();
    let jsonl_lines: Vec<String> = fs::read_to_string(out_dir.join(format!("{stem}.jsonl")))
        .unwrap()
        .lines()
        .map(str::to_string)
        .collect();
    assert_eq!(jsonl_lines.len(), n);
    for line in &jsonl_lines {
        let obj: serde_json::Map<String, serde_json::Value> = serde_json::from_str(line).unwrap();
        let keys: HashSet<&str> = obj.keys().map(String::as_str).collect();
        assert_eq!(keys, expected_keys);
    }

    // .png (PNG signature + IHDR >= 200x200)
    let png = fs::read(out_dir.join(format!("{stem}.png"))).unwrap();
    assert!(png.len() > 32);
    assert_eq!(&png[0..8], &[137, 80, 78, 71, 13, 10, 26, 10]);
    assert_eq!(&png[12..16], b"IHDR");
    let w = u32::from_be_bytes(png[16..20].try_into().unwrap());
    let h = u32::from_be_bytes(png[20..24].try_into().unwrap());
    assert!(w >= 200 && h >= 200);

    ff
}

// ============================================================================
// Non-Brittle `GECKO_PROFILE=1` Stderr Parser
// ============================================================================

#[derive(Debug, Clone)]
struct ProfileLine {
    stage: String,
    kv: BTreeMap<String, String>,
}

fn parse_profile_lines(stderr: &str) -> Vec<ProfileLine> {
    let mut out = Vec::new();
    for line in stderr.lines() {
        let Some((_, rest)) = line.split_once("[GECKO_PROFILE]") else {
            continue;
        };
        let mut tokens = rest.split_whitespace();
        let Some(stage) = tokens.next() else {
            continue;
        };
        let mut kv = BTreeMap::new();
        for tok in tokens {
            if let Some((k, v)) = tok.split_once('=') {
                kv.insert(k.to_string(), v.to_string());
            }
        }
        out.push(ProfileLine {
            stage: stage.to_string(),
            kv,
        });
    }
    out
}

fn used_minimizer(lines: &[ProfileLine]) -> bool {
    lines.iter().any(|l| {
        l.stage == "minimizer_query"
            || l.stage == "minimizer_build"
            || l.stage == "minimizer_multislab"
            || (l.stage == "all_vs_all_index"
                && l.kv.get("mode").map(String::as_str) == Some("minimizer"))
    })
}

fn used_exact_stream(lines: &[ProfileLine]) -> bool {
    lines.iter().any(|l| {
        l.stage == "stream_query"
            || l.stage == "index_build"
            || l.stage == "index_multislab"
            || (l.stage == "all_vs_all_index"
                && l.kv.get("mode").map(String::as_str) == Some("exact"))
    })
}

fn sum_u64_key(lines: &[ProfileLine], stage: &str, key: &str) -> Option<u64> {
    let mut found = false;
    let mut sum = 0u64;
    for l in lines {
        if l.stage == stage
            && let Some(val) = l.kv.get(key)
            && let Ok(n) = val.parse::<u64>()
        {
            found = true;
            sum = sum.saturating_add(n);
        }
    }
    found.then_some(sum)
}

fn first_str_key<'a>(lines: &'a [ProfileLine], stage: &str, key: &str) -> Option<&'a str> {
    lines
        .iter()
        .find(|l| l.stage == stage)
        .and_then(|l| l.kv.get(key).map(String::as_str))
}

// ============================================================================
// Shared CUDA Static Verification Harness (OnceLock Cargo JSON + cuobjdump)
// ============================================================================

const CUDA_REQUIRED_ARCHS: [&str; 4] = ["sm_75", "sm_80", "sm_86", "sm_89"];

struct CudaCheckResult {
    success: bool,
    exit_code: Option<i32>,
    out_dir: Option<PathBuf>,
    cfgs: Vec<String>,
    linked_libs: Vec<String>,
    warnings: Vec<String>,
    errors: Vec<String>,
    stderr: String,
}

#[derive(Debug, Clone)]
struct CudaResourceInfo {
    member: String,
    fn_name: String,
    reg: u32,
    stack: u64,
    local: u64,
}

#[derive(Debug, Clone, Default)]
struct CudaSassStats {
    member: String,
    fn_name: String,
    total_insns: usize,
    stl_count: usize,
    ldl_count: usize,
    atoms_count: usize,
    atomg_count: usize,
}

fn cached_cuda_cargo_check() -> &'static CudaCheckResult {
    static RESULT: OnceLock<CudaCheckResult> = OnceLock::new();
    RESULT.get_or_init(|| {
        let out = Command::new("cargo")
            .args([
                "check",
                "-p",
                "gecko-cli",
                "--features",
                "cuda",
                "--message-format=json",
            ])
            .current_dir(repo_root())
            .env_remove("GECKO_CUDA_ARCH")
            .output()
            .expect("failed to run cargo check --features cuda");

        let stdout = String::from_utf8_lossy(&out.stdout);
        let stderr = String::from_utf8_lossy(&out.stderr).into_owned();
        let mut out_dir = None;
        let mut cfgs = Vec::new();
        let mut linked_libs = Vec::new();
        let mut warnings = Vec::new();
        let mut errors = Vec::new();

        for line in stdout.lines() {
            let Ok(v) = serde_json::from_str::<serde_json::Value>(line) else {
                continue;
            };
            let reason = v["reason"].as_str().unwrap_or("");
            let pkg = v["package_id"].as_str().unwrap_or("");
            if reason == "build-script-executed" && pkg.contains("gecko-cuda") {
                if let Some(d) = v["out_dir"].as_str() {
                    out_dir = Some(PathBuf::from(d));
                }
                if let Some(arr) = v["cfgs"].as_array() {
                    cfgs = arr
                        .iter()
                        .filter_map(|x| x.as_str().map(str::to_string))
                        .collect();
                }
                if let Some(arr) = v["linked_libs"].as_array() {
                    linked_libs = arr
                        .iter()
                        .filter_map(|x| x.as_str().map(str::to_string))
                        .collect();
                }
            } else if reason == "compiler-message"
                && (pkg.contains("gecko-cuda") || pkg.contains("gecko-cli"))
            {
                let level = v["message"]["level"].as_str().unwrap_or("");
                let rendered = v["message"]["rendered"].as_str().unwrap_or("").to_string();
                if level == "warning" {
                    warnings.push(rendered);
                } else if level == "error" {
                    errors.push(rendered);
                }
            }
        }

        CudaCheckResult {
            success: out.status.success(),
            exit_code: out.status.code(),
            out_dir,
            cfgs,
            linked_libs,
            warnings,
            errors,
            stderr,
        }
    })
}

fn cuda_archive_path() -> PathBuf {
    let res = cached_cuda_cargo_check();
    assert!(
        res.success,
        "cargo check -p gecko-cli --features cuda failed (exit {:?}):\nerrors: {:?}\nstderr:\n{}",
        res.exit_code, res.errors, res.stderr
    );
    let out_dir = res
        .out_dir
        .as_ref()
        .expect("gecko-cuda build script out_dir must be reported by cargo check");
    let archive = out_dir.join("libgecko_cuda_kernels.a");
    assert!(
        archive.exists(),
        "missing static CUDA kernel archive {}",
        archive.display()
    );
    archive
}

fn is_third_party_cuda_kernel(mangled: &str) -> bool {
    const THIRD_PARTY_PREFIXES: &[&str] = &[
        "_ZN3cub",                 // cub::* and CCCL 2.x/3.x cub::CUB_*_NS::*
        "_ZN6thrust",              // thrust::* and thrust::THRUST_*_NS::*
        "_ZN4cuda",                // cuda::* and libcu++ cuda::std::*
        "_ZN4cccl",                // cccl::*
        "_ZN3std",                 // std::*
        "_ZN6nvidia",              // nvidia::*
        "_ZN18cooperative_groups", // cooperative_groups::*
    ];
    THIRD_PARTY_PREFIXES.iter().any(|p| mangled.starts_with(p))
}

fn extract_member_name(trimmed: &str) -> Option<String> {
    let rest = trimmed.strip_prefix("member ")?;
    let without_colon = rest.strip_suffix(':').unwrap_or(rest);
    let member_name = without_colon
        .rsplit_once(':')
        .map(|(_, obj)| obj)
        .unwrap_or(without_colon);
    Some(member_name.to_string())
}

fn parse_cuobjdump_archive_members(
    stdout: &str,
    entry_marker: &str,
) -> BTreeMap<String, Vec<String>> {
    let mut by_member: BTreeMap<String, Vec<String>> = BTreeMap::new();
    let mut current_member: Option<String> = None;
    for line in stdout.lines() {
        let trimmed = line.trim();
        if let Some(member_name) = extract_member_name(trimmed) {
            by_member.entry(member_name.clone()).or_default();
            current_member = Some(member_name);
        } else if trimmed.contains(entry_marker)
            && let Some(member) = current_member.as_ref()
        {
            by_member
                .entry(member.clone())
                .or_default()
                .push(trimmed.to_string());
        }
    }
    by_member
}

fn load_cuda_resource_usage(arch: &str) -> Vec<CudaResourceInfo> {
    let archive = cuda_archive_path();
    let out = Command::new("cuobjdump")
        .args([
            "-arch",
            arch,
            "--dump-resource-usage",
            archive.to_str().unwrap(),
        ])
        .output()
        .expect("failed to execute cuobjdump --dump-resource-usage");
    assert!(
        out.status.success(),
        "cuobjdump -arch {arch} --dump-resource-usage failed:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8_lossy(&out.stdout);
    let mut entries = Vec::new();
    let mut current_member = String::new();
    let mut pending_fn: Option<String> = None;

    for line in stdout.lines() {
        let trimmed = line.trim();
        if let Some(member) = extract_member_name(trimmed) {
            current_member = member;
            pending_fn = None;
        } else if let Some(rest) = trimmed.strip_prefix("Function ")
            && let Some(name) = rest.strip_suffix(':')
        {
            pending_fn = Some(name.trim().to_string());
        } else if trimmed.contains("REG:")
            && let Some(fn_name) = pending_fn.take()
        {
            let mut reg = 0u32;
            let mut stack = 0u64;
            let mut local = 0u64;
            for tok in trimmed.split_whitespace() {
                if let Some(v) = tok.strip_prefix("REG:") {
                    reg = v.parse().expect("valid REG integer");
                } else if let Some(v) = tok.strip_prefix("STACK:") {
                    stack = v.parse().expect("valid STACK integer");
                } else if let Some(v) = tok.strip_prefix("LOCAL:") {
                    local = v.parse().expect("valid LOCAL integer");
                }
            }
            entries.push(CudaResourceInfo {
                member: current_member.clone(),
                fn_name,
                reg,
                stack,
                local,
            });
        }
    }
    entries
}

fn cached_cuda_resource_usage(arch: &str) -> &'static Vec<CudaResourceInfo> {
    static RES_75: OnceLock<Vec<CudaResourceInfo>> = OnceLock::new();
    static RES_80: OnceLock<Vec<CudaResourceInfo>> = OnceLock::new();
    static RES_86: OnceLock<Vec<CudaResourceInfo>> = OnceLock::new();
    static RES_89: OnceLock<Vec<CudaResourceInfo>> = OnceLock::new();
    match arch {
        "sm_75" => RES_75.get_or_init(|| load_cuda_resource_usage("sm_75")),
        "sm_80" => RES_80.get_or_init(|| load_cuda_resource_usage("sm_80")),
        "sm_86" => RES_86.get_or_init(|| load_cuda_resource_usage("sm_86")),
        "sm_89" => RES_89.get_or_init(|| load_cuda_resource_usage("sm_89")),
        other => panic!("unsupported arch cache key: {other}"),
    }
}

fn load_cuda_sass_stats(arch: &str) -> Vec<CudaSassStats> {
    let archive = cuda_archive_path();
    let out = Command::new("cuobjdump")
        .args(["-arch", arch, "--dump-sass", archive.to_str().unwrap()])
        .output()
        .expect("failed to execute cuobjdump --dump-sass");
    assert!(
        out.status.success(),
        "cuobjdump -arch {arch} --dump-sass failed:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8_lossy(&out.stdout);
    let mut stats = Vec::new();
    let mut current_member = String::new();
    let mut active: Option<CudaSassStats> = None;

    for line in stdout.lines() {
        let trimmed = line.trim();
        if let Some(member) = extract_member_name(trimmed) {
            if let Some(cur) = active.take() {
                stats.push(cur);
            }
            current_member = member;
        } else if let Some(rest) = trimmed.strip_prefix("Function :") {
            if let Some(cur) = active.take() {
                stats.push(cur);
            }
            active = Some(CudaSassStats {
                member: current_member.clone(),
                fn_name: rest.trim().to_string(),
                ..CudaSassStats::default()
            });
        } else if let Some(cur) = active.as_mut()
            && trimmed.starts_with("/*")
            && let Some((_, after_pc)) = trimmed.split_once("*/")
            && after_pc.contains(';')
        {
            let body = after_pc.split(';').next().unwrap_or("");
            let mut toks = body.split_whitespace();
            let Some(first) = toks.next() else {
                continue;
            };
            let opcode = if first.starts_with('@') {
                let Some(next) = toks.next() else {
                    continue;
                };
                next
            } else {
                first
            };
            if opcode.starts_with('.') {
                continue;
            }
            let base_op = opcode.split('.').next().unwrap_or(opcode);
            cur.total_insns += 1;
            match base_op {
                "STL" => cur.stl_count += 1,
                "LDL" => cur.ldl_count += 1,
                "ATOMS" | "REDS" => cur.atoms_count += 1,
                "ATOMG" => cur.atomg_count += 1,
                _ => {}
            }
        }
    }
    if let Some(cur) = active.take() {
        stats.push(cur);
    }
    stats
}

fn cached_cuda_sass_stats(arch: &str) -> &'static Vec<CudaSassStats> {
    static SASS_75: OnceLock<Vec<CudaSassStats>> = OnceLock::new();
    static SASS_80: OnceLock<Vec<CudaSassStats>> = OnceLock::new();
    static SASS_86: OnceLock<Vec<CudaSassStats>> = OnceLock::new();
    static SASS_89: OnceLock<Vec<CudaSassStats>> = OnceLock::new();
    match arch {
        "sm_75" => SASS_75.get_or_init(|| load_cuda_sass_stats("sm_75")),
        "sm_80" => SASS_80.get_or_init(|| load_cuda_sass_stats("sm_80")),
        "sm_86" => SASS_86.get_or_init(|| load_cuda_sass_stats("sm_86")),
        "sm_89" => SASS_89.get_or_init(|| load_cuda_sass_stats("sm_89")),
        other => panic!("unsupported arch cache key: {other}"),
    }
}

// ============================================================================
// Shared Synthetic Fixture Builders for Wave 3 Workloads
// ============================================================================

/// Workload 3 (Explorer 1): 4-record multi-FASTA testing `update_noise` across `SEQ_SEP`
/// for both Case 3A (long backward-flank rejected HSP, `length > min_len`) and Case 3B
/// (short noise hit `length == 32 <= min_len`).
fn write_seq_sep_update_noise_fixture(dir: &Path) -> (PathBuf, PathBuf) {
    let r0_q = pseudo_dna(0xBEEF_0001, 400);
    let r0_r = format!("{}{}", pattern_mutate(&r0_q[..368], 3, 2), &r0_q[368..400]);

    let r1_q = pseudo_dna(0xBEEF_0002, 282);
    let r1_r = format!("{}{}", &r1_q[..32], pattern_mutate(&r1_q[32..282], 15, 1));

    let r2_q = pseudo_dna(0xBEEF_0003, 32);
    let r2_r = r2_q.clone();

    let r3_q = pseudo_dna(0xBEEF_0004, 282);
    let r3_r = format!("{}{}", &r3_q[..32], pattern_mutate(&r3_q[32..282], 15, 1));

    let q_path = dir.join("q.fa");
    let r_path = dir.join("r.fa");
    write_fasta(
        &q_path,
        &[("q0", &r0_q), ("q1", &r1_q), ("q2", &r2_q), ("q3", &r3_q)],
    );
    write_fasta(
        &r_path,
        &[("r0", &r0_r), ("r1", &r1_r), ("r2", &r2_r), ("r3", &r3_r)],
    );
    (q_path, r_path)
}

/// Workload 4 (Explorer 1): 50 immediate-death (64-bp) + 20 medium-death (128-bp) isolated
/// 32-mer noise hits + 1 real 300-bp HSP.
fn write_k32_fast_reject_fixture(dir: &Path) -> (PathBuf, PathBuf) {
    let mut q_seq = String::with_capacity(20_000);
    let mut r_seq = String::with_capacity(20_000);

    // 50 immediate-death 32-mer noise hits (try_fast_short_reject_64)
    for i in 0..50u64 {
        let blk = pseudo_dna(0xF00D_0000 + i, 232);
        q_seq.push_str(&blk);
        r_seq.push_str(&all_mismatch(&blk[..100]));
        r_seq.push_str(&blk[100..132]);
        r_seq.push_str(&all_mismatch(&blk[132..]));
    }

    // 20 medium-death 32-mer noise hits (try_fast_short_reject_128, min_len = 200 >= 160)
    for i in 0..20u64 {
        let blk = pseudo_dna(0xF00D_1000 + i, 292);
        q_seq.push_str(&blk);
        r_seq.push_str(&all_mismatch(&blk[..50]));
        r_seq.push_str(&pattern_mutate(&blk[50..130], 1, 3));
        r_seq.push_str(&blk[130..162]);
        r_seq.push_str(&pattern_mutate(&blk[162..242], 1, 3));
        r_seq.push_str(&all_mismatch(&blk[242..]));
    }

    // 1 real 300-bp accepted HSP
    let pad_l = pseudo_dna(0xF00D_9001, 100);
    let hsp = pseudo_dna(0xF00D_9002, 300);
    let pad_r = pseudo_dna(0xF00D_9003, 100);
    q_seq.push_str(&pad_l);
    q_seq.push_str(&hsp);
    q_seq.push_str(&pad_r);
    r_seq.push_str(&all_mismatch(&pad_l));
    r_seq.push_str(&hsp);
    r_seq.push_str(&all_mismatch(&pad_r));

    let q_path = dir.join("q.fa");
    let r_path = dir.join("r.fa");
    write_fasta(&q_path, &[("q", &q_seq)]);
    write_fasta(&r_path, &[("r", &r_seq)]);
    (q_path, r_path)
}

// ============================================================================
// TIER 1 — CPU Feature Coverage (12 Tests)
// ============================================================================

/// T1-CPU-1 (W3-F10): Default `--seeding` (omitted) MUST select `Auto` -> `Minimizer`
/// for contiguous `k >= 18` in both `compare` and `all-vs-all`, with 100% SHA-256 parity
/// against explicit `--seeding exact`.
#[test]
fn tier1_cpu_default_seeding_auto_selects_minimizer_for_k_ge_18() {
    let tmp = tempfile::tempdir().unwrap();

    // Synthetic pair with forward + reverse HSPs
    let f1 = pseudo_dna(0x1101, 220);
    let f2 = pseudo_dna(0x1102, 200);
    let q_seq = format!(
        "{}{}{}{}{}",
        pseudo_dna(0x1110, 400),
        f1,
        pseudo_dna(0x1111, 400),
        f2,
        pseudo_dna(0x1112, 400)
    );
    let r_seq = format!(
        "{}{}{}{}{}",
        pseudo_dna(0x1120, 450),
        f1,
        pseudo_dna(0x1121, 450),
        revcomp(&f2),
        pseudo_dna(0x1122, 450)
    );
    let q = tmp.path().join("q.fa");
    let r = tmp.path().join("r.fa");
    write_fasta(&q, &[("q", &q_seq)]);
    write_fasta(&r, &[("r", &r_seq)]);

    for k in ["18", "24", "32"] {
        let out_def = tmp.path().join(format!("def_k{k}"));
        let out_ex = tmp.path().join(format!("ex_k{k}"));

        let res_def = run_gecko_with_env(
            &[
                "compare",
                q.to_str().unwrap(),
                r.to_str().unwrap(),
                "-k",
                k,
                "-l",
                "100",
                "-s",
                "75",
                "-f",
                "frags",
                "-o",
                out_def.to_str().unwrap(),
                "--backend",
                "cpu",
            ],
            &[("GECKO_PROFILE", Some("1"))],
        );
        assert_success(&res_def, &format!("default seeding compare k={k}"));

        let stderr_def = String::from_utf8_lossy(&res_def.stderr);
        let lines_def = parse_profile_lines(&stderr_def);
        assert!(
            used_minimizer(&lines_def) && !used_exact_stream(&lines_def),
            "default seeding (omitted) for k={k} must auto-select minimizer (not exact stream).\nstderr:\n{stderr_def}"
        );

        let res_ex = run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            k,
            "-l",
            "100",
            "-s",
            "75",
            "--seeding",
            "exact",
            "-f",
            "frags",
            "-o",
            out_ex.to_str().unwrap(),
            "--backend",
            "cpu",
        ]);
        assert_success(&res_ex, &format!("exact seeding compare k={k}"));

        let parsed = parse_frags_file(&out_def.join("q-r.frags"));
        assert_eq!(parsed.fragments.len(), 2, "expected 2 HSPs at k={k}");
        assert_frags_sha256_identical(
            &out_ex.join("q-r.frags"),
            &out_def.join("q-r.frags"),
            &format!("compare default vs exact k={k}"),
        );
    }

    // Also verify `all-vs-all` without `--seeding` at k=32
    let gdir = tmp.path().join("genomes");
    fs::create_dir_all(&gdir).unwrap();
    let h12 = pseudo_dna(0x1151, 200);
    let h13 = pseudo_dna(0x1152, 200);
    let h23 = pseudo_dna(0x1153, 200);
    write_fasta(
        &gdir.join("g1.fa"),
        &[(
            "g1",
            &format!(
                "{}{}{}{}{}",
                pseudo_dna(0x1160, 300),
                h12,
                pseudo_dna(0x1161, 300),
                h13,
                pseudo_dna(0x1162, 300)
            ),
        )],
    );
    write_fasta(
        &gdir.join("g2.fa"),
        &[(
            "g2",
            &format!(
                "{}{}{}{}{}",
                pseudo_dna(0x1170, 400),
                revcomp(&h12),
                pseudo_dna(0x1171, 400),
                h23,
                pseudo_dna(0x1172, 400)
            ),
        )],
    );
    write_fasta(
        &gdir.join("g3.fa"),
        &[(
            "g3",
            &format!(
                "{}{}{}{}{}",
                pseudo_dna(0x1180, 500),
                h13,
                pseudo_dna(0x1181, 500),
                revcomp(&h23),
                pseudo_dna(0x1182, 500)
            ),
        )],
    );

    let ava_def = tmp.path().join("ava_def");
    let ava_ex = tmp.path().join("ava_ex");
    let res_ava_def = run_gecko_with_env(
        &[
            "all-vs-all",
            gdir.to_str().unwrap(),
            "-k",
            "32",
            "-l",
            "100",
            "-s",
            "75",
            "-f",
            "frags",
            "-o",
            ava_def.to_str().unwrap(),
            "--backend",
            "cpu",
        ],
        &[("GECKO_PROFILE", Some("1"))],
    );
    assert_success(&res_ava_def, "all-vs-all default seeding k=32");
    let stderr_ava = String::from_utf8_lossy(&res_ava_def.stderr);
    let lines_ava = parse_profile_lines(&stderr_ava);
    assert!(
        lines_ava.iter().any(|l| {
            l.stage == "all_vs_all_index"
                && l.kv.get("mode").map(String::as_str) == Some("minimizer")
        }),
        "all-vs-all default seeding at k=32 must emit all_vs_all_index mode=minimizer.\nstderr:\n{stderr_ava}"
    );
    assert!(
        lines_ava.iter().any(|l| l.stage == "minimizer_query"),
        "all-vs-all default seeding at k=32 must emit minimizer_query.\nstderr:\n{stderr_ava}"
    );

    assert_success(
        &run_gecko(&[
            "all-vs-all",
            gdir.to_str().unwrap(),
            "-k",
            "32",
            "-l",
            "100",
            "-s",
            "75",
            "--seeding",
            "exact",
            "-f",
            "frags",
            "-o",
            ava_ex.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "all-vs-all exact seeding k=32",
    );

    for pair in ["g1-g2.frags", "g1-g3.frags", "g2-g3.frags"] {
        let def_frags = ava_def.join(pair);
        assert_frags_sha256_identical(
            &ava_ex.join(pair),
            &def_frags,
            &format!("all-vs-all default vs exact {pair}"),
        );
        assert!(
            !parse_frags_file(&def_frags).fragments.is_empty(),
            "all-vs-all default seeding {pair} must emit at least one fragment"
        );
    }
}

/// T1-CPU-2 (Baseline): `--seeding auto` falls back cleanly to exact seeding for ineligible
/// modes (`k=17 < 18`, `--spaced-seed w14s20`, `--two-hit`, and `--cache-dir`).
#[test]
fn tier1_cpu_seeding_auto_fallbacks_short_k_spaced_two_hit_cache_dir() {
    let tmp = tempfile::tempdir().unwrap();
    let hsp_f = pseudo_dna(0x1201, 220);
    let hsp_r = pseudo_dna(0x1202, 200);
    let q_seq = format!(
        "{}{}{}{}{}",
        pseudo_dna(0x1210, 300),
        hsp_f,
        pseudo_dna(0x1211, 300),
        hsp_r,
        pseudo_dna(0x1212, 300)
    );
    let r_seq = format!(
        "{}{}{}{}{}",
        pseudo_dna(0x1220, 320),
        hsp_f,
        pseudo_dna(0x1221, 320),
        revcomp(&hsp_r),
        pseudo_dna(0x1222, 320)
    );
    let q = tmp.path().join("q.fa");
    let r = tmp.path().join("r.fa");
    write_fasta(&q, &[("q", &q_seq)]);
    write_fasta(&r, &[("r", &r_seq)]);

    // (a) k = 17 (< 18)
    // (b) --spaced-seed w14s20
    // (c) -k 32 --two-hit (placed at end after positionals)
    let cases: [(&str, Vec<&str>); 3] = [
        ("k17", vec!["-k", "17", "-l", "100", "-s", "75"]),
        (
            "spaced",
            vec!["--spaced-seed", "w14s20", "-l", "100", "-s", "75"],
        ),
        ("two_hit", vec!["-k", "32", "-l", "100", "-s", "75"]),
    ];

    for (label, extra) in cases {
        let out_auto = tmp.path().join(format!("auto_{label}"));
        let out_ex = tmp.path().join(format!("exact_{label}"));

        let mut args_auto = vec!["compare", q.to_str().unwrap(), r.to_str().unwrap()];
        args_auto.extend_from_slice(&extra);
        args_auto.extend_from_slice(&[
            "--seeding",
            "auto",
            "-f",
            "frags",
            "-o",
            out_auto.to_str().unwrap(),
            "--backend",
            "cpu",
        ]);
        if label == "two_hit" {
            args_auto.push("--two-hit");
        }
        let res_auto = run_gecko_with_env(&args_auto, &[("GECKO_PROFILE", Some("1"))]);
        assert_success(&res_auto, &format!("auto fallback {label}"));
        let stderr_auto = String::from_utf8_lossy(&res_auto.stderr);
        let lines_auto = parse_profile_lines(&stderr_auto);
        assert!(
            used_exact_stream(&lines_auto) && !used_minimizer(&lines_auto),
            "auto fallback ({label}) must use exact stream, not minimizer.\nstderr:\n{stderr_auto}"
        );

        let mut args_ex = vec!["compare", q.to_str().unwrap(), r.to_str().unwrap()];
        args_ex.extend_from_slice(&extra);
        args_ex.extend_from_slice(&[
            "--seeding",
            "exact",
            "-f",
            "frags",
            "-o",
            out_ex.to_str().unwrap(),
            "--backend",
            "cpu",
        ]);
        if label == "two_hit" {
            args_ex.push("--two-hit");
        }
        assert_success(&run_gecko(&args_ex), &format!("exact control {label}"));

        let parsed = parse_frags_file(&out_auto.join("q-r.frags"));
        assert_eq!(parsed.fragments.len(), 2, "expected 2 HSPs for {label}");
        assert_frags_sha256_identical(
            &out_ex.join("q-r.frags"),
            &out_auto.join("q-r.frags"),
            &format!("auto fallback vs exact ({label})"),
        );
    }

    // (d) -k 32 --cache-dir <dir>
    let cache_dir = tmp.path().join("gidx_cache");
    let out_cache_auto = tmp.path().join("auto_cache");
    let out_cache_ex = tmp.path().join("exact_cache");
    let res_cache = run_gecko_with_env(
        &[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "32",
            "-l",
            "100",
            "-s",
            "75",
            "--seeding",
            "auto",
            "--cache-dir",
            cache_dir.to_str().unwrap(),
            "-f",
            "frags",
            "-o",
            out_cache_auto.to_str().unwrap(),
            "--backend",
            "cpu",
        ],
        &[("GECKO_PROFILE", Some("1"))],
    );
    assert_success(&res_cache, "auto fallback cache-dir");
    let stderr_cache = String::from_utf8_lossy(&res_cache.stderr);
    let lines_cache = parse_profile_lines(&stderr_cache);
    assert!(
        !used_minimizer(&lines_cache),
        "--cache-dir must disable minimizer under --seeding auto.\nstderr:\n{stderr_cache}"
    );
    let gidx_count = fs::read_dir(&cache_dir)
        .expect("read cache_dir")
        .filter_map(Result::ok)
        .filter(|e| e.path().extension().and_then(|s| s.to_str()) == Some("gidx"))
        .count();
    assert!(
        gidx_count >= 1,
        "expected at least one .gidx file in {}",
        cache_dir.display()
    );

    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "32",
            "-l",
            "100",
            "-s",
            "75",
            "--seeding",
            "exact",
            "-f",
            "frags",
            "-o",
            out_cache_ex.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "exact control for cache-dir",
    );
    let auto_frags = out_cache_auto.join("q-r.frags");
    assert_frags_sha256_identical(
        &out_cache_ex.join("q-r.frags"),
        &auto_frags,
        "auto --cache-dir vs exact",
    );
    assert!(
        !parse_frags_file(&auto_frags).fragments.is_empty(),
        "auto --cache-dir must emit at least one fragment"
    );
}

/// T1-CPU-3 (W3-F1, W3-F9): Valley-trimming 74 bp MEM lossless recovery (`safe_b == 0`).
#[test]
fn tier1_cpu_valley_trimming_74bp_mem_lossless_recovery() {
    let tmp = tempfile::tempdir().unwrap();
    let qx = pseudo_dna(0xCAFE_0074, 2306);
    let rx = format!(
        "{}{}{}{}{}{}{}{}",
        all_mismatch(&qx[0..200]),
        pattern_mutate(&qx[200..1192], 15, 1),
        pattern_mutate(&qx[1192..1200], 7, 1),
        &qx[1200..1274],
        all_mismatch(&qx[1274..1309]),
        pattern_mutate(&qx[1309..2101], 5, 3),
        &qx[2101..2106],
        all_mismatch(&qx[2106..2306]),
    );
    assert_eq!(rx.len(), 2306);

    let q = tmp.path().join("q.fa");
    let r = tmp.path().join("r.fa");
    write_fasta(&q, &[("q", &qx)]);
    write_fasta(&r, &[("r", &rx)]);

    let expected_sim = (3792.0_f64 * 100.0 / (1074.0 * 4.0)) as f32;
    let modes: [(&str, Option<&str>); 3] = [
        ("exact", Some("exact")),
        ("minimizer", Some("minimizer")),
        ("default", None),
    ];
    let mut frags_paths = Vec::new();

    for (label, seeding) in modes {
        let out = tmp.path().join(label);
        let mut args = vec![
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "32",
            "-l",
            "1000",
            "-s",
            "75",
            "--freq-limit",
            "20",
            "-f",
            "frags",
            "-o",
            out.to_str().unwrap(),
            "--backend",
            "cpu",
        ];
        if let Some(m) = seeding {
            args.push("--seeding");
            args.push(m);
        }
        assert_success(
            &run_gecko(&args),
            &format!("74bp valley-trimming ({label})"),
        );
        let p = out.join("q-r.frags");
        let parsed = parse_frags_file(&p);
        assert_eq!(
            parsed.fragments.len(),
            1,
            "74bp valley-trimming ({label}): expected 1 recovered HSP, got {}",
            parsed.fragments.len()
        );
        let f = &parsed.fragments[0];
        assert_eq!(f.strand, 'f');
        assert_eq!(f.diag, 0);
        assert_eq!(f.x_start, 200);
        assert_eq!(f.x_end, 1273);
        assert_eq!(f.y_start, 200);
        assert_eq!(f.y_end, 1273);
        assert_eq!(f.length, 1074);
        assert_eq!(f.ident, 1011);
        assert_eq!(f.score, 3792);
        assert_eq!(
            f.similarity.to_bits(),
            expected_sim.to_bits(),
            "similarity bit mismatch: got {}, expected {}",
            f.similarity,
            expected_sim
        );
        frags_paths.push((label, p));
    }

    assert_frags_sha256_identical(&frags_paths[0].1, &frags_paths[1].1, "74bp exact vs min");
    assert_frags_sha256_identical(&frags_paths[0].1, &frags_paths[2].1, "74bp exact vs def");
}

/// T1-CPU-4 (W3-F1, W3-F9): Valley-trimming 44 bp short MEM (`span = 12 <= 16 = wl/2`).
#[test]
fn tier1_cpu_valley_trimming_44bp_short_mem_lossless_recovery() {
    let tmp = tempfile::tempdir().unwrap();
    let qx = pseudo_dna(0xCAFE_0044, 2276);
    let rx = format!(
        "{}{}{}{}{}{}{}{}",
        all_mismatch(&qx[0..200]),
        pattern_mutate(&qx[200..1192], 15, 1),
        pattern_mutate(&qx[1192..1200], 7, 1),
        &qx[1200..1244],
        all_mismatch(&qx[1244..1279]),
        pattern_mutate(&qx[1279..2071], 5, 3),
        &qx[2071..2076],
        all_mismatch(&qx[2076..2276]),
    );
    assert_eq!(rx.len(), 2276);

    let q = tmp.path().join("q.fa");
    let r = tmp.path().join("r.fa");
    write_fasta(&q, &[("q", &qx)]);
    write_fasta(&r, &[("r", &rx)]);

    let expected_sim = (3672.0_f64 * 100.0 / (1044.0 * 4.0)) as f32;
    let modes: [(&str, Option<&str>); 3] = [
        ("exact", Some("exact")),
        ("minimizer", Some("minimizer")),
        ("default", None),
    ];
    let mut frags_paths = Vec::new();

    for (label, seeding) in modes {
        let out = tmp.path().join(label);
        let mut args = vec![
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "32",
            "-l",
            "1000",
            "-s",
            "75",
            "--freq-limit",
            "20",
            "-f",
            "frags",
            "-o",
            out.to_str().unwrap(),
            "--backend",
            "cpu",
        ];
        if let Some(m) = seeding {
            args.push("--seeding");
            args.push(m);
        }
        assert_success(
            &run_gecko(&args),
            &format!("44bp valley-trimming ({label})"),
        );
        let p = out.join("q-r.frags");
        let parsed = parse_frags_file(&p);
        assert_eq!(
            parsed.fragments.len(),
            1,
            "44bp short-MEM valley-trimming ({label}): expected 1 recovered HSP, got {}",
            parsed.fragments.len()
        );
        let f = &parsed.fragments[0];
        assert_eq!(f.strand, 'f');
        assert_eq!(f.diag, 0);
        assert_eq!(f.x_start, 200);
        assert_eq!(f.x_end, 1243);
        assert_eq!(f.y_start, 200);
        assert_eq!(f.y_end, 1243);
        assert_eq!(f.length, 1044);
        assert_eq!(f.ident, 981);
        assert_eq!(f.score, 3672);
        assert_eq!(f.similarity.to_bits(), expected_sim.to_bits());
        frags_paths.push((label, p));
    }

    assert_frags_sha256_identical(&frags_paths[0].1, &frags_paths[1].1, "44bp exact vs min");
    assert_frags_sha256_identical(&frags_paths[0].1, &frags_paths[2].1, "44bp exact vs def");
}

/// T1-CPU-5a (W3-F1): Lossless `update_noise` & `SEQ_SEP` barrier under `--seeding exact`.
#[test]
fn tier1_cpu_update_noise_seq_sep_barrier_exact() {
    let tmp = tempfile::tempdir().unwrap();
    let (q, r) = write_seq_sep_update_noise_fixture(tmp.path());
    let out = tmp.path().join("exact");
    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "32",
            "-l",
            "200",
            "-s",
            "75",
            "--seeding",
            "exact",
            "-f",
            "frags",
            "-o",
            out.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "SEQ_SEP update_noise barrier (exact)",
    );

    let parsed = parse_frags_file(&out.join("q-r.frags"));
    assert_eq!(
        parsed.fragments.len(),
        2,
        "exact mode must recover both rec1 (after long bwd-flank rejected HSP) and rec3 (after short noise hit before SEQ_SEP), got {:?}",
        parsed.fragments
    );
    let mut coords: Vec<(u64, u64, u64, u64)> = parsed
        .fragments
        .iter()
        .map(|f| (f.seq_x, f.seq_y, f.x_start, f.x_end))
        .collect();
    coords.sort();
    assert_eq!(
        coords,
        vec![(1, 1, 400, 681), (3, 3, 714, 995)],
        "expected global .frags coordinates for rec1 and rec3"
    );
    for f in &parsed.fragments {
        assert_eq!(f.strand, 'f');
        assert_eq!(f.length, 282);
        assert_eq!(f.ident, 267);
        assert_eq!(f.score, 1008);
    }
}

/// T1-CPU-5b (W3-F9): Lossless `update_noise` & `SEQ_SEP` barrier under `--seeding minimizer`
/// and `--seeding auto`, with 100% SHA-256 parity vs `--seeding exact`.
#[test]
fn tier1_cpu_update_noise_seq_sep_barrier_minimizer() {
    let tmp = tempfile::tempdir().unwrap();
    let (q, r) = write_seq_sep_update_noise_fixture(tmp.path());

    for mode in ["exact", "minimizer", "auto"] {
        let out = tmp.path().join(mode);
        assert_success(
            &run_gecko(&[
                "compare",
                q.to_str().unwrap(),
                r.to_str().unwrap(),
                "-k",
                "32",
                "-l",
                "200",
                "-s",
                "75",
                "--seeding",
                mode,
                "-f",
                "frags",
                "-o",
                out.to_str().unwrap(),
                "--backend",
                "cpu",
            ]),
            &format!("SEQ_SEP update_noise barrier ({mode})"),
        );
    }

    let min_frags = parse_frags_file(&tmp.path().join("minimizer/q-r.frags"));
    assert_eq!(
        min_frags.fragments.len(),
        2,
        "minimizer mode must recover both rec1 and rec3 HSPs across SEQ_SEP, got {:?}",
        min_frags.fragments
    );
    for f in &min_frags.fragments {
        assert_eq!(f.length, 282);
        assert_eq!(f.ident, 267);
        assert_eq!(f.score, 1008);
    }
    assert_frags_sha256_identical(
        &tmp.path().join("exact/q-r.frags"),
        &tmp.path().join("minimizer/q-r.frags"),
        "SEQ_SEP barrier exact vs minimizer",
    );
    assert_frags_sha256_identical(
        &tmp.path().join("exact/q-r.frags"),
        &tmp.path().join("auto/q-r.frags"),
        "SEQ_SEP barrier exact vs auto",
    );
}

/// T1-CPU-6a (W3-F2): `k = 32` fast short reject (`0..=128` + 128-bp scalar reject) is active
/// on `--seeding exact` (`fast_reject >= 50` in `stream_query`).
#[test]
fn tier1_cpu_k32_fast_reject_active_exact() {
    let tmp = tempfile::tempdir().unwrap();
    let (q, r) = write_k32_fast_reject_fixture(tmp.path());
    let out = tmp.path().join("exact");
    let res = run_gecko_with_env(
        &[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "32",
            "-l",
            "200",
            "-s",
            "75",
            "--seeding",
            "exact",
            "-f",
            "frags",
            "-o",
            out.to_str().unwrap(),
            "--backend",
            "cpu",
        ],
        &[("GECKO_PROFILE", Some("1"))],
    );
    assert_success(&res, "k=32 fast_reject (exact)");

    let parsed = parse_frags_file(&out.join("q-r.frags"));
    assert_eq!(parsed.fragments.len(), 1, "expected 1 true 300-bp HSP");
    assert_eq!(parsed.fragments[0].length, 300);
    assert_eq!(parsed.fragments[0].score, 1200);

    let stderr = String::from_utf8_lossy(&res.stderr);
    let lines = parse_profile_lines(&stderr);
    let fast_rej = sum_u64_key(&lines, "stream_query", "fast_reject");
    assert!(
        fast_rej >= Some(50),
        "k=32 exact stream_query must activate fast_reject (>= 50), got {fast_rej:?}.\nstderr:\n{stderr}"
    );
}

/// T1-CPU-6b (W3-F9): `k = 32` fast short reject is active on `--seeding minimizer`
/// (`fast_reject >= 50` in `minimizer_query`) with 100% SHA-256 parity vs `exact`.
#[test]
fn tier1_cpu_k32_fast_reject_active_minimizer() {
    let tmp = tempfile::tempdir().unwrap();
    let (q, r) = write_k32_fast_reject_fixture(tmp.path());
    let out_ex = tmp.path().join("exact");
    let out_min = tmp.path().join("minimizer");

    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "32",
            "-l",
            "200",
            "-s",
            "75",
            "--seeding",
            "exact",
            "-f",
            "frags",
            "-o",
            out_ex.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "k=32 fast_reject control (exact)",
    );

    let res_min = run_gecko_with_env(
        &[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "32",
            "-l",
            "200",
            "-s",
            "75",
            "--seeding",
            "minimizer",
            "-f",
            "frags",
            "-o",
            out_min.to_str().unwrap(),
            "--backend",
            "cpu",
        ],
        &[("GECKO_PROFILE", Some("1"))],
    );
    assert_success(&res_min, "k=32 fast_reject (minimizer)");

    let min_frags = out_min.join("q-r.frags");
    assert_frags_sha256_identical(
        &out_ex.join("q-r.frags"),
        &min_frags,
        "k=32 fast_reject exact vs minimizer",
    );
    assert_eq!(
        parse_frags_file(&min_frags).fragments.len(),
        1,
        "k=32 fast_reject (minimizer) must emit the 1 true 300-bp HSP"
    );

    let stderr = String::from_utf8_lossy(&res_min.stderr);
    let lines = parse_profile_lines(&stderr);
    let fast_rej = sum_u64_key(&lines, "minimizer_query", "fast_reject");
    assert!(
        fast_rej >= Some(50),
        "k=32 minimizer_query must emit fast_reject >= 50 (None means field missing from minimizer_query telemetry), got {fast_rej:?}.\nstderr:\n{stderr}"
    );
}

/// T1-CPU-7 (W3-F6..F9): Hash-agnostic high-frequency repeat minimizer overflow routing (`M_CAP`).
/// Verifies 100% SHA-256 parity, exact `freq_limit=20` kept/dropped + reverse overflow HSPs,
/// and `exact_k_reject < 80_000` (vs `517,516` on `b3578ac`).
#[test]
fn tier1_cpu_minimizer_high_freq_repeat_overflow_parity_and_bounded_rejects() {
    let tmp = tempfile::tempdir().unwrap();
    let cores: [String; 3] = [
        pseudo_dna(0xABCD_1000, 27),
        pseudo_dna(0xABCD_1001, 27),
        pseudo_dna(0xABCD_1002, 27),
    ];

    let mut r_seq = String::with_capacity(96_000);
    let mut q_seq = String::with_capacity(96_000);

    // 3 * 500 disjoint repeat copies (ref: 0..500, query: 500..1000) -> 0 exact 32-mer hits.
    // Use different spacer lengths (20 bp in ref vs 23 bp in query) so the 27-bp cores never
    // align along a shared periodic diagonal corridor.
    for (c, core) in cores.iter().enumerate() {
        for i in 0..500usize {
            let fl_r = encode_base4_5mer(i);
            r_seq.push_str(&fl_r);
            r_seq.push_str(core);
            r_seq.push_str(&fl_r);
            r_seq.push_str(&pseudo_dna(0xABCD_2000 + (c * 1000 + i) as u64, 20));

            let fl_q = encode_base4_5mer(500 + i);
            q_seq.push_str(&fl_q);
            q_seq.push_str(core);
            q_seq.push_str(&fl_q);
            q_seq.push_str(&pseudo_dna(0xABCD_6000 + (c * 1000 + i) as u64, 23));
        }
    }

    // Helper: 237-bp matching block whose ONLY exact 32-mers lie strictly inside `fl + core0 + fl`
    // (mismatch immediately before `fl` at left[99] and immediately after `fl` at right[0]).
    let make_probe_pair = |seed: u64, fl: &str| -> (String, String) {
        let left_q = pseudo_dna(seed, 100);
        let right_q = pseudo_dna(seed + 1, 100);
        let left_r = format!(
            "{}{}",
            pattern_mutate(&left_q[..99], 15, 1),
            all_mismatch(&left_q[99..100])
        );
        let right_r = format!(
            "{}{}",
            all_mismatch(&right_q[..1]),
            pattern_mutate(&right_q[1..], 15, 1)
        );
        let mid = format!("{fl}{}{fl}", cores[0]);
        (
            format!("{left_q}{mid}{right_q}"),
            format!("{left_r}{mid}{right_r}"),
        )
    };

    // Probe 1: cnt == 20 <= freq_limit (KEPT, 1 'f' frag)
    let fl_kept = encode_base4_5mer(1000);
    let (p1_q, p1_r) = make_probe_pair(0xABCE_0100, &fl_kept);
    let g1l = pseudo_dna(0xABCE_0110, 60);
    let g1r = pseudo_dna(0xABCE_0111, 60);
    q_seq.push_str(&g1l);
    q_seq.push_str(&p1_q);
    q_seq.push_str(&g1r);
    r_seq.push_str(&all_mismatch(&g1l));
    r_seq.push_str(&p1_r);
    r_seq.push_str(&all_mismatch(&g1r));
    for extra in 0..19u64 {
        r_seq.push_str(&pseudo_dna(0xABCE_1000 + extra, 20));
        r_seq.push_str(&fl_kept);
        r_seq.push_str(&cores[0]);
        r_seq.push_str(&fl_kept);
        r_seq.push_str(&pseudo_dna(0xABCE_1100 + extra, 20));
    }

    // Probe 2: cnt == 21 > freq_limit (DROPPED, 0 frags)
    let fl_drop = encode_base4_5mer(1001);
    let (p2_q, p2_r) = make_probe_pair(0xABCE_0200, &fl_drop);
    let g2l = pseudo_dna(0xABCE_0210, 60);
    let g2r = pseudo_dna(0xABCE_0211, 60);
    q_seq.push_str(&g2l);
    q_seq.push_str(&p2_q);
    q_seq.push_str(&g2r);
    r_seq.push_str(&all_mismatch(&g2l));
    r_seq.push_str(&p2_r);
    r_seq.push_str(&all_mismatch(&g2r));
    for extra in 0..20u64 {
        r_seq.push_str(&pseudo_dna(0xABCE_2000 + extra, 20));
        r_seq.push_str(&fl_drop);
        r_seq.push_str(&cores[0]);
        r_seq.push_str(&fl_drop);
        r_seq.push_str(&pseudo_dna(0xABCE_2100 + extra, 20));
    }

    // Probe 3: cnt == 1 reverse-strand overflow HSP (KEPT, 1 'r' frag)
    let fl_rev = encode_base4_5mer(1002);
    let (p3_q_fwd, p3_r) = make_probe_pair(0xABCE_0300, &fl_rev);
    let g3l = pseudo_dna(0xABCE_0310, 60);
    let g3r = pseudo_dna(0xABCE_0311, 60);
    q_seq.push_str(&g3l);
    q_seq.push_str(&revcomp(&p3_q_fwd));
    q_seq.push_str(&g3r);
    r_seq.push_str(&all_mismatch(&revcomp(&g3r)));
    r_seq.push_str(&p3_r);
    r_seq.push_str(&all_mismatch(&revcomp(&g3l)));

    let q = tmp.path().join("q.fa");
    let r = tmp.path().join("r.fa");
    write_fasta(&q, &[("q", &q_seq)]);
    write_fasta(&r, &[("r", &r_seq)]);

    let out_ex = tmp.path().join("exact");
    let out_min = tmp.path().join("minimizer");
    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "32",
            "-l",
            "150",
            "-s",
            "75",
            "--freq-limit",
            "20",
            "--no-soft-mask",
            "--seeding",
            "exact",
            "-f",
            "frags",
            "-o",
            out_ex.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "HF repeat overflow (exact)",
    );

    let res_min = run_gecko_with_env(
        &[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "32",
            "-l",
            "150",
            "-s",
            "75",
            "--freq-limit",
            "20",
            "--no-soft-mask",
            "--seeding",
            "minimizer",
            "-f",
            "frags",
            "-o",
            out_min.to_str().unwrap(),
            "--backend",
            "cpu",
        ],
        &[("GECKO_PROFILE", Some("1"))],
    );
    assert_success(&res_min, "HF repeat overflow (minimizer)");

    assert_frags_sha256_identical(
        &out_ex.join("q-r.frags"),
        &out_min.join("q-r.frags"),
        "HF repeat overflow exact vs minimizer",
    );

    let parsed = parse_frags_file(&out_min.join("q-r.frags"));
    let fwd = parsed.fragments.iter().filter(|f| f.strand == 'f').count();
    let rev = parsed.fragments.iter().filter(|f| f.strand == 'r').count();
    assert_eq!(
        (parsed.fragments.len(), fwd, rev),
        (2, 1, 1),
        "expected exactly 2 frags (1 fwd cnt==20, 1 rev overflow; 0 from cnt==21), got {:?}",
        parsed.fragments
    );

    let stderr = String::from_utf8_lossy(&res_min.stderr);
    let lines = parse_profile_lines(&stderr);
    let exact_k_rej = sum_u64_key(&lines, "minimizer_query", "exact_k_reject")
        .expect("minimizer_query must report exact_k_reject");
    assert!(
        exact_k_rej < 80_000,
        "HF repeat overflow routing must bound exact_k_reject < 80_000, got {exact_k_rej}.\nstderr:\n{stderr}"
    );
}

/// T1-CPU-8 (W3-F10): `all-vs-all` supports `--no-soft-mask` and `--forward-only` with 100%
/// SHA-256 parity against pairwise `compare`.
#[test]
fn tier1_cpu_all_vs_all_no_soft_mask_and_forward_only_parity() {
    let tmp = tempfile::tempdir().unwrap();
    let gdir = tmp.path().join("genomes");
    fs::create_dir_all(&gdir).unwrap();

    let h_up_f = pseudo_dna(0x1801, 200);
    let h_up_r = pseudo_dna(0x1802, 200);
    let h_lo_f = pseudo_dna(0x1803, 200).to_ascii_lowercase();
    let h_lo_r = pseudo_dna(0x1804, 200).to_ascii_lowercase();

    let g0_seq = format!(
        "{}{}{}{}{}{}{}{}{}",
        pseudo_dna(0x1810, 200),
        h_up_f,
        pseudo_dna(0x1811, 200),
        h_up_r,
        pseudo_dna(0x1812, 200),
        h_lo_f,
        pseudo_dna(0x1813, 200),
        h_lo_r,
        pseudo_dna(0x1814, 200),
    );
    let g1_seq = format!(
        "{}{}{}{}{}{}{}{}{}",
        pseudo_dna(0x1820, 300),
        h_up_f,
        pseudo_dna(0x1821, 300),
        revcomp(&h_up_r),
        pseudo_dna(0x1822, 300),
        h_lo_f,
        pseudo_dna(0x1823, 300),
        revcomp(&h_lo_r),
        pseudo_dna(0x1824, 300),
    );
    let g2_seq = format!(
        "{}{}{}{}{}{}{}{}{}",
        pseudo_dna(0x1830, 450),
        h_up_f,
        pseudo_dna(0x1831, 450),
        revcomp(&h_up_r),
        pseudo_dna(0x1832, 450),
        h_lo_f,
        pseudo_dna(0x1833, 450),
        revcomp(&h_lo_r),
        pseudo_dna(0x1834, 450),
    );

    let p0 = gdir.join("g0.fa");
    let p1 = gdir.join("g1.fa");
    let p2 = gdir.join("g2.fa");
    write_fasta(&p0, &[("g0", &g0_seq)]);
    write_fasta(&p1, &[("g1", &g1_seq)]);
    write_fasta(&p2, &[("g2", &g2_seq)]);

    let pairs: [(&str, &Path, &Path); 3] = [
        ("g0-g1", &p0, &p1),
        ("g0-g2", &p0, &p2),
        ("g1-g2", &p1, &p2),
    ];

    // Soft-masked baseline (for strictly-more-frags check under --no-soft-mask)
    let out_base = tmp.path().join("ava_base");
    assert_success(
        &run_gecko(&[
            "all-vs-all",
            gdir.to_str().unwrap(),
            "-k",
            "32",
            "-l",
            "100",
            "-s",
            "75",
            "-f",
            "frags",
            "-o",
            out_base.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "all-vs-all soft-masked baseline",
    );

    let combos: [(&str, bool, bool); 3] = [
        ("no_soft_mask", true, false),
        ("forward_only", false, true),
        ("both", true, true),
    ];

    for (label, no_soft_mask, fwd_only) in combos {
        for seeding in [None, Some("exact")] {
            let tag = format!("{label}_{}", seeding.unwrap_or("default"));
            let ava_out = tmp.path().join(format!("ava_{tag}"));
            let mut ava_args = vec![
                "all-vs-all",
                gdir.to_str().unwrap(),
                "-k",
                "32",
                "-l",
                "100",
                "-s",
                "75",
                "-f",
                "frags",
                "-o",
                ava_out.to_str().unwrap(),
                "--backend",
                "cpu",
            ];
            if no_soft_mask {
                ava_args.push("--no-soft-mask");
            }
            if fwd_only {
                ava_args.push("--forward-only");
            }
            if let Some(s) = seeding {
                ava_args.push("--seeding");
                ava_args.push(s);
            }
            assert_success(&run_gecko(&ava_args), &format!("all-vs-all {tag}"));

            for (stem, qa, rb) in &pairs {
                let pair_out = tmp.path().join(format!("pair_{tag}_{stem}"));
                let mut cmp_args = vec![
                    "compare",
                    qa.to_str().unwrap(),
                    rb.to_str().unwrap(),
                    "-k",
                    "32",
                    "-l",
                    "100",
                    "-s",
                    "75",
                    "--seeding",
                    "exact",
                    "-f",
                    "frags",
                    "-o",
                    pair_out.to_str().unwrap(),
                    "--backend",
                    "cpu",
                ];
                if no_soft_mask {
                    cmp_args.push("--no-soft-mask");
                }
                if fwd_only {
                    cmp_args.push("--forward-only");
                }
                assert_success(&run_gecko(&cmp_args), &format!("pairwise {tag} {stem}"));

                let ava_frags_path = ava_out.join(format!("{stem}.frags"));
                let pair_frags_path = pair_out.join(format!("{stem}.frags"));
                assert_frags_sha256_identical(
                    &pair_frags_path,
                    &ava_frags_path,
                    &format!("all-vs-all vs pairwise {tag} {stem}"),
                );

                let parsed = parse_frags_file(&ava_frags_path);
                if fwd_only {
                    assert!(
                        !parsed.fragments.is_empty(),
                        "{tag} {stem}: expected at least one forward fragment under --forward-only"
                    );
                    assert!(
                        parsed.fragments.iter().all(|f| f.strand == 'f'),
                        "{tag} {stem}: expected zero 'r' fragments under --forward-only"
                    );
                }
                if no_soft_mask && !fwd_only {
                    let base_n = parse_frags_file(&out_base.join(format!("{stem}.frags")))
                        .fragments
                        .len();
                    assert!(
                        parsed.fragments.len() > base_n,
                        "{tag} {stem}: --no-soft-mask must recover more frags ({}) than soft-masked ({base_n})",
                        parsed.fragments.len()
                    );
                }
            }
        }
    }
}

/// T1-CPU-9 (Baseline + W3-M1 Regression Guard: W3-F3, W3-F4, W3-F5):
/// (a) `k=12` / `k=16` / `w14s20` with & without `--two-hit` (`two_hit_drop > 0`, dual-seed kept,
/// single-seed suppressed); (b) O(1) `'N'`-run bulk decay exact `fscore == 0` vs `fscore == -1`
/// boundary in BOTH forward and backward directions.
#[test]
fn tier1_cpu_exact_diagpair_two_hit_and_n_run_decay_boundary() {
    let tmp = tempfile::tempdir().unwrap();

    // Part (a): DirectCsr/HashCsr/Spaced two-hit dual-seed vs single-seed + noise drop
    let make_single_and_dual = |seed_base: u64, k_len: usize| -> (String, String) {
        let f0 = pseudo_dna(seed_base, 120);
        // Single-seed HSP: seed(k_len) + non-seeding flank (90 bp, mismatch every 4th base)
        let s1 = pseudo_dna(seed_base + 1, k_len);
        let t1_q = pseudo_dna(seed_base + 2, 90);
        let t1_r = format!(
            "{}{}",
            all_mismatch(&t1_q[..1]),
            mutate_every(&t1_q[1..], 4, 3)
        );
        let f1 = pseudo_dna(seed_base + 3, 120);

        // Dual-seed HSP: seed_a(k_len) + 12bp broken gap + seed_b(k_len) + 70bp broken tail
        let sa = pseudo_dna(seed_base + 4, k_len);
        let gap_q = pseudo_dna(seed_base + 5, 12);
        let gap_r = mutate_positions(&gap_q, &[0, 6, 11]);
        let sb = pseudo_dna(seed_base + 6, k_len);
        let t2_q = pseudo_dna(seed_base + 7, 70);
        let t2_r = format!(
            "{}{}",
            all_mismatch(&t2_q[..1]),
            mutate_every(&t2_q[1..], 5, 4)
        );
        let f2 = pseudo_dna(seed_base + 8, 120);

        // 20 isolated single-seed noise hits on DISTINCT diagonals to guarantee `two_hit_drop >= 20`
        let mut q_noise = String::new();
        let mut r_noise = String::new();
        for i in 0..20u64 {
            let blk = pseudo_dna(seed_base + 100 + i, k_len + 40);
            q_noise.push_str(&blk);
            r_noise.push_str(&all_mismatch(&blk[..20]));
            r_noise.push_str(&blk[20..20 + k_len]);
            r_noise.push_str(&all_mismatch(&blk[20 + k_len..]));
            r_noise.push_str(&pseudo_dna(seed_base + 300 + i, 5));
        }

        let q = format!("{f0}{s1}{t1_q}{f1}{sa}{gap_q}{sb}{t2_q}{f2}{q_noise}");
        let r = format!(
            "{}{s1}{t1_r}{}{sa}{gap_r}{sb}{t2_r}{}{r_noise}",
            all_mismatch(&f0),
            all_mismatch(&f1),
            all_mismatch(&f2)
        );
        (q, r)
    };

    let configs: [(&str, usize, Vec<&str>); 3] = [
        ("k10_direct", 10, vec!["-k", "10", "-l", "60", "-s", "45"]),
        ("k16_hash", 16, vec!["-k", "16", "-l", "60", "-s", "45"]),
        (
            "w14s20",
            20,
            vec!["--spaced-seed", "w14s20", "-l", "60", "-s", "45"],
        ),
    ];
    for (idx, (label, k_len, flags)) in configs.iter().enumerate() {
        let (qs, rs) = make_single_and_dual(0x1900 + (idx as u64) * 500, *k_len);
        let q = tmp.path().join(format!("q_{label}.fa"));
        let r = tmp.path().join(format!("r_{label}.fa"));
        write_fasta(&q, &[("q", &qs)]);
        write_fasta(&r, &[("r", &rs)]);

        let out_1 = tmp.path().join(format!("single_{label}"));
        let out_2 = tmp.path().join(format!("two_{label}"));

        let mut args1 = vec!["compare", q.to_str().unwrap(), r.to_str().unwrap()];
        args1.extend_from_slice(flags);
        args1.extend_from_slice(&[
            "-f",
            "frags",
            "-o",
            out_1.to_str().unwrap(),
            "--backend",
            "cpu",
        ]);
        assert_success(&run_gecko(&args1), &format!("single-hit {label}"));

        let mut args2 = vec!["compare", q.to_str().unwrap(), r.to_str().unwrap()];
        args2.extend_from_slice(flags);
        args2.extend_from_slice(&[
            "-f",
            "frags",
            "-o",
            out_2.to_str().unwrap(),
            "--backend",
            "cpu",
            "--two-hit",
        ]);
        let res2 = run_gecko_with_env(&args2, &[("GECKO_PROFILE", Some("1"))]);
        assert_success(&res2, &format!("two-hit {label}"));

        let f1 = parse_frags_file(&out_1.join(format!("q_{label}-r_{label}.frags"))).fragments;
        let f2 = parse_frags_file(&out_2.join(format!("q_{label}-r_{label}.frags"))).fragments;
        assert_eq!(f1.len(), 2, "{label}: single-hit must find both HSPs");
        assert_eq!(
            f2.len(),
            1,
            "{label}: --two-hit must suppress single-seed HSP and keep dual-seed HSP"
        );
        let stderr2 = String::from_utf8_lossy(&res2.stderr);
        let lines2 = parse_profile_lines(&stderr2);
        let dropped = sum_u64_key(&lines2, "stream_query", "two_hit_drop").unwrap_or(0);
        assert!(
            dropped >= 20,
            "{label}: expected two_hit_drop >= 20, got {dropped}.\nstderr:\n{stderr2}"
        );
    }

    // Part (b): O(1) 'N'-run bulk decay exact boundary (forward & backward)
    // 60M (+240) + 240 'N' (-240 -> fscore==0 survives!) + 80bp(15M,1X) -> len=379, ident=135, score=284
    // 60M (+240) + 241 'N' (-241 -> fscore==-1 dies!)   + 80bp(15M,1X) -> 0 frags at -l 350
    let pad_l = pseudo_dna(0x1990, 60);
    let pad_r = pseudo_dna(0x1991, 60);
    let exact60 = pseudo_dna(0x1992, 60);
    let flank80_q = pseudo_dna(0x1993, 80);
    let flank80_fwd_r = pattern_mutate(&flank80_q, 15, 1);
    let flank80_bwd_r = pattern_mutate_rev(&flank80_q, 1, 15);

    for (dir_label, is_bwd) in [("fwd", false), ("bwd", true)] {
        for (n_cnt, expect_hit) in [(240usize, true), (241usize, false)] {
            let n_run = "N".repeat(n_cnt);
            let (qs, rs) = if !is_bwd {
                (
                    format!("{pad_l}{exact60}{n_run}{flank80_q}{pad_r}"),
                    format!(
                        "{}{exact60}{n_run}{flank80_fwd_r}{}",
                        all_mismatch(&pad_l),
                        all_mismatch(&pad_r)
                    ),
                )
            } else {
                (
                    format!("{pad_l}{flank80_q}{n_run}{exact60}{pad_r}"),
                    format!(
                        "{}{flank80_bwd_r}{n_run}{exact60}{}",
                        all_mismatch(&pad_l),
                        all_mismatch(&pad_r)
                    ),
                )
            };
            let q = tmp.path().join(format!("qn_{dir_label}_{n_cnt}.fa"));
            let r = tmp.path().join(format!("rn_{dir_label}_{n_cnt}.fa"));
            write_fasta(&q, &[("q", &qs)]);
            write_fasta(&r, &[("r", &rs)]);

            let out = tmp.path().join(format!("out_{dir_label}_{n_cnt}"));
            assert_success(
                &run_gecko(&[
                    "compare",
                    q.to_str().unwrap(),
                    r.to_str().unwrap(),
                    "-k",
                    "32",
                    "-l",
                    "350",
                    "-s",
                    "15",
                    "--n-penalty",
                    "1",
                    "--seeding",
                    "exact",
                    "-f",
                    "frags",
                    "-o",
                    out.to_str().unwrap(),
                    "--backend",
                    "cpu",
                ]),
                &format!("N-run bulk decay {dir_label} n={n_cnt}"),
            );
            let stem = format!("qn_{dir_label}_{n_cnt}-rn_{dir_label}_{n_cnt}.frags");
            let parsed = parse_frags_file(&out.join(&stem));
            if expect_hit {
                assert_eq!(
                    parsed.fragments.len(),
                    1,
                    "{dir_label} n={n_cnt}: fscore==0 must survive and emit 1 merged frag"
                );
                let f = &parsed.fragments[0];
                assert_eq!(f.length, 379);
                assert_eq!(f.ident, 135);
                assert_eq!(f.score, 284);
            } else {
                assert_eq!(
                    parsed.fragments.len(),
                    0,
                    "{dir_label} n={n_cnt}: fscore==-1 must terminate and emit 0 frags"
                );
            }
        }
    }
}

/// T1-CPU-10 (Baseline + W3-F7/F8 Regression Guard): Sweep `k in {18, 20, 24, 28, 32}` on
/// multi-FASTA containing homopolymer/dinucleotide tie-stress HSPs, palindromic 32-mer HSPs,
/// forward + reverse HSPs, and soft-masked regions under both default and `--no-soft-mask`.
#[test]
fn tier1_cpu_minimizer_k18_to_k32_homopolymer_ties_and_softmask_parity() {
    let tmp = tempfile::tempdir().unwrap();
    let poly_a = "A".repeat(40);
    let dinuc_ac = "AC".repeat(25);
    let half_pal = pseudo_dna(0x1A01, 60);
    let pal120 = format!("{half_pal}{}", revcomp(&half_pal));
    assert_eq!(pal120, revcomp(&pal120));

    let h_tie_f = format!("{poly_a}{}{dinuc_ac}", pseudo_dna(0x1A02, 100));
    let h_rev = pseudo_dna(0x1A03, 180);
    let h_soft = pseudo_dna(0x1A04, 160).to_ascii_lowercase();
    let h_rec2 = pseudo_dna(0x1A05, 160);

    let q0 = format!(
        "{}{}{}{}{}{}{}",
        pseudo_dna(0x1A10, 150),
        h_tie_f,
        pseudo_dna(0x1A11, 150),
        pal120,
        pseudo_dna(0x1A12, 150),
        h_soft,
        pseudo_dna(0x1A13, 150)
    );
    let r0 = format!(
        "{}{}{}{}{}{}{}",
        pseudo_dna(0x1A20, 180),
        h_tie_f,
        pseudo_dna(0x1A21, 180),
        pal120,
        pseudo_dna(0x1A22, 180),
        h_soft,
        pseudo_dna(0x1A23, 180)
    );
    let q1 = format!(
        "{}{}{}{}{}",
        pseudo_dna(0x1A30, 150),
        h_rev,
        pseudo_dna(0x1A31, 150),
        h_rec2,
        pseudo_dna(0x1A32, 150)
    );
    let r1 = format!(
        "{}{}{}{}{}",
        pseudo_dna(0x1A40, 180),
        revcomp(&h_rev),
        pseudo_dna(0x1A41, 180),
        h_rec2,
        pseudo_dna(0x1A42, 180)
    );

    let q = tmp.path().join("q.fa");
    let r = tmp.path().join("r.fa");
    write_fasta(&q, &[("q0", &q0), ("q1", &q1)]);
    write_fasta(&r, &[("r0", &r0), ("r1", &r1)]);

    for k in ["18", "20", "24", "28", "32"] {
        for no_soft in [false, true] {
            let tag = format!("k{k}_{}", if no_soft { "unmasked" } else { "masked" });
            let out_ex = tmp.path().join(format!("ex_{tag}"));
            let out_min = tmp.path().join(format!("min_{tag}"));
            for (mode, out) in [("exact", &out_ex), ("minimizer", &out_min)] {
                let mut args = vec![
                    "compare",
                    q.to_str().unwrap(),
                    r.to_str().unwrap(),
                    "-k",
                    k,
                    "-l",
                    "90",
                    "-s",
                    "75",
                    "--seeding",
                    mode,
                    "-f",
                    "frags",
                    "-o",
                    out.to_str().unwrap(),
                    "--backend",
                    "cpu",
                ];
                if no_soft {
                    args.push("--no-soft-mask");
                }
                assert_success(&run_gecko(&args), &format!("sweep {tag} {mode}"));
            }
            assert_frags_sha256_identical(
                &out_ex.join("q-r.frags"),
                &out_min.join("q-r.frags"),
                &format!("minimizer vs exact sweep {tag}"),
            );
            let parsed = parse_frags_file(&out_min.join("q-r.frags"));
            let min_expected = if no_soft { 6 } else { 5 };
            assert!(
                parsed.fragments.len() >= min_expected,
                "{tag}: expected >= {min_expected} frags, got {}",
                parsed.fragments.len()
            );
        }
    }
}

// ============================================================================
// TIER 1 — CUDA Static SASS/PTX/Build Contracts (6 Tests)
// ============================================================================

/// T1-CUDA-1 (Baseline): `cargo check -p gecko-cli --features cuda` compiles cleanly with 0
/// warnings/errors, sets `has_cuda_kernels`, links `static=gecko_cuda_kernels`, and produces
/// `libgecko_cuda_kernels.a`.
#[test]
fn tier1_w3_cuda_cargo_check_features_cuda_clean() {
    let res = cached_cuda_cargo_check();
    assert!(
        res.success,
        "cargo check -p gecko-cli --features cuda failed (exit {:?}):\nerrors: {:?}\nstderr:\n{}",
        res.exit_code, res.errors, res.stderr
    );
    assert!(
        res.cfgs.iter().any(|c| c == "has_cuda_kernels"),
        "expected has_cuda_kernels cfg, got {:?}",
        res.cfgs
    );
    assert!(
        res.linked_libs
            .iter()
            .any(|l| l == "static=gecko_cuda_kernels"),
        "expected static=gecko_cuda_kernels in linked_libs, got {:?}",
        res.linked_libs
    );
    assert!(
        res.errors.is_empty(),
        "compiler errors in gecko-cuda/gecko-cli: {:?}",
        res.errors
    );
    assert!(
        res.warnings.is_empty(),
        "compiler warnings in gecko-cuda/gecko-cli: {:?}",
        res.warnings
    );
    let out_dir = res.out_dir.as_ref().expect("gecko-cuda out_dir");
    assert!(
        out_dir.join("libgecko_cuda_kernels.a").exists(),
        "missing libgecko_cuda_kernels.a in {}",
        out_dir.display()
    );
}

/// T1-CUDA-2 (Baseline + W3-M4 Regression Guard): EVERY `.o` member in `libgecko_cuda_kernels.a`
/// (including any new `cuda_minimizer.o`) embeds CUBINs for `sm_75`, `sm_80`, `sm_86`, `sm_89`
/// and PTX for `compute_89` (`.sm_89.ptx`).
#[test]
fn tier1_w3_cuda_all_archive_members_multi_arch_cubins_and_ptx() {
    let archive = cuda_archive_path();

    let ar_out = Command::new("ar")
        .args(["t", archive.to_str().unwrap()])
        .output()
        .expect("failed to execute ar t");
    assert!(
        ar_out.status.success(),
        "ar t failed:\n{}",
        String::from_utf8_lossy(&ar_out.stderr)
    );
    let ar_members: Vec<String> = String::from_utf8_lossy(&ar_out.stdout)
        .lines()
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(str::to_string)
        .collect();

    let elf_out = Command::new("cuobjdump")
        .args(["--list-elf", archive.to_str().unwrap()])
        .output()
        .expect("cuobjdump --list-elf");
    assert!(
        elf_out.status.success(),
        "cuobjdump --list-elf failed:\n{}",
        String::from_utf8_lossy(&elf_out.stderr)
    );
    let elf_by_member =
        parse_cuobjdump_archive_members(&String::from_utf8_lossy(&elf_out.stdout), "ELF file");

    let ptx_out = Command::new("cuobjdump")
        .args(["--list-ptx", archive.to_str().unwrap()])
        .output()
        .expect("cuobjdump --list-ptx");
    assert!(
        ptx_out.status.success(),
        "cuobjdump --list-ptx failed:\n{}",
        String::from_utf8_lossy(&ptx_out.stderr)
    );
    let ptx_by_member =
        parse_cuobjdump_archive_members(&String::from_utf8_lossy(&ptx_out.stdout), "PTX file");

    let mut members: BTreeSet<String> = [
        "cuda_index.o",
        "cuda_query.o",
        "cuda_extend.o",
        "cuda_bridge.o",
    ]
    .into_iter()
    .map(str::to_string)
    .collect();
    members.extend(ar_members);
    members.extend(elf_by_member.keys().cloned());
    members.extend(ptx_by_member.keys().cloned());

    let mut missing = Vec::new();
    for member in &members {
        let elfs = elf_by_member.get(member).map(Vec::as_slice).unwrap_or(&[]);
        for arch in CUDA_REQUIRED_ARCHS {
            let suffix = format!(".{arch}.cubin");
            if !elfs.iter().any(|l| l.ends_with(&suffix)) {
                missing.push(format!("{member}: missing {suffix}"));
            }
        }
        let ptxs = ptx_by_member.get(member).map(Vec::as_slice).unwrap_or(&[]);
        if !ptxs.iter().any(|l| l.ends_with(".sm_89.ptx")) {
            missing.push(format!("{member}: missing .sm_89.ptx"));
        }
    }

    assert!(
        missing.is_empty(),
        "CUDA archive multi-arch verification failed across {} members:\n  {}",
        members.len(),
        missing.join("\n  ")
    );
}

/// T1-CUDA-3 (Baseline + W3-M3/M4 Regression Guard): 100% of custom GECKO CUDA kernels across
/// `sm_75`, `sm_80`, `sm_86`, `sm_89` have `STACK:0`, `LOCAL:0`, `0 STL`, and `0 LDL`.
#[test]
fn tier1_w3_cuda_all_custom_kernels_zero_stack_and_no_spills_sm75_sm80_sm86_sm89() {
    // Warm all 4 arch caches in parallel
    std::thread::scope(|s| {
        for arch in CUDA_REQUIRED_ARCHS {
            s.spawn(move || {
                let _ = cached_cuda_resource_usage(arch);
                let _ = cached_cuda_sass_stats(arch);
            });
        }
    });

    let mut violations = Vec::new();
    for arch in CUDA_REQUIRED_ARCHS {
        let res_fns = cached_cuda_resource_usage(arch);
        let sass_fns = cached_cuda_sass_stats(arch);

        let custom_res: Vec<&CudaResourceInfo> = res_fns
            .iter()
            .filter(|f| !is_third_party_cuda_kernel(&f.fn_name))
            .collect();
        let custom_sass: Vec<&CudaSassStats> = sass_fns
            .iter()
            .filter(|f| !is_third_party_cuda_kernel(&f.fn_name))
            .collect();

        assert!(
            custom_res.len() >= 10,
            "{arch}: expected >= 10 custom CUDA kernels, got {}",
            custom_res.len()
        );
        for anchor in [
            "extend_short_kernel",
            "extend_long_warp_kernel",
            "l2_scatter_query_parts_kernel",
        ] {
            assert!(
                custom_res.iter().any(|f| f.fn_name.contains(anchor)),
                "{arch}: missing sanity anchor kernel '{anchor}' in custom resource usage"
            );
        }

        let sass_map: BTreeMap<(&str, &str), &CudaSassStats> = custom_sass
            .iter()
            .map(|s| ((s.member.as_str(), s.fn_name.as_str()), *s))
            .collect();
        assert_eq!(
            custom_sass.len(),
            custom_res.len(),
            "{arch}: custom kernel count mismatch between resource-usage ({}) and SASS ({})",
            custom_res.len(),
            custom_sass.len()
        );

        for r in &custom_res {
            if r.stack != 0 {
                violations.push(format!(
                    "{arch} {}:{}: STACK:{} (expected 0)",
                    r.member, r.fn_name, r.stack
                ));
            }
            if r.local != 0 {
                violations.push(format!(
                    "{arch} {}:{}: LOCAL:{} (expected 0)",
                    r.member, r.fn_name, r.local
                ));
            }
            match sass_map.get(&(r.member.as_str(), r.fn_name.as_str())) {
                None => violations.push(format!(
                    "{arch} {}:{}: present in resource-usage but missing in SASS",
                    r.member, r.fn_name
                )),
                Some(s) => {
                    if s.total_insns == 0 {
                        violations.push(format!(
                            "{arch} {}:{}: parsed 0 SASS instructions",
                            r.member, r.fn_name
                        ));
                    }
                    if s.stl_count != 0 || s.ldl_count != 0 {
                        violations.push(format!(
                            "{arch} {}:{}: {} STL, {} LDL (expected 0/0; total={})",
                            r.member, r.fn_name, s.stl_count, s.ldl_count, s.total_insns
                        ));
                    }
                }
            }
        }
    }

    assert!(
        violations.is_empty(),
        "Custom CUDA kernels violated 0-stack / 0-spill invariant:\n  {}",
        violations.join("\n  ")
    );
}

/// T1-CUDA-4 (W3-F13 Contract): `extend_long_warp_kernel` uses `<= 64` registers and `STACK:0`
/// on `sm_75` (unlocking 100% T4 occupancy).
#[test]
fn tier1_w3_cuda_extend_long_warp_max_64_regs_sm75() {
    let res_fns = cached_cuda_resource_usage("sm_75");
    let matches: Vec<&CudaResourceInfo> = res_fns
        .iter()
        .filter(|f| {
            !is_third_party_cuda_kernel(&f.fn_name) && f.fn_name.contains("extend_long_warp_kernel")
        })
        .collect();
    assert!(
        !matches.is_empty(),
        "sm_75 resource usage must contain extend_long_warp_kernel"
    );
    for entry in matches {
        assert!(
            entry.reg <= 64 && entry.stack == 0,
            "sm_75 {}:{}: expected REG<=64 and STACK:0, got REG:{} STACK:{}",
            entry.member,
            entry.fn_name,
            entry.reg,
            entry.stack
        );
    }
}

/// T1-CUDA-5 (W3-F13 Contract): `l2_scatter_query_parts_kernel` uses shared-memory atomics
/// (`ATOMS` / `REDS` > 0) on `sm_75` via CTA privatization.
#[test]
fn tier1_w3_cuda_l2_scatter_uses_shared_atomics_sm75() {
    let sass_fns = cached_cuda_sass_stats("sm_75");
    let matches: Vec<&CudaSassStats> = sass_fns
        .iter()
        .filter(|f| {
            !is_third_party_cuda_kernel(&f.fn_name)
                && f.fn_name.contains("l2_scatter_query_parts_kernel")
        })
        .collect();
    assert!(
        !matches.is_empty(),
        "sm_75 SASS must contain l2_scatter_query_parts_kernel"
    );
    for entry in matches {
        assert!(
            entry.total_insns > 0 && entry.atoms_count > 0,
            "sm_75 {}:{}: expected shared atomics (ATOMS/REDS > 0), got ATOMS:{} ATOMG:{} total_insns:{}",
            entry.member,
            entry.fn_name,
            entry.atoms_count,
            entry.atomg_count,
            entry.total_insns
        );
    }
}

/// T1-CUDA-6 (W3-F11 Contract): Unsafe `filter_sorted_hits_kernel` (`min_gap = wl/2`) is
/// completely removed from `libgecko_cuda_kernels.a` (both ELF symbols and all 4 archs).
#[test]
fn tier1_w3_cuda_filter_sorted_hits_kernel_removed() {
    let archive = cuda_archive_path();
    let sym_out = Command::new("cuobjdump")
        .args(["--dump-elf-symbols", archive.to_str().unwrap()])
        .output()
        .expect("failed to execute cuobjdump --dump-elf-symbols");
    assert!(
        sym_out.status.success(),
        "cuobjdump --dump-elf-symbols failed:\n{}",
        String::from_utf8_lossy(&sym_out.stderr)
    );
    let sym_stdout = String::from_utf8_lossy(&sym_out.stdout);
    assert!(
        sym_stdout.contains("cuda_query.o") && sym_stdout.contains("query_left_maximal_kernel"),
        "ELF symbol dump sanity check failed: expected cuda_query.o and query_left_maximal_kernel"
    );

    let bad_sym_lines: Vec<&str> = sym_stdout
        .lines()
        .filter(|l| l.contains("filter_sorted_hits_kernel"))
        .collect();
    assert!(
        bad_sym_lines.is_empty(),
        "filter_sorted_hits_kernel must be removed from ELF symbols, found:\n  {}",
        bad_sym_lines.join("\n  ")
    );

    let mut bad_res = Vec::new();
    for arch in CUDA_REQUIRED_ARCHS {
        for f in cached_cuda_resource_usage(arch) {
            if f.fn_name.contains("filter_sorted_hits_kernel") {
                bad_res.push(format!("{arch} {}:{}", f.member, f.fn_name));
            }
        }
    }
    assert!(
        bad_res.is_empty(),
        "filter_sorted_hits_kernel must be removed from all archs, found:\n  {}",
        bad_res.join("\n  ")
    );
}

// ============================================================================
// TIER 2 — Boundary, Corner & Negative Cases (7 Tests)
// ============================================================================

/// T2-1 (Baseline + W3-F6/F7 Regression Guard): Homopolymer, periodic & tie-heavy minimizer
/// SHA-256 parity across `k in {18, 20, 32}` and `forward_only in {false, true}`.
#[test]
fn tier2_w3_homopolymer_and_tie_heavy_minimizer_parity() {
    let tmp = tempfile::tempdir().unwrap();
    let poly_a = "A".repeat(120);
    let poly_t = "T".repeat(120);
    let dinuc_at = "AT".repeat(60);
    let tetranuc = "ACGT".repeat(30);
    let per7 = "ACGTTCA".repeat(18);
    let rep_hf = pseudo_dna(0x2100, 90);

    let hsp_f1 = format!(
        "{}{}{}",
        &poly_a[..40],
        pseudo_dna(30001, 140),
        &dinuc_at[..40]
    );
    let base_f2 = pseudo_dna(30002, 130);
    let hsp_f2_q = format!("{}{}{}", &tetranuc[..35], base_f2, &per7[..35]);
    let hsp_f2_r = format!(
        "{}{}{}",
        &tetranuc[..35],
        mutate_every(&base_f2, 35, 10),
        &per7[..35]
    );
    let hsp_r1 = format!(
        "{}{}{}",
        &poly_t[..40],
        pseudo_dna(30003, 140),
        &poly_a[..40]
    );

    // Separate each HSP onto distinct diagonals so adjacent blocks never merge
    let q_seq = format!(
        "{}{hsp_f1}{}{hsp_f2_q}{}{hsp_r1}{}{poly_a}{}{dinuc_at}{}",
        pseudo_dna(0x2110, 120),
        pseudo_dna(0x2111, 140),
        pseudo_dna(0x2112, 160),
        pseudo_dna(0x2113, 180),
        pseudo_dna(0x2114, 200),
        pseudo_dna(0x2115, 120),
    );
    let mut r_seq = format!(
        "{}{hsp_f1}{}{hsp_f2_r}{}{}{}{poly_a}{}{dinuc_at}{}",
        pseudo_dna(0x2120, 220),
        pseudo_dna(0x2121, 260),
        pseudo_dna(0x2122, 300),
        revcomp(&hsp_r1),
        pseudo_dna(0x2123, 340),
        pseudo_dna(0x2124, 380),
        pseudo_dna(0x2125, 220),
    );
    // Insert 70 copies of rep_hf in ref to exercise HF overflow alongside tie-heavy blocks
    for i in 0..70u64 {
        r_seq.push_str(&pseudo_dna(0x2150 + i, 25));
        r_seq.push_str(&rep_hf);
    }

    let q = tmp.path().join("q.fa");
    let r = tmp.path().join("r.fa");
    write_fasta(&q, &[("q", &q_seq)]);
    write_fasta(&r, &[("r", &r_seq)]);

    for k in ["18", "20", "32"] {
        for fwd_only in [false, true] {
            let tag = format!("k{k}_{}", if fwd_only { "fwd" } else { "both" });
            let out_ex = tmp.path().join(format!("ex_{tag}"));
            let out_min = tmp.path().join(format!("min_{tag}"));
            for (mode, out) in [("exact", &out_ex), ("minimizer", &out_min)] {
                let mut args = vec![
                    "compare",
                    q.to_str().unwrap(),
                    r.to_str().unwrap(),
                    "-k",
                    k,
                    "-l",
                    "80",
                    "-s",
                    "75",
                    "--freq-limit",
                    "20",
                    "--seeding",
                    mode,
                    "-f",
                    "frags",
                    "-o",
                    out.to_str().unwrap(),
                    "--backend",
                    "cpu",
                ];
                if fwd_only {
                    args.push("--forward-only");
                }
                assert_success(&run_gecko(&args), &format!("tie-heavy {tag} {mode}"));
            }
            assert_frags_sha256_identical(
                &out_ex.join("q-r.frags"),
                &out_min.join("q-r.frags"),
                &format!("tie-heavy exact vs minimizer {tag}"),
            );
            let parsed = parse_frags_file(&out_min.join("q-r.frags"));
            let min_frags = if fwd_only { 2 } else { 3 };
            assert!(
                parsed.fragments.len() >= min_frags,
                "{tag}: expected >= {min_frags} frags, got {}",
                parsed.fragments.len()
            );
        }
    }
}

/// T2-2 (Baseline + W3-F6/F7/F8 Regression Guard): Palindromic minimizers & exact `freq_limit`
/// boundaries (`FL=0`, `FL=1`, `FL=20`; `cnt == FL` kept vs `cnt == FL+1` dropped).
#[test]
fn tier2_w3_palindromic_minimizers_and_freq_limit_boundaries() {
    let tmp = tempfile::tempdir().unwrap();
    let half_a = pseudo_dna(30101, 60);
    let pal_a = format!("{half_a}{}", revcomp(&half_a));
    let half_b = pseudo_dna(30102, 60);
    let pal_b = format!("{half_b}{}", revcomp(&half_b));
    assert_eq!(pal_a, revcomp(&pal_a));
    assert_eq!(pal_b, revcomp(&pal_b));

    let np_keep = pseudo_dna(30103, 120);
    let np_drop = pseudo_dna(30104, 120);

    let build_pair = |copies_keep: usize, copies_drop: usize, seed_off: u64| -> (String, String) {
        let qf = [
            pseudo_dna(seed_off, 160),
            pseudo_dna(seed_off + 1, 160),
            pseudo_dna(seed_off + 2, 160),
            pseudo_dna(seed_off + 3, 160),
            pseudo_dna(seed_off + 4, 160),
        ];
        let q_seq = format!(
            "{}{pal_a}{}{pal_b}{}{np_keep}{}{np_drop}{}",
            qf[0], qf[1], qf[2], qf[3], qf[4]
        );

        let avoid_two = |a: u8, b: u8| -> char {
            *b"ACGT"
                .iter()
                .find(|&&c| c != a && c != b)
                .expect("at least 2 bases remain") as char
        };
        // 70-bp dual-strand guards: 2 * 70 * (-4) = -560 < -480, making cross-block diagonal
        // bleed mathematically impossible on both forward and reverse strands.
        let make_dual_strand_guards = |q_left: &str, q_right: &str| -> (String, String) {
            const G: usize = 70;
            let lt = &q_left.as_bytes()[q_left.len() - G..];
            let rh = &q_right.as_bytes()[..G];
            let gl: String = (0..G)
                .map(|k| avoid_two(lt[k], complement_base(rh[G - 1 - k])))
                .collect();
            let gr: String = (0..G)
                .map(|k| avoid_two(rh[k], complement_base(lt[G - 1 - k])))
                .collect();
            (gl, gr)
        };

        let mut r_seq = String::new();
        let mut idx = 0u64;
        let mut push_copies = |blk: &str, n: usize, q_left: &str, q_right: &str| {
            let (gl, gr) = make_dual_strand_guards(q_left, q_right);
            for _ in 0..n {
                idx += 1;
                r_seq.push_str(&gl);
                r_seq.push_str(blk);
                r_seq.push_str(&gr);
                r_seq.push_str(&pseudo_dna(seed_off + 900 + idx, 63));
            }
        };
        push_copies(&pal_a, copies_keep, &qf[0], &qf[1]);
        push_copies(&pal_b, copies_drop, &qf[1], &qf[2]);
        push_copies(&np_keep, copies_keep, &qf[2], &qf[3]);
        push_copies(&np_drop, copies_drop, &qf[3], &qf[4]);
        (q_seq, r_seq)
    };

    // Sub-fixture 1: FL=0 and FL=1
    let (q1s, r1s) = build_pair(1, 2, 0x2200);
    let q1 = tmp.path().join("q1.fa");
    let r1 = tmp.path().join("r1.fa");
    write_fasta(&q1, &[("q1", &q1s)]);
    write_fasta(&r1, &[("r1", &r1s)]);

    for k in ["20", "32"] {
        // FL = 0 -> 0 fragments (16-byte header only)
        let out_ex0 = tmp.path().join(format!("fl0_ex_k{k}"));
        let out_min0 = tmp.path().join(format!("fl0_min_k{k}"));
        for (mode, out) in [("exact", &out_ex0), ("minimizer", &out_min0)] {
            assert_success(
                &run_gecko(&[
                    "compare",
                    q1.to_str().unwrap(),
                    r1.to_str().unwrap(),
                    "-k",
                    k,
                    "-l",
                    "80",
                    "-s",
                    "80",
                    "--freq-limit",
                    "0",
                    "--seeding",
                    mode,
                    "-f",
                    "frags",
                    "-o",
                    out.to_str().unwrap(),
                    "--backend",
                    "cpu",
                ]),
                &format!("FL=0 k={k} {mode}"),
            );
        }
        assert_frags_sha256_identical(
            &out_ex0.join("q1-r1.frags"),
            &out_min0.join("q1-r1.frags"),
            &format!("FL=0 exact vs minimizer k={k}"),
        );
        assert_eq!(fs::read(out_min0.join("q1-r1.frags")).unwrap().len(), 16);

        // FL = 1 -> np_keep (1 'f') + pal_a (1 'f' + 1 'r') = 3 frags; np_drop & pal_b = 0
        let out_ex1 = tmp.path().join(format!("fl1_ex_k{k}"));
        let out_min1 = tmp.path().join(format!("fl1_min_k{k}"));
        for (mode, out) in [("exact", &out_ex1), ("minimizer", &out_min1)] {
            assert_success(
                &run_gecko(&[
                    "compare",
                    q1.to_str().unwrap(),
                    r1.to_str().unwrap(),
                    "-k",
                    k,
                    "-l",
                    "80",
                    "-s",
                    "80",
                    "--freq-limit",
                    "1",
                    "--seeding",
                    mode,
                    "-f",
                    "frags",
                    "-o",
                    out.to_str().unwrap(),
                    "--backend",
                    "cpu",
                ]),
                &format!("FL=1 k={k} {mode}"),
            );
        }
        assert_frags_sha256_identical(
            &out_ex1.join("q1-r1.frags"),
            &out_min1.join("q1-r1.frags"),
            &format!("FL=1 exact vs minimizer k={k}"),
        );
        let f1 = parse_frags_file(&out_min1.join("q1-r1.frags")).fragments;
        let f_cnt = f1.iter().filter(|f| f.strand == 'f').count();
        let r_cnt = f1.iter().filter(|f| f.strand == 'r').count();
        assert_eq!(
            (f1.len(), f_cnt, r_cnt),
            (3, 2, 1),
            "FL=1 k={k}: expected 3 frags (2 fwd, 1 rev), got {f1:?}"
        );
    }

    // Sub-fixture 2: FL=20 -> np_keep (20 'f') + pal_a (20 'f' + 20 'r') = 60 frags
    let (q2s, r2s) = build_pair(20, 21, 0x2280);
    let q2 = tmp.path().join("q2.fa");
    let r2 = tmp.path().join("r2.fa");
    write_fasta(&q2, &[("q2", &q2s)]);
    write_fasta(&r2, &[("r2", &r2s)]);

    for k in ["20", "32"] {
        let out_ex20 = tmp.path().join(format!("fl20_ex_k{k}"));
        let out_min20 = tmp.path().join(format!("fl20_min_k{k}"));
        for (mode, out) in [("exact", &out_ex20), ("minimizer", &out_min20)] {
            assert_success(
                &run_gecko(&[
                    "compare",
                    q2.to_str().unwrap(),
                    r2.to_str().unwrap(),
                    "-k",
                    k,
                    "-l",
                    "80",
                    "-s",
                    "80",
                    "--freq-limit",
                    "20",
                    "--seeding",
                    mode,
                    "-f",
                    "frags",
                    "-o",
                    out.to_str().unwrap(),
                    "--backend",
                    "cpu",
                ]),
                &format!("FL=20 k={k} {mode}"),
            );
        }
        assert_frags_sha256_identical(
            &out_ex20.join("q2-r2.frags"),
            &out_min20.join("q2-r2.frags"),
            &format!("FL=20 exact vs minimizer k={k}"),
        );
        let f20 = parse_frags_file(&out_min20.join("q2-r2.frags")).fragments;
        let f_cnt = f20.iter().filter(|f| f.strand == 'f').count();
        let r_cnt = f20.iter().filter(|f| f.strand == 'r').count();
        assert_eq!(
            (f20.len(), f_cnt, r_cnt),
            (60, 40, 20),
            "FL=20 k={k}: expected 60 frags (40 fwd, 20 rev), got len={}",
            f20.len()
        );
    }
}

/// T2-3 (W3-M1 + W3-M2 Contract): Backward-flank rejected HSP & short noise hit right before
/// `SEQ_SEP` (forward & reverse strands, exact & minimizer) must NEVER hide downstream HSPs.
#[test]
fn tier2_w3_backward_flank_rejected_hsp_seq_sep_barrier() {
    let tmp = tempfile::tempdir().unwrap();
    let q0_back = pseudo_dna(8101, 260);
    let r0_back = mutate_every(&q0_back, 5, 2);
    let seed0 = pseudo_dna(8102, 32);
    let seed1 = pseudo_dna(8103, 32);
    let q1_tail = pseudo_dna(8104, 188);
    let r1_tail = mutate_every(&q1_tail, 20, 10);

    let q0_long = format!("{q0_back}{seed0}");
    let r0_long = format!("{r0_back}{seed0}");
    let q1_rec = format!("{seed1}{q1_tail}");
    let r1_rec = format!("{seed1}{r1_tail}");
    let r0_rev_rec = revcomp(&r1_rec);
    let r1_rev_rec = revcomp(&r0_long);

    let subcases = [
        (
            "long_fwd",
            [("q0", q0_long.as_str()), ("q1", q1_rec.as_str())],
            [("r0", r0_long.as_str()), ("r1", r1_rec.as_str())],
            true,
            1u64,
            'f',
        ),
        (
            "short_fwd",
            [("q0", seed0.as_str()), ("q1", q1_rec.as_str())],
            [("r0", seed0.as_str()), ("r1", r1_rec.as_str())],
            true,
            1u64,
            'f',
        ),
        (
            "long_rev",
            [("q0", q0_long.as_str()), ("q1", q1_rec.as_str())],
            [("r0", r0_rev_rec.as_str()), ("r1", r1_rev_rec.as_str())],
            false,
            0u64,
            'r',
        ),
    ];

    for (label, q_recs, r_recs, fwd_only, exp_seq_y, exp_strand) in subcases {
        let q = tmp.path().join(format!("q_{label}.fa"));
        let r = tmp.path().join(format!("r_{label}.fa"));
        write_fasta(&q, &q_recs);
        write_fasta(&r, &r_recs);

        let out_ex = tmp.path().join(format!("ex_{label}"));
        let out_min = tmp.path().join(format!("min_{label}"));
        for (mode, out) in [("exact", &out_ex), ("minimizer", &out_min)] {
            let mut args = vec![
                "compare",
                q.to_str().unwrap(),
                r.to_str().unwrap(),
                "-k",
                "32",
                "-l",
                "200",
                "-s",
                "75",
                "--seeding",
                mode,
                "-f",
                "frags",
                "-o",
                out.to_str().unwrap(),
                "--backend",
                "cpu",
            ];
            if fwd_only {
                args.push("--forward-only");
            }
            assert_success(&run_gecko(&args), &format!("seq_sep {label} {mode}"));

            let stem = format!("q_{label}-r_{label}.frags");
            let parsed = parse_frags_file(&out.join(&stem));
            assert_eq!(
                parsed.fragments.len(),
                1,
                "{label} ({mode}): expected 1 recovered HSP across SEQ_SEP, got {:?}",
                parsed.fragments
            );
            let f = &parsed.fragments[0];
            assert_eq!(f.length, 220);
            assert_eq!(f.score, 808);
            assert_eq!(f.seq_x, 1);
            assert_eq!(f.seq_y, exp_seq_y);
            assert_eq!(f.strand, exp_strand);
        }

        let stem = format!("q_{label}-r_{label}.frags");
        assert_frags_sha256_identical(
            &out_ex.join(&stem),
            &out_min.join(&stem),
            &format!("seq_sep {label} exact vs minimizer"),
        );
    }
}

/// T2-4 (Baseline + W3-F5/F6..F9 Regression Guard): Dense `'N'`-runs bulk score decay,
/// `--n-penalty` sensitivity (1 vs 4), exact `128` vs `129` `'N'` drop-to-zero boundary in
/// forward AND backward directions, and `exact == minimizer` SHA-256 parity.
#[test]
fn tier2_w3_dense_n_runs_bulk_decay_and_n_penalty_parity() {
    let tmp = tempfile::tempdir().unwrap();

    // Subcase A: Interior 120-'N' bridge + flanking 200-'N' blocks
    let blk_l = pseudo_dna(30401, 100);
    let blk_r = pseudo_dna(30402, 100);
    let blk_adj = pseudo_dna(30403, 140);
    let fl1 = pseudo_dna(30404, 60);
    let fl2 = pseudo_dna(30405, 60);
    let fl3 = pseudo_dna(30406, 60);
    let n120 = "N".repeat(120);
    let n200 = "N".repeat(200);

    // Put `blk_l + N*120 + blk_r` and `blk_adj` on distinct diagonals so they never bleed
    let sep_diag = pseudo_dna(30407, 300);
    let qa = format!("{fl1}{blk_l}{n120}{blk_r}{fl2}{n200}{blk_adj}{n200}{fl3}");
    let ra = format!(
        "{}{blk_l}{n120}{blk_r}{}{sep_diag}{n200}{blk_adj}{n200}{}",
        all_mismatch(&fl1),
        all_mismatch(&fl2),
        all_mismatch(&fl3)
    );
    let qa_p = tmp.path().join("qa.fa");
    let ra_p = tmp.path().join("ra.fa");
    write_fasta(&qa_p, &[("qa", &qa)]);
    write_fasta(&ra_p, &[("ra", &ra)]);

    for (n_pen, exp_lens_scores) in [
        ("1", vec![(320u64, 680u64), (140, 560)]),
        ("4", vec![(100, 400), (100, 400), (140, 560)]),
    ] {
        let out_ex = tmp.path().join(format!("subA_ex_p{n_pen}"));
        let out_min = tmp.path().join(format!("subA_min_p{n_pen}"));
        for (mode, out) in [("exact", &out_ex), ("minimizer", &out_min)] {
            assert_success(
                &run_gecko(&[
                    "compare",
                    qa_p.to_str().unwrap(),
                    ra_p.to_str().unwrap(),
                    "-k",
                    "32",
                    "-l",
                    "80",
                    "-s",
                    "50",
                    "--n-penalty",
                    n_pen,
                    "--seeding",
                    mode,
                    "-f",
                    "frags",
                    "-o",
                    out.to_str().unwrap(),
                    "--backend",
                    "cpu",
                ]),
                &format!("dense N subA n_penalty={n_pen} {mode}"),
            );
        }
        assert_frags_sha256_identical(
            &out_ex.join("qa-ra.frags"),
            &out_min.join("qa-ra.frags"),
            &format!("dense N subA n_penalty={n_pen} exact vs minimizer"),
        );
        let mut got: Vec<(u64, u64)> = parse_frags_file(&out_min.join("qa-ra.frags"))
            .fragments
            .into_iter()
            .map(|f| (f.length, f.score))
            .collect();
        got.sort();
        let mut expected = exp_lens_scores;
        expected.sort();
        assert_eq!(
            got, expected,
            "subA n_penalty={n_pen} (len, score) mismatch"
        );
    }

    // Subcase B: Exact fscore==0 (128 'N') vs fscore==-1 (129 'N') in forward AND backward
    let seed = pseudo_dna(30410, 32);
    let tail_q = pseudo_dna(30411, 80);
    // Mismatches at 10, 30, 50, 70 -> both index 0 and index 79 are matches, max run 19 < 32
    let tail_r = mutate_every(&tail_q, 20, 10);
    let fl_l = pseudo_dna(30412, 60);
    let fl_r = pseudo_dna(30413, 60);

    for (dir_label, is_bwd) in [("fwd", false), ("bwd", true)] {
        for (n_cnt, expect_hit) in [(128usize, true), (129usize, false)] {
            let n_run = "N".repeat(n_cnt);
            let (qs, rs) = if !is_bwd {
                (
                    format!("{fl_l}{seed}{n_run}{tail_q}{fl_r}"),
                    format!(
                        "{}{seed}{n_run}{tail_r}{}",
                        all_mismatch(&fl_l),
                        all_mismatch(&fl_r)
                    ),
                )
            } else {
                (
                    format!("{fl_l}{tail_q}{n_run}{seed}{fl_r}"),
                    format!(
                        "{}{tail_r}{n_run}{seed}{}",
                        all_mismatch(&fl_l),
                        all_mismatch(&fl_r)
                    ),
                )
            };
            let q = tmp.path().join(format!("qb_{dir_label}_{n_cnt}.fa"));
            let r = tmp.path().join(format!("rb_{dir_label}_{n_cnt}.fa"));
            write_fasta(&q, &[("q", &qs)]);
            write_fasta(&r, &[("r", &rs)]);

            let out_ex = tmp.path().join(format!("subB_ex_{dir_label}_{n_cnt}"));
            let out_min = tmp.path().join(format!("subB_min_{dir_label}_{n_cnt}"));
            for (mode, out) in [("exact", &out_ex), ("minimizer", &out_min)] {
                assert_success(
                    &run_gecko(&[
                        "compare",
                        q.to_str().unwrap(),
                        r.to_str().unwrap(),
                        "-k",
                        "32",
                        "-l",
                        "100",
                        "-s",
                        "25",
                        "--n-penalty",
                        "1",
                        "--seeding",
                        mode,
                        "-f",
                        "frags",
                        "-o",
                        out.to_str().unwrap(),
                        "--backend",
                        "cpu",
                    ]),
                    &format!("dense N subB {dir_label} n={n_cnt} {mode}"),
                );
            }
            let stem = format!("qb_{dir_label}_{n_cnt}-rb_{dir_label}_{n_cnt}.frags");
            assert_frags_sha256_identical(
                &out_ex.join(&stem),
                &out_min.join(&stem),
                &format!("dense N subB {dir_label} n={n_cnt} exact vs minimizer"),
            );
            let parsed = parse_frags_file(&out_min.join(&stem));
            if expect_hit {
                assert_eq!(parsed.fragments.len(), 1);
                let f = &parsed.fragments[0];
                assert_eq!(f.length, 240);
                assert_eq!(f.ident, 108);
                assert_eq!(f.score, 288);
            } else {
                assert_eq!(parsed.fragments.len(), 0);
            }
        }
    }
}

/// T2-5 (Baseline + W3-F8 Regression Guard): Balanced multi-slab minimizer under `-m 256K` and
/// `-m 64K` on a 160 Kb reference (guaranteed `slabs >= 2` even after W3-F8 shrinks bytes/entry),
/// asserting 100% SHA-256 parity vs single-slab `--seeding exact -m 2G`.
#[test]
fn tier2_w3_balanced_multislab_minimizer_64k_256k_parity() {
    let tmp = tempfile::tempdir().unwrap();
    let hsp_f = pseudo_dna(0x2501, 240);
    let hsp_r = pseudo_dna(0x2502, 220);
    let hsp_lo = pseudo_dna(0x2503, 200).to_ascii_lowercase();
    let rep_keep = pseudo_dna(0x2504, 120);
    let rep_drop = pseudo_dna(0x2505, 120);

    let qf3 = pseudo_dna(0x2513, 32_000);
    let qf4 = pseudo_dna(0x2514, 16_000);
    let qf5 = pseudo_dna(0x2515, 16_000);
    let q_seq = format!(
        "{}{hsp_f}{}{hsp_r}{}{hsp_lo}{qf3}{rep_keep}{qf4}{rep_drop}{qf5}",
        pseudo_dna(0x2510, 32_000),
        pseudo_dna(0x2511, 32_000),
        pseudo_dna(0x2512, 32_000),
    );

    let keep_gl = all_mismatch(&qf3[qf3.len() - 40..]);
    let keep_gr = all_mismatch(&qf4[..40]);
    let drop_gl = all_mismatch(&qf4[qf4.len() - 40..]);
    let drop_gr = all_mismatch(&qf5[..40]);

    let mut r_seq = format!(
        "{}{hsp_f}{}{}{}{hsp_lo}{}",
        pseudo_dna(0x2520, 40_000),
        pseudo_dna(0x2521, 40_000),
        revcomp(&hsp_r),
        pseudo_dna(0x2522, 40_000),
        pseudo_dna(0x2523, 38_000),
    );
    for i in 0..3u64 {
        r_seq.push_str(&pseudo_dna(0x2530 + i, 80));
        r_seq.push_str(&keep_gl);
        r_seq.push_str(&rep_keep);
        r_seq.push_str(&keep_gr);
    }
    for i in 0..4u64 {
        r_seq.push_str(&pseudo_dna(0x2540 + i, 80));
        r_seq.push_str(&drop_gl);
        r_seq.push_str(&rep_drop);
        r_seq.push_str(&drop_gr);
    }
    r_seq.push_str(&pseudo_dna(0x2550, 2_000));

    let q = tmp.path().join("q.fa");
    let r = tmp.path().join("r.fa");
    write_fasta(&q, &[("q", &q_seq)]);
    write_fasta(&r, &[("r", &r_seq)]);

    let run_cfg = |name: &str, seeding: &str, mem: &str, no_soft: bool| -> (PathBuf, String) {
        let out = tmp.path().join(name);
        let mut args = vec![
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "32",
            "-l",
            "100",
            "-s",
            "80",
            "--freq-limit",
            "3",
            "--seeding",
            seeding,
            "-m",
            mem,
            "-f",
            "frags",
            "-o",
            out.to_str().unwrap(),
            "--backend",
            "cpu",
        ];
        if no_soft {
            args.push("--no-soft-mask");
        }
        let res = run_gecko_with_env(&args, &[("GECKO_PROFILE", Some("1"))]);
        assert_success(&res, name);
        (
            out.join("q-r.frags"),
            String::from_utf8_lossy(&res.stderr).into_owned(),
        )
    };

    let (p_ex_2g, _) = run_cfg("ex_2g", "exact", "2G", false);
    let (p_min_2g, _) = run_cfg("min_2g", "minimizer", "2G", false);
    let (p_min_256k, stderr_256k) = run_cfg("min_256k", "minimizer", "256K", false);
    let (p_min_64k, stderr_64k) = run_cfg("min_64k", "minimizer", "64K", false);

    assert_frags_sha256_identical(&p_ex_2g, &p_min_2g, "exact 2G vs minimizer 2G");
    assert_frags_sha256_identical(&p_ex_2g, &p_min_256k, "exact 2G vs minimizer 256K");
    assert_frags_sha256_identical(&p_ex_2g, &p_min_64k, "exact 2G vs minimizer 64K");
    assert_eq!(
        parse_frags_file(&p_ex_2g).fragments.len(),
        5,
        "soft-masked multislab expected 5 frags (hsp_f, hsp_r, 3x rep_keep)"
    );

    let slabs_256k = sum_u64_key(
        &parse_profile_lines(&stderr_256k),
        "minimizer_multislab",
        "slabs",
    )
    .unwrap_or_else(|| {
        panic!("expected minimizer_multislab slabs= under -m 256K.\nstderr:\n{stderr_256k}")
    });
    let slabs_64k = sum_u64_key(
        &parse_profile_lines(&stderr_64k),
        "minimizer_multislab",
        "slabs",
    )
    .unwrap_or_else(|| {
        panic!("expected minimizer_multislab slabs= under -m 64K.\nstderr:\n{stderr_64k}")
    });
    assert!(
        slabs_256k >= 2 && slabs_64k > slabs_256k,
        "expected slabs(64K)={slabs_64k} > slabs(256K)={slabs_256k} >= 2"
    );

    let (p_unmask_ex, _) = run_cfg("unmask_ex_2g", "exact", "2G", true);
    let (p_unmask_min, _) = run_cfg("unmask_min_256k", "minimizer", "256K", true);
    assert_frags_sha256_identical(
        &p_unmask_ex,
        &p_unmask_min,
        "unmasked exact 2G vs minimizer 256K",
    );
    assert_eq!(
        parse_frags_file(&p_unmask_min).fragments.len(),
        6,
        "unmasked multislab expected 6 frags (+ hsp_lo)"
    );
}

/// T2-6 (Baseline + W3-F4 Regression Guard): Exact inclusive `--two-hit` boundary (`dist - wl <= W`)
/// across DirectCsr (`-k 10`), HashCsr (`-k 12`), and Spaced HashCsr (`--spaced-seed w14s20`),
/// verifying `W=31` rejects, `W=32` accepts, and bare `--two-hit` == `--two-hit=40` SHA-256.
#[test]
fn tier2_w3_two_hit_boundaries_direct_csr_and_hash_csr() {
    let tmp = tempfile::tempdir().unwrap();

    // Build 3 exact-gap=32 fixtures where ONLY `pos1` and `pos2 = pos1 + span + 32` match
    let build_gap32 = |seed_base: u64,
                       span: usize,
                       gap_muts: &[usize],
                       tail_len: usize,
                       tail_mut_step: usize|
     -> (String, String) {
        let f0 = pseudo_dna(seed_base, 120);
        let sa = pseudo_dna(seed_base + 1, span);
        let gap_q = pseudo_dna(seed_base + 2, 32);
        let gap_r = mutate_positions(&gap_q, gap_muts);
        let sb = pseudo_dna(seed_base + 3, span);
        let tail_q = pseudo_dna(seed_base + 4, tail_len);
        let tail_r = mutate_every(&tail_q, tail_mut_step, 0);
        let f1 = pseudo_dna(seed_base + 5, 120);

        let mut q_noise = String::new();
        let mut r_noise = String::new();
        for i in 0..30u64 {
            let blk = pseudo_dna(seed_base + 100 + i, span + 40);
            q_noise.push_str(&blk);
            r_noise.push_str(&all_mismatch(&blk[..20]));
            r_noise.push_str(&blk[20..20 + span]);
            r_noise.push_str(&all_mismatch(&blk[20 + span..]));
        }

        let q = format!("{f0}{sa}{gap_q}{sb}{tail_q}{f1}{q_noise}");
        let r = format!(
            "{}{sa}{gap_r}{sb}{tail_r}{}{r_noise}",
            all_mismatch(&f0),
            all_mismatch(&f1)
        );
        (q, r)
    };

    let cases = [
        (
            "direct_k10",
            "DirectCsr",
            vec!["-k", "10", "-l", "80", "-s", "40"],
            build_gap32(0x2600, 10, &[0, 6, 12, 18, 24, 31], 70, 7),
        ),
        (
            "hash_k12",
            "HashCsr",
            vec!["-k", "12", "-l", "80", "-s", "40"],
            build_gap32(0x2650, 12, &[0, 8, 16, 24, 31], 70, 9),
        ),
        (
            "spaced_w14s20",
            "HashCsr",
            vec!["--spaced-seed", "w14s20", "-l", "80", "-s", "30"],
            build_gap32(
                0x26A0,
                20,
                &[0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 31],
                80,
                3,
            ),
        ),
    ];

    for (label, expected_kind, flags, (qs, rs)) in cases {
        let q = tmp.path().join(format!("q_{label}.fa"));
        let r = tmp.path().join(format!("r_{label}.fa"));
        write_fasta(&q, &[("q", &qs)]);
        write_fasta(&r, &[("r", &rs)]);

        let run_th = |sub: &str, th_arg: &str| -> (PathBuf, String) {
            let out = tmp.path().join(format!("{label}_{sub}"));
            let mut args = vec!["compare", q.to_str().unwrap(), r.to_str().unwrap()];
            args.extend_from_slice(&flags);
            args.extend_from_slice(&[
                "-f",
                "frags",
                "-o",
                out.to_str().unwrap(),
                "--backend",
                "cpu",
                th_arg,
            ]);
            let res = run_gecko_with_env(&args, &[("GECKO_PROFILE", Some("1"))]);
            assert_success(&res, &format!("{label} {th_arg}"));
            (
                out.join(format!("q_{label}-r_{label}.frags")),
                String::from_utf8_lossy(&res.stderr).into_owned(),
            )
        };

        let (p31, stderr31) = run_th("w31", "--two-hit=31");
        let (p32, _) = run_th("w32", "--two-hit=32");
        let (p40, _) = run_th("w40", "--two-hit=40");
        let (pdef, _) = run_th("wdef", "--two-hit");

        let lines31 = parse_profile_lines(&stderr31);
        assert_eq!(
            first_str_key(&lines31, "stream_query", "kind"),
            Some(expected_kind),
            "{label}: expected stream_query kind={expected_kind}.\nstderr:\n{stderr31}"
        );
        let dropped = sum_u64_key(&lines31, "stream_query", "two_hit_drop").unwrap_or(0);
        assert!(
            dropped >= 30,
            "{label}: expected two_hit_drop >= 30, got {dropped}"
        );

        assert_eq!(
            parse_frags_file(&p31).fragments.len(),
            0,
            "{label}: --two-hit=31 must strictly reject gap=32"
        );
        assert_eq!(
            parse_frags_file(&p32).fragments.len(),
            1,
            "{label}: --two-hit=32 must accept gap=32"
        );
        assert_eq!(
            parse_frags_file(&p40).fragments.len(),
            1,
            "{label}: --two-hit=40 must accept gap=32"
        );
        assert_frags_sha256_identical(
            &p40,
            &pdef,
            &format!("{label}: bare --two-hit must equal --two-hit=40"),
        );
    }
}

/// T2-7 (Baseline): Clean non-zero exit with zero panics (`exit != 101`) for invalid flags,
/// ranges, memory specs, and missing/empty inputs across `compare` and `all-vs-all`.
#[test]
fn tier2_w3_negative_cli_flags_compare_and_all_vs_all() {
    let tmp = tempfile::tempdir().unwrap();
    let q = tmp.path().join("q.fa");
    let r = tmp.path().join("r.fa");
    let empty = tmp.path().join("empty.fa");
    write_fasta(&q, &[("q", &pseudo_dna(0x2701, 200))]);
    write_fasta(&r, &[("r", &pseudo_dna(0x2702, 200))]);
    fs::write(&empty, b"").unwrap();

    let single_dir = tmp.path().join("one_genome");
    fs::create_dir_all(&single_dir).unwrap();
    write_fasta(
        &single_dir.join("only.fa"),
        &[("only", &pseudo_dna(0x2703, 200))],
    );

    let two_dir = tmp.path().join("two_genomes");
    fs::create_dir_all(&two_dir).unwrap();
    write_fasta(&two_dir.join("a.fa"), &[("a", &pseudo_dna(0x2704, 200))]);
    write_fasta(&two_dir.join("b.fa"), &[("b", &pseudo_dna(0x2705, 200))]);

    let missing = tmp.path().join("does_not_exist.fa");
    let missing_dir = tmp.path().join("no_such_dir");
    let out = tmp.path().join("out");

    let qs = q.to_str().unwrap();
    let rs = r.to_str().unwrap();
    let es = empty.to_str().unwrap();
    let ms = missing.to_str().unwrap();
    let os = out.to_str().unwrap();
    let td = two_dir.to_str().unwrap();
    let sd = single_dir.to_str().unwrap();
    let md = missing_dir.to_str().unwrap();

    let compare_cases: Vec<(&str, Vec<&str>)> = vec![
        (
            "seeding bogus",
            vec![
                "compare",
                qs,
                rs,
                "--seeding",
                "bogus",
                "-o",
                os,
                "--backend",
                "cpu",
            ],
        ),
        (
            "seeding empty",
            vec![
                "compare",
                qs,
                rs,
                "--seeding",
                "",
                "-o",
                os,
                "--backend",
                "cpu",
            ],
        ),
        (
            "two-hit=0",
            vec![
                "compare",
                qs,
                rs,
                "--two-hit=0",
                "-o",
                os,
                "--backend",
                "cpu",
            ],
        ),
        (
            "two-hit=-1",
            vec![
                "compare",
                qs,
                rs,
                "--two-hit=-1",
                "-o",
                os,
                "--backend",
                "cpu",
            ],
        ),
        (
            "two-hit=abc",
            vec![
                "compare",
                qs,
                rs,
                "--two-hit=abc",
                "-o",
                os,
                "--backend",
                "cpu",
            ],
        ),
        (
            "k=0",
            vec!["compare", qs, rs, "-k", "0", "-o", os, "--backend", "cpu"],
        ),
        (
            "k=3",
            vec!["compare", qs, rs, "-k", "3", "-o", os, "--backend", "cpu"],
        ),
        (
            "k=33",
            vec!["compare", qs, rs, "-k", "33", "-o", os, "--backend", "cpu"],
        ),
        (
            "m=0Z",
            vec!["compare", qs, rs, "-m", "0Z", "-o", os, "--backend", "cpu"],
        ),
        (
            "m empty",
            vec!["compare", qs, rs, "-m", "", "-o", os, "--backend", "cpu"],
        ),
        (
            "freq-limit -1",
            vec![
                "compare",
                qs,
                rs,
                "--freq-limit",
                "-1",
                "-o",
                os,
                "--backend",
                "cpu",
            ],
        ),
        (
            "freq-limit nope",
            vec![
                "compare",
                qs,
                rs,
                "--freq-limit",
                "nope",
                "-o",
                os,
                "--backend",
                "cpu",
            ],
        ),
        (
            "missing query",
            vec!["compare", ms, rs, "-o", os, "--backend", "cpu"],
        ),
        (
            "missing ref",
            vec!["compare", qs, ms, "-o", os, "--backend", "cpu"],
        ),
        (
            "empty fasta",
            vec!["compare", es, rs, "-o", os, "--backend", "cpu"],
        ),
    ];

    for (ctx, args) in compare_cases {
        assert_clean_failure(&run_gecko(&args), &format!("compare {ctx}"));
    }

    let ava_cases: Vec<(&str, Vec<&str>)> = vec![
        (
            "seeding bogus",
            vec![
                "all-vs-all",
                td,
                "--seeding",
                "bogus",
                "-o",
                os,
                "--backend",
                "cpu",
            ],
        ),
        (
            "two-hit=0",
            vec![
                "all-vs-all",
                td,
                "--two-hit=0",
                "-o",
                os,
                "--backend",
                "cpu",
            ],
        ),
        (
            "two-hit=xyz",
            vec![
                "all-vs-all",
                td,
                "--two-hit=xyz",
                "-o",
                os,
                "--backend",
                "cpu",
            ],
        ),
        (
            "unknown flag",
            vec![
                "all-vs-all",
                td,
                "--unknown-flag",
                "-o",
                os,
                "--backend",
                "cpu",
            ],
        ),
        (
            "k=3",
            vec!["all-vs-all", td, "-k", "3", "-o", os, "--backend", "cpu"],
        ),
        (
            "k=33",
            vec!["all-vs-all", td, "-k", "33", "-o", os, "--backend", "cpu"],
        ),
        (
            "s=105",
            vec!["all-vs-all", td, "-s", "105", "-o", os, "--backend", "cpu"],
        ),
        (
            "m=99Q",
            vec!["all-vs-all", td, "-m", "99Q", "-o", os, "--backend", "cpu"],
        ),
        (
            "non-pal spaced",
            vec![
                "all-vs-all",
                td,
                "--spaced-seed",
                "111011",
                "-o",
                os,
                "--backend",
                "cpu",
            ],
        ),
        (
            "missing dir",
            vec!["all-vs-all", md, "-o", os, "--backend", "cpu"],
        ),
        (
            "single genome dir",
            vec!["all-vs-all", sd, "-o", os, "--backend", "cpu"],
        ),
    ];

    for (ctx, args) in ava_cases {
        assert_clean_failure(&run_gecko(&args), &format!("all-vs-all {ctx}"));
    }
}

// ============================================================================
// TIER 3 — Cross-Feature Combinations (3 Tests)
// ============================================================================

/// T3-1 (Baseline + W3-F6..F9 Regression Guard): `--no-soft-mask` + `--seeding minimizer` +
/// `--freq-limit 20` + `--forward-only` + 3-record multi-FASTA + multi-slab `-m 256K` +
/// HF repeat overflow + `--format all` (all 9 formats) + `GECKO_PROFILE=1` + `--json-stats`.
#[test]
fn tier3_w3_no_soft_mask_minimizer_fl20_fwd_multifasta_multislab_all_formats_profile_json() {
    let tmp = tempfile::tempdir().unwrap();
    let h_up = pseudo_dna(0x3101, 200);
    let h_lo = pseudo_dna(0x3102, 200).to_ascii_lowercase();
    let h_rev = pseudo_dna(0x3103, 200);
    let rep_20 = pseudo_dna(0x3104, 90);
    let rep_21 = pseudo_dna(0x3105, 90);
    let hf_core = pseudo_dna(31050, 27);
    let h_rec2 = pseudo_dna(0x3106, 180);

    // Record 0: uppercase fwd, lowercase fwd, uppercase rev (~55 Kb filler)
    let q0 = format!(
        "{}{h_up}{}{h_lo}{}{h_rev}{}",
        pseudo_dna(0x3110, 18_000),
        pseudo_dna(0x3111, 18_000),
        pseudo_dna(0x3112, 18_000),
        pseudo_dna(0x3113, 2_000)
    );
    let r0 = format!(
        "{}{h_up}{}{h_lo}{}{}{}",
        pseudo_dna(0x3120, 18_000),
        pseudo_dna(0x3121, 18_000),
        pseudo_dna(0x3122, 18_000),
        revcomp(&h_rev),
        pseudo_dna(0x3123, 2_000)
    );

    // Record 1: rep_20 (20x in ref == FL) + rep_21 (21x in ref > FL) + 75 HF overflow distinct 32-mers
    let mut q1 = format!(
        "{}{rep_20}{}{rep_21}{}",
        pseudo_dna(0x3130, 25_000),
        pseudo_dna(0x3131, 25_000),
        pseudo_dna(0x3132, 4_000)
    );
    let mut r1 = pseudo_dna(0x3140, 25_000);
    for i in 0..20u64 {
        let g = pseudo_dna(0x3141 + i, 20);
        r1.push_str(&all_mismatch(&g));
        r1.push_str(&rep_20);
        r1.push_str(&pseudo_dna(0x3160 + i, 40));
    }
    for i in 0..21u64 {
        let g = pseudo_dna(0x3180 + i, 20);
        r1.push_str(&all_mismatch(&g));
        r1.push_str(&rep_21);
        r1.push_str(&pseudo_dna(0x31A0 + i, 40));
    }
    for i in 0..75usize {
        let fl_r = encode_base4_5mer(i);
        r1.push_str(&fl_r);
        r1.push_str(&hf_core);
        r1.push_str(&fl_r);
        r1.push_str(&pseudo_dna(0x31C0 + i as u64, 20));

        let fl_q = encode_base4_5mer(100 + i);
        q1.push_str(&fl_q);
        q1.push_str(&hf_core);
        q1.push_str(&fl_q);
        q1.push_str(&pseudo_dna(0x32C0 + i as u64, 20));
    }
    r1.push_str(&pseudo_dna(0x31F0, 25_000));

    // Record 2: cross-record 180-bp HSP (~55 Kb filler)
    let q2 = format!(
        "{}{h_rec2}{}",
        pseudo_dna(0x3200, 28_000),
        pseudo_dna(0x3201, 28_000)
    );
    let r2 = format!(
        "{}{h_rec2}{}",
        pseudo_dna(0x3210, 28_000),
        pseudo_dna(0x3211, 28_000)
    );

    let q = tmp.path().join("q.fa");
    let r = tmp.path().join("r.fa");
    write_fasta(&q, &[("q0", &q0), ("q1", &q1), ("q2", &q2)]);
    write_fasta(&r, &[("r0", &r0), ("r1", &r1), ("r2", &r2)]);

    let out_exact = tmp.path().join("out_exact");
    let out_combo = tmp.path().join("out_combo");

    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "32",
            "-l",
            "80",
            "-s",
            "80",
            "--no-soft-mask",
            "--freq-limit",
            "20",
            "--forward-only",
            "--seeding",
            "exact",
            "-m",
            "2G",
            "-f",
            "frags",
            "-o",
            out_exact.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "T3-1 exact oracle",
    );

    let res_combo = run_gecko_with_env(
        &[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "32",
            "-l",
            "80",
            "-s",
            "80",
            "--no-soft-mask",
            "--freq-limit",
            "20",
            "--forward-only",
            "--seeding",
            "minimizer",
            "-m",
            "256K",
            "-f",
            "all",
            "--json-stats",
            "-o",
            out_combo.to_str().unwrap(),
            "--backend",
            "cpu",
        ],
        &[("GECKO_PROFILE", Some("1"))],
    );
    assert_success(&res_combo, "T3-1 combo");

    assert_frags_sha256_identical(
        &out_exact.join("q-r.frags"),
        &out_combo.join("q-r.frags"),
        "T3-1 exact 2G vs minimizer 256K combo",
    );

    let stdout_str = String::from_utf8_lossy(&res_combo.stdout);
    assert!(
        !stdout_str.contains("[GECKO_PROFILE]"),
        "GECKO_PROFILE must never pollute --json-stats stdout"
    );
    let json: serde_json::Value =
        serde_json::from_slice(&res_combo.stdout).expect("valid --json-stats stdout");

    let ff = validate_all_9_formats(&out_combo, "q-r", true);
    assert_eq!(json["reverse_fragments"].as_u64(), Some(0));
    assert_eq!(json["fragments"].as_u64(), Some(ff.fragments.len() as u64));
    assert!(
        ff.fragments.len() >= 23,
        "expected >= 23 fragments (h_up + h_lo + 20x rep_20 + q2/r2), got {}",
        ff.fragments.len()
    );
    assert!(
        ff.fragments.iter().any(|f| f.seq_x == 2 && f.seq_y == 2),
        "expected (seq_x=2, seq_y=2) fragment in multi-FASTA output"
    );

    let stderr_str = String::from_utf8_lossy(&res_combo.stderr);
    let lines = parse_profile_lines(&stderr_str);
    let slabs = sum_u64_key(&lines, "minimizer_multislab", "slabs")
        .expect("expected minimizer_multislab slabs= on stderr");
    assert!(
        slabs >= 2,
        "expected slabs >= 2 under -m 256K, got {slabs}.\nstderr:\n{stderr_str}"
    );
}

/// T3-2 (W3-M2 Contract: W3-F10): `all-vs-all --no-soft-mask --forward-only --seeding auto`
/// across 3 unequal-length genomes matches pairwise `compare --seeding exact` SHA-256 bit-for-bit.
#[test]
fn tier3_w3_all_vs_all_no_soft_mask_forward_only_seeding_auto_vs_pairwise_exact() {
    let tmp = tempfile::tempdir().unwrap();
    let gdir = tmp.path().join("gdir");
    fs::create_dir_all(&gdir).unwrap();

    let h_up = pseudo_dna(0x3201, 200);
    let h_lo = pseudo_dna(0x3202, 180).to_ascii_lowercase();
    let h_rev = pseudo_dna(0x3203, 180);

    let s1 = format!(
        "{}{h_up}{}{h_lo}{}{h_rev}{}",
        pseudo_dna(0x3210, 400),
        pseudo_dna(0x3211, 400),
        pseudo_dna(0x3212, 400),
        pseudo_dna(0x3213, 300)
    );
    let s2 = format!(
        "{}{h_up}{}{h_lo}{}{}{}",
        pseudo_dna(0x3220, 900),
        pseudo_dna(0x3221, 900),
        pseudo_dna(0x3222, 900),
        revcomp(&h_rev),
        pseudo_dna(0x3223, 800)
    );
    let s3 = format!(
        "{}{h_up}{}{h_lo}{}{}{}",
        pseudo_dna(0x3230, 1400),
        pseudo_dna(0x3231, 1400),
        pseudo_dna(0x3232, 1400),
        revcomp(&h_rev),
        pseudo_dna(0x3233, 1300)
    );

    let p1 = gdir.join("g1.fa");
    let p2 = gdir.join("g2.fa");
    let p3 = gdir.join("g3.fa");
    write_fasta(&p1, &[("g1", &s1)]);
    write_fasta(&p2, &[("g2", &s2)]);
    write_fasta(&p3, &[("g3", &s3)]);

    let out_ava = tmp.path().join("out_ava");
    assert_success(
        &run_gecko(&[
            "all-vs-all",
            gdir.to_str().unwrap(),
            "-k",
            "24",
            "-l",
            "100",
            "-s",
            "80",
            "--freq-limit",
            "20",
            "--no-soft-mask",
            "--forward-only",
            "--seeding",
            "auto",
            "-f",
            "frags,paf",
            "-o",
            out_ava.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "T3-2 all-vs-all --no-soft-mask --forward-only",
    );

    let out_masked = tmp.path().join("out_ava_masked");
    assert_success(
        &run_gecko(&[
            "all-vs-all",
            gdir.to_str().unwrap(),
            "-k",
            "24",
            "-l",
            "100",
            "-s",
            "80",
            "--freq-limit",
            "20",
            "--forward-only",
            "--seeding",
            "auto",
            "-f",
            "frags",
            "-o",
            out_masked.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "T3-2 all-vs-all soft-masked control",
    );

    let pairs: [(&str, &Path, &Path, usize, usize); 3] = [
        ("g1-g2", &p1, &p2, 2, 1),
        ("g1-g3", &p1, &p3, 2, 1),
        ("g2-g3", &p2, &p3, 3, 2),
    ];

    for (stem, qa, rb, exp_unmasked, exp_masked) in pairs {
        let out_pair = tmp.path().join(format!("pair_{stem}"));
        assert_success(
            &run_gecko(&[
                "compare",
                qa.to_str().unwrap(),
                rb.to_str().unwrap(),
                "-k",
                "24",
                "-l",
                "100",
                "-s",
                "80",
                "--freq-limit",
                "20",
                "--no-soft-mask",
                "--forward-only",
                "--seeding",
                "exact",
                "-f",
                "frags",
                "-o",
                out_pair.to_str().unwrap(),
                "--backend",
                "cpu",
            ]),
            &format!("T3-2 pairwise exact {stem}"),
        );

        let ava_p = out_ava.join(format!("{stem}.frags"));
        let pair_p = out_pair.join(format!("{stem}.frags"));
        assert_frags_sha256_identical(&pair_p, &ava_p, &format!("T3-2 {stem} ava vs pairwise"));

        let f_unmask = parse_frags_file(&ava_p).fragments;
        let f_mask = parse_frags_file(&out_masked.join(format!("{stem}.frags"))).fragments;
        assert!(f_unmask.iter().all(|f| f.strand == 'f'));
        assert_eq!(f_unmask.len(), exp_unmasked, "T3-2 {stem} unmasked count");
        assert_eq!(f_mask.len(), exp_masked, "T3-2 {stem} masked count");
        assert!(out_ava.join(format!("{stem}.paf")).exists());
    }
}

/// T3-3 (Baseline + W3-F3/F4/F5 Regression Guard): `--spaced-seed w14s20` + `--two-hit=40` +
/// `GECKO_PROFILE=1` + `--format all` + `--json-stats` on divergent multi-FASTA homologs.
#[test]
fn tier3_w3_spaced_w14s20_two_hit_profile_format_all() {
    let tmp = tempfile::tempdir().unwrap();

    // w14s20 = "11110101100110101111" -> don't-care at 4, 6, 9, 10, 13, 15.
    // Mutating at (4, 10, 15) every 20 bp yields 85% identity, 0 contiguous 14-mers,
    // and multiple spaced w14s20 hits every 20 bp (<= 40 bp two-hit window)!
    let mutate_w14s20_dc = |s: &str| -> String {
        s.bytes()
            .enumerate()
            .map(|(i, b)| {
                let m = i % 20;
                if m == 4 || m == 10 || m == 15 {
                    alt_base(b)
                } else {
                    b
                }
            })
            .map(|b| b as char)
            .collect()
    };

    let h_fwd_q = pseudo_dna(0x3301, 240);
    let h_fwd_r = mutate_w14s20_dc(&h_fwd_q);
    let h_rev_q = pseudo_dna(0x3302, 240);
    let h_rev_r = revcomp(&mutate_w14s20_dc(&h_rev_q));

    // Inject 40 isolated 20-bp single-hit noise seeds on DISTINCT diagonals -> two_hit_drop >= 40
    let mut q_noise = String::new();
    let mut r_noise = String::new();
    for i in 0..40u64 {
        let blk = pseudo_dna(0x3350 + i, 60);
        q_noise.push_str(&blk);
        r_noise.push_str(&all_mismatch(&blk[..20]));
        r_noise.push_str(&blk[20..40]);
        r_noise.push_str(&all_mismatch(&blk[40..]));
        r_noise.push_str(&pseudo_dna(0x33A0 + i, 5));
    }

    let q0 = format!(
        "{}{h_fwd_q}{}{q_noise}",
        pseudo_dna(0x3310, 3_000),
        pseudo_dna(0x3311, 3_000)
    );
    let r0 = format!(
        "{}{h_fwd_r}{}{r_noise}",
        pseudo_dna(0x3320, 3_000),
        pseudo_dna(0x3321, 3_000)
    );
    let q1 = format!(
        "{}{h_rev_q}{}",
        pseudo_dna(0x3330, 3_000),
        pseudo_dna(0x3331, 3_000)
    );
    let r1 = format!(
        "{}{h_rev_r}{}",
        pseudo_dna(0x3340, 3_000),
        pseudo_dna(0x3341, 3_000)
    );

    let q = tmp.path().join("q.fa");
    let r = tmp.path().join("r.fa");
    write_fasta(&q, &[("q0", &q0), ("q1", &q1)]);
    write_fasta(&r, &[("r0", &r0), ("r1", &r1)]);

    let out_single = tmp.path().join("out_single");
    let out_th = tmp.path().join("out_twohit");

    let res_single = run_gecko(&[
        "compare",
        q.to_str().unwrap(),
        r.to_str().unwrap(),
        "--spaced-seed",
        "w14s20",
        "-l",
        "80",
        "-s",
        "50",
        "-f",
        "frags",
        "--json-stats",
        "-o",
        out_single.to_str().unwrap(),
        "--backend",
        "cpu",
    ]);
    assert_success(&res_single, "T3-3 single-hit");
    let j_single: serde_json::Value = serde_json::from_slice(&res_single.stdout).unwrap();

    let res_th = run_gecko_with_env(
        &[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "--spaced-seed",
            "w14s20",
            "--two-hit=40",
            "-l",
            "80",
            "-s",
            "50",
            "-f",
            "all",
            "--json-stats",
            "-o",
            out_th.to_str().unwrap(),
            "--backend",
            "cpu",
        ],
        &[("GECKO_PROFILE", Some("1"))],
    );
    assert_success(&res_th, "T3-3 two-hit combo");

    let stdout_th = String::from_utf8_lossy(&res_th.stdout);
    assert!(!stdout_th.contains("[GECKO_PROFILE]"));
    let j_th: serde_json::Value = serde_json::from_slice(&res_th.stdout).unwrap();

    let hits_single = j_single["hits_used"].as_u64().unwrap();
    let hits_th = j_th["hits_used"].as_u64().unwrap();
    assert!(
        hits_th < hits_single,
        "two-hit must reduce hits_used ({hits_th} vs single-hit {hits_single})"
    );

    let stderr_th = String::from_utf8_lossy(&res_th.stderr);
    let lines_th = parse_profile_lines(&stderr_th);
    assert_eq!(
        first_str_key(&lines_th, "stream_query", "kind"),
        Some("HashCsr")
    );
    let drop_cnt = sum_u64_key(&lines_th, "stream_query", "two_hit_drop").unwrap_or(0);
    assert!(
        drop_cnt >= 40,
        "expected two_hit_drop >= 40, got {drop_cnt}.\nstderr:\n{stderr_th}"
    );

    let ff = validate_all_9_formats(&out_th, "q-r", false);
    assert_eq!(
        ff.fragments.len(),
        2,
        "expected 2 divergent homologs (1f, 1r)"
    );
    assert_eq!(ff.fragments.iter().filter(|f| f.strand == 'f').count(), 1);
    assert_eq!(ff.fragments.iter().filter(|f| f.strand == 'r').count(), 1);
}

// ============================================================================
// TIER 4 — Real-World Application Scenarios (5 Tests)
// ============================================================================

/// T4-1 (W3-M1 + W3-M2 Contract): Real unmasked primate chromosome comparison at Colab T4 params
/// (`chimp_chr21[25M..27M]` x `gorilla_chr21[29.65M..31.75M]`, `--no-soft-mask`, `WL=32, L=1000,
/// S=75, FL=20`). Verifies 100% SHA-256 parity across `exact`, `minimizer`, and default CLI,
/// `fast_reject > 0`, default auto-minimizer selection, and bounded `exact_k_reject < 360_000`.
#[test]
fn tier4_w3_real_chimp_gorilla_unmasked_colab_params_parity_and_hf_overflow() {
    let chimp = repo_root().join("chimp_chr21.fa");
    let gorilla = repo_root().join("gorilla_chr21.fa");
    if !chimp.exists() || !gorilla.exists() {
        eprintln!(
            "SKIPPED tier4_w3_real_chimp_gorilla_unmasked_colab_params_parity_and_hf_overflow: real FASTA absent"
        );
        return;
    }

    let tmp = tempfile::tempdir().unwrap();
    let q = tmp.path().join("chimp_2m.fa");
    let r = tmp.path().join("gorilla_2m1.fa");
    let q_len = slice_real_fasta_to_temp(&chimp, &q, "chimp_25M_27M", 25_000_000, 27_000_000);
    let r_len =
        slice_real_fasta_to_temp(&gorilla, &r, "gorilla_2965M_3175M", 29_650_000, 31_750_000);
    assert_eq!(q_len, 2_000_000);
    assert_eq!(r_len, 2_100_000);

    let run_mode = |label: &str, seeding: Option<&str>| -> (PathBuf, String) {
        let out = tmp.path().join(label);
        let mut args = vec![
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "32",
            "-l",
            "1000",
            "-s",
            "75",
            "--freq-limit",
            "20",
            "--no-soft-mask",
            "-f",
            "frags",
            "-o",
            out.to_str().unwrap(),
            "--backend",
            "cpu",
        ];
        if let Some(m) = seeding {
            args.push("--seeding");
            args.push(m);
        }
        let res = run_gecko_with_env(&args, &[("GECKO_PROFILE", Some("1"))]);
        assert_success(&res, &format!("T4-1 ({label})"));
        (
            out.join("chimp_2m-gorilla_2m1.frags"),
            String::from_utf8_lossy(&res.stderr).into_owned(),
        )
    };

    let (p_ex, stderr_ex) = run_mode("exact", Some("exact"));
    let (p_min, stderr_min) = run_mode("minimizer", Some("minimizer"));
    let (p_def, stderr_def) = run_mode("default", None);

    assert_frags_sha256_identical(&p_ex, &p_min, "T4-1 exact vs minimizer");
    assert_frags_sha256_identical(&p_ex, &p_def, "T4-1 exact vs default");

    let frags = parse_frags_file(&p_ex).fragments;
    assert!(
        frags.len() >= 500,
        "T4-1 expected >= 500 unmasked chimp x gorilla fragments, got {}",
        frags.len()
    );

    let ex_lines = parse_profile_lines(&stderr_ex);
    let fast_rej = sum_u64_key(&ex_lines, "stream_query", "fast_reject");
    assert!(
        fast_rej > Some(0),
        "T4-1 exact k=32 must have fast_reject > 0 (W3-F2), got {fast_rej:?}.\nstderr:\n{stderr_ex}"
    );

    let def_lines = parse_profile_lines(&stderr_def);
    assert!(
        used_minimizer(&def_lines) && !used_exact_stream(&def_lines),
        "T4-1 default CLI (no --seeding) must use minimizer_query (W3-F10).\nstderr:\n{stderr_def}"
    );

    let min_lines = parse_profile_lines(&stderr_min);
    let exact_k_rej = sum_u64_key(&min_lines, "minimizer_query", "exact_k_reject")
        .expect("minimizer_query must report exact_k_reject");
    assert!(
        exact_k_rej < 360_000,
        "T4-1 unmasked HF repeat overflow must reduce exact_k_reject < 360_000 (b3578ac baseline=428,847), got {exact_k_rej}.\nstderr:\n{stderr_min}"
    );
}

/// T4-2 (W3-M2 Contract): Real unmasked distant mammalian synteny (`chimp_chr21[19.4M..20.6M]` x
/// `mouse_chr16[90.95M..92.05M]`, 0% N, `--no-soft-mask`, `WL=32, L=1000, S=75, FL=20`).
/// Verifies 100% SHA-256 parity across `exact`, `minimizer`, and default CLI, default auto-minimizer
/// selection, and `exact_k_reject < 50_000` (collapsing `260,982`).
#[test]
fn tier4_w3_real_chimp_mouse_unmasked_syntenic_parity_and_hf_overflow() {
    let chimp = repo_root().join("chimp_chr21.fa");
    let mouse = repo_root().join("mouse_chr16.fa");
    if !chimp.exists() || !mouse.exists() {
        eprintln!(
            "SKIPPED tier4_w3_real_chimp_mouse_unmasked_syntenic_parity_and_hf_overflow: real FASTA absent"
        );
        return;
    }

    let tmp = tempfile::tempdir().unwrap();
    let q = tmp.path().join("chimp_1m2.fa");
    let r = tmp.path().join("mouse_1m1.fa");
    let q_len = slice_real_fasta_to_temp(&chimp, &q, "chimp_194M_206M", 19_400_000, 20_600_000);
    let r_len = slice_real_fasta_to_temp(&mouse, &r, "mouse_9095M_9205M", 90_950_000, 92_050_000);
    assert_eq!(q_len, 1_200_000);
    assert_eq!(r_len, 1_100_000);

    let run_mode = |label: &str, seeding: Option<&str>| -> (PathBuf, String) {
        let out = tmp.path().join(label);
        let mut args = vec![
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "32",
            "-l",
            "1000",
            "-s",
            "75",
            "--freq-limit",
            "20",
            "--no-soft-mask",
            "-f",
            "frags",
            "-o",
            out.to_str().unwrap(),
            "--backend",
            "cpu",
        ];
        if let Some(m) = seeding {
            args.push("--seeding");
            args.push(m);
        }
        let res = run_gecko_with_env(&args, &[("GECKO_PROFILE", Some("1"))]);
        assert_success(&res, &format!("T4-2 ({label})"));
        (
            out.join("chimp_1m2-mouse_1m1.frags"),
            String::from_utf8_lossy(&res.stderr).into_owned(),
        )
    };

    let (p_ex, _) = run_mode("exact", Some("exact"));
    let (p_min, stderr_min) = run_mode("minimizer", Some("minimizer"));
    let (p_def, stderr_def) = run_mode("default", None);

    assert_frags_sha256_identical(&p_ex, &p_min, "T4-2 exact vs minimizer");
    assert_frags_sha256_identical(&p_ex, &p_def, "T4-2 exact vs default");

    let frags = parse_frags_file(&p_ex).fragments;
    assert!(
        frags.len() >= 3,
        "T4-2 expected >= 3 real chimp x mouse syntenic HSPs, got {}",
        frags.len()
    );

    let def_lines = parse_profile_lines(&stderr_def);
    assert!(
        used_minimizer(&def_lines) && !used_exact_stream(&def_lines),
        "T4-2 default CLI (no --seeding) must use minimizer_query (W3-F10).\nstderr:\n{stderr_def}"
    );

    let min_lines = parse_profile_lines(&stderr_min);
    let exact_k_rej = sum_u64_key(&min_lines, "minimizer_query", "exact_k_reject")
        .expect("minimizer_query must report exact_k_reject");
    assert!(
        exact_k_rej < 50_000,
        "T4-2 unmasked chimp x mouse HF overflow must collapse exact_k_reject < 50_000 (b3578ac baseline=260,982), got {exact_k_rej}.\nstderr:\n{stderr_min}"
    );
}

/// T4-3 (Baseline): Soft-masked real primate chromosome comparison (`chimp_chr21[25M..26M]` x
/// `gorilla_chr21[29.65M..30.70M]`, `WL=32, L=1000, S=75, FL=20`) — 100% SHA-256 `.frags` parity
/// across `exact`, `minimizer`, and `auto`.
#[test]
fn tier4_w3_real_chimp_gorilla_softmasked_colab_params_sha256_parity() {
    let chimp = repo_root().join("chimp_chr21.fa");
    let gorilla = repo_root().join("gorilla_chr21.fa");
    if !chimp.exists() || !gorilla.exists() {
        eprintln!(
            "SKIPPED tier4_w3_real_chimp_gorilla_softmasked_colab_params_sha256_parity: real FASTA absent"
        );
        return;
    }

    let tmp = tempfile::tempdir().unwrap();
    let q = tmp.path().join("chimp_1m.fa");
    let r = tmp.path().join("gorilla_1m05.fa");
    assert_eq!(
        slice_real_fasta_to_temp(&chimp, &q, "chimp_25M_26M", 25_000_000, 26_000_000),
        1_000_000
    );
    assert_eq!(
        slice_real_fasta_to_temp(&gorilla, &r, "gorilla_2965M_3070M", 29_650_000, 30_700_000),
        1_050_000
    );

    let mut paths = Vec::new();
    for mode in ["exact", "minimizer", "auto"] {
        let out = tmp.path().join(mode);
        assert_success(
            &run_gecko(&[
                "compare",
                q.to_str().unwrap(),
                r.to_str().unwrap(),
                "-k",
                "32",
                "-l",
                "1000",
                "-s",
                "75",
                "--freq-limit",
                "20",
                "--seeding",
                mode,
                "-f",
                "frags",
                "-o",
                out.to_str().unwrap(),
                "--backend",
                "cpu",
            ]),
            &format!("T4-3 ({mode})"),
        );
        paths.push(out.join("chimp_1m-gorilla_1m05.frags"));
    }

    let frags = parse_frags_file(&paths[0]).fragments;
    assert!(
        frags.len() >= 200,
        "T4-3 expected >= 200 soft-masked chimp x gorilla fragments, got {}",
        frags.len()
    );
    assert_frags_sha256_identical(&paths[0], &paths[1], "T4-3 exact vs minimizer");
    assert_frags_sha256_identical(&paths[0], &paths[2], "T4-3 exact vs auto");
}

/// T4-4 (Baseline + W3-F4 DirectCsr Early Two-Hit Stress): Real `k=12` `DirectCsr` run
/// (`chimp_chr21[25M..25.15M]` x `gorilla_chr21[22M..30.5M]` 8.50 Mb ref triggers `DirectCsr`!).
#[test]
fn tier4_w3_real_chimp_gorilla_k12_direct_csr_two_hit_recall_and_noise_reduction() {
    let chimp = repo_root().join("chimp_chr21.fa");
    let gorilla = repo_root().join("gorilla_chr21.fa");
    if !chimp.exists() || !gorilla.exists() {
        eprintln!(
            "SKIPPED tier4_w3_real_chimp_gorilla_k12_direct_csr_two_hit_recall_and_noise_reduction: real FASTA absent"
        );
        return;
    }

    let tmp = tempfile::tempdir().unwrap();
    let q = tmp.path().join("chimp_150k.fa");
    let r = tmp.path().join("gorilla_8m5.fa");
    assert_eq!(
        slice_real_fasta_to_temp(&chimp, &q, "chimp_25M_2515M", 25_000_000, 25_150_000),
        150_000
    );
    assert_eq!(
        slice_real_fasta_to_temp(&gorilla, &r, "gorilla_22M_305M", 22_000_000, 30_500_000),
        8_500_000
    );

    let out_single = tmp.path().join("single");
    let out_th = tmp.path().join("two_hit");

    let res_single = run_gecko(&[
        "compare",
        q.to_str().unwrap(),
        r.to_str().unwrap(),
        "-k",
        "12",
        "-l",
        "200",
        "-s",
        "75",
        "--freq-limit",
        "20",
        "-f",
        "frags",
        "--json-stats",
        "-o",
        out_single.to_str().unwrap(),
        "--backend",
        "cpu",
    ]);
    assert_success(&res_single, "T4-4 single-hit");
    let j_single: serde_json::Value = serde_json::from_slice(&res_single.stdout).unwrap();

    let res_th = run_gecko_with_env(
        &[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "12",
            "-l",
            "200",
            "-s",
            "75",
            "--freq-limit",
            "20",
            "--two-hit=40",
            "-f",
            "frags",
            "--json-stats",
            "-o",
            out_th.to_str().unwrap(),
            "--backend",
            "cpu",
        ],
        &[("GECKO_PROFILE", Some("1"))],
    );
    assert_success(&res_th, "T4-4 two-hit");
    let j_th: serde_json::Value = serde_json::from_slice(&res_th.stdout).unwrap();

    let stderr_th = String::from_utf8_lossy(&res_th.stderr);
    let lines_th = parse_profile_lines(&stderr_th);
    assert_eq!(
        first_str_key(&lines_th, "stream_query", "kind"),
        Some("DirectCsr"),
        "T4-4 must exercise DirectCsr on 8.50 Mb ref at k=12.\nstderr:\n{stderr_th}"
    );
    let dropped = sum_u64_key(&lines_th, "stream_query", "two_hit_drop").unwrap_or(0);
    assert!(
        dropped > 50_000,
        "T4-4 expected two_hit_drop > 50_000, got {dropped}"
    );

    let frags_single = j_single["fragments"].as_u64().unwrap();
    let frags_th = j_th["fragments"].as_u64().unwrap();
    let hits_single = j_single["hits_used"].as_u64().unwrap();
    let hits_th = j_th["hits_used"].as_u64().unwrap();
    assert!(
        frags_th >= 100 && frags_th * 100 >= frags_single * 85,
        "T4-4 expected >= 100 frags and >= 85% recall (single={frags_single}, two_hit={frags_th})"
    );
    assert!(
        hits_th * 50 < hits_single,
        "T4-4 expected >50x hits_used reduction (single={hits_single}, two_hit={hits_th})"
    );
}

/// T4-5 (Baseline + W3-F4 HashCsr Early Two-Hit Stress): Real divergent mammalian synteny
/// (`chimp_chr21[19.9M..20.2M]` x `mouse_chr16[91.35M..91.65M]`, 0% N) with `--spaced-seed w14s20`
/// and `--two-hit=40`.
#[test]
fn tier4_w3_real_chimp_mouse_w14s20_syntenic_two_hit_recall_and_noise_reduction() {
    let chimp = repo_root().join("chimp_chr21.fa");
    let mouse = repo_root().join("mouse_chr16.fa");
    if !chimp.exists() || !mouse.exists() {
        eprintln!(
            "SKIPPED tier4_w3_real_chimp_mouse_w14s20_syntenic_two_hit_recall_and_noise_reduction: real FASTA absent"
        );
        return;
    }

    let tmp = tempfile::tempdir().unwrap();
    let q = tmp.path().join("chimp_300k.fa");
    let r = tmp.path().join("mouse_300k.fa");
    assert_eq!(
        slice_real_fasta_to_temp(&chimp, &q, "chimp_199M_202M", 19_900_000, 20_200_000),
        300_000
    );
    assert_eq!(
        slice_real_fasta_to_temp(&mouse, &r, "mouse_9135M_9165M", 91_350_000, 91_650_000),
        300_000
    );

    let out_single = tmp.path().join("single");
    let out_th = tmp.path().join("two_hit");

    let res_single = run_gecko(&[
        "compare",
        q.to_str().unwrap(),
        r.to_str().unwrap(),
        "--spaced-seed",
        "w14s20",
        "-l",
        "100",
        "-s",
        "60",
        "--freq-limit",
        "50",
        "-f",
        "frags",
        "--json-stats",
        "-o",
        out_single.to_str().unwrap(),
        "--backend",
        "cpu",
    ]);
    assert_success(&res_single, "T4-5 single-hit");
    let j_single: serde_json::Value = serde_json::from_slice(&res_single.stdout).unwrap();

    let res_th = run_gecko_with_env(
        &[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "--spaced-seed",
            "w14s20",
            "--two-hit=40",
            "-l",
            "100",
            "-s",
            "60",
            "--freq-limit",
            "50",
            "-f",
            "frags",
            "--json-stats",
            "-o",
            out_th.to_str().unwrap(),
            "--backend",
            "cpu",
        ],
        &[("GECKO_PROFILE", Some("1"))],
    );
    assert_success(&res_th, "T4-5 two-hit");
    let j_th: serde_json::Value = serde_json::from_slice(&res_th.stdout).unwrap();

    let stderr_th = String::from_utf8_lossy(&res_th.stderr);
    let lines_th = parse_profile_lines(&stderr_th);
    assert_eq!(
        first_str_key(&lines_th, "stream_query", "kind"),
        Some("HashCsr"),
        "T4-5 must exercise HashCsr for w14s20.\nstderr:\n{stderr_th}"
    );
    let dropped = sum_u64_key(&lines_th, "stream_query", "two_hit_drop").unwrap_or(0);
    assert!(
        dropped > 500,
        "T4-5 expected two_hit_drop > 500, got {dropped}"
    );

    let frags_single = j_single["fragments"].as_u64().unwrap();
    let frags_th = j_th["fragments"].as_u64().unwrap();
    let hits_single = j_single["hits_used"].as_u64().unwrap();
    let hits_th = j_th["hits_used"].as_u64().unwrap();
    assert!(
        frags_single >= 50 && frags_th * 100 >= frags_single * 85,
        "T4-5 expected frags_single >= 50 and >= 85% recall (single={frags_single}, two_hit={frags_th})"
    );
    assert!(
        hits_th * 4 < hits_single,
        "T4-5 expected >4x hits_used reduction (single={hits_single}, two_hit={hits_th})"
    );
}
