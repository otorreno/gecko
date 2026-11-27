//! Opaque-box 4-tier E2E test suite for GECKO v2 (Wave 2: CPU & CUDA optimization +
//! comparison mechanism challenge).
//!
//! Design rules:
//! - Strictly opaque-box: invokes the compiled `gecko` binary (and `cargo`/`nvcc`/`cuobjdump`
//!   for CUDA static verification) with zero dependency on `gecko_core` internals.
//! - Deterministic on CPU-only VMs: all comparison tests pass `--backend cpu` unless
//!   explicitly testing backend dispatch or CUDA static compilation/SASS.
//! - Isolated hygiene: every test writes solely inside its own `tempfile::TempDir` and
//!   never modifies repository files, `.agents/teamwork/`, or untracked benchmark FASTAs.
//! - 4-tier naming:
//!   - `tier1_*`: Feature coverage (contiguous k-mers, spaced seeds, strands/soft-mask/gzip,
//!     `freq_limit`/memory budgets, multi-FASTA/`N`/IUPAC, output formats/`all-vs-all`,
//!     M3 `--seeding`, M3 `--two-hit`, M3 `GECKO_PROFILE=1`, M4 CUDA static checks).
//!   - `tier2_*`: Boundary, corner, and negative CLI/input cases.
//!   - `tier3_*`: Cross-feature combinations.
//!   - `tier4_*`: Real-world synthetic and real-chromosome scenarios.

use std::cmp::Reverse;
use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::OnceLock;

use flate2::Compression;
use flate2::write::GzEncoder;

// ============================================================================
// Common Opaque-Box Test Helpers & Parsers
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

/// Build a sequence of length `len` whose bases differ from `orig` at every position.
fn all_mismatch(orig: &str) -> String {
    orig.bytes().map(alt_base).map(|b| b as char).collect()
}

/// Build a non-seeding homologous flank of length `len`: matches `orig` except a mismatch
/// every `period` bases (so no contiguous run of length `>= period` matches).
fn broken_seed_match(orig: &str, period: usize) -> String {
    assert!(period >= 2);
    orig.bytes()
        .enumerate()
        .map(|(i, b)| {
            if i % period == period - 1 {
                alt_base(b)
            } else {
                b
            }
        })
        .map(|b| b as char)
        .collect()
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

fn write_fasta_gz(path: &Path, records: &[(&str, &str)]) {
    let file = File::create(path).expect("create gz fasta");
    let mut enc = GzEncoder::new(file, Compression::default());
    for (header, seq) in records {
        writeln!(enc, ">{header}").unwrap();
        for chunk in seq.as_bytes().chunks(80) {
            enc.write_all(chunk).unwrap();
            enc.write_all(b"\n").unwrap();
        }
    }
    enc.finish().expect("finish gzip");
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
    let reader = BufReader::with_capacity(256 * 1024, file);
    let mut out = File::create(dst_path).expect("create sliced fasta");
    writeln!(out, ">{header}").unwrap();

    let mut pos = 0usize;
    let mut written = 0usize;
    let mut line = String::new();
    let mut buf = Vec::with_capacity(60);
    let mut reader = reader;

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
// Opaque-Box Binary `.frags` Parser (Big-Endian: 16-byte header + 93 B / frag)
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
    fn canonical_key(&self) -> (char, u64, u64, i64, u64, u64) {
        (
            self.strand,
            self.seq_x,
            self.seq_y,
            self.diag,
            self.x_start,
            self.x_end,
        )
    }

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

fn canonicalize_frags(mut frags: Vec<ParsedFrag>) -> Vec<ParsedFrag> {
    frags.sort_by_key(ParsedFrag::canonical_key);
    frags
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

fn assert_frags_bytes_identical(a: &Path, b: &Path, ctx: &str) {
    let ba = fs::read(a).unwrap_or_else(|e| panic!("{ctx}: read {}: {e}", a.display()));
    let bb = fs::read(b).unwrap_or_else(|e| panic!("{ctx}: read {}: {e}", b.display()));
    if ba != bb {
        let fa = parse_frags_file(a);
        let fb = parse_frags_file(b);
        panic!(
            "{ctx}: .frags files differ ({} bytes / {} frags vs {} bytes / {} frags)",
            ba.len(),
            fa.fragments.len(),
            bb.len(),
            fb.fragments.len()
        );
    }
}

// ============================================================================
// CUDA Static Verification Shared Helpers (OnceLock Cargo JSON + nvcc/cuobjdump)
// ============================================================================

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

fn verify_cuda_extend_zero_stack_for_arch(arch: &str) {
    let tmp = tempfile::tempdir().unwrap();
    let cuda_dir = repo_root().join("crates/gecko-cuda/cuda");
    let extend_cu = cuda_dir.join("cuda_extend.cu");
    let obj_path = tmp.path().join(format!("cuda_extend_{arch}.o"));

    let nvcc_out = Command::new("nvcc")
        .args([
            "-O3",
            "-std=c++17",
            &format!("-I{}", cuda_dir.display()),
            &format!("-arch={arch}"),
            "--ptxas-options=-v",
            "-c",
            extend_cu.to_str().unwrap(),
            "-o",
            obj_path.to_str().unwrap(),
        ])
        .output()
        .expect("failed to execute nvcc");

    assert!(
        nvcc_out.status.success(),
        "nvcc failed for {arch}:\nstderr:\n{}",
        String::from_utf8_lossy(&nvcc_out.stderr)
    );

    // 1. Parse `ptxas` verbose stderr for every entry function
    let ptxas_stderr = String::from_utf8_lossy(&nvcc_out.stderr);
    let mut ptxas_funcs: Vec<(String, u64, u64, u64)> = Vec::new();
    let mut cur_fn: Option<String> = None;
    for line in ptxas_stderr.lines() {
        if let Some(rest) = line.split("Function properties for ").nth(1) {
            cur_fn = Some(rest.trim().to_string());
        } else if line.contains("bytes stack frame") {
            let fn_name = cur_fn
                .take()
                .expect("bytes stack frame line without preceding Function properties line");
            let parts: Vec<&str> = line.split(',').collect();
            assert_eq!(parts.len(), 3, "unexpected ptxas properties line: {line}");
            let parse_leading_u64 = |s: &str| -> u64 {
                s.split_whitespace()
                    .next()
                    .expect("number token")
                    .parse::<u64>()
                    .expect("parse u64")
            };
            let stack_bytes = parse_leading_u64(parts[0]);
            let spill_stores = parse_leading_u64(parts[1]);
            let spill_loads = parse_leading_u64(parts[2]);
            ptxas_funcs.push((fn_name, stack_bytes, spill_stores, spill_loads));
        }
    }
    assert!(
        ptxas_funcs.len() >= 2,
        "{arch}: expected >= 2 entry kernels in ptxas stderr, parsed {}:\n{ptxas_stderr}",
        ptxas_funcs.len()
    );
    assert!(
        ptxas_funcs
            .iter()
            .any(|(n, _, _, _)| n.contains("extend_short_kernel") || n.contains("extend")),
        "{arch}: no extension kernel found in ptxas functions: {ptxas_funcs:?}"
    );
    let mut violations: Vec<String> = Vec::new();
    for (fn_name, stack, stores, loads) in &ptxas_funcs {
        if *stack != 0 {
            violations.push(format!(
                "{arch}: ptxas reported {stack} bytes stack frame for {fn_name} (expected 0)"
            ));
        }
        if *stores != 0 {
            violations.push(format!(
                "{arch}: ptxas reported {stores} bytes spill stores for {fn_name} (expected 0)"
            ));
        }
        if *loads != 0 {
            violations.push(format!(
                "{arch}: ptxas reported {loads} bytes spill loads for {fn_name} (expected 0)"
            ));
        }
    }

    // 2. Parse `cuobjdump --dump-resource-usage`
    let res_out = Command::new("cuobjdump")
        .args(["--dump-resource-usage", obj_path.to_str().unwrap()])
        .output()
        .expect("failed to run cuobjdump --dump-resource-usage");
    assert!(
        res_out.status.success(),
        "cuobjdump --dump-resource-usage failed for {arch}:\n{}",
        String::from_utf8_lossy(&res_out.stderr)
    );
    let res_stdout = String::from_utf8_lossy(&res_out.stdout);
    let mut res_funcs: Vec<(String, u64)> = Vec::new();
    let res_lines: Vec<&str> = res_stdout.lines().collect();
    for i in 0..res_lines.len() {
        let trimmed = res_lines[i].trim();
        if let Some(rest) = trimmed.strip_prefix("Function ")
            && let Some(fn_name) = rest.strip_suffix(':')
        {
            let next_line = res_lines.get(i + 1).expect("missing resource line");
            let stack_tok = next_line
                .split_whitespace()
                .find_map(|tok| tok.strip_prefix("STACK:"))
                .unwrap_or_else(|| panic!("missing STACK token in line: {next_line}"));
            let stack_val: u64 = stack_tok.parse().expect("parse STACK value");
            res_funcs.push((fn_name.to_string(), stack_val));
        }
    }
    assert!(
        res_funcs.len() >= 2,
        "{arch}: expected >= 2 functions in cuobjdump resource usage, got {}:\n{res_stdout}",
        res_funcs.len()
    );
    for (fn_name, stack) in &res_funcs {
        if *stack != 0 {
            violations.push(format!(
                "{arch}: cuobjdump --dump-resource-usage reported STACK:{stack} for {fn_name} (expected STACK:0)"
            ));
        }
    }

    // 3. Parse `cuobjdump --dump-sass` and count STL/LDL opcodes per function
    let sass_out = Command::new("cuobjdump")
        .args(["--dump-sass", obj_path.to_str().unwrap()])
        .output()
        .expect("failed to run cuobjdump --dump-sass");
    assert!(
        sass_out.status.success(),
        "cuobjdump --dump-sass failed for {arch}:\n{}",
        String::from_utf8_lossy(&sass_out.stderr)
    );
    let sass_stdout = String::from_utf8_lossy(&sass_out.stdout);
    // Map function_name -> (total_insns, stl_count, ldl_count)
    let mut sass_stats: BTreeMap<String, (usize, usize, usize)> = BTreeMap::new();
    let mut active_fn: Option<String> = None;
    for line in sass_stdout.lines() {
        let trimmed = line.trim();
        if let Some(rest) = trimmed.strip_prefix("Function :") {
            let name = rest.trim().to_string();
            sass_stats.entry(name.clone()).or_insert((0, 0, 0));
            active_fn = Some(name);
            continue;
        }
        if let Some(ref fn_name) = active_fn
            && trimmed.starts_with("/*")
            && let Some(after_pc) = trimmed.split_once("*/").map(|(_, r)| r)
        {
            let body = after_pc.split(';').next().unwrap_or(after_pc);
            let mut tokens = body.split_whitespace();
            let first = match tokens.next() {
                Some(t) => t,
                None => continue,
            };
            let opcode = if first.starts_with('@') {
                match tokens.next() {
                    Some(t) => t,
                    None => continue,
                }
            } else {
                first
            };
            let base_op = opcode.split('.').next().unwrap_or(opcode);
            let entry = sass_stats.get_mut(fn_name).unwrap();
            entry.0 += 1;
            if base_op == "STL" {
                entry.1 += 1;
            } else if base_op == "LDL" {
                entry.2 += 1;
            }
        }
    }
    assert!(
        sass_stats.len() >= 2,
        "{arch}: expected >= 2 functions in cuobjdump --dump-sass, got {}",
        sass_stats.len()
    );
    for (fn_name, (total, stl, ldl)) in &sass_stats {
        assert!(
            *total > 0,
            "{arch}: parsed 0 SASS instructions for {fn_name}"
        );
        if *stl != 0 || *ldl != 0 {
            violations.push(format!(
                "{arch}: SASS for {fn_name} contains {stl} STL and {ldl} LDL instructions (expected 0/0; total={total})"
            ));
        }
    }

    assert!(
        violations.is_empty(),
        "CUDA zero-stack static verification failed for {arch}:\n{}",
        violations.join("\n")
    );
}

/// Build a synthetic pair of genomes containing:
/// - 2 exact forward HSPs (one 220 bp, one 160 bp with 5% SNPs)
/// - 1 exact inverted HSP (200 bp)
/// - 1 short 3-copy repeat block (80 bp)
fn build_multi_hsp_pair(dir: &Path) -> (PathBuf, PathBuf) {
    let hsp_f1 = pseudo_dna(1001, 220);
    let hsp_f2 = pseudo_dna(1002, 160);
    let hsp_f2_mut = mutate_every(&hsp_f2, 40, 12);
    let hsp_r1 = pseudo_dna(1003, 200);
    let hsp_r1_rc = revcomp(&hsp_r1);
    let rep = pseudo_dna(1004, 80);

    let q_seq = format!(
        "{}{}{}{}{}{}{}{}{}",
        pseudo_dna(1101, 300),
        hsp_f1,
        pseudo_dna(1102, 300),
        hsp_f2,
        pseudo_dna(1103, 300),
        hsp_r1,
        pseudo_dna(1104, 300),
        rep,
        pseudo_dna(1105, 300),
    );
    let r_seq = format!(
        "{}{}{}{}{}{}{}{}{}{}{}{}{}",
        pseudo_dna(1201, 250),
        hsp_f1,
        pseudo_dna(1202, 250),
        rep,
        pseudo_dna(1203, 250),
        hsp_f2_mut,
        pseudo_dna(1204, 250),
        rep,
        pseudo_dna(1205, 250),
        hsp_r1_rc,
        pseudo_dna(1206, 250),
        rep,
        pseudo_dna(1207, 250),
    );

    let q = dir.join("q.fa");
    let r = dir.join("r.fa");
    write_fasta(&q, &[("q1", &q_seq)]);
    write_fasta(&r, &[("r1", &r_seq)]);
    (q, r)
}

// ============================================================================
// TIER 1 — Feature Coverage
// ============================================================================

// --- 1A. Contiguous k-mers (5 tests) ---

#[test]
fn tier1_contiguous_k8_direct_csr() {
    let tmp = tempfile::tempdir().unwrap();
    let (q, r) = build_multi_hsp_pair(tmp.path());
    let out = tmp.path().join("out");

    let res = run_gecko(&[
        "compare",
        q.to_str().unwrap(),
        r.to_str().unwrap(),
        "-k",
        "8",
        "-l",
        "100",
        "-s",
        "70",
        "--freq-limit",
        "2",
        "-f",
        "frags",
        "-o",
        out.to_str().unwrap(),
        "--backend",
        "cpu",
    ]);
    assert_success(&res, "k=8 DirectCsr compare");

    let ff = parse_frags_file(&out.join("q-r.frags"));
    assert_eq!(ff.fragments.len(), 3, "expected 2 fwd + 1 rev HSP");
    assert_eq!(ff.fragments.iter().filter(|f| f.strand == 'f').count(), 2);
    assert_eq!(ff.fragments.iter().filter(|f| f.strand == 'r').count(), 1);
    assert!(
        ff.fragments
            .iter()
            .any(|f| f.strand == 'f' && f.length >= 220)
    );
    assert!(
        ff.fragments
            .iter()
            .any(|f| f.strand == 'r' && f.length >= 200)
    );
}

#[test]
fn tier1_contiguous_k12_default_word_length() {
    let tmp = tempfile::tempdir().unwrap();
    let (q, r) = build_multi_hsp_pair(tmp.path());
    let out_def = tmp.path().join("out_def");
    let out_exp = tmp.path().join("out_exp");

    let res_def = run_gecko(&[
        "compare",
        q.to_str().unwrap(),
        r.to_str().unwrap(),
        "-l",
        "100",
        "-s",
        "70",
        "-f",
        "frags",
        "-o",
        out_def.to_str().unwrap(),
        "--backend",
        "cpu",
    ]);
    assert_success(&res_def, "default k=12 compare");

    let res_exp = run_gecko(&[
        "compare",
        q.to_str().unwrap(),
        r.to_str().unwrap(),
        "-k",
        "12",
        "-l",
        "100",
        "-s",
        "70",
        "-f",
        "frags",
        "-o",
        out_exp.to_str().unwrap(),
        "--backend",
        "cpu",
    ]);
    assert_success(&res_exp, "explicit k=12 compare");

    assert_frags_bytes_identical(
        &out_def.join("q-r.frags"),
        &out_exp.join("q-r.frags"),
        "default -k vs explicit -k 12",
    );
    let ff = parse_frags_file(&out_def.join("q-r.frags"));
    assert_eq!(ff.fragments.len(), 3);
}

#[test]
fn tier1_contiguous_k16_hash_csr() {
    let tmp = tempfile::tempdir().unwrap();
    let (q, r) = build_multi_hsp_pair(tmp.path());
    let out = tmp.path().join("out");

    let res = run_gecko(&[
        "compare",
        q.to_str().unwrap(),
        r.to_str().unwrap(),
        "-k",
        "16",
        "-l",
        "100",
        "-s",
        "70",
        "--freq-limit",
        "2",
        "-f",
        "frags",
        "-o",
        out.to_str().unwrap(),
        "--backend",
        "cpu",
    ]);
    assert_success(&res, "k=16 HashCsr compare");

    let ff = parse_frags_file(&out.join("q-r.frags"));
    assert_eq!(ff.fragments.len(), 3);
    for f in &ff.fragments {
        assert_eq!(f.x_end - f.x_start + 1, f.length);
        assert_eq!(f.y_end - f.y_start + 1, f.length);
        assert!(f.length > 100);
        assert!(f.similarity > 70.0);
    }
}

#[test]
fn tier1_contiguous_k20_longer_seed() {
    let tmp = tempfile::tempdir().unwrap();
    // Block 1: 100 bp with mismatches every 19 bp -> max exact match is 18 bp (< 20), NOT seeded!
    // Block 2: 100 bp exact match -> seeded at k=20!
    let b1 = pseudo_dna(2001, 100);
    let b1_broken = broken_seed_match(&b1, 19);
    let b2 = pseudo_dna(2002, 100);

    let f_l = pseudo_dna(2010, 120);
    let f_m = pseudo_dna(2011, 120);
    let f_r = pseudo_dna(2012, 120);
    let q = tmp.path().join("q.fa");
    let r = tmp.path().join("r.fa");
    write_fasta(&q, &[("q", &format!("{f_l}{b1}{f_m}{b2}{f_r}"))]);
    // Place b1_broken and b2 on distinct diagonals with all_mismatch boundary guards around b2
    write_fasta(
        &r,
        &[(
            "r",
            &format!(
                "{}{}{}{}{}{}",
                pseudo_dna(2020, 120),
                b1_broken,
                pseudo_dna(2021, 200),
                all_mismatch(&f_m),
                b2,
                all_mismatch(&f_r)
            ),
        )],
    );

    let out = tmp.path().join("out");
    let res = run_gecko(&[
        "compare",
        q.to_str().unwrap(),
        r.to_str().unwrap(),
        "-k",
        "20",
        "-l",
        "50",
        "-s",
        "60",
        "-f",
        "frags",
        "-o",
        out.to_str().unwrap(),
        "--backend",
        "cpu",
    ]);
    assert_success(&res, "k=20 compare");
    let ff = parse_frags_file(&out.join("q-r.frags"));
    assert_eq!(
        ff.fragments.len(),
        1,
        "only the 20+-mer exact block should seed"
    );
    assert_eq!(ff.fragments[0].x_start, 120 + 100 + 120);
    assert_eq!(ff.fragments[0].length, 100);
}

#[test]
fn tier1_contiguous_k32_max_word_length() {
    let tmp = tempfile::tempdir().unwrap();
    let b_short = pseudo_dna(3001, 120);
    let b_short_broken = broken_seed_match(&b_short, 31); // max run 30 < 32
    let b_exact = pseudo_dna(3002, 120);

    let f_l = pseudo_dna(3010, 150);
    let f_m = pseudo_dna(3011, 150);
    let f_r = pseudo_dna(3012, 150);
    let q = tmp.path().join("q.fa");
    let r = tmp.path().join("r.fa");
    write_fasta(&q, &[("q", &format!("{f_l}{b_short}{f_m}{b_exact}{f_r}"))]);
    write_fasta(
        &r,
        &[(
            "r",
            &format!(
                "{}{}{}{}{}{}",
                pseudo_dna(3020, 150),
                b_short_broken,
                pseudo_dna(3021, 220),
                all_mismatch(&f_m),
                b_exact,
                all_mismatch(&f_r)
            ),
        )],
    );

    let out = tmp.path().join("out");
    let res = run_gecko(&[
        "compare",
        q.to_str().unwrap(),
        r.to_str().unwrap(),
        "-k",
        "32",
        "-l",
        "60",
        "-s",
        "60",
        "-f",
        "frags",
        "-o",
        out.to_str().unwrap(),
        "--backend",
        "cpu",
    ]);
    assert_success(&res, "k=32 compare");
    let ff = parse_frags_file(&out.join("q-r.frags"));
    assert_eq!(ff.fragments.len(), 1, "only >=32bp exact block should seed");
    assert_eq!(ff.fragments[0].length, 120);
    assert_eq!(ff.fragments[0].ident, 120);
}

// --- 1B. Spaced seeds (5 tests) ---

#[test]
fn tier1_spaced_seed_preset_w12s18() {
    let tmp = tempfile::tempdir().unwrap();
    // Divergent HSP: mutate every 9th base -> max contiguous match is 8 bp (<12),
    // so contiguous -k 12 misses it completely, while w12s18 ("111010110011010111") hits it!
    // Put mutations specifically at don't-care offsets of w12s18 (indices 3, 8, 9, 14) repeated
    // every 18 bp: max contiguous match is 8 bp, identity is 14/18 = 77.8%, sim = 55.6% (>40%).
    let hsp = pseudo_dna(4001, 180);
    let mut mut_pos = Vec::new();
    for blk in 0..10 {
        for &dc in &[3usize, 9, 14] {
            mut_pos.push(blk * 18 + dc);
        }
    }
    let hsp_div = mutate_positions(&hsp, &mut_pos);

    let q = tmp.path().join("q.fa");
    let r = tmp.path().join("r.fa");
    write_fasta(
        &q,
        &[(
            "q",
            &format!("{}{}{}", pseudo_dna(4010, 150), hsp, pseudo_dna(4011, 150)),
        )],
    );
    write_fasta(
        &r,
        &[(
            "r",
            &format!(
                "{}{}{}",
                pseudo_dna(4020, 150),
                hsp_div,
                pseudo_dna(4021, 150)
            ),
        )],
    );

    let out_contig = tmp.path().join("contig");
    let out_spaced = tmp.path().join("spaced");
    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "12",
            "-l",
            "80",
            "-s",
            "40",
            "-f",
            "frags",
            "-o",
            out_contig.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "contiguous k=12 on divergent HSP",
    );
    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "--spaced-seed",
            "w12s18",
            "-l",
            "80",
            "-s",
            "40",
            "-f",
            "frags",
            "-o",
            out_spaced.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "spaced w12s18 on divergent HSP",
    );

    assert_eq!(
        parse_frags_file(&out_contig.join("q-r.frags"))
            .fragments
            .len(),
        0
    );
    let ff = parse_frags_file(&out_spaced.join("q-r.frags"));
    assert_eq!(ff.fragments.len(), 1);
    assert!(ff.fragments[0].length >= 160);
}

#[test]
fn tier1_spaced_seed_preset_w14s20_case_insensitive() {
    let tmp = tempfile::tempdir().unwrap();
    // w14s20 = "11110101100110101111" -> don't-care at 4, 6, 9, 10, 13, 15
    let hsp = pseudo_dna(4101, 200);
    let mut mut_pos = Vec::new();
    for blk in 0..10 {
        for &dc in &[4usize, 10, 15] {
            mut_pos.push(blk * 20 + dc);
        }
    }
    let hsp_div = mutate_positions(&hsp, &mut_pos);

    let q = tmp.path().join("q.fa");
    let r = tmp.path().join("r.fa");
    write_fasta(
        &q,
        &[(
            "q",
            &format!("{}{}{}", pseudo_dna(4110, 120), hsp, pseudo_dna(4111, 120)),
        )],
    );
    write_fasta(
        &r,
        &[(
            "r",
            &format!(
                "{}{}{}",
                pseudo_dna(4120, 120),
                hsp_div,
                pseudo_dna(4121, 120)
            ),
        )],
    );

    let out_lo = tmp.path().join("lo");
    let out_hi = tmp.path().join("hi");
    assert_success(
        &run_gecko(&[
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
            "-o",
            out_lo.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "w14s20 lowercase",
    );
    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "--spaced-seed",
            "W14S20",
            "-l",
            "80",
            "-s",
            "50",
            "-f",
            "frags",
            "-o",
            out_hi.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "W14S20 uppercase",
    );

    assert_frags_bytes_identical(
        &out_lo.join("q-r.frags"),
        &out_hi.join("q-r.frags"),
        "w14s20 vs W14S20",
    );
    let ff = parse_frags_file(&out_lo.join("q-r.frags"));
    assert_eq!(ff.fragments.len(), 1);
    assert!(ff.fragments[0].length >= 180);
}

#[test]
fn tier1_spaced_seed_preset_w16s22() {
    let tmp = tempfile::tempdir().unwrap();
    // w16s22 = "1111011010110101101111" -> don't-care at 4, 7, 9, 12, 14, 17
    let hsp = pseudo_dna(4201, 220);
    let mut mut_pos = Vec::new();
    for blk in 0..10 {
        for &dc in &[4usize, 12, 17] {
            mut_pos.push(blk * 22 + dc);
        }
    }
    let hsp_div = mutate_positions(&hsp, &mut_pos);

    let q = tmp.path().join("q.fa");
    let r = tmp.path().join("r.fa");
    write_fasta(
        &q,
        &[(
            "q",
            &format!("{}{}{}", pseudo_dna(4210, 120), hsp, pseudo_dna(4211, 120)),
        )],
    );
    write_fasta(
        &r,
        &[(
            "r",
            &format!(
                "{}{}{}",
                pseudo_dna(4220, 120),
                hsp_div,
                pseudo_dna(4221, 120)
            ),
        )],
    );

    let out = tmp.path().join("out");
    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "--spaced-seed",
            "w16s22",
            "-l",
            "80",
            "-s",
            "50",
            "-f",
            "frags",
            "-o",
            out.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "w16s22 compare",
    );
    let ff = parse_frags_file(&out.join("q-r.frags"));
    assert_eq!(ff.fragments.len(), 1);
    assert!(ff.fragments[0].length >= 200);
}

#[test]
fn tier1_spaced_seed_custom_palindromic_syntax_parity() {
    let tmp = tempfile::tempdir().unwrap();
    let (q, r) = build_multi_hsp_pair(tmp.path());

    let out_preset = tmp.path().join("preset");
    let out_01 = tmp.path().join("bin01");
    let out_sym = tmp.path().join("sym");

    for (spec, outdir) in [
        ("w12s18", &out_preset),
        ("111010110011010111", &out_01),
        ("###_#_##--##_#_###", &out_sym),
    ] {
        assert_success(
            &run_gecko(&[
                "compare",
                q.to_str().unwrap(),
                r.to_str().unwrap(),
                "--spaced-seed",
                spec,
                "-l",
                "100",
                "-s",
                "70",
                "-f",
                "frags",
                "-o",
                outdir.to_str().unwrap(),
                "--backend",
                "cpu",
            ]),
            spec,
        );
    }

    assert_frags_bytes_identical(
        &out_preset.join("q-r.frags"),
        &out_01.join("q-r.frags"),
        "w12s18 vs 0/1 mask",
    );
    assert_frags_bytes_identical(
        &out_preset.join("q-r.frags"),
        &out_sym.join("q-r.frags"),
        "w12s18 vs #/_/- mask",
    );
    assert_eq!(
        parse_frags_file(&out_preset.join("q-r.frags"))
            .fragments
            .len(),
        3
    );
}

#[test]
fn tier1_spaced_seed_all_ones_matches_contiguous_and_overrides_k() {
    let tmp = tempfile::tempdir().unwrap();
    let (q, r) = build_multi_hsp_pair(tmp.path());
    let out_contig = tmp.path().join("contig");
    let out_all1 = tmp.path().join("all1");

    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "12",
            "-l",
            "100",
            "-s",
            "70",
            "-f",
            "frags",
            "-o",
            out_contig.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "contiguous k=12",
    );
    // Explicitly pass -k 8 alongside 12-bp all-'1' mask to verify --spaced-seed overrides -k
    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "8",
            "--spaced-seed",
            "111111111111",
            "-l",
            "100",
            "-s",
            "70",
            "-f",
            "frags",
            "-o",
            out_all1.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "all-ones spaced seed overriding -k 8",
    );

    assert_frags_bytes_identical(
        &out_contig.join("q-r.frags"),
        &out_all1.join("q-r.frags"),
        "all-1s spaced seed vs contiguous -k 12",
    );
}

// --- 1C. Forward-only / strands / soft-masking / gzip (5 tests) ---

#[test]
fn tier1_forward_only_bit_exact_subset_and_zero_reverse() {
    let tmp = tempfile::tempdir().unwrap();
    let (q, r) = build_multi_hsp_pair(tmp.path());
    let out_both = tmp.path().join("both");
    let out_fwd = tmp.path().join("fwd");

    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "12",
            "-l",
            "100",
            "-s",
            "70",
            "-f",
            "frags",
            "-o",
            out_both.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "dual strand",
    );
    let res_fwd = run_gecko(&[
        "compare",
        q.to_str().unwrap(),
        r.to_str().unwrap(),
        "-k",
        "12",
        "-l",
        "100",
        "-s",
        "70",
        "--forward-only",
        "--json-stats",
        "-f",
        "frags",
        "-o",
        out_fwd.to_str().unwrap(),
        "--backend",
        "cpu",
    ]);
    assert_success(&res_fwd, "forward-only");

    let stats: serde_json::Value = serde_json::from_slice(&res_fwd.stdout).expect("valid json");
    assert_eq!(stats["reverse_fragments"], 0);
    assert_eq!(stats["forward_fragments"], 2);

    let both_f: Vec<ParsedFrag> = parse_frags_file(&out_both.join("q-r.frags"))
        .fragments
        .into_iter()
        .filter(|f| f.strand == 'f')
        .collect();
    let fwd_all = parse_frags_file(&out_fwd.join("q-r.frags")).fragments;
    assert_eq!(
        fwd_all, both_f,
        "--forward-only must equal exact 'f' subset"
    );
}

#[test]
fn tier1_both_strands_forward_and_inverted_hsps() {
    let tmp = tempfile::tempdir().unwrap();
    let (q, r) = build_multi_hsp_pair(tmp.path());
    let out = tmp.path().join("out");

    let res = run_gecko(&[
        "compare",
        q.to_str().unwrap(),
        r.to_str().unwrap(),
        "-k",
        "12",
        "-l",
        "100",
        "-s",
        "70",
        "--json-stats",
        "-f",
        "frags",
        "-o",
        out.to_str().unwrap(),
        "--backend",
        "cpu",
    ]);
    assert_success(&res, "both strands");
    let stats: serde_json::Value = serde_json::from_slice(&res.stdout).unwrap();
    assert_eq!(stats["fragments"], 3);
    assert_eq!(stats["forward_fragments"], 2);
    assert_eq!(stats["reverse_fragments"], 1);
    assert!(
        stats["backend"].as_str().unwrap().starts_with("cpu"),
        "expected cpu backend, got {}",
        stats["backend"]
    );

    let ff = parse_frags_file(&out.join("q-r.frags"));
    assert_eq!(ff.fragments.len(), 3);
    // Forward precedes reverse
    assert_eq!(ff.fragments[0].strand, 'f');
    assert_eq!(ff.fragments[1].strand, 'f');
    assert_eq!(ff.fragments[2].strand, 'r');
    // Reverse fragment still has y_start <= y_end in RC coordinate space
    assert!(ff.fragments[2].y_start <= ff.fragments[2].y_end);
}

#[test]
fn tier1_soft_mask_blocks_seeding_but_allows_extension() {
    let tmp = tempfile::tempdir().unwrap();
    // HSP1: 120 bp completely lowercase -> blocked from seeding under default soft_mask.
    // HSP2: 50 bp lowercase + 20 bp UPPERCASE anchor + 50 bp lowercase -> seeds on anchor
    //       and extends cleanly across BOTH lowercase flanks to full 120 bp!
    let hsp1 = pseudo_dna(5001, 120).to_ascii_lowercase();
    let hsp2_left = pseudo_dna(5002, 50).to_ascii_lowercase();
    let hsp2_mid = pseudo_dna(5003, 20); // uppercase anchor >= k=12
    let hsp2_right = pseudo_dna(5004, 50).to_ascii_lowercase();
    let hsp2 = format!("{hsp2_left}{hsp2_mid}{hsp2_right}");

    let f_l = pseudo_dna(5010, 100);
    let f_m = pseudo_dna(5011, 100);
    let f_r = pseudo_dna(5012, 100);
    let q = tmp.path().join("q.fa");
    let r = tmp.path().join("r.fa");
    write_fasta(&q, &[("q", &format!("{f_l}{hsp1}{f_m}{hsp2}{f_r}"))]);
    // Place hsp1 and hsp2 on distinct diagonals with all_mismatch guards around hsp2
    write_fasta(
        &r,
        &[(
            "r",
            &format!(
                "{}{}{}{}{}{}",
                pseudo_dna(5020, 100),
                hsp1,
                pseudo_dna(5021, 200),
                all_mismatch(&f_m),
                hsp2,
                all_mismatch(&f_r)
            ),
        )],
    );

    let out = tmp.path().join("out");
    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "12",
            "-l",
            "80",
            "-s",
            "90",
            "-f",
            "frags",
            "-o",
            out.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "soft-mask extension",
    );

    let ff = parse_frags_file(&out.join("q-r.frags"));
    assert_eq!(ff.fragments.len(), 1, "only anchored HSP2 should seed");
    assert!(
        ff.fragments[0].length >= 120,
        "extension must cross lowercase flanks (got len {})",
        ff.fragments[0].length
    );
    assert_eq!(ff.fragments[0].ident, ff.fragments[0].length);
}

#[test]
fn tier1_no_soft_mask_matches_uppercase_bit_exact() {
    let tmp = tempfile::tempdir().unwrap();
    let hsp = pseudo_dna(5101, 160);
    let q_up = format!("{}{}{}", pseudo_dna(5110, 120), hsp, pseudo_dna(5111, 120));
    let r_up = format!("{}{}{}", pseudo_dna(5120, 120), hsp, pseudo_dna(5121, 120));

    let q_up_p = tmp.path().join("q_up.fa");
    let r_up_p = tmp.path().join("r_up.fa");
    let q_lo_p = tmp.path().join("q_lo.fa");
    let r_lo_p = tmp.path().join("r_lo.fa");
    write_fasta(&q_up_p, &[("q", &q_up)]);
    write_fasta(&r_up_p, &[("r", &r_up)]);
    write_fasta(&q_lo_p, &[("q", &q_up.to_ascii_lowercase())]);
    write_fasta(&r_lo_p, &[("r", &r_up.to_ascii_lowercase())]);

    let out_up = tmp.path().join("up");
    let out_lo = tmp.path().join("lo");
    assert_success(
        &run_gecko(&[
            "compare",
            q_up_p.to_str().unwrap(),
            r_up_p.to_str().unwrap(),
            "-k",
            "12",
            "-l",
            "80",
            "-s",
            "80",
            "--prefix",
            "cmp",
            "-f",
            "frags",
            "-o",
            out_up.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "uppercase",
    );
    assert_success(
        &run_gecko(&[
            "compare",
            q_lo_p.to_str().unwrap(),
            r_lo_p.to_str().unwrap(),
            "-k",
            "12",
            "-l",
            "80",
            "-s",
            "80",
            "--no-soft-mask",
            "--prefix",
            "cmp",
            "-f",
            "frags",
            "-o",
            out_lo.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "lowercase --no-soft-mask",
    );

    assert_frags_bytes_identical(
        &out_up.join("cmp.frags"),
        &out_lo.join("cmp.frags"),
        "--no-soft-mask lowercase vs uppercase",
    );
    assert_eq!(
        parse_frags_file(&out_up.join("cmp.frags")).fragments.len(),
        1
    );
}

#[test]
fn tier1_gzip_input_parity() {
    let tmp = tempfile::tempdir().unwrap();
    let hsp_f = pseudo_dna(5201, 180);
    let hsp_r = pseudo_dna(5202, 160);
    let q_seq = format!(
        "{}{}{}{}{}",
        pseudo_dna(5210, 120),
        hsp_f,
        pseudo_dna(5211, 120),
        hsp_r,
        pseudo_dna(5212, 120)
    );
    let r_seq = format!(
        "{}{}{}{}{}",
        pseudo_dna(5220, 120),
        hsp_f,
        pseudo_dna(5221, 120),
        revcomp(&hsp_r),
        pseudo_dna(5222, 120)
    );

    let q_fa = tmp.path().join("q.fa");
    let r_fa = tmp.path().join("r.fa");
    let q_gz = tmp.path().join("q.fa.gz");
    let r_gz = tmp.path().join("r.fa.gz");
    write_fasta(&q_fa, &[("q", &q_seq)]);
    write_fasta(&r_fa, &[("r", &r_seq)]);
    write_fasta_gz(&q_gz, &[("q", &q_seq)]);
    write_fasta_gz(&r_gz, &[("r", &r_seq)]);

    let out_plain = tmp.path().join("plain");
    let out_gz = tmp.path().join("gz");
    assert_success(
        &run_gecko(&[
            "compare",
            q_fa.to_str().unwrap(),
            r_fa.to_str().unwrap(),
            "-k",
            "12",
            "-l",
            "80",
            "-s",
            "80",
            "-f",
            "frags",
            "-o",
            out_plain.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "plain fasta",
    );
    assert_success(
        &run_gecko(&[
            "compare",
            q_gz.to_str().unwrap(),
            r_gz.to_str().unwrap(),
            "-k",
            "12",
            "-l",
            "80",
            "-s",
            "80",
            "-f",
            "frags",
            "-o",
            out_gz.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "gzipped fasta",
    );

    assert_frags_bytes_identical(
        &out_plain.join("q-r.frags"),
        &out_gz.join("q-r.frags"),
        "plain vs .fa.gz",
    );
    assert_eq!(
        parse_frags_file(&out_gz.join("q-r.frags")).fragments.len(),
        2
    );
}

// --- 1D. `freq_limit` & multi-slab / `-m` budget parity (5 tests) ---

#[test]
fn tier1_freq_limit_per_strand_orientation() {
    let tmp = tempfile::tempdir().unwrap();
    // Non-palindromic block W appears 4x forward in ref, while rc(W) appears 1x in ref.
    // Query contains W once.
    // With --freq-limit 2: all 4 forward copies of W in ref have cnt_fwd=4 > 2 (DROPPED),
    // while the 1 copy of rc(W) has cnt_rev=1 <= 2 (KEPT -> produces 1 reverse fragment!).
    let w = pseudo_dna(6001, 100);
    let w_rc = revcomp(&w);
    assert_ne!(w, w_rc);

    let ql = pseudo_dna(6010, 150);
    let qr = pseudo_dna(6011, 150);
    let ql_mm = all_mismatch(&ql[ql.len() - 40..]);
    let qr_mm = all_mismatch(&qr[..40]);
    let q_seq = format!("{ql}{w}{qr}");
    // Guard every copy of w with all_mismatch of q's adjacent flanks so no boundary k-mer
    // straddles outside w, and guard w_rc with all_mismatch of revcomp(qr/ql) so length == 100.
    let rc_left_mm = all_mismatch(&revcomp(&qr[..40]));
    let rc_right_mm = all_mismatch(&revcomp(&ql[ql.len() - 40..]));
    let r_seq = format!(
        "{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}",
        pseudo_dna(6020, 80),
        ql_mm,
        w,
        qr_mm,
        pseudo_dna(6021, 80),
        ql_mm,
        w,
        qr_mm,
        pseudo_dna(6022, 80),
        ql_mm,
        w,
        qr_mm,
        pseudo_dna(6023, 80),
        ql_mm,
        w,
        qr_mm,
        pseudo_dna(6024, 80),
        format_args!("{rc_left_mm}{w_rc}{rc_right_mm}"),
        pseudo_dna(6025, 80),
    );

    let q = tmp.path().join("q.fa");
    let r = tmp.path().join("r.fa");
    write_fasta(&q, &[("q", &q_seq)]);
    write_fasta(&r, &[("r", &r_seq)]);

    let out = tmp.path().join("out");
    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "12",
            "-l",
            "60",
            "-s",
            "80",
            "--freq-limit",
            "2",
            "-f",
            "frags",
            "-o",
            out.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "per-strand freq_limit",
    );

    let ff = parse_frags_file(&out.join("q-r.frags"));
    assert_eq!(
        ff.fragments.len(),
        1,
        "expected 4 forward copies dropped and 1 reverse copy kept, got {:?}",
        ff.fragments
    );
    assert_eq!(ff.fragments[0].strand, 'r');
    assert_eq!(ff.fragments[0].length, 100);
}

#[test]
fn tier1_freq_limit_repeat_filtering_preserves_unique_hsps() {
    let tmp = tempfile::tempdir().unwrap();
    let rep = pseudo_dna(6101, 90);
    let uniq = pseudo_dna(6102, 140);

    let ql = pseudo_dna(6110, 120);
    let qm = pseudo_dna(6111, 120);
    let qr = pseudo_dna(6112, 120);
    let rep_l_mm = all_mismatch(&ql[ql.len() - 40..]);
    let rep_r_mm = all_mismatch(&qm[..40]);
    let uniq_l_mm = all_mismatch(&qm[qm.len() - 40..]);
    let uniq_r_mm = all_mismatch(&qr[..40]);

    let q_seq = format!("{ql}{rep}{qm}{uniq}{qr}");
    let mut r_seq = pseudo_dna(6120, 120);
    for i in 0..5 {
        r_seq.push_str(&rep_l_mm);
        r_seq.push_str(&rep);
        r_seq.push_str(&rep_r_mm);
        // Distinct spacer lengths so no repeat copy shares a diagonal with uniq
        r_seq.push_str(&pseudo_dna(6121 + i, 90 + (i as usize) * 17));
    }
    r_seq.push_str(&uniq_l_mm);
    r_seq.push_str(&uniq);
    r_seq.push_str(&uniq_r_mm);
    r_seq.push_str(&pseudo_dna(6130, 120));

    let q = tmp.path().join("q.fa");
    let r = tmp.path().join("r.fa");
    write_fasta(&q, &[("q", &q_seq)]);
    write_fasta(&r, &[("r", &r_seq)]);

    let out_hi = tmp.path().join("hi");
    let out_lo = tmp.path().join("lo");
    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "12",
            "-l",
            "60",
            "-s",
            "80",
            "--freq-limit",
            "10",
            "-f",
            "frags",
            "-o",
            out_hi.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "freq-limit 10",
    );
    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "12",
            "-l",
            "60",
            "-s",
            "80",
            "--freq-limit",
            "2",
            "-f",
            "frags",
            "-o",
            out_lo.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "freq-limit 2",
    );

    let f_hi = parse_frags_file(&out_hi.join("q-r.frags")).fragments;
    let f_lo = parse_frags_file(&out_lo.join("q-r.frags")).fragments;
    assert_eq!(f_hi.len(), 6, "5 repeat copies + 1 unique");
    assert_eq!(f_lo.len(), 1, "only unique HSP preserved");
    assert_eq!(f_lo[0].length, 140);
    assert!(f_hi.contains(&f_lo[0]), "unique HSP must be bit-identical");
}

#[test]
fn tier1_tiny_memory_ondisk_spill_canonical_parity() {
    let tmp = tempfile::tempdir().unwrap();
    // Embed a 6,000-bp forward HSP + 250-bp reverse HSP so raw_hits > 4,096 (the minimum
    // in-memory Hit buffer in hits.rs), guaranteeing spilled_runs > 0 under -m 64K!
    let hsp_f = pseudo_dna(6201, 6_000);
    let hsp_r = pseudo_dna(6202, 250);
    let q_seq = format!(
        "{}{}{}{}{}",
        pseudo_dna(6210, 8_000),
        hsp_f,
        pseudo_dna(6211, 8_000),
        hsp_r,
        pseudo_dna(6212, 8_000)
    );
    let r_seq = format!(
        "{}{}{}{}{}",
        pseudo_dna(6220, 8_000),
        hsp_f,
        pseudo_dna(6221, 9_500),
        revcomp(&hsp_r),
        pseudo_dna(6222, 8_000)
    );

    let q = tmp.path().join("q.fa");
    let r = tmp.path().join("r.fa");
    write_fasta(&q, &[("q", &q_seq)]);
    write_fasta(&r, &[("r", &r_seq)]);

    let out_2g = tmp.path().join("m2g");
    let out_64k = tmp.path().join("m64k");

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
            "80",
            "-m",
            "2G",
            "-f",
            "frags",
            "-o",
            out_2g.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "-m 2G baseline",
    );
    let res_64k = run_gecko(&[
        "compare",
        q.to_str().unwrap(),
        r.to_str().unwrap(),
        "-k",
        "32",
        "--seeding",
        "exact",
        "-l",
        "100",
        "-s",
        "80",
        "-m",
        "64K",
        "--json-stats",
        "-f",
        "frags",
        "-o",
        out_64k.to_str().unwrap(),
        "--backend",
        "cpu",
    ]);
    assert_success(&res_64k, "-m 64K on-disk spill");

    let stats: serde_json::Value = serde_json::from_slice(&res_64k.stdout).unwrap();
    let spilled = stats["spilled_runs"].as_u64().unwrap_or(0);
    assert!(
        spilled > 0,
        "expected spilled_runs > 0 under -m 64K, got {spilled} (stats={stats})"
    );

    let f_2g = canonicalize_frags(parse_frags_file(&out_2g.join("q-r.frags")).fragments);
    let f_64k = canonicalize_frags(parse_frags_file(&out_64k.join("q-r.frags")).fragments);
    assert_eq!(f_2g.len(), 2);
    assert_eq!(f_2g, f_64k, "OnDisk spill must match single-slab streaming");
}

#[test]
fn tier1_cache_dir_gidx_canonical_parity() {
    let tmp = tempfile::tempdir().unwrap();
    let hsp_f = pseudo_dna(6301, 200);
    let hsp_r = pseudo_dna(6302, 180);
    let q_seq = format!(
        "{}{}{}{}{}",
        pseudo_dna(6310, 400),
        hsp_f,
        pseudo_dna(6311, 400),
        hsp_r,
        pseudo_dna(6312, 400)
    );
    let r_seq = format!(
        "{}{}{}{}{}",
        pseudo_dna(6320, 400),
        hsp_f,
        pseudo_dna(6321, 400),
        revcomp(&hsp_r),
        pseudo_dna(6322, 400)
    );
    let q = tmp.path().join("q.fa");
    let r = tmp.path().join("r.fa");
    write_fasta(&q, &[("q", &q_seq)]);
    write_fasta(&r, &[("r", &r_seq)]);

    let cache = tmp.path().join("cache");
    let out_stream = tmp.path().join("stream");
    let out_c1 = tmp.path().join("c1");
    let out_c2 = tmp.path().join("c2");

    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "16",
            "-l",
            "80",
            "-s",
            "80",
            "-f",
            "frags",
            "-o",
            out_stream.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "streaming",
    );
    for (out, lbl) in [(&out_c1, "cache run 1"), (&out_c2, "cache run 2")] {
        assert_success(
            &run_gecko(&[
                "compare",
                q.to_str().unwrap(),
                r.to_str().unwrap(),
                "-k",
                "16",
                "-l",
                "80",
                "-s",
                "80",
                "--cache-dir",
                cache.to_str().unwrap(),
                "-f",
                "frags",
                "-o",
                out.to_str().unwrap(),
                "--backend",
                "cpu",
            ]),
            lbl,
        );
    }

    // Verify at least one .gidx file exists with magic b"GECKOIDX"
    let gidx_files: Vec<PathBuf> = fs::read_dir(&cache)
        .unwrap()
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| p.extension().and_then(|s| s.to_str()) == Some("gidx"))
        .collect();
    assert!(!gidx_files.is_empty(), "expected .gidx files in cache dir");
    for gf in &gidx_files {
        let head = &fs::read(gf).unwrap()[..8];
        assert_eq!(head, b"GECKOIDX", "invalid .gidx magic in {}", gf.display());
    }

    let f_s = canonicalize_frags(parse_frags_file(&out_stream.join("q-r.frags")).fragments);
    let f_c1 = canonicalize_frags(parse_frags_file(&out_c1.join("q-r.frags")).fragments);
    let f_c2 = canonicalize_frags(parse_frags_file(&out_c2.join("q-r.frags")).fragments);
    assert_eq!(f_s.len(), 2);
    assert_eq!(f_s, f_c1);
    assert_eq!(f_c1, f_c2);
}

#[test]
fn tier1_spaced_seed_tight_memory_multi_slab_parity() {
    let tmp = tempfile::tempdir().unwrap();
    // ~4.10 Mb synthetic reference (> 4_000_000 k-mers) + tiny 3 Kb query.
    // With `--spaced-seed w14s20 -m 1M`, Branch 1 & 2 fail and Branch 3 runs multi-slab
    // partition streaming with max_kmers = 4_000_000 (2 slabs: [0..128) and [128..256)).
    let hsp_f = pseudo_dna(6401, 200);
    let hsp_r = pseudo_dna(6402, 180);
    let q_seq = format!(
        "{}{}{}{}{}",
        pseudo_dna(6410, 1_000),
        hsp_f,
        pseudo_dna(6411, 1_000),
        hsp_r,
        pseudo_dna(6412, 1_000)
    );
    let mut r_seq = String::with_capacity(4_105_000);
    r_seq.push_str(&pseudo_dna(6420, 1_500_000));
    r_seq.push_str(&hsp_f);
    r_seq.push_str(&pseudo_dna(6421, 1_500_000));
    r_seq.push_str(&revcomp(&hsp_r));
    r_seq.push_str(&pseudo_dna(6422, 1_100_000));

    let q = tmp.path().join("q.fa");
    let r = tmp.path().join("r.fa");
    write_fasta(&q, &[("q", &q_seq)]);
    write_fasta(&r, &[("r", &r_seq)]);

    let out_2g = tmp.path().join("m2g");
    let out_1m = tmp.path().join("m1m");

    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "--spaced-seed",
            "w14s20",
            "-l",
            "100",
            "-s",
            "80",
            "-m",
            "2G",
            "-f",
            "frags",
            "-o",
            out_2g.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "spaced w14s20 -m 2G",
    );
    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "--spaced-seed",
            "w14s20",
            "-l",
            "100",
            "-s",
            "80",
            "-m",
            "1M",
            "-f",
            "frags",
            "-o",
            out_1m.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "spaced w14s20 -m 1M multi-slab",
    );

    assert_frags_bytes_identical(
        &out_2g.join("q-r.frags"),
        &out_1m.join("q-r.frags"),
        "spaced-seed 2-slab (-m 1M) vs single-slab (-m 2G)",
    );
    assert_eq!(
        parse_frags_file(&out_1m.join("q-r.frags")).fragments.len(),
        2
    );
}

// --- 1E. Multi-FASTA & `N`/IUPAC (5 tests) ---

#[test]
fn tier1_multifasta_seq_indices_and_local_coords_across_formats() {
    let tmp = tempfile::tempdir().unwrap();
    let hsp01 = pseudo_dna(7001, 140); // q0 -> r1 (forward)
    let hsp20 = pseudo_dna(7002, 150); // q2 -> r0 (forward)

    let q0_l = pseudo_dna(7010, 80);
    let q0_r = pseudo_dna(7011, 80);
    let q2_l = pseudo_dna(7013, 60);
    let q2_r = pseudo_dna(7014, 90);
    let q0 = format!("{q0_l}{hsp01}{q0_r}");
    let q1 = pseudo_dna(7012, 200);
    let q2 = format!("{q2_l}{hsp20}{q2_r}");

    // Guard immediate flanks with all_mismatch so coordinates are exact without random edge bleed
    let r0 = format!(
        "{}{}{}{}{}",
        pseudo_dna(7020, 30),
        all_mismatch(&q2_l[q2_l.len() - 40..]),
        hsp20,
        all_mismatch(&q2_r[..40]),
        pseudo_dna(7021, 40)
    );
    let r1 = format!(
        "{}{}{}{}{}",
        pseudo_dna(7022, 10),
        all_mismatch(&q0_l[q0_l.len() - 40..]),
        hsp01,
        all_mismatch(&q0_r[..40]),
        pseudo_dna(7023, 60)
    );
    let r2 = pseudo_dna(7024, 180);

    let q = tmp.path().join("q.fa");
    let r = tmp.path().join("r.fa");
    write_fasta(
        &q,
        &[("q_chr0 desc", &q0), ("q_chr1", &q1), ("q_chr2", &q2)],
    );
    write_fasta(&r, &[("r_chr0", &r0), ("r_chr1", &r1), ("r_chr2", &r2)]);

    let out = tmp.path().join("out");
    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "12",
            "-l",
            "100",
            "-s",
            "80",
            "-f",
            "frags,paf,jsonl,bedpe,blast",
            "-o",
            out.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "multi-FASTA compare",
    );

    let ff = parse_frags_file(&out.join("q-r.frags"));
    assert_eq!(ff.x_total, (q0.len() + q1.len() + q2.len()) as u64);
    assert_eq!(ff.y_total, (r0.len() + r1.len() + r2.len()) as u64);
    assert_eq!(ff.fragments.len(), 2);

    // Check separator-space diagonal invariant: diag == (x_start + seq_x) - (y_start + seq_y)
    for f in &ff.fragments {
        let expected_diag = (f.x_start + f.seq_x) as i64 - (f.y_start + f.seq_y) as i64;
        assert_eq!(
            f.diag, expected_diag,
            "multi-FASTA separator-space diag mismatch"
        );
    }

    // Check JSONL local coordinates
    let jsonl = fs::read_to_string(out.join("q-r.jsonl")).unwrap();
    let recs: Vec<serde_json::Value> = jsonl
        .lines()
        .map(|l| serde_json::from_str(l).unwrap())
        .collect();
    assert_eq!(recs.len(), 2);
    let rec_01 = recs.iter().find(|v| v["query_id"] == "q_chr0").unwrap();
    assert_eq!(rec_01["target_id"], "r_chr1");
    assert_eq!(rec_01["seq_x"], 0);
    assert_eq!(rec_01["seq_y"], 1);
    assert_eq!(rec_01["query_start"], 80);
    assert_eq!(rec_01["query_end"], 80 + 140);
    assert_eq!(rec_01["target_start"], 50);
    assert_eq!(rec_01["target_end"], 50 + 140);

    let rec_20 = recs.iter().find(|v| v["query_id"] == "q_chr2").unwrap();
    assert_eq!(rec_20["target_id"], "r_chr0");
    assert_eq!(rec_20["seq_x"], 2);
    assert_eq!(rec_20["seq_y"], 0);
    assert_eq!(rec_20["query_start"], 60);
    assert_eq!(rec_20["query_end"], 60 + 150);
    assert_eq!(rec_20["target_start"], 70);
    assert_eq!(rec_20["target_end"], 70 + 150);
}

#[test]
fn tier1_multifasta_reverse_strand_seq_y_forward_index_mapping() {
    let tmp = tempfile::tempdir().unwrap();
    // 3 ref records: place inverted HSP on r_chr0 (seq_y=0) and another on r_chr2 (seq_y=2)
    let h0 = pseudo_dna(7101, 140);
    let h2 = pseudo_dna(7102, 160);

    let q_l = pseudo_dna(7110, 80);
    let q_m = pseudo_dna(7111, 80);
    let q_r = pseudo_dna(7112, 80);
    let q0 = format!("{q_l}{h0}{q_m}{h2}{q_r}");

    // Guard inverted flanks with all_mismatch(revcomp(...)) so exact local coords match
    let r0 = format!(
        "{}{}{}{}",
        all_mismatch(&revcomp(&q_m[..40])),
        revcomp(&h0),
        all_mismatch(&revcomp(&q_l[q_l.len() - 40..])),
        pseudo_dna(7121, 20)
    );
    let r1 = pseudo_dna(7122, 200);
    let r2 = format!(
        "{}{}{}{}{}",
        pseudo_dna(7123, 15),
        all_mismatch(&revcomp(&q_r[..40])),
        revcomp(&h2),
        all_mismatch(&revcomp(&q_m[q_m.len() - 40..])),
        pseudo_dna(7124, 5)
    );

    let q = tmp.path().join("q.fa");
    let r = tmp.path().join("r.fa");
    write_fasta(&q, &[("q0", &q0)]);
    write_fasta(&r, &[("r0", &r0), ("r1", &r1), ("r2", &r2)]);

    let out = tmp.path().join("out");
    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "12",
            "-l",
            "100",
            "-s",
            "80",
            "-f",
            "frags,paf,jsonl",
            "-o",
            out.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "multi-FASTA reverse seq_y mapping",
    );

    let ff = parse_frags_file(&out.join("q-r.frags"));
    assert_eq!(ff.fragments.len(), 2);
    assert!(ff.fragments.iter().all(|f| f.strand == 'r'));
    assert!(
        ff.fragments.iter().any(|f| f.seq_y == 0 && f.length == 140),
        "r0 inverted HSP must have original forward seq_y == 0"
    );
    assert!(
        ff.fragments.iter().any(|f| f.seq_y == 2 && f.length == 160),
        "r2 inverted HSP must have original forward seq_y == 2"
    );

    let jsonl = fs::read_to_string(out.join("q-r.jsonl")).unwrap();
    let recs: Vec<serde_json::Value> = jsonl
        .lines()
        .map(|l| serde_json::from_str(l).unwrap())
        .collect();
    let j0 = recs.iter().find(|v| v["target_id"] == "r0").unwrap();
    assert_eq!(j0["strand"], "-");
    assert_eq!(j0["seq_y"], 0);
    assert_eq!(j0["target_start"], 40);
    assert_eq!(j0["target_end"], 40 + 140);

    let j2 = recs.iter().find(|v| v["target_id"] == "r2").unwrap();
    assert_eq!(j2["strand"], "-");
    assert_eq!(j2["seq_y"], 2);
    assert_eq!(j2["target_start"], 55);
    assert_eq!(j2["target_end"], 55 + 160);
}

#[test]
fn tier1_n_bases_break_seeding_and_penalize_in_extension() {
    let tmp = tempfile::tempdir().unwrap();
    // Block 1: N every 10 bp -> no valid 12-mer seed -> 0 hits.
    // Block 2: 30 bp exact + 5 Ns (N vs N) + 30 bp exact -> seeds on 30-bp flanks,
    //          extends through all 5 Ns! Total len = 65, ident = 60 (N vs N NOT a match!),
    //          score = 60*4 - 5*1 = 235.
    let mut b1 = pseudo_dna(7201, 100).into_bytes();
    for i in (0..b1.len()).step_by(10) {
        b1[i] = b'N';
    }
    let b1_str = String::from_utf8(b1).unwrap();
    let b2 = format!("{}NNNNN{}", pseudo_dna(7202, 30), pseudo_dna(7203, 30));

    let f_l = pseudo_dna(7210, 100);
    let f_m = pseudo_dna(7211, 100);
    let f_r = pseudo_dna(7212, 100);
    let q = tmp.path().join("q.fa");
    let r = tmp.path().join("r.fa");
    write_fasta(&q, &[("q", &format!("{f_l}{b1_str}{f_m}{b2}{f_r}"))]);
    // Distinct diagonals + all_mismatch guards around b2 so length/ident/score are exact
    write_fasta(
        &r,
        &[(
            "r",
            &format!(
                "{}{}{}{}{}{}",
                pseudo_dna(7220, 100),
                b1_str,
                pseudo_dna(7221, 200),
                all_mismatch(&f_m),
                b2,
                all_mismatch(&f_r)
            ),
        )],
    );

    let out = tmp.path().join("out");
    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "12",
            "-l",
            "50",
            "-s",
            "60",
            "-f",
            "frags",
            "-o",
            out.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "N behavior",
    );

    let ff = parse_frags_file(&out.join("q-r.frags"));
    assert_eq!(
        ff.fragments.len(),
        1,
        "only anchored Block 2 should produce an HSP"
    );
    let f = &ff.fragments[0];
    assert_eq!(f.length, 65);
    assert_eq!(f.ident, 60, "N vs N must NOT count as identity match");
    assert_eq!(f.score, 60 * 4 - 5, "default n_penalty is -1 per N");
}

#[test]
fn tier1_iupac_breaks_seeding_and_scores_in_extension() {
    let tmp = tempfile::tempdir().unwrap();
    // Block 1: 'R' every 10 bp -> no valid 12-mer seed.
    // Block 2 (R vs R): 30 bp exact + 4 'R's + 30 bp exact -> len 64, ident 64, score 256.
    // Block 3 (R vs Y): 30 bp exact + 4 'R' vs 4 'Y' + 30 bp exact -> len 64, ident 60, score 60*4 - 4*4 = 224.
    let mut b1 = pseudo_dna(7301, 100).into_bytes();
    for i in (0..b1.len()).step_by(10) {
        b1[i] = b'R';
    }
    let b1_str = String::from_utf8(b1).unwrap();
    let left2 = pseudo_dna(7302, 30);
    let right2 = pseudo_dna(7303, 30);
    let left3 = pseudo_dna(7304, 30);
    let right3 = pseudo_dna(7305, 30);

    let f0 = pseudo_dna(7310, 100);
    let f1 = pseudo_dna(7311, 100);
    let f2 = pseudo_dna(7312, 100);
    let f3 = pseudo_dna(7313, 100);
    let q_seq = format!("{f0}{b1_str}{f1}{left2}RRRR{right2}{f2}{left3}RRRR{right3}{f3}");
    // Place each block on a distinct diagonal with all_mismatch boundary guards
    let r_seq = format!(
        "{}{b1_str}{}{}{left2}RRRR{right2}{}{}{left3}YYYY{right3}{}",
        pseudo_dna(7320, 100),
        pseudo_dna(7321, 180),
        all_mismatch(&f1),
        all_mismatch(&f2[..50]),
        pseudo_dna(7322, 180),
        all_mismatch(&f3),
    );

    let q = tmp.path().join("q.fa");
    let r = tmp.path().join("r.fa");
    write_fasta(&q, &[("q", &q_seq)]);
    write_fasta(&r, &[("r", &r_seq)]);

    let out = tmp.path().join("out");
    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "12",
            "-l",
            "50",
            "-s",
            "60",
            "-f",
            "frags",
            "-o",
            out.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "IUPAC behavior",
    );

    let ff = parse_frags_file(&out.join("q-r.frags"));
    assert_eq!(ff.fragments.len(), 2);
    let f_match = ff
        .fragments
        .iter()
        .find(|f| f.ident == 64)
        .expect("R vs R match");
    assert_eq!(f_match.length, 64);
    assert_eq!(f_match.score, 64 * 4);
    let f_mm = ff
        .fragments
        .iter()
        .find(|f| f.ident == 60)
        .expect("R vs Y mismatch");
    assert_eq!(f_mm.length, 64);
    assert_eq!(f_mm.score, 60 * 4 - 4 * 4);
}

#[test]
fn tier1_rna_u_treated_identically_to_t() {
    let tmp = tempfile::tempdir().unwrap();
    let hsp_f = pseudo_dna(7401, 150);
    let hsp_r = pseudo_dna(7402, 140);
    let q_dna = format!(
        "{}{}{}{}{}",
        pseudo_dna(7410, 100),
        hsp_f,
        pseudo_dna(7411, 100).to_ascii_lowercase(),
        hsp_r,
        pseudo_dna(7412, 100)
    );
    let r_dna = format!(
        "{}{}{}{}{}",
        pseudo_dna(7420, 100),
        hsp_f,
        pseudo_dna(7421, 100).to_ascii_lowercase(),
        revcomp(&hsp_r),
        pseudo_dna(7422, 100)
    );
    let q_rna = q_dna.replace('T', "U").replace('t', "u");
    let r_rna = r_dna.replace('T', "U").replace('t', "u");

    let qd = tmp.path().join("qd.fa");
    let rd = tmp.path().join("rd.fa");
    let qr = tmp.path().join("qr.fa");
    let rr = tmp.path().join("rr.fa");
    write_fasta(&qd, &[("q", &q_dna)]);
    write_fasta(&rd, &[("r", &r_dna)]);
    write_fasta(&qr, &[("q", &q_rna)]);
    write_fasta(&rr, &[("r", &r_rna)]);

    let out_dna = tmp.path().join("dna");
    let out_rna = tmp.path().join("rna");
    assert_success(
        &run_gecko(&[
            "compare",
            qd.to_str().unwrap(),
            rd.to_str().unwrap(),
            "-k",
            "12",
            "-l",
            "80",
            "-s",
            "80",
            "--prefix",
            "cmp",
            "-f",
            "frags",
            "-o",
            out_dna.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "DNA compare",
    );
    assert_success(
        &run_gecko(&[
            "compare",
            qr.to_str().unwrap(),
            rr.to_str().unwrap(),
            "-k",
            "12",
            "-l",
            "80",
            "-s",
            "80",
            "--prefix",
            "cmp",
            "-f",
            "frags",
            "-o",
            out_rna.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "RNA compare",
    );

    assert_frags_bytes_identical(
        &out_dna.join("cmp.frags"),
        &out_rna.join("cmp.frags"),
        "RNA U/u vs DNA T/t",
    );
    assert_eq!(
        parse_frags_file(&out_rna.join("cmp.frags")).fragments.len(),
        2
    );
}

// --- 1F. Output formats & `all-vs-all` (5 tests) ---

#[test]
fn tier1_output_frags_wire_layout_and_sorting() {
    let tmp = tempfile::tempdir().unwrap();
    let (q, r) = build_multi_hsp_pair(tmp.path());
    let out = tmp.path().join("out");

    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "12",
            "-l",
            "60",
            "-s",
            "60",
            "-f",
            "frags",
            "-o",
            out.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "frags wire layout",
    );

    let frags_path = out.join("q-r.frags");
    let raw = fs::read(&frags_path).unwrap();
    let ff = parse_frags_file(&frags_path);
    assert_eq!(raw.len(), 16 + FRAG_WIRE_SIZE * ff.fragments.len());
    assert!(ff.fragments.len() >= 5);

    // Check exact dedup_fragments sort order:
    // (strand, diag, seq_x, seq_y, x_start, Reverse(x_end), Reverse(score))
    for w in ff.fragments.windows(2) {
        let k0 = (
            w[0].strand,
            w[0].diag,
            w[0].seq_x,
            w[0].seq_y,
            w[0].x_start,
            Reverse(w[0].x_end),
            Reverse(w[0].score),
        );
        let k1 = (
            w[1].strand,
            w[1].diag,
            w[1].seq_x,
            w[1].seq_y,
            w[1].x_start,
            Reverse(w[1].x_end),
            Reverse(w[1].score),
        );
        assert!(k0 <= k1, "frags not sorted: {:?} > {:?}", w[0], w[1]);
    }
    for f in &ff.fragments {
        assert_eq!(f.block, 0);
        assert_eq!(f.x_end - f.x_start + 1, f.length);
        assert_eq!(f.y_end - f.y_start + 1, f.length);
        assert!(f.ident <= f.length);
    }
}

#[test]
fn tier1_output_csv_and_inf_exact_schema() {
    let tmp = tempfile::tempdir().unwrap();
    let (q, r) = build_multi_hsp_pair(tmp.path());
    let out = tmp.path().join("out");

    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "12",
            "-l",
            "100",
            "-s",
            "70",
            "-f",
            "frags,csv,inf",
            "-o",
            out.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "csv + inf schema",
    );

    let ff = parse_frags_file(&out.join("q-r.frags"));
    assert!(
        !ff.fragments.is_empty(),
        "csv + inf schema fixture must emit at least one fragment"
    );
    let inf_text = fs::read_to_string(out.join("q-r.frags.INF")).unwrap();
    let inf_lines: Vec<&str> = inf_text.lines().collect();
    assert_eq!(inf_lines.len(), 14, ".frags.INF must have exactly 14 lines");
    assert_eq!(inf_lines[13], "=".repeat(56));
    assert!(inf_lines[2].ends_with("q.fa"));
    assert!(inf_lines[3].ends_with("r.fa"));
    assert!(inf_lines[6].ends_with(&ff.x_total.to_string()));
    assert!(inf_lines[7].ends_with(&ff.y_total.to_string()));
    assert!(inf_lines[12].ends_with(&ff.fragments.len().to_string()));

    let csv_text = fs::read_to_string(out.join("q-r.csv")).unwrap();
    let csv_lines: Vec<&str> = csv_text.lines().collect();
    assert_eq!(csv_lines.len(), 17 + ff.fragments.len());
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
        assert_eq!(cols[1].parse::<u64>().unwrap(), f.x_start);
        assert_eq!(cols[3].parse::<u64>().unwrap(), f.x_end);
        let ys: u64 = cols[2].parse().unwrap();
        let ye: u64 = cols[4].parse().unwrap();
        if f.strand == 'f' {
            assert_eq!(ys, f.y_start);
            assert_eq!(ye, f.y_end);
            assert!(ys <= ye);
        } else {
            assert_eq!(ys, ff.y_total - f.y_start - 1);
            assert_eq!(ye, ff.y_total - f.y_end - 1);
            assert!(ys >= ye, "reverse CSV must have yStart >= yEnd");
        }
        assert_eq!(cols[5], f.strand.to_string());
        assert_eq!(cols[7].parse::<u64>().unwrap(), f.length);
        assert_eq!(cols[8].parse::<u64>().unwrap(), f.score);
        assert_eq!(cols[9].parse::<u64>().unwrap(), f.ident);
        assert!(cols[10].contains('.') && cols[10].split('.').nth(1).unwrap().len() == 2);
        assert!(cols[11].contains('.') && cols[11].split('.').nth(1).unwrap().len() == 2);
    }
}

#[test]
fn tier1_output_mat_exact_dimensions_and_lazy_omission_parity() {
    let tmp = tempfile::tempdir().unwrap();
    let (q, r) = build_multi_hsp_pair(tmp.path());
    let out_leg = tmp.path().join("legacy");
    let out_frags = tmp.path().join("frags_only");

    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "12",
            "-l",
            "100",
            "-s",
            "70",
            "-f",
            "legacy",
            "-o",
            out_leg.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "--format legacy",
    );
    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "12",
            "-l",
            "100",
            "-s",
            "70",
            "-f",
            "frags",
            "-o",
            out_frags.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "--format frags",
    );

    let ff = parse_frags_file(&out_leg.join("q-r.frags"));
    let mat = parse_mat_file(&out_leg.join("q-r.frags.MAT"));
    assert_mat_covers_frags(&mat, &ff.fragments);

    assert!(
        !out_frags.join("q-r.frags.MAT").exists(),
        "--format frags must NOT emit .frags.MAT"
    );
    assert_frags_bytes_identical(
        &out_leg.join("q-r.frags"),
        &out_frags.join("q-r.frags"),
        "--format legacy vs --format frags",
    );
}

#[test]
fn tier1_output_modern_formats_paf_blast_bedpe_jsonl_png() {
    let tmp = tempfile::tempdir().unwrap();
    let (q, r) = build_multi_hsp_pair(tmp.path());
    let out = tmp.path().join("all");

    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "12",
            "-l",
            "100",
            "-s",
            "70",
            "-f",
            "all",
            "-o",
            out.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "--format all",
    );

    for ext in [
        "q-r.frags",
        "q-r.csv",
        "q-r.frags.INF",
        "q-r.frags.MAT",
        "q-r.paf",
        "q-r.blast.tsv",
        "q-r.bedpe",
        "q-r.jsonl",
        "q-r.png",
    ] {
        assert!(
            out.join(ext).exists(),
            "missing expected --format all file {ext}"
        );
    }

    let ff = parse_frags_file(&out.join("q-r.frags"));
    let n = ff.fragments.len();
    assert_eq!(n, 3);

    // PAF: 15 columns
    let paf_lines: Vec<String> = fs::read_to_string(out.join("q-r.paf"))
        .unwrap()
        .lines()
        .map(str::to_string)
        .collect();
    assert_eq!(paf_lines.len(), n);
    for (f, line) in ff.fragments.iter().zip(&paf_lines) {
        let cols: Vec<&str> = line.split('\t').collect();
        assert_eq!(cols.len(), 15, "PAF must have 15 cols: {line}");
        assert_eq!(cols[4], if f.strand == 'f' { "+" } else { "-" });
        assert_eq!(cols[9].parse::<u64>().unwrap(), f.ident);
        assert_eq!(cols[10].parse::<u64>().unwrap(), f.length);
        assert_eq!(cols[11], "255");
        assert_eq!(cols[12], format!("AS:i:{}", f.score));
        assert!(cols[13].starts_with("si:f:"));
        assert_eq!(cols[14], format!("cg:Z:{}M", f.length));
    }

    // BLAST tabular: 12 columns + reverse sstart >= send + Karlin-Altschul
    let blast_lines: Vec<String> = fs::read_to_string(out.join("q-r.blast.tsv"))
        .unwrap()
        .lines()
        .map(str::to_string)
        .collect();
    assert_eq!(blast_lines.len(), n);
    let lambda = 3.0_f64.ln() / 4.0;
    let k_ka = 0.1_f64;
    let mn = (ff.x_total as f64) * (ff.y_total as f64);
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
            assert!(sstart >= send, "reverse BLAST must have sstart >= send");
        }
        let expected_e = (k_ka * mn * (-lambda * (f.score as f64)).exp()).max(0.0);
        let expected_bits = ((lambda * (f.score as f64) - k_ka.ln()) / 2.0_f64.ln()).max(0.0);
        let got_e: f64 = cols[10].parse().unwrap();
        let got_bits: f64 = cols[11].parse().unwrap();
        assert!((got_bits - expected_bits).abs() <= 0.15);
        assert!(
            (got_e - expected_e).abs() <= expected_e.abs() * 0.05 + 1e-12,
            "evalue mismatch: got {got_e}, expected {expected_e}"
        );
    }

    // BEDPE: 10 columns + frag_<i>
    let bedpe_lines: Vec<String> = fs::read_to_string(out.join("q-r.bedpe"))
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

    // JSONL: exact 16 keys per object
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
    let jsonl_lines: Vec<String> = fs::read_to_string(out.join("q-r.jsonl"))
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

    // PNG: magic + IHDR dimensions >= 200x200
    let png = fs::read(out.join("q-r.png")).unwrap();
    assert!(png.len() > 32);
    assert_eq!(&png[0..8], &[137, 80, 78, 71, 13, 10, 26, 10]);
    assert_eq!(&png[12..16], b"IHDR");
    let w = u32::from_be_bytes(png[16..20].try_into().unwrap());
    let h = u32::from_be_bytes(png[20..24].try_into().unwrap());
    assert!(w >= 200 && h >= 200);
}

#[test]
fn tier1_all_vs_all_unequal_lengths_pairwise_parity_and_skip_existing() {
    let tmp = tempfile::tempdir().unwrap();
    let gdir = tmp.path().join("genomes");
    fs::create_dir_all(&gdir).unwrap();

    // 3 unequal-length genomes: g1 (short ~900 bp), g2 (medium ~1900 bp), g3 (long ~3200 bp)
    let h12 = pseudo_dna(7501, 200);
    let h13 = pseudo_dna(7502, 180);
    let h23 = pseudo_dna(7503, 220);

    let s1 = format!(
        "{}{}{}{}{}",
        pseudo_dna(7510, 200),
        h12,
        pseudo_dna(7511, 150),
        h13,
        pseudo_dna(7512, 170)
    );
    let s2 = format!(
        "{}{}{}{}{}",
        pseudo_dna(7520, 600),
        revcomp(&h12),
        pseudo_dna(7521, 500),
        h23,
        pseudo_dna(7522, 400)
    );
    let s3 = format!(
        "{}{}{}{}{}",
        pseudo_dna(7530, 1000),
        h13,
        pseudo_dna(7531, 900),
        revcomp(&h23),
        pseudo_dna(7532, 900)
    );

    let p1 = gdir.join("g1.fa");
    let p2 = gdir.join("g2.fa");
    let p3 = gdir.join("g3.fa");
    write_fasta(&p1, &[("g1", &s1)]);
    write_fasta(&p2, &[("g2", &s2)]);
    write_fasta(&p3, &[("g3", &s3)]);

    let ava_out = tmp.path().join("ava");
    assert_success(
        &run_gecko(&[
            "all-vs-all",
            gdir.to_str().unwrap(),
            "-k",
            "12",
            "-l",
            "100",
            "-s",
            "70",
            "-f",
            "frags",
            "-o",
            ava_out.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "all-vs-all unequal lengths",
    );

    // Every upper-triangle pair must be 100% bit-identical to standalone `gecko compare`
    let pairs = [
        (&p1, &p2, "g1-g2"),
        (&p1, &p3, "g1-g3"),
        (&p2, &p3, "g2-g3"),
    ];
    for (qa, rb, name) in pairs {
        let single_out = tmp.path().join(format!("single_{name}"));
        assert_success(
            &run_gecko(&[
                "compare",
                qa.to_str().unwrap(),
                rb.to_str().unwrap(),
                "-k",
                "12",
                "-l",
                "100",
                "-s",
                "70",
                "-f",
                "frags",
                "-o",
                single_out.to_str().unwrap(),
                "--backend",
                "cpu",
            ]),
            &format!("standalone compare {name}"),
        );
        assert_frags_bytes_identical(
            &ava_out.join(format!("{name}.frags")),
            &single_out.join(format!("{name}.frags")),
            &format!("all-vs-all vs compare for {name}"),
        );
        assert_eq!(
            parse_frags_file(&ava_out.join(format!("{name}.frags")))
                .fragments
                .len(),
            1
        );
    }

    // Verify --skip-existing preserves sentinel content in an existing .frags file
    let sentinel_path = ava_out.join("g1-g2.frags");
    let sentinel_bytes = b"SENTINEL_EXISTING_FRAGS";
    fs::write(&sentinel_path, sentinel_bytes).unwrap();
    fs::remove_file(ava_out.join("g2-g3.frags")).unwrap();
    assert_success(
        &run_gecko(&[
            "all-vs-all",
            gdir.to_str().unwrap(),
            "-k",
            "12",
            "-l",
            "100",
            "-s",
            "70",
            "--skip-existing",
            "-f",
            "frags",
            "-o",
            ava_out.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "all-vs-all --skip-existing",
    );
    assert_eq!(
        fs::read(&sentinel_path).unwrap(),
        sentinel_bytes,
        "--skip-existing must not overwrite existing .frags"
    );
    assert!(
        ava_out.join("g2-g3.frags").exists(),
        "--skip-existing must still generate missing pairs"
    );
}

// --- 1G. New M3 `--seeding` (5 tests; expected FAIL on HEAD 5237ca9, PASS after M3) ---

fn assert_minimizer_vs_exact_parity(k: usize) {
    let tmp = tempfile::tempdir().unwrap();
    let (q, r) = build_multi_hsp_pair(tmp.path());
    let out_exact = tmp.path().join("exact");
    let out_min = tmp.path().join("min");
    let k_str = k.to_string();

    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            &k_str,
            "-l",
            "60",
            "-s",
            "60",
            "--freq-limit",
            "2",
            "--seeding",
            "exact",
            "-f",
            "frags",
            "-o",
            out_exact.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        &format!("--seeding exact k={k}"),
    );
    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            &k_str,
            "-l",
            "60",
            "-s",
            "60",
            "--freq-limit",
            "2",
            "--seeding",
            "minimizer",
            "-f",
            "frags",
            "-o",
            out_min.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        &format!("--seeding minimizer k={k}"),
    );

    assert_frags_bytes_identical(
        &out_exact.join("q-r.frags"),
        &out_min.join("q-r.frags"),
        &format!("--seeding minimizer vs exact at k={k}"),
    );
    let ff = parse_frags_file(&out_min.join("q-r.frags"));
    assert!(
        ff.fragments.len() >= 3,
        "expected >= 3 HSPs at k={k}, got {}",
        ff.fragments.len()
    );
}

#[test]
fn tier1_seeding_minimizer_k18_exact_parity() {
    assert_minimizer_vs_exact_parity(18);
}

#[test]
fn tier1_seeding_minimizer_k20_exact_parity() {
    assert_minimizer_vs_exact_parity(20);
}

#[test]
fn tier1_seeding_minimizer_k24_exact_parity() {
    assert_minimizer_vs_exact_parity(24);
}

#[test]
fn tier1_seeding_minimizer_k32_exact_parity() {
    assert_minimizer_vs_exact_parity(32);
}

#[test]
fn tier1_seeding_auto_parity_and_short_k_fallback() {
    let tmp = tempfile::tempdir().unwrap();
    let (q, r) = build_multi_hsp_pair(tmp.path());

    for k_str in ["12", "24"] {
        let out_exact = tmp.path().join(format!("exact_{k_str}"));
        let out_auto = tmp.path().join(format!("auto_{k_str}"));
        assert_success(
            &run_gecko(&[
                "compare",
                q.to_str().unwrap(),
                r.to_str().unwrap(),
                "-k",
                k_str,
                "-l",
                "60",
                "-s",
                "60",
                "--seeding",
                "exact",
                "-f",
                "frags",
                "-o",
                out_exact.to_str().unwrap(),
                "--backend",
                "cpu",
            ]),
            &format!("--seeding exact k={k_str}"),
        );
        assert_success(
            &run_gecko(&[
                "compare",
                q.to_str().unwrap(),
                r.to_str().unwrap(),
                "-k",
                k_str,
                "-l",
                "60",
                "-s",
                "60",
                "--seeding",
                "auto",
                "-f",
                "frags",
                "-o",
                out_auto.to_str().unwrap(),
                "--backend",
                "cpu",
            ]),
            &format!("--seeding auto k={k_str}"),
        );
        assert_frags_bytes_identical(
            &out_exact.join("q-r.frags"),
            &out_auto.join("q-r.frags"),
            &format!("--seeding auto vs exact at k={k_str}"),
        );
    }

    // Also confirm --seeding minimizer gracefully falls back on short k=12
    let out_min12 = tmp.path().join("min_12");
    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "12",
            "-l",
            "60",
            "-s",
            "60",
            "--seeding",
            "minimizer",
            "-f",
            "frags",
            "-o",
            out_min12.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "--seeding minimizer k=12 fallback",
    );
    assert_frags_bytes_identical(
        &tmp.path().join("exact_12/q-r.frags"),
        &out_min12.join("q-r.frags"),
        "--seeding minimizer k=12 fallback vs exact",
    );
}

// --- 1H. New M3 `--two-hit` (5 tests; expected FAIL on HEAD 5237ca9, PASS after M3) ---

/// Construct an HSP of ~56 bp that has EXACTLY ONE 12-mer match (center 12 bp),
/// flanked by 22 bp on each side with a mismatch every 11 bp (no other 12-mer seeds,
/// yet 92.8% identity / 85.7% similarity so single-hit strongly accepts it).
fn make_single_seed_k12_hsp(seed: u64) -> (String, String) {
    let left = pseudo_dna(seed, 22);
    let mid = pseudo_dna(seed + 1, 12);
    let right = pseudo_dna(seed + 2, 22);
    let q = format!("{left}{mid}{right}");
    // Mismatch at positions 0 and 11 in left (so runs are 10 bp and 10 bp), and at 0 and 11 in right
    let left_m = mutate_positions(&left, &[10, 21]);
    let right_m = mutate_positions(&right, &[0, 11]);
    let r = format!("{left_m}{mid}{right_m}");
    (q, r)
}

/// Construct a dual-seed k=12 HSP with two 12-mer seeds separated by `gap_bp` bases
/// (where every `10` bp in the gap has a mismatch so NO intermediate 12-mer exists).
fn make_dual_seed_k12_hsp(seed: u64, gap_bp: usize) -> (String, String) {
    let flank_l = pseudo_dna(seed, 16);
    let s1 = pseudo_dna(seed + 1, 12);
    let gap = pseudo_dna(seed + 2, gap_bp);
    let s2 = pseudo_dna(seed + 3, 12);
    let flank_r = pseudo_dna(seed + 4, 16);

    let q = format!("{flank_l}{s1}{gap}{s2}{flank_r}");
    let flank_l_m = mutate_positions(&flank_l, &[5, 15]);
    let mut gap_mut_pos = vec![0];
    let mut p = 10;
    while p < gap_bp {
        gap_mut_pos.push(p);
        p += 10;
    }
    if gap_bp > 1 {
        gap_mut_pos.push(gap_bp - 1);
    }
    let gap_m = mutate_positions(&gap, &gap_mut_pos);
    let flank_r_m = mutate_positions(&flank_r, &[0, 10]);
    let r = format!("{flank_l_m}{s1}{gap_m}{s2}{flank_r_m}");
    (q, r)
}

#[test]
fn tier1_two_hit_default_window_preserves_multi_seed_hsp() {
    let tmp = tempfile::tempdir().unwrap();
    let (q, r) = build_multi_hsp_pair(tmp.path());
    let out = tmp.path().join("twohit");

    // Portable Clap rule: place bare `--two-hit` at the END after positional args
    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "12",
            "-l",
            "100",
            "-s",
            "70",
            "-f",
            "frags",
            "-o",
            out.to_str().unwrap(),
            "--backend",
            "cpu",
            "--two-hit",
        ]),
        "bare --two-hit (default window 40)",
    );

    let ff = parse_frags_file(&out.join("q-r.frags"));
    assert_eq!(
        ff.fragments.len(),
        3,
        "all 3 multi-seed HSPs must be preserved"
    );
    assert_eq!(ff.fragments.iter().filter(|f| f.strand == 'f').count(), 2);
    assert_eq!(ff.fragments.iter().filter(|f| f.strand == 'r').count(), 1);
}

#[test]
fn tier1_two_hit_suppresses_single_seed_keeps_dual_seed() {
    let tmp = tempfile::tempdir().unwrap();
    let (single_q, single_r) = make_single_seed_k12_hsp(8001);
    let (dual_q, dual_r) = make_dual_seed_k12_hsp(8010, 5); // gap=5, start-to-start=17 <= 40

    let f0 = pseudo_dna(8020, 150);
    let f1 = pseudo_dna(8021, 150);
    let f2 = pseudo_dna(8022, 150);
    let q_seq = format!("{f0}{single_q}{f1}{dual_q}{f2}");
    // Put single and dual on distinct diagonals with all_mismatch guards
    let r_seq = format!(
        "{}{single_r}{}{}{}{dual_r}{}",
        all_mismatch(&f0),
        all_mismatch(&f1[..50]),
        pseudo_dna(8031, 200),
        all_mismatch(&f1[f1.len() - 50..]),
        all_mismatch(&f2)
    );

    let q = tmp.path().join("q.fa");
    let r = tmp.path().join("r.fa");
    write_fasta(&q, &[("q", &q_seq)]);
    write_fasta(&r, &[("r", &r_seq)]);

    let out_single = tmp.path().join("single_hit");
    let out_two = tmp.path().join("two_hit");

    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "12",
            "-l",
            "40",
            "-s",
            "60",
            "-f",
            "frags",
            "-o",
            out_single.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "single-hit baseline",
    );
    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "12",
            "-l",
            "40",
            "-s",
            "60",
            "--two-hit=40",
            "-f",
            "frags",
            "-o",
            out_two.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "--two-hit=40",
    );

    let f1 = parse_frags_file(&out_single.join("q-r.frags")).fragments;
    let f2 = parse_frags_file(&out_two.join("q-r.frags")).fragments;
    assert_eq!(
        f1.len(),
        2,
        "single-hit must recover both single-seed and dual-seed HSPs"
    );
    assert_eq!(
        f2.len(),
        1,
        "--two-hit=40 must suppress single-seed HSP and keep dual-seed HSP"
    );
    // Surviving fragment must be the dual-seed HSP (at x_start >= 150 + single_q.len() + 150)
    assert!(
        f2[0].x_start >= (150 + single_q.len() + 150) as u64,
        "surviving fragment must be the dual-seed HSP, got x_start={}",
        f2[0].x_start
    );
}

#[test]
fn tier1_two_hit_explicit_window_boundary() {
    let tmp = tempfile::tempdir().unwrap();
    // Dual-seed HSP with gap=43 bp (start-to-start=55 bp):
    //   --two-hit=30 -> BOTH gap(43) > 30 and start-to-start(55) > 30 -> REJECTED!
    //   --two-hit=80 -> BOTH gap(43) <= 80 and start-to-start(55) <= 80 -> ACCEPTED!
    // Plus a close dual-seed HSP (gap=5 bp, start-to-start=17 bp) accepted by both.
    let (wide_q, wide_r) = make_dual_seed_k12_hsp(8101, 43);
    let (close_q, close_r) = make_dual_seed_k12_hsp(8110, 5);

    let f0 = pseudo_dna(8120, 150);
    let f1 = pseudo_dna(8121, 150);
    let f2 = pseudo_dna(8122, 150);
    let q_seq = format!("{f0}{wide_q}{f1}{close_q}{f2}");
    // Put wide and close on distinct diagonals with all_mismatch guards
    let r_seq = format!(
        "{}{wide_r}{}{}{}{close_r}{}",
        all_mismatch(&f0),
        all_mismatch(&f1[..50]),
        pseudo_dna(8131, 200),
        all_mismatch(&f1[f1.len() - 50..]),
        all_mismatch(&f2)
    );

    let q = tmp.path().join("q.fa");
    let r = tmp.path().join("r.fa");
    write_fasta(&q, &[("q", &q_seq)]);
    write_fasta(&r, &[("r", &r_seq)]);

    let out_w30 = tmp.path().join("w30");
    let out_w80 = tmp.path().join("w80");
    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "12",
            "-l",
            "40",
            "-s",
            "60",
            "--two-hit=30",
            "-f",
            "frags",
            "-o",
            out_w30.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "--two-hit=30",
    );
    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "12",
            "-l",
            "40",
            "-s",
            "60",
            "--two-hit=80",
            "-f",
            "frags",
            "-o",
            out_w80.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "--two-hit=80",
    );

    let f30 = parse_frags_file(&out_w30.join("q-r.frags")).fragments;
    let f80 = parse_frags_file(&out_w80.join("q-r.frags")).fragments;
    assert_eq!(
        f30.len(),
        1,
        "--two-hit=30 must keep only close dual-seed HSP"
    );
    assert_eq!(
        f80.len(),
        2,
        "--two-hit=80 must keep both close and wide dual-seed HSPs"
    );
}

#[test]
fn tier1_two_hit_reduces_noise_hits_on_divergent_background() {
    let tmp = tempfile::tempdir().unwrap();
    // 25 Kb x 25 Kb random background at k=8 produces ~9,500 random single-hit extensions
    // (25000^2 / 65536 * 2 strands ~ 19,000), plus 1 embedded 200-bp true HSP.
    let hsp = pseudo_dna(8201, 200);
    let q_seq = format!(
        "{}{}{}",
        pseudo_dna(8210, 12_500),
        hsp,
        pseudo_dna(8211, 12_500)
    );
    let r_seq = format!(
        "{}{}{}",
        pseudo_dna(8220, 12_500),
        hsp,
        pseudo_dna(8221, 12_500)
    );

    let q = tmp.path().join("q.fa");
    let r = tmp.path().join("r.fa");
    write_fasta(&q, &[("q", &q_seq)]);
    write_fasta(&r, &[("r", &r_seq)]);

    let out_base = tmp.path().join("base");
    let out_two = tmp.path().join("two");

    let res_base = run_gecko(&[
        "compare",
        q.to_str().unwrap(),
        r.to_str().unwrap(),
        "-k",
        "8",
        "-l",
        "100",
        "-s",
        "80",
        "--json-stats",
        "-f",
        "frags",
        "-o",
        out_base.to_str().unwrap(),
        "--backend",
        "cpu",
    ]);
    assert_success(&res_base, "noisy background single-hit");

    let res_two = run_gecko(&[
        "compare",
        q.to_str().unwrap(),
        r.to_str().unwrap(),
        "-k",
        "8",
        "-l",
        "100",
        "-s",
        "80",
        "--two-hit=40",
        "--json-stats",
        "-f",
        "frags",
        "-o",
        out_two.to_str().unwrap(),
        "--backend",
        "cpu",
    ]);
    assert_success(&res_two, "noisy background --two-hit=40");

    let s_base: serde_json::Value = serde_json::from_slice(&res_base.stdout).unwrap();
    let s_two: serde_json::Value = serde_json::from_slice(&res_two.stdout).unwrap();

    let hits_base = s_base["hits_used"].as_u64().unwrap();
    let hits_two = s_two["hits_used"].as_u64().unwrap();
    assert!(
        hits_base > 1_000,
        "expected substantial random noise hits, got {hits_base}"
    );
    assert!(
        hits_two * 10 < hits_base,
        "expected --two-hit=40 to cut noise hits by >10x (base={hits_base}, two_hit={hits_two})"
    );
    assert_eq!(
        parse_frags_file(&out_two.join("q-r.frags")).fragments.len(),
        1
    );
}

#[test]
fn tier1_two_hit_with_spaced_seed_w12s18() {
    let tmp = tempfile::tempdir().unwrap();
    // Divergent multi-seed HSP matching w12s18 across 180 bp
    let hsp = pseudo_dna(8301, 180);
    let mut mut_pos = Vec::new();
    for blk in 0..10 {
        for &dc in &[3usize, 9, 14] {
            mut_pos.push(blk * 18 + dc);
        }
    }
    let hsp_div = mutate_positions(&hsp, &mut_pos);

    let q = tmp.path().join("q.fa");
    let r = tmp.path().join("r.fa");
    write_fasta(
        &q,
        &[(
            "q",
            &format!("{}{}{}", pseudo_dna(8310, 200), hsp, pseudo_dna(8311, 200)),
        )],
    );
    write_fasta(
        &r,
        &[(
            "r",
            &format!(
                "{}{}{}",
                pseudo_dna(8320, 200),
                hsp_div,
                pseudo_dna(8321, 200)
            ),
        )],
    );

    let out = tmp.path().join("out");
    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "--spaced-seed",
            "w12s18",
            "--two-hit=40",
            "-l",
            "80",
            "-s",
            "40",
            "-f",
            "frags",
            "-o",
            out.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "--spaced-seed w12s18 --two-hit=40",
    );
    let ff = parse_frags_file(&out.join("q-r.frags"));
    assert_eq!(ff.fragments.len(), 1);
    assert!(ff.fragments[0].length >= 160);
}

// --- 1I. New M3 `GECKO_PROFILE=1` (2 tests; expected FAIL on HEAD 5237ca9, PASS after M3) ---

#[test]
fn tier1_gecko_profile_emits_stderr_telemetry_and_exact_frags() {
    let tmp = tempfile::tempdir().unwrap();
    let (q, r) = build_multi_hsp_pair(tmp.path());
    let out_def = tmp.path().join("def");
    let out_prof = tmp.path().join("prof");

    let res_def = run_gecko_with_env(
        &[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "16",
            "-l",
            "80",
            "-s",
            "70",
            "-f",
            "frags",
            "-o",
            out_def.to_str().unwrap(),
            "--backend",
            "cpu",
        ],
        &[("GECKO_PROFILE", None)],
    );
    assert_success(&res_def, "default without GECKO_PROFILE");

    let res_prof = run_gecko_with_env(
        &[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "16",
            "-l",
            "80",
            "-s",
            "70",
            "-f",
            "frags",
            "-o",
            out_prof.to_str().unwrap(),
            "--backend",
            "cpu",
        ],
        &[("GECKO_PROFILE", Some("1"))],
    );
    assert_success(&res_prof, "with GECKO_PROFILE=1");

    assert_frags_bytes_identical(
        &out_def.join("q-r.frags"),
        &out_prof.join("q-r.frags"),
        "GECKO_PROFILE=1 must not change .frags output",
    );

    let stderr_def = String::from_utf8_lossy(&res_def.stderr);
    let stderr_prof = String::from_utf8_lossy(&res_prof.stderr);
    assert!(
        stderr_prof.lines().count() > stderr_def.lines().count(),
        "GECKO_PROFILE=1 must emit additional profiling lines to stderr!\ndefault ({} lines):\n{stderr_def}\nprofile ({} lines):\n{stderr_prof}",
        stderr_def.lines().count(),
        stderr_prof.lines().count()
    );
    let lower = stderr_prof.to_ascii_lowercase();
    assert!(
        lower.contains("profile")
            || lower.contains("index")
            || lower.contains("bloom")
            || lower.contains("extend")
            || lower.contains("query")
            || lower.contains("bypass"),
        "GECKO_PROFILE=1 stderr should contain sub-phase/counter telemetry:\n{stderr_prof}"
    );
}

#[test]
fn tier1_gecko_profile_preserves_clean_json_stats_stdout() {
    let tmp = tempfile::tempdir().unwrap();
    let (q, r) = build_multi_hsp_pair(tmp.path());
    let out = tmp.path().join("prof_json");

    let res = run_gecko_with_env(
        &[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "16",
            "-l",
            "80",
            "-s",
            "70",
            "--json-stats",
            "-f",
            "frags",
            "-o",
            out.to_str().unwrap(),
            "--backend",
            "cpu",
        ],
        &[("GECKO_PROFILE", Some("1"))],
    );
    assert_success(&res, "GECKO_PROFILE=1 --json-stats");

    let json: serde_json::Value = serde_json::from_slice(&res.stdout).unwrap_or_else(|e| {
        panic!(
            "stdout must remain valid JSON under GECKO_PROFILE=1: {e}\nstdout:\n{}",
            String::from_utf8_lossy(&res.stdout)
        )
    });
    assert!(
        json["backend"].as_str().unwrap().starts_with("cpu"),
        "expected cpu backend, got {}",
        json["backend"]
    );
    assert!(json["fragments"].as_u64().unwrap() >= 3);

    let stderr = String::from_utf8_lossy(&res.stderr);
    assert!(
        !stderr.trim().is_empty(),
        "GECKO_PROFILE=1 with --json-stats must emit profiling telemetry to stderr (not stdout)"
    );
}

// --- 1J. New M4 CUDA Static Verification (5 tests: 1 baseline PASS + 4 M4 contract) ---

#[test]
fn tier1_cuda_cargo_check_features_cuda_clean() {
    let res = cached_cuda_cargo_check();
    assert!(
        res.success,
        "cargo check -p gecko-cli --features cuda failed (exit {:?}):\nstderr:\n{}",
        res.exit_code, res.stderr
    );
    assert!(
        res.cfgs.iter().any(|c| c == "has_cuda_kernels"),
        "gecko-cuda build.rs did not emit has_cuda_kernels cfg: {:?}",
        res.cfgs
    );
    assert!(
        res.linked_libs
            .iter()
            .any(|l| l == "static=gecko_cuda_kernels"),
        "gecko-cuda build.rs did not link static=gecko_cuda_kernels: {:?}",
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

fn parse_cuobjdump_archive_members(
    stdout: &str,
    entry_marker: &str,
) -> BTreeMap<String, Vec<String>> {
    let mut by_member: BTreeMap<String, Vec<String>> = BTreeMap::new();
    let mut current_member: Option<String> = None;
    for line in stdout.lines() {
        let trimmed = line.trim();
        if let Some(rest) = trimmed.strip_prefix("member ") {
            let without_colon = rest.strip_suffix(':').unwrap_or(rest);
            let member_name = without_colon
                .rsplit_once(':')
                .map(|(_, obj)| obj)
                .unwrap_or(without_colon)
                .to_string();
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

#[test]
fn tier1_cuda_build_targets_sm75_sm80_sm86_sm89_and_compute89_ptx() {
    let res = cached_cuda_cargo_check();
    assert!(res.success, "cargo check --features cuda must succeed");
    let out_dir = res.out_dir.as_ref().expect("gecko-cuda out_dir");
    let archive = out_dir.join("libgecko_cuda_kernels.a");
    assert!(archive.exists(), "missing {}", archive.display());

    // Inspect embedded CUBINs and PTX across EVERY archive member in libgecko_cuda_kernels.a
    let elf_out = Command::new("cuobjdump")
        .args(["--list-elf", archive.to_str().unwrap()])
        .output()
        .expect("cuobjdump --list-elf");
    assert!(
        elf_out.status.success(),
        "cuobjdump --list-elf failed:\n{}",
        String::from_utf8_lossy(&elf_out.stderr)
    );
    let elf_stdout = String::from_utf8_lossy(&elf_out.stdout);
    let elf_by_member = parse_cuobjdump_archive_members(&elf_stdout, "ELF file");

    let ptx_out = Command::new("cuobjdump")
        .args(["--list-ptx", archive.to_str().unwrap()])
        .output()
        .expect("cuobjdump --list-ptx");
    assert!(
        ptx_out.status.success(),
        "cuobjdump --list-ptx failed:\n{}",
        String::from_utf8_lossy(&ptx_out.stderr)
    );
    let ptx_stdout = String::from_utf8_lossy(&ptx_out.stdout);
    let ptx_by_member = parse_cuobjdump_archive_members(&ptx_stdout, "PTX file");

    // Verify all 4 standard CUDA translation units plus any additional archive members
    let mut members: BTreeSet<String> = [
        "cuda_index.o",
        "cuda_query.o",
        "cuda_extend.o",
        "cuda_bridge.o",
    ]
    .into_iter()
    .map(String::from)
    .collect();
    members.extend(elf_by_member.keys().cloned());
    members.extend(ptx_by_member.keys().cloned());

    let mut missing = Vec::new();
    for member in &members {
        let elf_lines = elf_by_member.get(member);
        for req_arch in [
            ".sm_75.cubin",
            ".sm_80.cubin",
            ".sm_86.cubin",
            ".sm_89.cubin",
        ] {
            let has_arch = elf_lines.is_some_and(|ls| ls.iter().any(|l| l.ends_with(req_arch)));
            if !has_arch {
                missing.push(format!("{member}: missing {req_arch}"));
            }
        }
        let ptx_lines = ptx_by_member.get(member);
        let has_ptx89 = ptx_lines.is_some_and(|ls| ls.iter().any(|l| l.ends_with(".sm_89.ptx")));
        if !has_ptx89 {
            missing.push(format!("{member}: missing .sm_89.ptx"));
        }
    }

    assert!(
        missing.is_empty(),
        "libgecko_cuda_kernels.a is missing required per-member CUDA targets:\n  {}\n\ncuobjdump --list-elf:\n{elf_stdout}\ncuobjdump --list-ptx:\n{ptx_stdout}",
        missing.join("\n  ")
    );
}

#[test]
fn tier1_cuda_extend_kernels_zero_stack_sm75() {
    verify_cuda_extend_zero_stack_for_arch("sm_75");
}

#[test]
fn tier1_cuda_extend_kernels_zero_stack_sm80() {
    verify_cuda_extend_zero_stack_for_arch("sm_80");
}

#[test]
fn tier1_cuda_extend_kernels_zero_stack_sm86() {
    verify_cuda_extend_zero_stack_for_arch("sm_86");
}

// ============================================================================
// TIER 2 — Boundary, Corner & Negative Cases (15 tests; ALL baseline PASS)
// ============================================================================

#[test]
fn tier2_palindromic_kmer_and_palindromic_hsp() {
    let tmp = tempfile::tempdir().unwrap();
    // Construct an exact 120-bp DNA palindrome: half + revcomp(half)
    let half = pseudo_dna(9001, 60);
    let pal = format!("{half}{}", revcomp(&half));
    assert_eq!(pal, revcomp(&pal));

    let q_seq = format!("{}{}{}", pseudo_dna(9010, 100), pal, pseudo_dna(9011, 100));
    let r_seq = format!("{}{}{}", pseudo_dna(9020, 100), pal, pseudo_dna(9021, 100));
    let q = tmp.path().join("q.fa");
    let r = tmp.path().join("r.fa");
    write_fasta(&q, &[("q", &q_seq)]);
    write_fasta(&r, &[("r", &r_seq)]);

    let out_both = tmp.path().join("both");
    let out_fwd = tmp.path().join("fwd");
    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "12",
            "-l",
            "80",
            "-s",
            "90",
            "-f",
            "frags",
            "-o",
            out_both.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "palindromic HSP dual-strand",
    );
    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "12",
            "-l",
            "80",
            "-s",
            "90",
            "--forward-only",
            "-f",
            "frags",
            "-o",
            out_fwd.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "palindromic HSP forward-only",
    );

    let f_both = parse_frags_file(&out_both.join("q-r.frags")).fragments;
    assert_eq!(
        f_both.len(),
        2,
        "palindromic HSP must produce both 'f' and 'r'"
    );
    assert_eq!(f_both[0].strand, 'f');
    assert_eq!(f_both[1].strand, 'r');

    let f_fwd = parse_frags_file(&out_fwd.join("q-r.frags")).fragments;
    assert_eq!(f_fwd.len(), 1);
    assert_eq!(f_fwd[0].strand, 'f');
}

#[test]
fn tier2_freq_limit_exact_boundary_kept_vs_dropped() {
    let tmp = tempfile::tempdir().unwrap();
    // Ref contains block A exactly 2 times, block B exactly 3 times.
    // With --freq-limit 2: cnt(A)==2 <= 2 KEPT (2 frags), cnt(B)==3 > 2 DROPPED (0 frags).
    let a = pseudo_dna(9101, 100);
    let b = pseudo_dna(9102, 100);

    let ql = pseudo_dna(9110, 120);
    let qm = pseudo_dna(9111, 120);
    let qr = pseudo_dna(9112, 120);
    let a_l_mm = all_mismatch(&ql[ql.len() - 40..]);
    let a_r_mm = all_mismatch(&qm[..40]);
    let b_l_mm = all_mismatch(&qm[qm.len() - 40..]);
    let b_r_mm = all_mismatch(&qr[..40]);

    let q_seq = format!("{ql}{a}{qm}{b}{qr}");
    let r_seq = format!(
        "{}{}{a}{}{}{}{a}{}{}{}{b}{}{}{}{b}{}{}{}{b}{}{}",
        pseudo_dna(9120, 80),
        a_l_mm,
        a_r_mm,
        pseudo_dna(9121, 90),
        a_l_mm,
        a_r_mm,
        pseudo_dna(9122, 110),
        b_l_mm,
        b_r_mm,
        pseudo_dna(9123, 85),
        b_l_mm,
        b_r_mm,
        pseudo_dna(9124, 95),
        b_l_mm,
        b_r_mm,
        pseudo_dna(9125, 80),
    );

    let q = tmp.path().join("q.fa");
    let r = tmp.path().join("r.fa");
    write_fasta(&q, &[("q", &q_seq)]);
    write_fasta(&r, &[("r", &r_seq)]);

    let out = tmp.path().join("out");
    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "12",
            "-l",
            "60",
            "-s",
            "80",
            "--freq-limit",
            "2",
            "-f",
            "frags",
            "-o",
            out.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "freq_limit exact boundary",
    );

    let ff = parse_frags_file(&out.join("q-r.frags"));
    assert_eq!(
        ff.fragments.len(),
        2,
        "cnt==FL kept (2x A), cnt==FL+1 dropped (3x B)"
    );
    assert!(
        ff.fragments
            .iter()
            .all(|f| f.x_start == 120 && f.length == 100)
    );
}

#[test]
fn tier2_freq_limit_zero_drops_all_kmers() {
    let tmp = tempfile::tempdir().unwrap();
    let (q, r) = build_multi_hsp_pair(tmp.path());
    let out = tmp.path().join("out");

    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "12",
            "-l",
            "40",
            "-s",
            "50",
            "--freq-limit",
            "0",
            "-f",
            "frags",
            "-o",
            out.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "--freq-limit 0",
    );
    let ff = parse_frags_file(&out.join("q-r.frags"));
    assert_eq!(ff.fragments.len(), 0);
    assert_eq!(fs::read(out.join("q-r.frags")).unwrap().len(), 16);
}

#[test]
fn tier2_sequence_shorter_than_k() {
    let tmp = tempfile::tempdir().unwrap();
    let q = tmp.path().join("q.fa");
    let r = tmp.path().join("r.fa");
    write_fasta(&q, &[("q", "ACGTACGTAC")]); // 10 bp < k=12
    write_fasta(&r, &[("r", "ACGTACGTAC")]);

    let out = tmp.path().join("out");
    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "12",
            "-l",
            "5",
            "-s",
            "50",
            "-f",
            "frags",
            "-o",
            out.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "seq shorter than k",
    );
    let ff = parse_frags_file(&out.join("q-r.frags"));
    assert_eq!(ff.x_total, 10);
    assert_eq!(ff.y_total, 10);
    assert_eq!(ff.fragments.len(), 0);
}

#[test]
fn tier2_all_n_and_all_lowercase_sequences() {
    let tmp = tempfile::tempdir().unwrap();
    let q_n = tmp.path().join("qn.fa");
    let r_n = tmp.path().join("rn.fa");
    write_fasta(&q_n, &[("qn", &"N".repeat(200))]);
    write_fasta(&r_n, &[("rn", &"N".repeat(200))]);

    let out_n = tmp.path().join("out_n");
    assert_success(
        &run_gecko(&[
            "compare",
            q_n.to_str().unwrap(),
            r_n.to_str().unwrap(),
            "-k",
            "12",
            "-l",
            "40",
            "-s",
            "50",
            "-f",
            "frags",
            "-o",
            out_n.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "all-N sequences",
    );
    assert_eq!(
        parse_frags_file(&out_n.join("qn-rn.frags")).fragments.len(),
        0
    );

    let seq_lo = pseudo_dna(9201, 200).to_ascii_lowercase();
    let q_lo = tmp.path().join("qlo.fa");
    let r_lo = tmp.path().join("rlo.fa");
    write_fasta(&q_lo, &[("qlo", &seq_lo)]);
    write_fasta(&r_lo, &[("rlo", &seq_lo)]);

    let out_masked = tmp.path().join("masked");
    let out_unmasked = tmp.path().join("unmasked");
    assert_success(
        &run_gecko(&[
            "compare",
            q_lo.to_str().unwrap(),
            r_lo.to_str().unwrap(),
            "-k",
            "12",
            "-l",
            "40",
            "-s",
            "80",
            "-f",
            "frags",
            "-o",
            out_masked.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "all-lowercase default soft-mask",
    );
    assert_eq!(
        parse_frags_file(&out_masked.join("qlo-rlo.frags"))
            .fragments
            .len(),
        0
    );

    assert_success(
        &run_gecko(&[
            "compare",
            q_lo.to_str().unwrap(),
            r_lo.to_str().unwrap(),
            "-k",
            "12",
            "-l",
            "40",
            "-s",
            "80",
            "--no-soft-mask",
            "-f",
            "frags",
            "-o",
            out_unmasked.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "all-lowercase --no-soft-mask",
    );
    let ff = parse_frags_file(&out_unmasked.join("qlo-rlo.frags"));
    assert_eq!(ff.fragments.len(), 1);
    assert_eq!(ff.fragments[0].length, 200);
}

#[test]
fn tier2_strict_min_length_boundary() {
    let tmp = tempfile::tempdir().unwrap();
    // Exact 60-bp HSP bounded by all-mismatch flanks so length is deterministically 60
    let flank_q1 = pseudo_dna(9301, 80);
    let flank_q2 = pseudo_dna(9302, 80);
    let hsp = pseudo_dna(9303, 60);

    let q_seq = format!("{flank_q1}{hsp}{flank_q2}");
    let r_seq = format!(
        "{}{}{}",
        all_mismatch(&flank_q1),
        hsp,
        all_mismatch(&flank_q2)
    );
    let q = tmp.path().join("q.fa");
    let r = tmp.path().join("r.fa");
    write_fasta(&q, &[("q", &q_seq)]);
    write_fasta(&r, &[("r", &r_seq)]);

    let out_eq = tmp.path().join("l60");
    let out_lt = tmp.path().join("l59");
    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "12",
            "-l",
            "60",
            "-s",
            "80",
            "-f",
            "frags",
            "-o",
            out_eq.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "-l 60 strict > boundary",
    );
    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "12",
            "-l",
            "59",
            "-s",
            "80",
            "-f",
            "frags",
            "-o",
            out_lt.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "-l 59 strict > boundary",
    );

    assert_eq!(
        parse_frags_file(&out_eq.join("q-r.frags")).fragments.len(),
        0,
        "length == min_length (60) must be REJECTED (strict >)"
    );
    let ff = parse_frags_file(&out_lt.join("q-r.frags"));
    assert_eq!(ff.fragments.len(), 1);
    assert_eq!(ff.fragments[0].length, 60);
}

#[test]
fn tier2_strict_min_similarity_boundary_and_score_formula() {
    let tmp = tempfile::tempdir().unwrap();
    // Exact 100-bp HSP with 10 interior mismatches (90% identity -> score=320, similarity=80.0%)
    let flank_q1 = pseudo_dna(9401, 80);
    let flank_q2 = pseudo_dna(9402, 80);
    let hsp = pseudo_dna(9403, 100);
    // 10 interior mismatches spaced cleanly so first 15 bp & last 15 bp are exact matches
    let mut_pos = [16, 23, 30, 37, 44, 51, 58, 65, 72, 79];
    let hsp_mut = mutate_positions(&hsp, &mut_pos);

    let q_seq = format!("{flank_q1}{hsp}{flank_q2}");
    let r_seq = format!(
        "{}{}{}",
        all_mismatch(&flank_q1),
        hsp_mut,
        all_mismatch(&flank_q2)
    );
    let q = tmp.path().join("q.fa");
    let r = tmp.path().join("r.fa");
    write_fasta(&q, &[("q", &q_seq)]);
    write_fasta(&r, &[("r", &r_seq)]);

    let out_eq = tmp.path().join("s80");
    let out_lt = tmp.path().join("s799");
    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "12",
            "-l",
            "50",
            "-s",
            "80",
            "-f",
            "frags",
            "-o",
            out_eq.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "-s 80 strict >",
    );
    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "12",
            "-l",
            "50",
            "-s",
            "79.9",
            "-f",
            "frags",
            "-o",
            out_lt.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "-s 79.9 strict >",
    );

    assert_eq!(
        parse_frags_file(&out_eq.join("q-r.frags")).fragments.len(),
        0,
        "similarity == min_similarity (80.0) must be REJECTED (strict >)"
    );
    let ff = parse_frags_file(&out_lt.join("q-r.frags"));
    assert_eq!(ff.fragments.len(), 1);
    let f = &ff.fragments[0];
    assert_eq!(f.length, 100);
    assert_eq!(f.ident, 90);
    assert_eq!(f.score, 320);
    assert!(
        (f.similarity - 80.0).abs() < 1e-4,
        "expected 80.0, got {}",
        f.similarity
    );
}

#[test]
fn tier2_extreme_word_lengths_k4_and_k32() {
    let tmp = tempfile::tempdir().unwrap();
    let hsp = pseudo_dna(9501, 140);
    let q = tmp.path().join("q.fa");
    let r = tmp.path().join("r.fa");
    write_fasta(
        &q,
        &[(
            "q",
            &format!("{}{}{}", pseudo_dna(9510, 120), hsp, pseudo_dna(9511, 120)),
        )],
    );
    write_fasta(
        &r,
        &[(
            "r",
            &format!("{}{}{}", pseudo_dna(9520, 120), hsp, pseudo_dna(9521, 120)),
        )],
    );

    for k_str in ["4", "32"] {
        let out = tmp.path().join(format!("k{k_str}"));
        assert_success(
            &run_gecko(&[
                "compare",
                q.to_str().unwrap(),
                r.to_str().unwrap(),
                "-k",
                k_str,
                "-l",
                "100",
                "-s",
                "90",
                "-f",
                "frags",
                "-o",
                out.to_str().unwrap(),
                "--backend",
                "cpu",
            ]),
            &format!("extreme word_length k={k_str}"),
        );
        let ff = parse_frags_file(&out.join("q-r.frags"));
        assert_eq!(ff.fragments.len(), 1);
        assert!(ff.fragments[0].length >= 140);
    }
}

#[test]
fn tier2_multifasta_seq_sep_prevents_cross_record_bleed() {
    let tmp = tempfile::tempdir().unwrap();
    let a = pseudo_dna(9601, 80);
    let b = pseudo_dna(9602, 80);
    let ab_concat = format!("{a}{b}");
    let ab_spaced = format!("{a}C{b}");
    let a_rc = revcomp(&a);
    let b_rc = revcomp(&b);

    let run_case =
        |sub: &str, q_recs: &[(&str, &str)], r_recs: &[(&str, &str)]| -> ParsedFragsFile {
            let dir = tmp.path().join(sub);
            fs::create_dir_all(&dir).unwrap();
            let q = dir.join("q.fa");
            let r = dir.join("r.fa");
            let out = dir.join("out");
            write_fasta(&q, q_recs);
            write_fasta(&r, r_recs);
            assert_success(
                &run_gecko(&[
                    "compare",
                    q.to_str().unwrap(),
                    r.to_str().unwrap(),
                    "-k",
                    "12",
                    "-l",
                    "50",
                    "-s",
                    "90",
                    "-f",
                    "frags",
                    "-o",
                    out.to_str().unwrap(),
                    "--backend",
                    "cpu",
                ]),
                &format!("SEQ_SEP boundary {sub}"),
            );
            parse_frags_file(&out.join("q-r.frags"))
        };

    // Case A (0-width concatenation check): q = [a+b] (160 bp) vs r = [a, b] (SEQ_SEP at raw y=80).
    // Note: `.frags` stores base-space coordinates (excluding SEQ_SEP) and separator-space `diag`.
    let ff_a = run_case("case_a", &[("q0", &ab_concat)], &[("r0", &a), ("r1", &b)]);
    assert_eq!(ff_a.x_total, 160);
    assert_eq!(ff_a.y_total, 160);
    assert_eq!(
        ff_a.fragments.len(),
        2,
        "case A: must emit 2 separate 80-bp frags, never 1 merged 160-bp frag"
    );
    assert!(
        ff_a.fragments
            .iter()
            .all(|f| f.length == 80 && f.strand == 'f')
    );
    assert!(ff_a.fragments.iter().any(|f| {
        f.seq_x == 0
            && f.seq_y == 0
            && f.x_start == 0
            && f.x_end == 79
            && f.y_start == 0
            && f.y_end == 79
            && f.diag == 0
    }));
    assert!(ff_a.fragments.iter().any(|f| {
        f.seq_x == 0
            && f.seq_y == 1
            && f.x_start == 80
            && f.x_end == 159
            && f.y_start == 80
            && f.y_end == 159
            && f.diag == -1
    }));

    // Case B (same-diagonal vy == SEQ_SEP on diag = 0):
    // q = [a + "C" + b] (161 bp single record) vs r = [a, b] (SEQ_SEP at raw y=80).
    // Both a and b sit on separator-space diag = 0; if extension fails to stop at vy == SEQ_SEP,
    // fscore crosses (80,80) and merges into one 161-bp fragment spanning x=80!
    let ff_b = run_case("case_b", &[("q0", &ab_spaced)], &[("r0", &a), ("r1", &b)]);
    assert_eq!(ff_b.x_total, 161);
    assert_eq!(ff_b.y_total, 160);
    assert_eq!(
        ff_b.fragments.len(),
        2,
        "case B (vy==SEQ_SEP on diag=0): must emit 2 separate 80-bp frags, got {:?}",
        ff_b.fragments
    );
    for f in &ff_b.fragments {
        assert_eq!(f.strand, 'f');
        assert_eq!(f.diag, 0, "case B: both halves must lie on diag=0: {f:?}");
        assert_eq!(
            f.length, 80,
            "case B: fragment must not bleed across SEQ_SEP: {f:?}"
        );
        assert!(
            !(f.x_start <= 80 && f.x_end >= 80),
            "case B: fragment spans spacer index x=80 (SEQ_SEP): {f:?}"
        );
    }
    assert!(ff_b.fragments.iter().any(|f| {
        f.seq_x == 0
            && f.seq_y == 0
            && f.x_start == 0
            && f.x_end == 79
            && f.y_start == 0
            && f.y_end == 79
    }));
    assert!(ff_b.fragments.iter().any(|f| {
        f.seq_x == 0
            && f.seq_y == 1
            && f.x_start == 81
            && f.x_end == 160
            && f.y_start == 80
            && f.y_end == 159
    }));

    // Case C (same-diagonal vx == SEQ_SEP on diag = 0):
    // Symmetric: q = [a, b] (SEQ_SEP at raw x=80) vs r = [a + "C" + b] (single record, empty special_runs).
    let ff_c = run_case("case_c", &[("q0", &a), ("q1", &b)], &[("r0", &ab_spaced)]);
    assert_eq!(ff_c.x_total, 160);
    assert_eq!(ff_c.y_total, 161);
    assert_eq!(
        ff_c.fragments.len(),
        2,
        "case C (vx==SEQ_SEP on diag=0): must emit 2 separate 80-bp frags, got {:?}",
        ff_c.fragments
    );
    for f in &ff_c.fragments {
        assert_eq!(f.strand, 'f');
        assert_eq!(f.diag, 0, "case C: both halves must lie on diag=0: {f:?}");
        assert_eq!(
            f.length, 80,
            "case C: fragment must not bleed across SEQ_SEP: {f:?}"
        );
        assert!(
            !(f.y_start <= 80 && f.y_end >= 80),
            "case C: fragment spans spacer index y=80 (SEQ_SEP): {f:?}"
        );
    }
    assert!(ff_c.fragments.iter().any(|f| {
        f.seq_x == 0
            && f.seq_y == 0
            && f.x_start == 0
            && f.x_end == 79
            && f.y_start == 0
            && f.y_end == 79
    }));
    assert!(ff_c.fragments.iter().any(|f| {
        f.seq_x == 1
            && f.seq_y == 0
            && f.x_start == 80
            && f.x_end == 159
            && f.y_start == 81
            && f.y_end == 160
    }));

    // Case D (same-diagonal vx == SEQ_SEP && vy == SEQ_SEP on diag = 0):
    // Catches broken 2-bit-only fast paths that compare SEQ_SEP (00) == SEQ_SEP (00) as a +4 match!
    let ff_d = run_case(
        "case_d",
        &[("q0", &a), ("q1", &b)],
        &[("r0", &a), ("r1", &b)],
    );
    assert_eq!(ff_d.x_total, 160);
    assert_eq!(ff_d.y_total, 160);
    assert_eq!(
        ff_d.fragments.len(),
        2,
        "case D (vx==vy==SEQ_SEP on diag=0): must emit 2 separate 80-bp frags, got {:?}",
        ff_d.fragments
    );
    for f in &ff_d.fragments {
        assert_eq!(f.strand, 'f');
        assert_eq!(f.diag, 0, "case D: both halves must lie on diag=0: {f:?}");
        assert_eq!(
            f.length, 80,
            "case D: fragment must not bleed across SEQ_SEP: {f:?}"
        );
    }
    assert!(ff_d.fragments.iter().any(|f| {
        f.seq_x == 0
            && f.seq_y == 0
            && f.x_start == 0
            && f.x_end == 79
            && f.y_start == 0
            && f.y_end == 79
    }));
    assert!(ff_d.fragments.iter().any(|f| {
        f.seq_x == 1
            && f.seq_y == 1
            && f.x_start == 80
            && f.x_end == 159
            && f.y_start == 80
            && f.y_end == 159
    }));

    // Case E (reverse-strand same-diagonal vy_rc == SEQ_SEP on diag = 0):
    // revcomp([r0=rc(b), r1=rc(a)]) == a + SEQ_SEP + b -> tests SEQ_SEP preservation in reverse_complement().
    let ff_e = run_case(
        "case_e",
        &[("q0", &ab_spaced)],
        &[("r0", &b_rc), ("r1", &a_rc)],
    );
    assert_eq!(ff_e.x_total, 161);
    assert_eq!(ff_e.y_total, 160);
    assert_eq!(
        ff_e.fragments.len(),
        2,
        "case E (reverse-strand vy_rc==SEQ_SEP on diag=0): must emit 2 separate 80-bp frags, got {:?}",
        ff_e.fragments
    );
    for f in &ff_e.fragments {
        assert_eq!(f.strand, 'r');
        assert_eq!(
            f.diag, 0,
            "case E: both reverse halves must lie on diag=0: {f:?}"
        );
        assert_eq!(
            f.length, 80,
            "case E: reverse fragment must not bleed across SEQ_SEP: {f:?}"
        );
        assert!(
            !(f.x_start <= 80 && f.x_end >= 80),
            "case E: reverse fragment spans spacer index x=80 (SEQ_SEP): {f:?}"
        );
    }
    assert!(ff_e.fragments.iter().any(|f| {
        f.seq_x == 0
            && f.seq_y == 1
            && f.x_start == 0
            && f.x_end == 79
            && f.y_start == 0
            && f.y_end == 79
    }));
    assert!(ff_e.fragments.iter().any(|f| {
        f.seq_x == 0
            && f.seq_y == 0
            && f.x_start == 81
            && f.x_end == 160
            && f.y_start == 80
            && f.y_end == 159
    }));
}

#[test]
fn tier2_multifasta_empty_middle_and_first_records() {
    let tmp = tempfile::tempdir().unwrap();
    let h1 = pseudo_dna(9701, 100);
    let h3 = pseudo_dna(9702, 110);

    let q = tmp.path().join("q.fa");
    let r = tmp.path().join("r.fa");
    write_fasta(&q, &[("s0", ""), ("s1", &h1), ("s2", ""), ("s3", &h3)]);
    write_fasta(&r, &[("t0", ""), ("t1", &h1), ("t2", ""), ("t3", &h3)]);

    let out = tmp.path().join("out");
    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "12",
            "-l",
            "60",
            "-s",
            "90",
            "-f",
            "frags",
            "-o",
            out.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "multi-FASTA empty records",
    );

    let ff = parse_frags_file(&out.join("q-r.frags"));
    assert_eq!(ff.x_total, 210);
    assert_eq!(ff.y_total, 210);
    assert_eq!(ff.fragments.len(), 2);
    assert!(
        ff.fragments
            .iter()
            .any(|f| f.seq_x == 1 && f.seq_y == 1 && f.length == 100)
    );
    assert!(
        ff.fragments
            .iter()
            .any(|f| f.seq_x == 3 && f.seq_y == 3 && f.length == 110)
    );
}

#[test]
fn tier2_extension_drop_to_zero_boundary() {
    let tmp = tempfile::tempdir().unwrap();
    // k=8, match=+4, mismatch=-4.
    // Head: 8 exact matches -> seed_score = +32.
    // Dip: 8 mismatches -> fscore = 0 (CONTINUES!) vs 9 mismatches -> fscore = -4 (STOPS!).
    // Tail: 27 bp of 3*(6 matches + 1 mismatch) + 6 matches -> NO 8-mer seeds in tail,
    //       net +84 over 27 bp -> peak reaches +84 > 32, total len = 43 bp, sim = 48.8% > 40%.
    let flank_l = pseudo_dna(9801, 40);
    let head = pseudo_dna(9802, 8);
    let dip8 = pseudo_dna(9803, 8);
    let dip9 = pseudo_dna(9804, 9);
    let tail = pseudo_dna(9805, 27);
    let flank_r = pseudo_dna(9806, 40);

    let tail_m = mutate_positions(&tail, &[6, 13, 20]); // max run = 6 < 8

    // Variant A: dip to exact 0
    let qa = tmp.path().join("qa.fa");
    let ra = tmp.path().join("ra.fa");
    write_fasta(
        &qa,
        &[("qa", &format!("{flank_l}{head}{dip8}{tail}{flank_r}"))],
    );
    write_fasta(
        &ra,
        &[(
            "ra",
            &format!(
                "{}{}{}{}{}",
                all_mismatch(&flank_l),
                head,
                all_mismatch(&dip8),
                tail_m,
                all_mismatch(&flank_r)
            ),
        )],
    );

    // Variant B: dip to -4
    let qb = tmp.path().join("qb.fa");
    let rb = tmp.path().join("rb.fa");
    write_fasta(
        &qb,
        &[("qb", &format!("{flank_l}{head}{dip9}{tail}{flank_r}"))],
    );
    write_fasta(
        &rb,
        &[(
            "rb",
            &format!(
                "{}{}{}{}{}",
                all_mismatch(&flank_l),
                head,
                all_mismatch(&dip9),
                tail_m,
                all_mismatch(&flank_r)
            ),
        )],
    );

    let out_a = tmp.path().join("out_a");
    let out_b = tmp.path().join("out_b");
    assert_success(
        &run_gecko(&[
            "compare",
            qa.to_str().unwrap(),
            ra.to_str().unwrap(),
            "-k",
            "8",
            "-l",
            "35",
            "-s",
            "40",
            "--forward-only",
            "-f",
            "frags",
            "-o",
            out_a.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "fscore == 0 boundary",
    );
    assert_success(
        &run_gecko(&[
            "compare",
            qb.to_str().unwrap(),
            rb.to_str().unwrap(),
            "-k",
            "8",
            "-l",
            "35",
            "-s",
            "40",
            "--forward-only",
            "-f",
            "frags",
            "-o",
            out_b.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "fscore < 0 boundary",
    );

    let fa = parse_frags_file(&out_a.join("qa-ra.frags")).fragments;
    let fb = parse_frags_file(&out_b.join("qb-rb.frags")).fragments;
    assert_eq!(
        fa.len(),
        1,
        "fscore == 0 must continue extending and reach peak 84"
    );
    assert_eq!(fa[0].length, 43);
    assert_eq!(fa[0].score, 84);
    assert_eq!(
        fb.len(),
        0,
        "fscore == -4 (< 0) must terminate extension immediately"
    );
}

#[test]
fn tier2_negative_cli_flags_formats_memory_missing_empty_fasta() {
    let tmp = tempfile::tempdir().unwrap();
    let valid_fa = tmp.path().join("valid.fa");
    let empty_fa = tmp.path().join("empty.fa");
    let missing_fa = tmp.path().join("does_not_exist.fa");
    write_fasta(&valid_fa, &[("v", &pseudo_dna(9901, 100))]);
    fs::write(&empty_fa, "").unwrap();
    let out = tmp.path().join("out");

    let cases: &[(&str, &[&str])] = &[
        ("unknown flag", &["compare", "--no-such-flag"]),
        (
            "missing ref positional",
            &["compare", valid_fa.to_str().unwrap()],
        ),
        (
            "invalid --format",
            &[
                "compare",
                valid_fa.to_str().unwrap(),
                valid_fa.to_str().unwrap(),
                "--format",
                "not_a_format",
                "-o",
                out.to_str().unwrap(),
                "--backend",
                "cpu",
            ],
        ),
        (
            "invalid -m",
            &[
                "compare",
                valid_fa.to_str().unwrap(),
                valid_fa.to_str().unwrap(),
                "-m",
                "10X",
                "-o",
                out.to_str().unwrap(),
                "--backend",
                "cpu",
            ],
        ),
        (
            "missing FASTA file",
            &[
                "compare",
                missing_fa.to_str().unwrap(),
                valid_fa.to_str().unwrap(),
                "-o",
                out.to_str().unwrap(),
                "--backend",
                "cpu",
            ],
        ),
        (
            "empty FASTA file",
            &[
                "compare",
                empty_fa.to_str().unwrap(),
                valid_fa.to_str().unwrap(),
                "-o",
                out.to_str().unwrap(),
                "--backend",
                "cpu",
            ],
        ),
    ];

    for (label, args) in cases {
        let res = run_gecko(args);
        assert_clean_failure(&res, label);
    }
}

#[test]
fn tier2_negative_invalid_spaced_seed_masks() {
    let tmp = tempfile::tempdir().unwrap();
    let valid_fa = tmp.path().join("v.fa");
    write_fasta(&valid_fa, &[("v", &pseudo_dna(9911, 100))]);
    let out = tmp.path().join("out");
    let too_long = "1".repeat(33);

    for (label, mask) in [
        ("asymmetric", "111011"),
        ("invalid char", "11211"),
        ("span < 4", "111"),
        ("span > 32", too_long.as_str()),
        ("leading/trailing 0", "011110"),
        ("weight < 4", "100001"),
    ] {
        let res = run_gecko(&[
            "compare",
            valid_fa.to_str().unwrap(),
            valid_fa.to_str().unwrap(),
            "--spaced-seed",
            mask,
            "-o",
            out.to_str().unwrap(),
            "--backend",
            "cpu",
        ]);
        assert_clean_failure(&res, &format!("invalid spaced seed ({label}: {mask})"));
    }
}

#[test]
fn tier2_negative_invalid_numeric_ranges_and_all_vs_all() {
    let tmp = tempfile::tempdir().unwrap();
    let valid_fa = tmp.path().join("v.fa");
    write_fasta(&valid_fa, &[("v", &pseudo_dna(9921, 100))]);
    let out = tmp.path().join("out");

    let one_dir = tmp.path().join("one_genome");
    fs::create_dir_all(&one_dir).unwrap();
    write_fasta(
        &one_dir.join("only.fa"),
        &[("only", &pseudo_dna(9922, 100))],
    );
    let missing_dir = tmp.path().join("no_such_dir");

    let cases: &[(&str, &[&str])] = &[
        (
            "-k 3 (< 4)",
            &[
                "compare",
                valid_fa.to_str().unwrap(),
                valid_fa.to_str().unwrap(),
                "-k",
                "3",
                "-o",
                out.to_str().unwrap(),
                "--backend",
                "cpu",
            ],
        ),
        (
            "-k 33 (> 32)",
            &[
                "compare",
                valid_fa.to_str().unwrap(),
                valid_fa.to_str().unwrap(),
                "-k",
                "33",
                "-o",
                out.to_str().unwrap(),
                "--backend",
                "cpu",
            ],
        ),
        (
            "-s 101 (> 100)",
            &[
                "compare",
                valid_fa.to_str().unwrap(),
                valid_fa.to_str().unwrap(),
                "-s",
                "101",
                "-o",
                out.to_str().unwrap(),
                "--backend",
                "cpu",
            ],
        ),
        (
            "-s=-5 (< 0)",
            &[
                "compare",
                valid_fa.to_str().unwrap(),
                valid_fa.to_str().unwrap(),
                "-s=-5",
                "-o",
                out.to_str().unwrap(),
                "--backend",
                "cpu",
            ],
        ),
        (
            "--match-score 0",
            &[
                "compare",
                valid_fa.to_str().unwrap(),
                valid_fa.to_str().unwrap(),
                "--match-score",
                "0",
                "-o",
                out.to_str().unwrap(),
                "--backend",
                "cpu",
            ],
        ),
        (
            "all-vs-all < 2 genomes",
            &[
                "all-vs-all",
                one_dir.to_str().unwrap(),
                "-o",
                out.to_str().unwrap(),
                "--backend",
                "cpu",
            ],
        ),
        (
            "all-vs-all missing dir",
            &[
                "all-vs-all",
                missing_dir.to_str().unwrap(),
                "-o",
                out.to_str().unwrap(),
                "--backend",
                "cpu",
            ],
        ),
    ];

    for (label, args) in cases {
        assert_clean_failure(&run_gecko(args), label);
    }
}

#[test]
fn tier2_negative_invalid_m3_flags() {
    let tmp = tempfile::tempdir().unwrap();
    let valid_fa = tmp.path().join("v.fa");
    write_fasta(&valid_fa, &[("v", &pseudo_dna(9931, 100))]);
    let out = tmp.path().join("out");

    for (label, flag_args) in [
        ("--seeding bogus", vec!["--seeding", "bogus"]),
        ("--two-hit=notanumber", vec!["--two-hit=notanumber"]),
    ] {
        let mut args = vec![
            "compare",
            valid_fa.to_str().unwrap(),
            valid_fa.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--backend",
            "cpu",
        ];
        args.extend(flag_args);
        assert_clean_failure(&run_gecko(&args), label);
    }
}

// ============================================================================
// TIER 3 — Cross-Feature Combinations (8 tests: 4 baseline + 4 M3 contract)
// ============================================================================

#[test]
fn tier3_minimizer_forward_only_freq_limit_multifasta_softmask_parity() {
    let tmp = tempfile::tempdir().unwrap();
    let h1 = pseudo_dna(10001, 160);
    let h2 = format!(
        "{}{}{}",
        pseudo_dna(10002, 40).to_ascii_lowercase(),
        pseudo_dna(10003, 40),
        pseudo_dna(10004, 40).to_ascii_lowercase()
    );
    let h_rev = pseudo_dna(10005, 150);
    let rep = pseudo_dna(10006, 90);

    let q = tmp.path().join("q.fa");
    let r = tmp.path().join("r.fa");
    write_fasta(
        &q,
        &[
            (
                "q0",
                &format!(
                    "{}{}{}{}",
                    pseudo_dna(10010, 100),
                    h1,
                    pseudo_dna(10011, 100),
                    rep
                ),
            ),
            (
                "q1",
                &format!(
                    "{}{}{}{}",
                    pseudo_dna(10012, 100),
                    h2,
                    pseudo_dna(10013, 100),
                    h_rev
                ),
            ),
        ],
    );
    write_fasta(
        &r,
        &[
            (
                "r0",
                &format!(
                    "{}{}{}{}{}{}",
                    pseudo_dna(10020, 100),
                    rep,
                    pseudo_dna(10021, 50),
                    rep,
                    pseudo_dna(10022, 50),
                    rep
                ),
            ),
            (
                "r1",
                &format!(
                    "{}{}{}{}",
                    pseudo_dna(10023, 100),
                    h1,
                    pseudo_dna(10024, 100),
                    revcomp(&h_rev)
                ),
            ),
            (
                "r2",
                &format!("{}{}{}", pseudo_dna(10025, 100), h2, pseudo_dna(10026, 100)),
            ),
        ],
    );

    let out_exact = tmp.path().join("exact");
    let out_min = tmp.path().join("min");
    for (mode, outdir) in [("exact", &out_exact), ("minimizer", &out_min)] {
        assert_success(
            &run_gecko(&[
                "compare",
                q.to_str().unwrap(),
                r.to_str().unwrap(),
                "-k",
                "20",
                "-l",
                "80",
                "-s",
                "80",
                "--forward-only",
                "--freq-limit",
                "2",
                "--seeding",
                mode,
                "-f",
                "frags",
                "-o",
                outdir.to_str().unwrap(),
                "--backend",
                "cpu",
            ]),
            &format!("tier3 minimizer combo ({mode})"),
        );
    }

    assert_frags_bytes_identical(
        &out_exact.join("q-r.frags"),
        &out_min.join("q-r.frags"),
        "minimizer + forward-only + freq-limit + multi-FASTA + soft-mask",
    );
    let ff = parse_frags_file(&out_min.join("q-r.frags"));
    assert_eq!(
        ff.fragments.len(),
        2,
        "expected h1 + h2 forward only (rep & rev excluded)"
    );
    assert!(ff.fragments.iter().all(|f| f.strand == 'f'));
}

#[test]
fn tier3_minimizer_constrained_memory_budget_parity() {
    let tmp = tempfile::tempdir().unwrap();
    let hsp_f = pseudo_dna(10101, 240);
    let hsp_r = pseudo_dna(10102, 220);
    let q_seq = format!(
        "{}{}{}{}{}",
        pseudo_dna(10110, 20_000),
        hsp_f,
        pseudo_dna(10111, 20_000),
        hsp_r,
        pseudo_dna(10112, 20_000)
    );
    let r_seq = format!(
        "{}{}{}{}{}",
        pseudo_dna(10120, 20_000),
        hsp_f,
        pseudo_dna(10121, 20_000),
        revcomp(&hsp_r),
        pseudo_dna(10122, 20_000)
    );

    let q = tmp.path().join("q.fa");
    let r = tmp.path().join("r.fa");
    write_fasta(&q, &[("q", &q_seq)]);
    write_fasta(&r, &[("r", &r_seq)]);

    let out_exact = tmp.path().join("exact_2g");
    let out_min = tmp.path().join("min_64k");
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
            "80",
            "--freq-limit",
            "5",
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
        "exact -m 2G",
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
            "80",
            "--freq-limit",
            "5",
            "--seeding",
            "minimizer",
            "-m",
            "64K",
            "-f",
            "frags",
            "-o",
            out_min.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "minimizer -m 64K",
    );

    let fe = canonicalize_frags(parse_frags_file(&out_exact.join("q-r.frags")).fragments);
    let fm = canonicalize_frags(parse_frags_file(&out_min.join("q-r.frags")).fragments);
    assert_eq!(fe.len(), 2);
    assert_eq!(fe, fm);
}

#[test]
fn tier3_two_hit_spaced_seed_format_all() {
    let tmp = tempfile::tempdir().unwrap();
    // w14s20 = "11110101100110101111"
    let hsp = pseudo_dna(10201, 200);
    let mut mut_pos = Vec::new();
    for blk in 0..10 {
        for &dc in &[4usize, 10, 15] {
            mut_pos.push(blk * 20 + dc);
        }
    }
    let hsp_div = mutate_positions(&hsp, &mut_pos);

    let q = tmp.path().join("q.fa");
    let r = tmp.path().join("r.fa");
    write_fasta(
        &q,
        &[(
            "q",
            &format!(
                "{}{}{}",
                pseudo_dna(10210, 150),
                hsp,
                pseudo_dna(10211, 150)
            ),
        )],
    );
    write_fasta(
        &r,
        &[(
            "r",
            &format!(
                "{}{}{}",
                pseudo_dna(10220, 150),
                hsp_div,
                pseudo_dna(10221, 150)
            ),
        )],
    );

    let out = tmp.path().join("all");
    assert_success(
        &run_gecko(&[
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
            "-o",
            out.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "--two-hit=40 --spaced-seed w14s20 --format all",
    );

    for ext in [
        "q-r.frags",
        "q-r.csv",
        "q-r.frags.INF",
        "q-r.frags.MAT",
        "q-r.paf",
        "q-r.blast.tsv",
        "q-r.bedpe",
        "q-r.jsonl",
        "q-r.png",
    ] {
        assert!(out.join(ext).exists(), "missing {ext}");
    }
    let ff = parse_frags_file(&out.join("q-r.frags"));
    assert_eq!(ff.fragments.len(), 1);
    let mat = parse_mat_file(&out.join("q-r.frags.MAT"));
    assert_mat_covers_frags(&mat, &ff.fragments);
    assert_eq!(
        fs::read_to_string(out.join("q-r.paf"))
            .unwrap()
            .lines()
            .count(),
        1
    );
    assert_eq!(
        fs::read_to_string(out.join("q-r.jsonl"))
            .unwrap()
            .lines()
            .count(),
        1
    );
}

#[test]
fn tier3_gecko_profile_all_vs_all_mat_coverage() {
    let tmp = tempfile::tempdir().unwrap();
    let gdir = tmp.path().join("genomes");
    fs::create_dir_all(&gdir).unwrap();
    let shared = pseudo_dna(10301, 240);
    for (name, seed) in [("a.fa", 10310), ("b.fa", 10320), ("c.fa", 10330)] {
        let seq = format!(
            "{}{}{}",
            pseudo_dna(seed, 200),
            shared,
            pseudo_dna(seed + 1, 200)
        );
        write_fasta(&gdir.join(name), &[(name, &seq)]);
    }

    let out_def = tmp.path().join("def");
    let out_prof = tmp.path().join("prof");
    let res_def = run_gecko_with_env(
        &[
            "all-vs-all",
            gdir.to_str().unwrap(),
            "-k",
            "12",
            "-l",
            "100",
            "-s",
            "80",
            "-f",
            "legacy",
            "-o",
            out_def.to_str().unwrap(),
            "--backend",
            "cpu",
        ],
        &[("GECKO_PROFILE", None)],
    );
    assert_success(&res_def, "all-vs-all default");

    let res_prof = run_gecko_with_env(
        &[
            "all-vs-all",
            gdir.to_str().unwrap(),
            "-k",
            "12",
            "-l",
            "100",
            "-s",
            "80",
            "-f",
            "legacy",
            "-o",
            out_prof.to_str().unwrap(),
            "--backend",
            "cpu",
        ],
        &[("GECKO_PROFILE", Some("1"))],
    );
    assert_success(&res_prof, "all-vs-all GECKO_PROFILE=1");

    let s_def = String::from_utf8_lossy(&res_def.stderr);
    let s_prof = String::from_utf8_lossy(&res_prof.stderr);
    assert!(
        s_prof.lines().count() > s_def.lines().count(),
        "GECKO_PROFILE=1 in all-vs-all must emit extra profiling lines on stderr:\ndef:\n{s_def}\nprof:\n{s_prof}"
    );

    for pair in ["a-b", "a-c", "b-c"] {
        assert_frags_bytes_identical(
            &out_def.join(format!("{pair}.frags")),
            &out_prof.join(format!("{pair}.frags")),
            &format!("all-vs-all GECKO_PROFILE=1 parity for {pair}"),
        );
        let ff = parse_frags_file(&out_prof.join(format!("{pair}.frags")));
        assert_eq!(ff.fragments.len(), 1);
        let mat = parse_mat_file(&out_prof.join(format!("{pair}.frags.MAT")));
        assert_mat_covers_frags(&mat, &ff.fragments);
    }
}

#[test]
fn tier3_direct_csr_k10_vs_hash_csr_k16_k32_mat_exactness() {
    let tmp = tempfile::tempdir().unwrap();
    let (q, r) = build_multi_hsp_pair(tmp.path());

    let mut frags_by_k = Vec::new();
    for k_str in ["10", "16", "32"] {
        let out_leg = tmp.path().join(format!("leg_{k_str}"));
        let out_f = tmp.path().join(format!("frags_{k_str}"));
        assert_success(
            &run_gecko(&[
                "compare",
                q.to_str().unwrap(),
                r.to_str().unwrap(),
                "-k",
                k_str,
                "-l",
                "100",
                "-s",
                "70",
                "--freq-limit",
                "2",
                "-f",
                "legacy",
                "-o",
                out_leg.to_str().unwrap(),
                "--backend",
                "cpu",
            ]),
            &format!("legacy k={k_str}"),
        );
        assert_success(
            &run_gecko(&[
                "compare",
                q.to_str().unwrap(),
                r.to_str().unwrap(),
                "-k",
                k_str,
                "-l",
                "100",
                "-s",
                "70",
                "--freq-limit",
                "2",
                "-f",
                "frags",
                "-o",
                out_f.to_str().unwrap(),
                "--backend",
                "cpu",
            ]),
            &format!("frags k={k_str}"),
        );
        assert!(!out_f.join("q-r.frags.MAT").exists());
        assert_frags_bytes_identical(
            &out_leg.join("q-r.frags"),
            &out_f.join("q-r.frags"),
            &format!("k={k_str} legacy vs frags"),
        );
        let ff = parse_frags_file(&out_leg.join("q-r.frags"));
        let mat = parse_mat_file(&out_leg.join("q-r.frags.MAT"));
        assert_mat_covers_frags(&mat, &ff.fragments);
        assert_eq!(ff.fragments.len(), 3);
        frags_by_k.push(ff.fragments);
    }
    // All 3 HSPs contain >= 32bp exact runs (hsp_f2_mut mutates every 40bp -> 39bp exact runs),
    // so DirectCsr k=10, HashCsr k=16, and HashCsr k=32 all recover the exact same 3 fragments.
    assert_eq!(
        frags_by_k[0], frags_by_k[1],
        "DirectCsr k=10 must match HashCsr k=16"
    );
    assert_eq!(
        frags_by_k[0], frags_by_k[2],
        "DirectCsr k=10 must match HashCsr k=32"
    );
}

#[test]
fn tier3_spaced_seed_soft_mask_forward_only_multifasta_all_formats() {
    let tmp = tempfile::tempdir().unwrap();
    // w14s20 = "11110101100110101111"
    let h_up = pseudo_dna(10501, 180);
    let h_lo = pseudo_dna(10502, 160).to_ascii_lowercase();
    let h_rev = pseudo_dna(10503, 160);

    let q = tmp.path().join("q.fa");
    let r = tmp.path().join("r.fa");
    write_fasta(
        &q,
        &[
            (
                "q0",
                &format!("{}{}{}", pseudo_dna(10510, 80), h_up, pseudo_dna(10511, 80)),
            ),
            (
                "q1",
                &format!(
                    "{}{}{}{}",
                    pseudo_dna(10512, 80),
                    h_lo,
                    pseudo_dna(10513, 80),
                    h_rev
                ),
            ),
        ],
    );
    write_fasta(
        &r,
        &[
            (
                "r0",
                &format!("{}{}{}", pseudo_dna(10520, 80), h_lo, pseudo_dna(10521, 80)),
            ),
            (
                "r1",
                &format!(
                    "{}{}{}{}",
                    pseudo_dna(10522, 80),
                    h_up,
                    pseudo_dna(10523, 80),
                    revcomp(&h_rev)
                ),
            ),
        ],
    );

    let out_masked = tmp.path().join("masked");
    let out_unmasked = tmp.path().join("unmasked");
    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "--spaced-seed",
            "w14s20",
            "--forward-only",
            "-l",
            "80",
            "-s",
            "80",
            "-f",
            "all",
            "-o",
            out_masked.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "spaced + soft-mask + forward-only + multi-FASTA + all",
    );
    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "--spaced-seed",
            "w14s20",
            "--forward-only",
            "--no-soft-mask",
            "-l",
            "80",
            "-s",
            "80",
            "-f",
            "frags",
            "-o",
            out_unmasked.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "spaced + --no-soft-mask + forward-only",
    );

    let f_m = parse_frags_file(&out_masked.join("q-r.frags")).fragments;
    let f_u = parse_frags_file(&out_unmasked.join("q-r.frags")).fragments;
    assert_eq!(
        f_m.len(),
        1,
        "soft-mask + forward-only keeps only uppercase forward h_up"
    );
    assert_eq!(f_m[0].seq_x, 0);
    assert_eq!(f_m[0].seq_y, 1);
    assert_eq!(f_m[0].strand, 'f');

    assert_eq!(
        f_u.len(),
        2,
        "--no-soft-mask + forward-only recovers both h_up and h_lo"
    );
    assert!(f_u.iter().all(|f| f.strand == 'f'));
    assert!(out_masked.join("q-r.png").exists());
    assert_eq!(
        fs::read_to_string(out_masked.join("q-r.paf"))
            .unwrap()
            .lines()
            .count(),
        1
    );
}

#[test]
fn tier3_percent_length_custom_scoring_n_penalty() {
    let tmp = tempfile::tempdir().unwrap();
    // Shorter sequence length = 500 bp; --percent-length -l 10 => min_len = 50 (strict > 50).
    // HSP1: 45 bp exact (< 50 -> rejected).
    // HSP2: 28 bp + 4 Ns + 28 bp = 60 bp (> 50 -> accepted!).
    // Custom scoring: match=5, mismatch=3, n_penalty=2 -> score = 56*5 - 4*2 = 272,
    // max_possible = 60*5 = 300, similarity = 272/300*100 = 90.6667%.
    let h1 = pseudo_dna(10601, 45);
    let h2 = format!("{}NNNN{}", pseudo_dna(10602, 28), pseudo_dna(10603, 28));
    let f0 = pseudo_dna(10610, 120);
    let f1 = pseudo_dna(10611, 120);
    let f2 = pseudo_dna(10612, 155);
    let q_seq = format!("{f0}{h1}{f1}{h2}{f2}");
    assert_eq!(q_seq.len(), 500);
    let r_seq = format!(
        "{}{}{h1}{}{}{}{h2}{}{}",
        pseudo_dna(10620, 160),
        all_mismatch(&f0[f0.len() - 40..]),
        all_mismatch(&f1[..40]),
        pseudo_dna(10621, 120),
        all_mismatch(&f1[f1.len() - 40..]),
        all_mismatch(&f2[..40]),
        pseudo_dna(10622, 160),
    );

    let q = tmp.path().join("q.fa");
    let r = tmp.path().join("r.fa");
    write_fasta(&q, &[("q", &q_seq)]);
    write_fasta(&r, &[("r", &r_seq)]);

    let out = tmp.path().join("out");
    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "12",
            "--percent-length",
            "-l",
            "10",
            "-s",
            "85",
            "--match-score",
            "5",
            "--mismatch-penalty",
            "3",
            "--n-penalty",
            "2",
            "-f",
            "frags",
            "-o",
            out.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "percent-length + custom scoring + N penalty",
    );

    let ff = parse_frags_file(&out.join("q-r.frags"));
    assert_eq!(ff.fragments.len(), 1);
    let f = &ff.fragments[0];
    assert_eq!(f.length, 60);
    assert_eq!(f.ident, 56);
    assert_eq!(f.score, 56 * 5 - 4 * 2);
    let expected_sim = ((272.0_f64) * 100.0 / 300.0_f64) as f32;
    assert!((f.similarity - expected_sim).abs() < 1e-3);
}

#[test]
fn tier3_all_vs_all_spaced_seed_gzip_mixed_inputs() {
    let tmp = tempfile::tempdir().unwrap();
    let gdir = tmp.path().join("genomes");
    fs::create_dir_all(&gdir).unwrap();

    let hsp = pseudo_dna(10701, 180);
    let mut mut_pos = Vec::new();
    for blk in 0..10 {
        for &dc in &[3usize, 9, 14] {
            mut_pos.push(blk * 18 + dc);
        }
    }
    let hsp_div = mutate_positions(&hsp, &mut_pos);

    let p_a = gdir.join("alpha.fa");
    let p_b = gdir.join("beta.fasta");
    let p_c = gdir.join("gamma.fa.gz");
    write_fasta(
        &p_a,
        &[(
            "alpha",
            &format!(
                "{}{}{}",
                pseudo_dna(10710, 200),
                hsp,
                pseudo_dna(10711, 200)
            ),
        )],
    );
    write_fasta(
        &p_b,
        &[(
            "beta",
            &format!(
                "{}{}{}",
                pseudo_dna(10720, 300),
                hsp_div,
                pseudo_dna(10721, 300)
            ),
        )],
    );
    write_fasta_gz(
        &p_c,
        &[(
            "gamma",
            &format!(
                "{}{}{}",
                pseudo_dna(10730, 400),
                revcomp(&hsp),
                pseudo_dna(10731, 400)
            ),
        )],
    );

    let out_ava = tmp.path().join("ava");
    assert_success(
        &run_gecko(&[
            "all-vs-all",
            gdir.to_str().unwrap(),
            "--spaced-seed",
            "w12s18",
            "--freq-limit",
            "50",
            "-l",
            "80",
            "-s",
            "40",
            "-f",
            "frags",
            "-o",
            out_ava.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "all-vs-all spaced seed mixed gzip inputs",
    );

    for (qa, rb, pair) in [
        (&p_a, &p_b, "alpha-beta"),
        (&p_a, &p_c, "alpha-gamma"),
        (&p_b, &p_c, "beta-gamma"),
    ] {
        let out_s = tmp.path().join(format!("s_{pair}"));
        assert_success(
            &run_gecko(&[
                "compare",
                qa.to_str().unwrap(),
                rb.to_str().unwrap(),
                "--spaced-seed",
                "w12s18",
                "--freq-limit",
                "50",
                "-l",
                "80",
                "-s",
                "40",
                "-f",
                "frags",
                "-o",
                out_s.to_str().unwrap(),
                "--backend",
                "cpu",
            ]),
            &format!("compare {pair}"),
        );
        assert_frags_bytes_identical(
            &out_ava.join(format!("{pair}.frags")),
            &out_s.join(format!("{pair}.frags")),
            &format!("all-vs-all mixed gzip parity for {pair}"),
        );
        assert_eq!(
            parse_frags_file(&out_ava.join(format!("{pair}.frags")))
                .fragments
                .len(),
            1
        );
    }
}

// ============================================================================
// TIER 4 — Real-World Application Scenarios (6 tests: 3 baseline + 3 M3)
// ============================================================================

#[test]
fn tier4_synthetic_multichrom_synteny_snps_inversions_repeats_softmask() {
    let tmp = tempfile::tempdir().unwrap();
    // Realistic 3-chr x 3-chr genome comparison (~25 Kb per genome):
    // - chr1 <-> chr1: 2 conserved syntenic blocks (4% SNPs), 1 inverted block (3% SNPs)
    // - chr2 <-> chr2: 1 conserved block flanked by soft-masked transposon repeats
    // - chr3 -> chr1: inter-chromosomal translocation block (5% SNPs)
    // - High-copy 150-bp repeat occurring 6x in ref (suppressed by --freq-limit 3)
    let syn1 = pseudo_dna(11001, 1200);
    let syn2 = pseudo_dna(11002, 900);
    let inv1 = pseudo_dna(11003, 1000);
    let syn_c2 = pseudo_dna(11004, 1100);
    let transloc = pseudo_dna(11005, 800);
    let te_soft = pseudo_dna(11006, 400).to_ascii_lowercase();
    let rep = pseudo_dna(11007, 150);

    let q_c1 = format!(
        "{}{}{}{}{}{}{}{}{}",
        pseudo_dna(11101, 2000),
        syn1,
        pseudo_dna(11102, 1500),
        inv1,
        pseudo_dna(11103, 1500),
        syn2,
        pseudo_dna(11104, 1000),
        rep,
        pseudo_dna(11105, 1000)
    );
    let q_c2 = format!(
        "{}{}{}{}{}",
        pseudo_dna(11106, 2000),
        te_soft,
        syn_c2,
        te_soft,
        pseudo_dna(11107, 2000)
    );
    let q_c3 = format!(
        "{}{}{}",
        pseudo_dna(11108, 2500),
        transloc,
        pseudo_dna(11109, 2500)
    );

    let r_c1 = format!(
        "{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}",
        pseudo_dna(11201, 1800),
        mutate_every(&syn1, 25, 10),
        pseudo_dna(11202, 1200),
        rep,
        pseudo_dna(11203, 1200),
        revcomp(&mutate_every(&inv1, 33, 15)),
        pseudo_dna(11204, 1000),
        rep,
        pseudo_dna(11205, 1000),
        mutate_every(&syn2, 25, 12),
        pseudo_dna(11206, 1000),
        mutate_every(&transloc, 20, 8),
        pseudo_dna(11207, 1000),
        rep,
        pseudo_dna(11208, 1000),
    );
    let r_c2 = format!(
        "{}{}{}{}{}{}{}{}{}",
        pseudo_dna(11209, 1500),
        rep,
        pseudo_dna(11210, 1000),
        te_soft,
        mutate_every(&syn_c2, 30, 10),
        te_soft,
        pseudo_dna(11211, 1000),
        rep,
        pseudo_dna(11212, 1500)
    );
    let r_c3 = format!(
        "{}{}{}",
        pseudo_dna(11213, 2500),
        rep,
        pseudo_dna(11214, 2500)
    );

    let q = tmp.path().join("spA.fa");
    let r = tmp.path().join("spB.fa");
    write_fasta(&q, &[("chr1", &q_c1), ("chr2", &q_c2), ("chr3", &q_c3)]);
    write_fasta(&r, &[("chr1", &r_c1), ("chr2", &r_c2), ("chr3", &r_c3)]);

    let out = tmp.path().join("out");
    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "16",
            "-l",
            "400",
            "-s",
            "70",
            "--freq-limit",
            "3",
            "-f",
            "frags,paf,jsonl",
            "-o",
            out.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "multi-chromosome synthetic synteny",
    );

    let ff = parse_frags_file(&out.join("spA-spB.frags"));
    assert_eq!(
        ff.fragments.len(),
        5,
        "expected 5 major syntenic/inverted/translocated blocks"
    );
    assert!(
        ff.fragments
            .iter()
            .any(|f| f.seq_x == 0 && f.seq_y == 0 && f.strand == 'r' && f.length >= 950)
    );
    assert!(
        ff.fragments
            .iter()
            .any(|f| f.seq_x == 1 && f.seq_y == 1 && f.strand == 'f' && f.length >= 1800)
    );
    assert!(
        ff.fragments
            .iter()
            .any(|f| f.seq_x == 2 && f.seq_y == 0 && f.strand == 'f' && f.length >= 750)
    );

    let paf_lines = fs::read_to_string(out.join("spA-spB.paf")).unwrap();
    assert_eq!(paf_lines.lines().count(), 5);
    assert!(paf_lines.contains("chr3\t") && paf_lines.contains("\t+\tchr1\t"));
}

#[test]
fn tier4_divergent_homolog_discovery_spaced_seeds_and_two_hit() {
    let tmp = tempfile::tempdir().unwrap();
    // Divergent homolog block (360 bp) with mutations at w12s18 don't-care offsets (every 18 bp)
    // so contiguous k=12 misses it, embedded in a 15 Kb x 15 Kb noisy genome.
    let hsp = pseudo_dna(12001, 360);
    let mut mut_pos = Vec::new();
    for blk in 0..20 {
        for &dc in &[3usize, 9, 14] {
            mut_pos.push(blk * 18 + dc);
        }
    }
    let hsp_div = mutate_positions(&hsp, &mut_pos);

    let q_seq = format!(
        "{}{}{}",
        pseudo_dna(12010, 7500),
        hsp,
        pseudo_dna(12011, 7500)
    );
    let r_seq = format!(
        "{}{}{}",
        pseudo_dna(12020, 7500),
        hsp_div,
        pseudo_dna(12021, 7500)
    );
    let q = tmp.path().join("q.fa");
    let r = tmp.path().join("r.fa");
    write_fasta(&q, &[("q", &q_seq)]);
    write_fasta(&r, &[("r", &r_seq)]);

    let out_contig = tmp.path().join("contig");
    let out_spaced = tmp.path().join("spaced");
    let out_twohit = tmp.path().join("twohit");

    assert_success(
        &run_gecko(&[
            "compare",
            q.to_str().unwrap(),
            r.to_str().unwrap(),
            "-k",
            "12",
            "-l",
            "100",
            "-s",
            "40",
            "-f",
            "frags",
            "-o",
            out_contig.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "divergent homolog contiguous k=12",
    );
    let res_sp = run_gecko(&[
        "compare",
        q.to_str().unwrap(),
        r.to_str().unwrap(),
        "--spaced-seed",
        "w12s18",
        "-l",
        "100",
        "-s",
        "40",
        "--json-stats",
        "-f",
        "frags",
        "-o",
        out_spaced.to_str().unwrap(),
        "--backend",
        "cpu",
    ]);
    assert_success(&res_sp, "divergent homolog w12s18");

    let res_th = run_gecko(&[
        "compare",
        q.to_str().unwrap(),
        r.to_str().unwrap(),
        "--spaced-seed",
        "w12s18",
        "--two-hit=40",
        "-l",
        "100",
        "-s",
        "40",
        "--json-stats",
        "-f",
        "frags",
        "-o",
        out_twohit.to_str().unwrap(),
        "--backend",
        "cpu",
    ]);
    assert_success(&res_th, "divergent homolog w12s18 --two-hit=40");

    assert_eq!(
        parse_frags_file(&out_contig.join("q-r.frags"))
            .fragments
            .len(),
        0
    );
    assert_eq!(
        parse_frags_file(&out_spaced.join("q-r.frags"))
            .fragments
            .len(),
        1
    );
    let f_th = parse_frags_file(&out_twohit.join("q-r.frags")).fragments;
    assert_eq!(
        f_th.len(),
        1,
        "two-hit must retain the true divergent homolog"
    );
    assert!(f_th[0].length >= 340);

    let j_sp: serde_json::Value = serde_json::from_slice(&res_sp.stdout).unwrap();
    let j_th: serde_json::Value = serde_json::from_slice(&res_th.stdout).unwrap();
    assert!(
        j_th["hits_used"].as_u64().unwrap() < j_sp["hits_used"].as_u64().unwrap(),
        "two-hit must reduce noise hits_used vs single-hit (spaced={}, two_hit={})",
        j_sp["hits_used"],
        j_th["hits_used"]
    );
}

#[test]
fn tier4_multi_genome_all_vs_all_matrix_and_symmetry() {
    let tmp = tempfile::tempdir().unwrap();
    let gdir = tmp.path().join("species");
    fs::create_dir_all(&gdir).unwrap();

    // 4 unequal-length genomes sharing a core conserved gene (300 bp) at varying divergence
    let core = pseudo_dna(13001, 300);
    let s1 = format!(
        "{}{}{}",
        pseudo_dna(13010, 600),
        core,
        pseudo_dna(13011, 600)
    );
    let s2 = format!(
        "{}{}{}",
        pseudo_dna(13020, 1100),
        mutate_every(&core, 50, 15),
        pseudo_dna(13021, 1100)
    );
    let s3 = format!(
        "{}{}{}",
        pseudo_dna(13030, 1800),
        mutate_every(&core, 30, 10),
        pseudo_dna(13031, 1900)
    );
    let s4 = format!(
        "{}{}{}",
        pseudo_dna(13040, 2800),
        revcomp(&core),
        pseudo_dna(13041, 2900)
    );

    let paths = [
        ("sp1.fa", &s1),
        ("sp2.fa", &s2),
        ("sp3.fa", &s3),
        ("sp4.fa", &s4),
    ];
    for (name, seq) in paths {
        write_fasta(&gdir.join(name), &[(name, seq)]);
    }

    let out_ava = tmp.path().join("ava");
    assert_success(
        &run_gecko(&[
            "all-vs-all",
            gdir.to_str().unwrap(),
            "-k",
            "12",
            "-l",
            "150",
            "-s",
            "70",
            "-f",
            "frags",
            "-o",
            out_ava.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "4-genome all-vs-all",
    );

    // All C(4,2) = 6 pairs must exist and contain the shared core HSP
    for pair in [
        "sp1-sp2", "sp1-sp3", "sp1-sp4", "sp2-sp3", "sp2-sp4", "sp3-sp4",
    ] {
        let p = out_ava.join(format!("{pair}.frags"));
        assert!(p.exists(), "missing pair {pair}");
        let ff = parse_frags_file(&p);
        assert_eq!(ff.fragments.len(), 1, "expected 1 core HSP for {pair}");
    }

    // Symmetry check: compare sp1 vs sp2 and sp2 vs sp1 (forward HSP swaps x <-> y)
    let out_ab = tmp.path().join("ab");
    let out_ba = tmp.path().join("ba");
    let p1 = gdir.join("sp1.fa");
    let p2 = gdir.join("sp2.fa");
    assert_success(
        &run_gecko(&[
            "compare",
            p1.to_str().unwrap(),
            p2.to_str().unwrap(),
            "-k",
            "12",
            "-l",
            "150",
            "-s",
            "70",
            "-f",
            "frags",
            "-o",
            out_ab.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "compare sp1 sp2",
    );
    assert_success(
        &run_gecko(&[
            "compare",
            p2.to_str().unwrap(),
            p1.to_str().unwrap(),
            "-k",
            "12",
            "-l",
            "150",
            "-s",
            "70",
            "-f",
            "frags",
            "-o",
            out_ba.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "compare sp2 sp1",
    );

    assert_frags_bytes_identical(
        &out_ava.join("sp1-sp2.frags"),
        &out_ab.join("sp1-sp2.frags"),
        "all-vs-all sp1-sp2 vs standalone compare",
    );
    let fab = &parse_frags_file(&out_ab.join("sp1-sp2.frags")).fragments[0];
    let fba = &parse_frags_file(&out_ba.join("sp2-sp1.frags")).fragments[0];
    assert_eq!(fab.strand, 'f');
    assert_eq!(fba.strand, 'f');
    assert_eq!(fab.x_start, fba.y_start);
    assert_eq!(fab.x_end, fba.y_end);
    assert_eq!(fab.y_start, fba.x_start);
    assert_eq!(fab.y_end, fba.x_end);
    assert_eq!(fab.length, fba.length);
    assert_eq!(fab.ident, fba.ident);
    assert_eq!(fab.score, fba.score);
}

#[test]
fn tier4_real_chromosome_chimp_gorilla_baseline_synteny_and_softmask() {
    let chimp_src = repo_root().join("chimp_chr21.fa");
    let gorilla_src = repo_root().join("gorilla_chr21.fa");
    if !chimp_src.exists() || !gorilla_src.exists() {
        eprintln!(
            "SKIPPED tier4_real_chromosome_chimp_gorilla_baseline_synteny_and_softmask: real FASTAs absent"
        );
        return;
    }

    let tmp = tempfile::tempdir().unwrap();
    let chimp_slice = tmp.path().join("chimp_25m.fa");
    let gorilla_slice = tmp.path().join("gorilla_29m.fa");
    let n_c = slice_real_fasta_to_temp(
        &chimp_src,
        &chimp_slice,
        "chimp_chr21_25m",
        25_000_000,
        25_200_000,
    );
    let n_g = slice_real_fasta_to_temp(
        &gorilla_src,
        &gorilla_slice,
        "gorilla_chr21_29m",
        29_700_000,
        29_950_000,
    );
    assert_eq!(n_c, 200_000);
    assert_eq!(n_g, 250_000);

    let out_def = tmp.path().join("def");
    let out_nosm = tmp.path().join("nosm");
    assert_success(
        &run_gecko(&[
            "compare",
            chimp_slice.to_str().unwrap(),
            gorilla_slice.to_str().unwrap(),
            "-k",
            "16",
            "-l",
            "100",
            "-s",
            "70",
            "--freq-limit",
            "50",
            "-f",
            "legacy,paf",
            "-o",
            out_def.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "real chimp x gorilla soft-masked",
    );
    assert_success(
        &run_gecko(&[
            "compare",
            chimp_slice.to_str().unwrap(),
            gorilla_slice.to_str().unwrap(),
            "-k",
            "16",
            "-l",
            "100",
            "-s",
            "70",
            "--freq-limit",
            "50",
            "--no-soft-mask",
            "-f",
            "frags",
            "-o",
            out_nosm.to_str().unwrap(),
            "--backend",
            "cpu",
        ]),
        "real chimp x gorilla --no-soft-mask",
    );

    let ff_def = parse_frags_file(&out_def.join("chimp_25m-gorilla_29m.frags"));
    let ff_nosm = parse_frags_file(&out_nosm.join("chimp_25m-gorilla_29m.frags"));
    assert!(
        ff_def.fragments.len() >= 10,
        "expected >= 10 real synteny HSPs on chimp[25M..25.2M] x gorilla[29.7M..29.95M], got {}",
        ff_def.fragments.len()
    );
    assert!(
        ff_nosm.fragments.len() >= ff_def.fragments.len(),
        "--no-soft-mask ({}) should discover >= soft-masked ({}) fragments on real primate DNA",
        ff_nosm.fragments.len(),
        ff_def.fragments.len()
    );
    let mat = parse_mat_file(&out_def.join("chimp_25m-gorilla_29m.frags.MAT"));
    assert_mat_covers_frags(&mat, &ff_def.fragments);
    let paf_count = fs::read_to_string(out_def.join("chimp_25m-gorilla_29m.paf"))
        .unwrap()
        .lines()
        .count();
    assert_eq!(paf_count, ff_def.fragments.len());
}

#[test]
fn tier4_real_chromosome_chimp_gorilla_minimizer_vs_exact_parity() {
    let chimp_src = repo_root().join("chimp_chr21.fa");
    let gorilla_src = repo_root().join("gorilla_chr21.fa");
    if !chimp_src.exists() || !gorilla_src.exists() {
        eprintln!(
            "SKIPPED tier4_real_chromosome_chimp_gorilla_minimizer_vs_exact_parity: real FASTAs absent"
        );
        return;
    }

    let tmp = tempfile::tempdir().unwrap();
    let chimp_slice = tmp.path().join("chimp_500k.fa");
    let gorilla_slice = tmp.path().join("gorilla_550k.fa");
    let n_c = slice_real_fasta_to_temp(
        &chimp_src,
        &chimp_slice,
        "chimp_chr21_500k",
        25_000_000,
        25_500_000,
    );
    let n_g = slice_real_fasta_to_temp(
        &gorilla_src,
        &gorilla_slice,
        "gorilla_chr21_550k",
        29_650_000,
        30_200_000,
    );
    assert_eq!(n_c, 500_000);
    assert_eq!(n_g, 550_000);

    for k_str in ["32", "20"] {
        let out_exact = tmp.path().join(format!("exact_{k_str}"));
        let out_min = tmp.path().join(format!("min_{k_str}"));
        assert_success(
            &run_gecko(&[
                "compare",
                chimp_slice.to_str().unwrap(),
                gorilla_slice.to_str().unwrap(),
                "-k",
                k_str,
                "-l",
                "100",
                "-s",
                "70",
                "--freq-limit",
                "50",
                "--seeding",
                "exact",
                "-f",
                "frags",
                "-o",
                out_exact.to_str().unwrap(),
                "--backend",
                "cpu",
            ]),
            &format!("real chimp x gorilla exact k={k_str}"),
        );
        assert_success(
            &run_gecko(&[
                "compare",
                chimp_slice.to_str().unwrap(),
                gorilla_slice.to_str().unwrap(),
                "-k",
                k_str,
                "-l",
                "100",
                "-s",
                "70",
                "--freq-limit",
                "50",
                "--seeding",
                "minimizer",
                "-f",
                "frags",
                "-o",
                out_min.to_str().unwrap(),
                "--backend",
                "cpu",
            ]),
            &format!("real chimp x gorilla minimizer k={k_str}"),
        );

        let f_ex = out_exact.join("chimp_500k-gorilla_550k.frags");
        let f_mn = out_min.join("chimp_500k-gorilla_550k.frags");
        assert_frags_bytes_identical(
            &f_ex,
            &f_mn,
            &format!("real chimp x gorilla minimizer vs exact k={k_str}"),
        );
        let parsed = parse_frags_file(&f_mn);
        assert!(
            parsed.fragments.len() >= 10,
            "expected >= 10 real synteny fragments at k={k_str}, got {}",
            parsed.fragments.len()
        );
    }
}

#[test]
fn tier4_real_chromosome_chimp_mouse_divergent_two_hit() {
    let chimp_src = repo_root().join("chimp_chr21.fa");
    let mouse_src = repo_root().join("mouse_chr16.fa");
    if !chimp_src.exists() || !mouse_src.exists() {
        eprintln!(
            "SKIPPED tier4_real_chromosome_chimp_mouse_divergent_two_hit: real FASTAs absent"
        );
        return;
    }

    let tmp = tempfile::tempdir().unwrap();
    let chimp_slice = tmp.path().join("chimp_20m.fa");
    let mouse_slice = tmp.path().join("mouse_5m.fa");
    // Mouse slice starts at 5_000_000 (past the 3.05 Mb leading N-telomere!)
    let n_c = slice_real_fasta_to_temp(
        &chimp_src,
        &chimp_slice,
        "chimp_chr21_20m",
        20_000_000,
        20_250_000,
    );
    let n_m = slice_real_fasta_to_temp(
        &mouse_src,
        &mouse_slice,
        "mouse_chr16_5m",
        5_000_000,
        5_250_000,
    );
    assert_eq!(n_c, 250_000);
    assert_eq!(n_m, 250_000);

    let out_base = tmp.path().join("base");
    let out_two = tmp.path().join("twohit");
    let res_base = run_gecko(&[
        "compare",
        chimp_slice.to_str().unwrap(),
        mouse_slice.to_str().unwrap(),
        "--spaced-seed",
        "w14s20",
        "-l",
        "40",
        "-s",
        "40",
        "--freq-limit",
        "50",
        "--json-stats",
        "-f",
        "frags",
        "-o",
        out_base.to_str().unwrap(),
        "--backend",
        "cpu",
    ]);
    assert_success(&res_base, "real chimp x mouse w14s20 single-hit");

    let res_two = run_gecko(&[
        "compare",
        chimp_slice.to_str().unwrap(),
        mouse_slice.to_str().unwrap(),
        "--spaced-seed",
        "w14s20",
        "--two-hit=40",
        "-l",
        "40",
        "-s",
        "40",
        "--freq-limit",
        "50",
        "--json-stats",
        "-f",
        "frags",
        "-o",
        out_two.to_str().unwrap(),
        "--backend",
        "cpu",
    ]);
    assert_success(&res_two, "real chimp x mouse w14s20 --two-hit=40");

    let j_base: serde_json::Value = serde_json::from_slice(&res_base.stdout).unwrap();
    let j_two: serde_json::Value = serde_json::from_slice(&res_two.stdout).unwrap();
    let hits_base = j_base["hits_used"].as_u64().unwrap();
    let hits_two = j_two["hits_used"].as_u64().unwrap();
    assert!(
        hits_base > 0,
        "expected non-zero w14s20 hits on real chimp[20M..20.25M] x mouse[5M..5.25M]"
    );
    assert!(
        hits_two < hits_base,
        "--two-hit=40 must reduce hits_used on real divergent mammalian DNA (base={hits_base}, two_hit={hits_two})"
    );
}
