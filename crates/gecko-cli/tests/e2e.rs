use std::fs;
use std::io::Write;
use std::path::PathBuf;
use std::process::Command;

use flate2::Compression;
use flate2::write::GzEncoder;
use gecko_core::formats::csv::read_csv_path;
use gecko_core::formats::frags::read_frags_path;

fn gecko_bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_gecko"))
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

fn revcomp(s: &str) -> String {
    s.bytes()
        .rev()
        .map(gecko_core::seq::complement_base)
        .map(|b| b as char)
        .collect()
}

fn mutate(s: &str, every: usize) -> String {
    s.bytes()
        .enumerate()
        .map(|(i, b)| {
            if i > 0 && i % every == 0 {
                match b {
                    b'A' => b'C',
                    b'C' => b'G',
                    b'G' => b'T',
                    _ => b'A',
                }
            } else {
                b
            }
        })
        .map(|b| b as char)
        .collect()
}

#[test]
fn e2e_compare_all_formats_gzip_align_borders_filter_stats() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();

    // Build two multi-FASTA genomes with known forward (with 5% divergence) and inverted HSPs.
    let hsp_fwd = pseudo_dna(11, 250);
    let hsp_fwd_mut = mutate(&hsp_fwd, 20); // ~95% identity, ~90% similarity
    let hsp_rev = pseudo_dna(22, 220);
    let hsp_rev_rc = revcomp(&hsp_rev);

    let q_fa = format!(
        ">q_chr1 query chromosome 1\n{}{}{}\n>q_chr2 query chromosome 2\n{}{}{}\n",
        pseudo_dna(101, 150),
        hsp_fwd,
        pseudo_dna(102, 150),
        pseudo_dna(103, 120),
        hsp_rev,
        pseudo_dna(104, 120),
    );
    let r_fa = format!(
        ">r_chr1 ref chromosome 1\n{}{}{}\n>r_chr2 ref chromosome 2\n{}{}{}\n",
        pseudo_dna(201, 180),
        hsp_fwd_mut,
        pseudo_dna(202, 180),
        pseudo_dna(203, 140),
        hsp_rev_rc,
        pseudo_dna(204, 140),
    );

    let q_path = dir.join("query.fasta");
    let r_path = dir.join("ref.fa");
    fs::write(&q_path, &q_fa).unwrap();
    fs::write(&r_path, &r_fa).unwrap();

    // Also create a gzipped version of reference to test transparent .gz input.
    let r_gz_path = dir.join("ref.fa.gz");
    {
        let mut enc = GzEncoder::new(
            fs::File::create(&r_gz_path).unwrap(),
            Compression::default(),
        );
        enc.write_all(r_fa.as_bytes()).unwrap();
        enc.finish().unwrap();
    }

    let outdir = dir.join("out");
    let cache_dir = dir.join("cache");
    let status = Command::new(gecko_bin())
        .args([
            "compare",
            q_path.to_str().unwrap(),
            r_path.to_str().unwrap(),
            "-l",
            "100",
            "-s",
            "60",
            "-k",
            "12",
            "--format",
            "all",
            "-o",
            outdir.to_str().unwrap(),
            "--cache-dir",
            cache_dir.to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(status.success());

    let frags_path = outdir.join("query-ref.frags");
    let csv_path = outdir.join("query-ref.csv");
    let inf_path = outdir.join("query-ref.frags.INF");
    let mat_path = outdir.join("query-ref.frags.MAT");
    let paf_path = outdir.join("query-ref.paf");
    let blast_path = outdir.join("query-ref.blast.tsv");
    let bedpe_path = outdir.join("query-ref.bedpe");
    let jsonl_path = outdir.join("query-ref.jsonl");

    for p in [
        &frags_path,
        &csv_path,
        &inf_path,
        &mat_path,
        &paf_path,
        &blast_path,
        &bedpe_path,
        &jsonl_path,
    ] {
        assert!(p.exists(), "missing expected output {}", p.display());
    }

    let ff = read_frags_path(&frags_path).unwrap();
    assert_eq!(ff.fragments.len(), 2);
    assert!(
        ff.fragments
            .iter()
            .any(|f| f.strand == 'f' && f.seq_x == 0 && f.seq_y == 0)
    );
    assert!(
        ff.fragments
            .iter()
            .any(|f| f.strand == 'r' && f.seq_x == 1 && f.seq_y == 1)
    );

    let csv_ff = read_csv_path(&csv_path).unwrap();
    assert_eq!(csv_ff.fragments.len(), 2);

    // Check PAF lines contain resolved sequence names and proper strands
    let paf = fs::read_to_string(&paf_path).unwrap();
    let paf_lines: Vec<&str> = paf.lines().collect();
    assert_eq!(paf_lines.len(), 2);
    assert!(
        paf_lines
            .iter()
            .any(|l| l.starts_with("q_chr1\t") && l.contains("\t+\tr_chr1\t"))
    );
    assert!(
        paf_lines
            .iter()
            .any(|l| l.starts_with("q_chr2\t") && l.contains("\t-\tr_chr2\t"))
    );

    // Check BLAST tabular, BEDPE, JSONL line counts
    assert_eq!(fs::read_to_string(&blast_path).unwrap().lines().count(), 2);
    assert_eq!(fs::read_to_string(&bedpe_path).unwrap().lines().count(), 2);
    assert_eq!(fs::read_to_string(&jsonl_path).unwrap().lines().count(), 2);

    // Compare using gzipped reference and confirm identical fragments
    let out_gz = dir.join("out_gz");
    let status = Command::new(gecko_bin())
        .args([
            "compare",
            q_path.to_str().unwrap(),
            r_gz_path.to_str().unwrap(),
            "-l",
            "100",
            "-s",
            "60",
            "-k",
            "12",
            "--format",
            "frags",
            "-o",
            out_gz.to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(status.success());
    let ff_gz = read_frags_path(out_gz.join("query-ref.frags")).unwrap();
    assert_eq!(ff.fragments, ff_gz.fragments);

    // Test `gecko align` from both .frags and .csv
    let align_frags = dir.join("align_from_frags.txt");
    let align_csv = dir.join("align_from_csv.txt");
    assert!(
        Command::new(gecko_bin())
            .args([
                "align",
                frags_path.to_str().unwrap(),
                q_path.to_str().unwrap(),
                r_path.to_str().unwrap(),
                "-o",
                align_frags.to_str().unwrap(),
            ])
            .status()
            .unwrap()
            .success()
    );
    assert!(
        Command::new(gecko_bin())
            .args([
                "align",
                csv_path.to_str().unwrap(),
                q_path.to_str().unwrap(),
                r_path.to_str().unwrap(),
                "-o",
                align_csv.to_str().unwrap(),
            ])
            .status()
            .unwrap()
            .success()
    );
    assert_eq!(
        fs::read_to_string(&align_frags).unwrap(),
        fs::read_to_string(&align_csv).unwrap()
    );

    // Test `gecko borders` with --emit-query
    let borders_out = dir.join("borders.txt");
    let query_borders_fa = dir.join("query-25.fasta");
    assert!(
        Command::new(gecko_bin())
            .args([
                "borders",
                csv_path.to_str().unwrap(),
                q_path.to_str().unwrap(),
                r_path.to_str().unwrap(),
                "-b",
                "25",
                "-o",
                borders_out.to_str().unwrap(),
                "--emit-query",
                query_borders_fa.to_str().unwrap(),
            ])
            .status()
            .unwrap()
            .success()
    );
    assert!(
        fs::read_to_string(&query_borders_fa)
            .unwrap()
            .contains(">S0-")
    );

    // Test `gecko filter` (verifying similarity threshold works properly on 0..=100 scale!)
    let filtered_csv = dir.join("filtered.csv");
    assert!(
        Command::new(gecko_bin())
            .args([
                "filter",
                frags_path.to_str().unwrap(),
                "-l",
                "100",
                "-s",
                "95",
                "-o",
                filtered_csv.to_str().unwrap(),
            ])
            .status()
            .unwrap()
            .success()
    );
    // Only the 100%-identical inverted HSP passes >= 95% similarity (the mutated one is ~80%).
    let filtered_ff = read_csv_path(&filtered_csv).unwrap();
    assert_eq!(filtered_ff.fragments.len(), 1);
    assert_eq!(filtered_ff.fragments[0].strand, 'r');

    // Verify --format all also emitted a valid PNG dotplot
    let png_path = outdir.join("query-ref.png");
    let png_bytes = fs::read(&png_path).unwrap();
    assert_eq!(&png_bytes[0..8], &[137, 80, 78, 71, 13, 10, 26, 10]);

    // Test standalone `gecko plot` (dark mode + similarity heatmap + FASTA contig guides)
    let custom_png = dir.join("custom_plot.png");
    assert!(
        Command::new(gecko_bin())
            .args([
                "plot",
                csv_path.to_str().unwrap(),
                "-o",
                custom_png.to_str().unwrap(),
                "--width",
                "600",
                "--height",
                "600",
                "--color",
                "similarity",
                "--origin",
                "bottom-left",
                "--dark",
                "--query",
                q_path.to_str().unwrap(),
                "--reference",
                r_path.to_str().unwrap(),
            ])
            .status()
            .unwrap()
            .success()
    );
    let custom_bytes = fs::read(&custom_png).unwrap();
    assert_eq!(&custom_bytes[0..8], &[137, 80, 78, 71, 13, 10, 26, 10]);

    // Test `gecko stats` on .frags and .gidx
    let stats_out = Command::new(gecko_bin())
        .args(["stats", frags_path.to_str().unwrap(), "--json"])
        .output()
        .unwrap();
    assert!(stats_out.status.success());
    let stats_json: serde_json::Value = serde_json::from_slice(&stats_out.stdout).unwrap();
    assert_eq!(stats_json["fragments"], 2);
}

#[test]
fn e2e_all_vs_all_compares_all_pairs() {
    let tmp = tempfile::tempdir().unwrap();
    let gdir = tmp.path().join("genomes");
    let odir = tmp.path().join("out");
    fs::create_dir_all(&gdir).unwrap();

    let shared = pseudo_dna(777, 300);
    for (name, seed) in [("g1.fa", 1), ("g2.fa", 2), ("g3.fa", 3)] {
        let content = format!(
            ">{name}\n{}{}{}\n",
            pseudo_dna(seed, 100),
            shared,
            pseudo_dna(seed + 10, 100)
        );
        fs::write(gdir.join(name), content).unwrap();
    }

    let status = Command::new(gecko_bin())
        .args([
            "all-vs-all",
            gdir.to_str().unwrap(),
            "-l",
            "100",
            "-s",
            "80",
            "-k",
            "12",
            "--format",
            "frags,paf",
            "-o",
            odir.to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(status.success());

    for pair in ["g1-g2", "g1-g3", "g2-g3"] {
        assert!(odir.join(format!("{pair}.frags")).exists());
        assert!(odir.join(format!("{pair}.paf")).exists());
        let ff = read_frags_path(odir.join(format!("{pair}.frags"))).unwrap();
        assert_eq!(ff.fragments.len(), 1);
    }
}
