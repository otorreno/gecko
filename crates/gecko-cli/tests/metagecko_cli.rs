use std::fs::{self, File};
use std::io::Write;
use std::process::Command;

use flate2::Compression;
use flate2::write::GzEncoder;
use gecko_core::seq::complement_base;
use tempfile::tempdir;

fn make_dna(seed: u64, len: usize) -> String {
    const ACGT: [u8; 4] = *b"ACGT";
    let mut x = seed;
    let mut s = String::with_capacity(len);
    for _ in 0..len {
        x = x
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        s.push(ACGT[((x >> 62) & 3) as usize] as char);
    }
    s
}

fn rc_str(s: &str) -> String {
    s.bytes()
        .rev()
        .map(complement_base)
        .map(|b| b as char)
        .collect()
}

fn write_fasta(path: &std::path::Path, records: &[(&str, &str)]) {
    let mut f = File::create(path).unwrap();
    for (id, seq) in records {
        writeln!(f, ">{id}").unwrap();
        writeln!(f, "{seq}").unwrap();
    }
}

fn format_fastq(reads: &[(String, String)]) -> Vec<u8> {
    let mut out = Vec::new();
    for (id, seq) in reads {
        writeln!(out, "@{id}").unwrap();
        writeln!(out, "{seq}").unwrap();
        writeln!(out, "+").unwrap();
        let qual = "I".repeat(seq.len());
        writeln!(out, "{qual}").unwrap();
    }
    out
}

#[test]
fn metagecko_multi_sample_fastq_gz_and_batch_streaming() {
    let tmp = tempdir().unwrap();
    let refs_dir = tmp.path().join("refs");
    fs::create_dir_all(&refs_dir).unwrap();

    // Create 3 bacterial reference genomes in a directory (A and B share a 200 bp conserved region)
    let shared = make_dna(555, 200);
    let mut ga = make_dna(101, 6000);
    let mut gb = make_dna(202, 6000);
    let gc = make_dna(303, 6000);
    ga.replace_range(1000..1200, &shared);
    gb.replace_range(2000..2200, &shared);

    write_fasta(&refs_dir.join("bac_A.fna"), &[("chrA", &ga)]);
    write_fasta(&refs_dir.join("bac_B.fna"), &[("chrB", &gb)]);
    write_fasta(&refs_dir.join("bac_C.fna"), &[("chrC", &gc)]);

    // Sample 1 (plain FASTQ): 10 reads from bac_A, 2 low-abundance reads from bac_C, 2 shared A/B reads, 3 noise
    let mut s1_reads = Vec::new();
    for i in 0..10 {
        let off = i * 400 + 50;
        let seq = if i % 2 == 0 {
            ga[off..off + 150].to_string()
        } else {
            rc_str(&ga[off..off + 150])
        };
        s1_reads.push((format!("s1_a_{i}"), seq));
    }
    s1_reads.push(("s1_c_0".into(), gc[200..350].to_string()));
    s1_reads.push(("s1_c_1".into(), rc_str(&gc[3000..3150])));
    s1_reads.push(("s1_shared_0".into(), shared[10..160].to_string()));
    s1_reads.push(("s1_shared_1".into(), rc_str(&shared[20..170])));
    for i in 0..3 {
        s1_reads.push((format!("s1_noise_{i}"), make_dna(9000 + i, 150)));
    }
    let s1_path = tmp.path().join("sample1.fastq");
    fs::write(&s1_path, format_fastq(&s1_reads)).unwrap();

    // Sample 2 (gzipped FASTQ): 12 reads from bac_B, 2 noise reads
    let mut s2_reads = Vec::new();
    for i in 0..12 {
        let off = i * 350 + 100;
        let seq = if i % 2 == 0 {
            gb[off..off + 150].to_string()
        } else {
            rc_str(&gb[off..off + 150])
        };
        s2_reads.push((format!("s2_b_{i}"), seq));
    }
    s2_reads.push(("s2_noise_0".into(), make_dna(7001, 150)));
    s2_reads.push(("s2_noise_1".into(), make_dna(7002, 150)));
    let s2_path = tmp.path().join("sample2.fastq.gz");
    {
        let f = File::create(&s2_path).unwrap();
        let mut gz = GzEncoder::new(f, Compression::fast());
        gz.write_all(&format_fastq(&s2_reads)).unwrap();
        gz.finish().unwrap();
    }

    let out_dir = tmp.path().join("out");
    let metagecko_bin = env!("CARGO_BIN_EXE_metagecko");

    // Force multi-batch streaming with tiny --batch-size 1K to exercise chunk boundaries!
    let status = Command::new(metagecko_bin)
        .args([
            s1_path.to_str().unwrap(),
            s2_path.to_str().unwrap(),
            "--refs",
            refs_dir.to_str().unwrap(),
            "-k",
            "16",
            "-l",
            "80",
            "-s",
            "80",
            "--min-read-cov",
            "70",
            "--batch-size",
            "1K",
            "--backend",
            "cpu",
            "-f",
            "all",
            "-o",
            out_dir.to_str().unwrap(),
        ])
        .status()
        .expect("failed to execute metagecko binary");
    assert!(status.success(), "metagecko exited with failure");

    // Verify all outputs exist
    assert!(out_dir.join("sample1.abundance.tsv").exists());
    assert!(out_dir.join("sample1.abundance.json").exists());
    assert!(out_dir.join("sample1.reads.tsv").exists());
    assert!(out_dir.join("sample1.csv").exists());
    assert!(out_dir.join("sample1.frags").exists());
    assert!(out_dir.join("sample2.abundance.tsv").exists());
    assert!(out_dir.join("abundance_matrix.rel.tsv").exists());
    assert!(out_dir.join("abundance_matrix.reads.tsv").exists());

    // Parse sample1 JSON and assert exact read accounting
    let s1_json: serde_json::Value =
        serde_json::from_str(&fs::read_to_string(out_dir.join("sample1.abundance.json")).unwrap())
            .unwrap();
    assert_eq!(s1_json["reads_analyzed"].as_u64(), Some(17));
    assert_eq!(s1_json["mapped_reads"].as_u64(), Some(14));
    assert_eq!(s1_json["unique_mapped_reads"].as_u64(), Some(12));
    assert_eq!(s1_json["multi_mapped_reads"].as_u64(), Some(2));
    assert!(s1_json["batches_processed"].as_u64().unwrap() >= 2);

    let profiles = s1_json["profiles"].as_array().unwrap();
    let find_prof = |name: &str| {
        profiles
            .iter()
            .find(|p| p["genome_name"].as_str() == Some(name))
            .unwrap()
    };
    let pa = find_prof("bac_A");
    let pb = find_prof("bac_B");
    let pc = find_prof("bac_C");

    // Note: bac_A read 2 (off=850..1000) & read 3 (off=1250..1400) do not overlap shared (1000..1200)
    assert_eq!(pa["unique_reads"].as_u64(), Some(10));
    assert_eq!(pa["multi_reads"].as_u64(), Some(2));
    assert!((pa["assigned_reads"].as_f64().unwrap() - 11.0).abs() < 1e-6);

    assert_eq!(pb["unique_reads"].as_u64(), Some(0));
    assert_eq!(pb["multi_reads"].as_u64(), Some(2));
    assert!((pb["assigned_reads"].as_f64().unwrap() - 1.0).abs() < 1e-6);

    assert_eq!(pc["unique_reads"].as_u64(), Some(2));
    assert_eq!(pc["covered_bases"].as_u64(), Some(300));

    // Parse sample2 JSON
    let s2_json: serde_json::Value =
        serde_json::from_str(&fs::read_to_string(out_dir.join("sample2.abundance.json")).unwrap())
            .unwrap();
    assert_eq!(s2_json["mapped_reads"].as_u64(), Some(12));
    let p2b = s2_json["profiles"]
        .as_array()
        .unwrap()
        .iter()
        .find(|p| p["genome_name"].as_str() == Some("bac_B"))
        .unwrap();
    assert_eq!(p2b["unique_reads"].as_u64(), Some(12));
}
