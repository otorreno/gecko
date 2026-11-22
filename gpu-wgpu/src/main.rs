//! gecko_wgpu CLI entry point — orchestrates seq -> index -> gpu -> output.

mod gpu;
mod index;
mod output;
mod seq;

use std::path::Path;
use std::process::ExitCode;
use std::time::Instant;

use gpu::WgpuContext;
use index::WgpuRefIndex;
use output::{finalize_and_dedup, write_outputs};
use seq::FastSeq;

fn print_usage(prog: &str) {
    eprintln!(
        "Usage:\n  {prog} --check [--allow-software]\n  \
         {prog} <query.fa> <ref.fa> <out_prefix> <length> <similarity> <WL> <fixedL> [threads] [freqLimit] [--allow-software]"
    );
}

fn main() -> ExitCode {
    let raw_args: Vec<String> = std::env::args().collect();
    let allow_software = raw_args.iter().any(|a| a == "--allow-software");
    let args: Vec<String> = raw_args
        .into_iter()
        .filter(|a| a != "--allow-software")
        .collect();

    if args.len() == 2 && args[1] == "--check" {
        match WgpuContext::init(allow_software) {
            Ok(ctx) => {
                println!(
                    "WGPU adapter: {} ({}, software={})",
                    ctx.adapter_name, ctx.adapter_backend, ctx.is_software
                );
                return ExitCode::SUCCESS;
            }
            Err(e) => {
                eprintln!("{e}");
                return ExitCode::from(2);
            }
        }
    }

    if args.len() < 8 || args.len() > 10 {
        print_usage(&args[0]);
        return ExitCode::FAILURE;
    }

    let seq_x_path = &args[1];
    let seq_y_path = &args[2];
    let out_prefix = &args[3];
    let l_min: u64 = args[4].parse().expect("invalid length");
    let sim_th: u64 = args[5].parse().expect("invalid similarity");
    let wl: usize = args[6].parse().expect("invalid WL");
    let fixed_l: i32 = args[7].parse().expect("invalid fixedL");
    let freq_limit: u32 = if args.len() >= 10 {
        args[9].parse().expect("invalid freqLimit")
    } else {
        1000
    };

    if !(8..=32).contains(&wl) {
        eprintln!("Error: WL must be between 8 and 32 (got {wl})");
        return ExitCode::FAILURE;
    }

    if args.len() >= 9 {
        if let Ok(t) = args[8].parse::<usize>() {
            if t > 0 {
                let _ = rayon::ThreadPoolBuilder::new().num_threads(t).build_global();
            }
        }
    }

    let ctx = match WgpuContext::init(allow_software) {
        Ok(c) => c,
        Err(e) => {
            eprintln!("Error: {e}");
            return ExitCode::from(2);
        }
    };

    let t0 = Instant::now();

    let (res_x, res_y) = std::thread::scope(|s| {
        let hx = s.spawn(|| FastSeq::load_fasta(Path::new(seq_x_path)));
        let ry = FastSeq::load_fasta(Path::new(seq_y_path));
        (hx.join().expect("load_fasta X thread panicked"), ry)
    });
    let s_x = match res_x {
        Ok(s) => s,
        Err(e) => {
            eprintln!("Error: {e}");
            return ExitCode::FAILURE;
        }
    };
    let s_y = match res_y {
        Ok(s) => s,
        Err(e) => {
            eprintln!("Error: {e}");
            return ExitCode::FAILURE;
        }
    };
    let s_y_rc = s_y.reverse_complement();
    let t_load = t0.elapsed().as_secs_f64();

    let min_len = if fixed_l != 0 {
        l_min
    } else {
        ((s_x.n_bases.min(s_y.n_bases) as f64) * (l_min as f64 / 100.0)) as u64
    };

    let idx = WgpuRefIndex::build(&s_y, wl, freq_limit);
    let t_index = t0.elapsed().as_secs_f64();

    let gpu_stats = match ctx.run_query_and_extend(
        &s_x,
        &s_y,
        &s_y_rc,
        &idx,
        wl as u32,
        min_len as u32,
        sim_th as u32,
    ) {
        Ok(st) => st,
        Err(e) => {
            eprintln!("GPU execution error: {e}");
            return ExitCode::FAILURE;
        }
    };
    let t_gpu = t0.elapsed().as_secs_f64();

    let (frags_f, frags_r, coverage) =
        finalize_and_dedup(&gpu_stats.raw_frags, &s_x, &s_y, &s_y_rc);

    if let Err(e) = write_outputs(
        out_prefix,
        seq_x_path,
        seq_y_path,
        &s_x,
        &s_y,
        l_min,
        sim_th,
        gpu_stats.leftmax_hits,
        gpu_stats.extended_hits,
        &frags_f,
        &frags_r,
        &coverage,
    ) {
        eprintln!("Output write error: {e}");
        return ExitCode::FAILURE;
    }
    let t_end = t0.elapsed().as_secs_f64();

    println!(
        "gecko_wgpu summary ({}, {}):\n\
         \x20 SeqX: {} bp ({} seqs) | SeqY: {} bp ({} seqs)\n\
         \x20 Index: WGPU-Hash-CSR (WL={}, FL={}, P={}, kept={}/{})\n\
         \x20 Left-max seeds: {} | Extended: {} | Fragments: {} (F:{} R:{})\n\
         \x20 Timing: load={:.3}s  index={:.3}s  wgpu_query+ext={:.3}s  write={:.3}s  TOTAL={:.3}s",
        gpu_stats.adapter_name,
        gpu_stats.adapter_backend,
        s_x.n_bases,
        s_x.n_seqs,
        s_y.n_bases,
        s_y.n_seqs,
        wl,
        freq_limit,
        idx.bucket_bits,
        idx.n_kept,
        idx.n_valid,
        gpu_stats.leftmax_hits,
        gpu_stats.extended_hits,
        frags_f.len() + frags_r.len(),
        frags_f.len(),
        frags_r.len(),
        t_load,
        t_index - t_load,
        t_gpu - t_index,
        t_end - t_gpu,
        t_end
    );

    ExitCode::SUCCESS
}
