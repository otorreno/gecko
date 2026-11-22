//! Standalone `metagecko` binary: Rust + GPU (CUDA / WGPU / AVX2 CPU) metagenomic
//! read-to-genome profiling and differential abundance engine built on GECKO v2.

#[allow(dead_code)]
mod backend;
mod meta_cmd;

use anyhow::{Context, Result};
use clap::Parser;

#[derive(Parser)]
#[command(
    name = "metagecko",
    version,
    about = "MetaGECKO v2 (Rust + GPU): Fine-grained metagenomic profiling, low-abundance detection, and comparative abundance analysis",
    long_about = None
)]
struct MetaGeckoCli {
    /// Number of worker threads (defaults to available logical CPUs).
    #[arg(short = 't', long)]
    threads: Option<usize>,

    #[command(flatten)]
    args: meta_cmd::MetaArgs,
}

fn main() -> Result<()> {
    let cli = MetaGeckoCli::parse();
    if let Some(threads) = cli.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .context("failed to configure rayon thread pool")?;
    }
    meta_cmd::run_meta(cli.args)
}
