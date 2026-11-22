//! Unified compute backend probing and dispatch (`auto | cuda | wgpu | cpu`).

use std::path::Path;
use std::time::Instant;

use anyhow::{Context, Result, bail};
use clap::ValueEnum;
use gecko_core::pipeline::{CompareConfig, CompareResult, compare_genomes};
use gecko_core::seq::Genome;

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum BackendArg {
    /// Prefer NVIDIA CUDA, then hardware WGPU (Vulkan/Metal/DX12), else fast AVX2 CPU.
    Auto,
    /// Native NVIDIA CUDA + CUB GPU backend (requires build with `--features cuda`).
    Cuda,
    /// Cross-vendor WGPU compute shader backend (Vulkan / Metal / DX12).
    Wgpu,
    /// Multi-threaded fast CPU backend (canonical Hash-CSR/Direct-CSR + AVX2 SIMD).
    Cpu,
}

impl std::fmt::Display for BackendArg {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            BackendArg::Auto => f.write_str("auto"),
            BackendArg::Cuda => f.write_str("cuda"),
            BackendArg::Wgpu => f.write_str("wgpu"),
            BackendArg::Cpu => f.write_str("cpu"),
        }
    }
}

pub struct BackendRunOutput {
    pub query: Genome,
    pub reference: Genome,
    pub result: CompareResult,
    pub backend_label: String,
    pub load_secs: f64,
    pub compute_secs: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SelectedBackend {
    Cuda,
    Wgpu,
    Cpu,
}

fn probe_cuda_info() -> std::result::Result<String, String> {
    #[cfg(feature = "cuda")]
    {
        gecko_cuda::probe_cuda().map(|info| {
            format!(
                "cuda ({}, SM {}.{}, {:.1} GB VRAM)",
                info.name,
                info.sm_major,
                info.sm_minor,
                info.vram_bytes as f64 / (1024.0 * 1024.0 * 1024.0)
            )
        })
    }
    #[cfg(not(feature = "cuda"))]
    {
        Err("not compiled in (rebuild with `cargo build --release --features cuda`)".to_string())
    }
}

fn probe_wgpu_info(allow_software: bool) -> std::result::Result<String, String> {
    #[cfg(feature = "wgpu")]
    {
        gecko_wgpu::probe_wgpu(allow_software).map(|info| {
            format!(
                "wgpu ({}, {}, software={})",
                info.name, info.backend, info.is_software
            )
        })
    }
    #[cfg(not(feature = "wgpu"))]
    {
        let _ = allow_software;
        Err("not compiled in (rebuild with `--features wgpu`)".to_string())
    }
}

fn cpu_label() -> String {
    let threads = rayon::current_num_threads();
    #[cfg(target_arch = "x86_64")]
    let simd = if is_x86_feature_detected!("avx2") {
        "AVX2"
    } else {
        "scalar"
    };
    #[cfg(not(target_arch = "x86_64"))]
    let simd = "portable";
    format!("cpu ({simd}, {threads} threads)")
}

/// Print all compiled and detected GECKO compute backends.
pub fn print_backends(allow_software_gpu: bool) {
    println!("GECKO compute backends:");
    match probe_cuda_info() {
        Ok(desc) => println!("  cuda : AVAILABLE - {desc}"),
        Err(msg) => println!("  cuda : unavailable - {msg}"),
    }
    match probe_wgpu_info(allow_software_gpu) {
        Ok(desc) => println!("  wgpu : AVAILABLE - {desc}"),
        Err(msg) => println!("  wgpu : unavailable - {msg}"),
    }
    println!("  cpu  : AVAILABLE - {}", cpu_label());
}

fn resolve_backend(
    requested: BackendArg,
    allow_software_gpu: bool,
    config: &CompareConfig,
    cache_dir: Option<&Path>,
) -> Result<SelectedBackend> {
    let gpu_eligible = cache_dir.is_none() && (8..=32).contains(&config.word_length);
    match requested {
        BackendArg::Cpu => Ok(SelectedBackend::Cpu),
        BackendArg::Cuda => {
            if let Err(e) = probe_cuda_info() {
                bail!("CUDA backend requested but unavailable: {e}");
            }
            if !(8..=32).contains(&config.word_length) {
                bail!(
                    "CUDA backend requires word_length in 8..=32 (got {}); use --backend cpu",
                    config.word_length
                );
            }
            Ok(SelectedBackend::Cuda)
        }
        BackendArg::Wgpu => {
            if let Err(e) = probe_wgpu_info(allow_software_gpu) {
                bail!("WGPU backend requested but unavailable: {e}");
            }
            if !(8..=32).contains(&config.word_length) {
                bail!(
                    "WGPU backend requires word_length in 8..=32 (got {}); use --backend cpu",
                    config.word_length
                );
            }
            Ok(SelectedBackend::Wgpu)
        }
        BackendArg::Auto => {
            if gpu_eligible {
                if probe_cuda_info().is_ok() {
                    return Ok(SelectedBackend::Cuda);
                }
                if probe_wgpu_info(allow_software_gpu).is_ok() {
                    return Ok(SelectedBackend::Wgpu);
                }
            }
            Ok(SelectedBackend::Cpu)
        }
    }
}

/// Load FASTAs and execute pairwise comparison on the resolved backend.
pub fn run_compare_from_paths(
    query_path: &Path,
    ref_path: &Path,
    config: &CompareConfig,
    cache_dir: Option<&Path>,
    requested: BackendArg,
    allow_software_gpu: bool,
) -> Result<BackendRunOutput> {
    let selected = resolve_backend(requested, allow_software_gpu, config, cache_dir)?;

    match selected {
        SelectedBackend::Cuda => {
            #[cfg(feature = "cuda")]
            {
                let t0 = Instant::now();
                let (query, reference, result, info) =
                    gecko_cuda::compare_cuda_paths_overlapped(query_path, ref_path, config)?;
                let total_compute = t0.elapsed().as_secs_f64();
                let label = format!(
                    "cuda ({}, SM {}.{})",
                    info.name, info.sm_major, info.sm_minor
                );
                Ok(BackendRunOutput {
                    query,
                    reference,
                    result,
                    backend_label: label,
                    load_secs: 0.0,
                    compute_secs: total_compute,
                })
            }
            #[cfg(not(feature = "cuda"))]
            {
                bail!("CUDA backend not compiled in (rebuild with `--features cuda`)");
            }
        }
        SelectedBackend::Wgpu => {
            #[cfg(feature = "wgpu")]
            {
                let t0 = Instant::now();
                let (query, reference) = Genome::load_pair_parallel(query_path, ref_path)
                    .with_context(|| {
                        format!(
                            "failed to load FASTA inputs {} / {}",
                            query_path.display(),
                            ref_path.display()
                        )
                    })?;
                let load_secs = t0.elapsed().as_secs_f64();
                let t1 = Instant::now();
                let (result, info) =
                    gecko_wgpu::compare_wgpu(&query, &reference, config, allow_software_gpu)?;
                let compute_secs = t1.elapsed().as_secs_f64();
                let label = format!("wgpu ({}, {})", info.name, info.backend);
                Ok(BackendRunOutput {
                    query,
                    reference,
                    result,
                    backend_label: label,
                    load_secs,
                    compute_secs,
                })
            }
            #[cfg(not(feature = "wgpu"))]
            {
                let _ = (query_path, ref_path, config, allow_software_gpu);
                bail!("WGPU backend not compiled in (rebuild with `--features wgpu`)");
            }
        }
        SelectedBackend::Cpu => {
            let t0 = Instant::now();
            let (query, reference) = Genome::load_pair_parallel(query_path, ref_path)
                .with_context(|| {
                    format!(
                        "failed to load FASTA inputs {} / {}",
                        query_path.display(),
                        ref_path.display()
                    )
                })?;
            let load_secs = t0.elapsed().as_secs_f64();
            let t1 = Instant::now();
            let result = compare_genomes(&query, &reference, config, cache_dir)?;
            let compute_secs = t1.elapsed().as_secs_f64();
            Ok(BackendRunOutput {
                query,
                reference,
                result,
                backend_label: cpu_label(),
                load_secs,
                compute_secs,
            })
        }
    }
}
