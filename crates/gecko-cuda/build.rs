//! Optional native CUDA build script.
//!
//! When built WITHOUT `--features cuda-native`, this script exits immediately so the
//! workspace compiles cleanly on macOS, Windows, and CPU-only Linux CI runners.
//! When built WITH `--features cuda-native`, it compiles `cuda/*.cu` with `nvcc` into
//! `libgecko_cuda_kernels.a` and links `cudart`.

use std::env;
use std::path::PathBuf;
use std::process::Command;

fn main() {
    println!("cargo:rustc-check-cfg=cfg(has_cuda_kernels)");
    println!("cargo:rerun-if-env-changed=GECKO_CUDA_ARCH");
    println!("cargo:rerun-if-env-changed=CUDA_HOME");
    println!("cargo:rerun-if-env-changed=CUDA_PATH");

    let cuda_sources = [
        "cuda/cuda_common.cuh",
        "cuda/cuda_index.cu",
        "cuda/cuda_query.cu",
        "cuda/cuda_extend.cu",
        "cuda/cuda_bridge.cu",
    ];
    for src in &cuda_sources {
        println!("cargo:rerun-if-changed={src}");
    }

    if env::var_os("CARGO_FEATURE_CUDA_NATIVE").is_none() {
        return;
    }

    let out_dir = PathBuf::from(env::var("OUT_DIR").expect("OUT_DIR not set"));
    let nvcc = env::var("NVCC").unwrap_or_else(|_| "nvcc".to_string());
    let arch_str = env::var("GECKO_CUDA_ARCH").unwrap_or_else(|_| {
        "-gencode arch=compute_75,code=sm_75 -gencode arch=compute_80,code=sm_80 -gencode arch=compute_89,code=sm_89".to_string()
    });
    let arch_flags: Vec<&str> = arch_str.split_whitespace().collect();

    let cu_files = [
        ("cuda/cuda_index.cu", "cuda_index.o"),
        ("cuda/cuda_query.cu", "cuda_query.o"),
        ("cuda/cuda_extend.cu", "cuda_extend.o"),
        ("cuda/cuda_bridge.cu", "cuda_bridge.o"),
    ];

    let mut obj_paths = Vec::with_capacity(cu_files.len());
    for (src, obj_name) in &cu_files {
        let obj_path = out_dir.join(obj_name);
        let mut cmd = Command::new(&nvcc);
        cmd.args(["-O3", "-std=c++17", "-Icuda", "-Xcompiler", "-fPIC", "-c"]);
        cmd.args(&arch_flags);
        cmd.arg(src);
        cmd.arg("-o");
        cmd.arg(&obj_path);

        let status = cmd.status().unwrap_or_else(|e| {
            panic!("Failed to execute `{nvcc}` (is CUDA Toolkit installed?): {e}");
        });
        assert!(
            status.success(),
            "nvcc failed compiling {src} with status {status}"
        );
        obj_paths.push(obj_path);
    }

    let lib_path = out_dir.join("libgecko_cuda_kernels.a");
    let mut ar = Command::new("ar");
    ar.arg("rcs").arg(&lib_path);
    for obj in &obj_paths {
        ar.arg(obj);
    }
    let status = ar.status().expect("Failed to execute `ar`");
    assert!(status.success(), "`ar` failed with status {status}");

    println!("cargo:rustc-link-search=native={}", out_dir.display());
    println!("cargo:rustc-link-lib=static=gecko_cuda_kernels");

    // Locate CUDA runtime library directory
    let cuda_home = env::var("CUDA_HOME")
        .or_else(|_| env::var("CUDA_PATH"))
        .unwrap_or_else(|_| "/usr/local/cuda".to_string());
    let candidate_lib_dirs = [
        format!("{cuda_home}/lib64"),
        format!("{cuda_home}/lib/x64"),
        "/usr/lib/x86_64-linux-gnu".to_string(),
        "/usr/local/cuda/lib64".to_string(),
    ];
    for dir in &candidate_lib_dirs {
        if PathBuf::from(dir).exists() {
            println!("cargo:rustc-link-search=native={dir}");
        }
    }
    println!("cargo:rustc-link-lib=dylib=cudart");
    println!("cargo:rustc-link-lib=dylib=stdc++");
    println!("cargo:rustc-cfg=has_cuda_kernels");
}
