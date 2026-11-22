#!/bin/bash

set -e

BACKEND="${GECKO_BACKEND:-auto}"
EXTRA_WGPU_FLAGS=""

usage() {
   echo " ==== GECKO FAST HELP ..."
   echo ""
   echo "   Unified drop-in replacement for workflow.sh with automatic backend selection:"
   echo "   ./gecko/bin/workflow_fast.sh [--backend auto|cuda|wgpu|cpu] [--allow-software] query reference length similarity wordLength 1 [threads] [freqLimit]"
   echo ""
   echo "   Backends:"
   echo "     auto  - Prefer NVIDIA CUDA (gecko_cuda), then hardware WGPU (gecko_wgpu), else AVX2 CPU (gecko_fast)"
   echo "     cuda  - Native CUDA + CUB GPU backend (bin/gecko_cuda)"
   echo "     wgpu  - Cross-vendor Rust WGPU backend: Vulkan/Metal/DX12 (bin/gecko_wgpu)"
   echo "     cpu   - Multi-threaded AVX2 CPU backend (bin/gecko_fast)"
   echo ""
   echo "   Outputs GECKO-compatible files under ./results/<query>-<reference>.{frags,frags.INF,frags.MAT,csv}"
   echo ""
   echo " ===="
}

POSITIONAL=()
while [ $# -gt 0 ]; do
    case "$1" in
        -h|--help)
            usage
            exit 0
            ;;
        --backend)
            BACKEND="$2"
            shift 2
            ;;
        --backend=*)
            BACKEND="${1#*=}"
            shift
            ;;
        --allow-software)
            EXTRA_WGPU_FLAGS="--allow-software"
            shift
            ;;
        *)
            POSITIONAL+=("$1")
            shift
            ;;
    esac
done
set -- "${POSITIONAL[@]}"

if [ $# -lt 6 ] || [ $# -gt 8 ]; then
   usage
   exit 1
fi

BINDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

queryPath="$1"
refPath="$2"
length="$3"
similarity="$4"
WL="$5"
fixedL="$6"
threads="${7:-}"
freqLimit="${8:-}"

seqXName=$(basename "$queryPath")
seqXName="${seqXName%.*}"
seqYName=$(basename "$refPath")
seqYName="${seqYName%.*}"

mkdir -p results
outPrefix="results/${seqXName}-${seqYName}"

resolve_backend() {
    case "$BACKEND" in
        cuda)
            if [ ! -x "${BINDIR}/gecko_cuda" ]; then
                echo "Error: CUDA backend requested but ${BINDIR}/gecko_cuda is not built (run: make -C src gecko_cuda)." >&2
                exit 1
            fi
            echo "cuda"
            ;;
        wgpu)
            if [ ! -x "${BINDIR}/gecko_wgpu" ]; then
                echo "Error: WGPU backend requested but ${BINDIR}/gecko_wgpu is not built (run: make -C src gecko_wgpu)." >&2
                exit 1
            fi
            echo "wgpu"
            ;;
        cpu)
            if [ ! -x "${BINDIR}/gecko_fast" ]; then
                echo "Error: CPU backend requested but ${BINDIR}/gecko_fast is not built (run: make -C src gecko_fast)." >&2
                exit 1
            fi
            echo "cpu"
            ;;
        auto)
            if [ -x "${BINDIR}/gecko_cuda" ] && "${BINDIR}/gecko_cuda" --check >/dev/null 2>&1; then
                echo "cuda"
            elif [ -x "${BINDIR}/gecko_wgpu" ] && "${BINDIR}/gecko_wgpu" --check ${EXTRA_WGPU_FLAGS} >/dev/null 2>&1; then
                echo "wgpu"
            elif [ -x "${BINDIR}/gecko_fast" ]; then
                echo "cpu"
            else
                echo "Error: no GECKO fast backend binary found in ${BINDIR}." >&2
                exit 1
            fi
            ;;
        *)
            echo "Error: unknown backend '$BACKEND' (expected auto|cuda|wgpu|cpu)." >&2
            exit 1
            ;;
    esac
}

SELECTED=$(resolve_backend)
echo "[workflow_fast] Selected backend: ${SELECTED}"

case "$SELECTED" in
    cuda) BIN="${BINDIR}/gecko_cuda"; EXTRA="" ;;
    wgpu) BIN="${BINDIR}/gecko_wgpu"; EXTRA="${EXTRA_WGPU_FLAGS}" ;;
    cpu)  BIN="${BINDIR}/gecko_fast"; EXTRA="" ;;
esac

if [ -n "$freqLimit" ]; then
    "${BIN}" "$queryPath" "$refPath" "$outPrefix" "$length" "$similarity" "$WL" "$fixedL" "$threads" "$freqLimit" ${EXTRA}
elif [ -n "$threads" ]; then
    "${BIN}" "$queryPath" "$refPath" "$outPrefix" "$length" "$similarity" "$WL" "$fixedL" "$threads" ${EXTRA}
else
    "${BIN}" "$queryPath" "$refPath" "$outPrefix" "$length" "$similarity" "$WL" "$fixedL" ${EXTRA}
fi
