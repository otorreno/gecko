#!/usr/bin/env python3
"""Generate publication-grade 300-DPI PDF & PNG figures for the BMC Bioinformatics GECKO 2 manuscript."""
import csv
import json
import re
import zipfile
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

PAPER_DIR = Path(__file__).resolve().parent
FIG_DIR = PAPER_DIR / "figures"
FIG_DIR.mkdir(exist_ok=True)

ZIP_PATH = Path("/home/antigravity/gecko2_paper_results.zip")
CHR24_MD = Path("/home/antigravity/human_vs_chimp_backend_timings_(1).md")

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 9.5,
    "axes.titlesize": 10.5,
    "axes.titleweight": "bold",
    "axes.labelsize": 9.5,
    "xtick.labelsize": 8.2,
    "ytick.labelsize": 8.5,
    "legend.fontsize": 8.0,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

def load_24chr_md(path: Path):
    rows = []
    for line in path.read_text().splitlines():
        if not line.startswith("| chr"):
            continue
        parts = [p.strip() for p in line.strip().strip("|").split("|")]
        ch = parts[0].split()[0]
        rows.append({
            "chr": ch,
            "cpu_wall": float(parts[3].rstrip("s")),
            "c_fast_wall": float(parts[7].rstrip("s")),
            "cuda_wall": float(parts[13].rstrip("s")),
        })
    return rows

def main():
    zf = zipfile.ZipFile(ZIP_PATH)
    e1 = json.loads(zf.read("exp1_original_paper_pairwise.json"))
    e2_rows = list(csv.DictReader(zf.read("exp2_all_vs_all_45pairs.csv").decode().splitlines()))
    chr24 = load_24chr_md(CHR24_MD)

    # =========================================================================
    # FIGURE 2: 4-Panel Primary Performance Overview
    # =========================================================================
    fig, axes = plt.subplots(2, 2, figsize=(13.2, 9.0))

    # --- Panel A: Exp 1 (2015 Table 2 Replication, 4 kb -> 249 Mb) ---
    ax = axes[0, 0]
    ids = ["tylcv", "buchnera", "ecoli", "drosophila", "hs_pt_chr1"]
    short_names = ["TYLCV\n(4 kb)", "Buchnera\n(0.64 Mb)", "E. coli\n(5 Mb)", "Drosophila\n(20 Mb)", "HS×PT chr1\n(249 Mb)"]
    x = np.arange(len(ids))

    pub_2015 = [0.36, 1.60, 17.20, 48.72, 11848.15]
    ax.plot(x, pub_2015, marker="D", linestyle=":", color="#7f8c8d", linewidth=1.6, markersize=5.5,
            label="GECKO v1 (2015 published)")

    for b, lbl, col, mk, ls in [
        ("legacy_v1", "Legacy GECKO v1 (Colab, disk)", "#95a5a6", "s", "--"),
        ("c_fast",    "C gecko_fast (-march=native)",  "#e67e22", "^", "-"),
        ("cpu",       "GECKO 2 Rust CPU (2 vCPUs)",    "#2980b9", "o", "-"),
        ("cuda",      "GECKO 2 CUDA (Tesla T4)",       "#27ae60", "*", "-"),
    ]:
        ys = [e1[i][b]["wall_s"] if (b in e1[i] and e1[i][b].get("ok")) else np.nan for i in ids]
        ms = 9.0 if mk == "*" else 6.0
        ax.plot(x, ys, marker=mk, linestyle=ls, color=col, linewidth=2.0, markersize=ms, label=lbl)

    # Annotate key speedups on Panel A
    ax.annotate("340× vs v1\n(0 B disk)", xy=(3, e1["drosophila"]["cuda"]["wall_s"]),
                xytext=(2.18, 0.08), fontsize=7.8, fontweight="bold", color="#1e8449",
                arrowprops=dict(arrowstyle="->", color="#1e8449", lw=1.0))
    ax.annotate("2,282× vs 2015 v1\n12.0× vs Rust CPU", xy=(4, e1["hs_pt_chr1"]["cuda"]["wall_s"]),
                xytext=(3.05, 0.14), fontsize=7.6, fontweight="bold", color="#1e8449",
                arrowprops=dict(arrowstyle="->", color="#1e8449", lw=1.0))

    ax.set_yscale("log")
    ax.set_ylim(0.005, 40000)
    ax.set_xticks(x); ax.set_xticklabels(short_names)
    ax.set_ylabel("Wall Time (s, log scale)")
    ax.set_title("(A) Original 2015 Table 2 Replication (4 kb → 249 Mb)")
    ax.grid(True, which="both", alpha=0.22)
    ax.legend(loc="upper left", frameon=True, framealpha=0.92)

    # --- Panel B: Exp 2 (45-Pair Mammalian Chr1 All-vs-All vs 2015 Table 4) ---
    ax = axes[0, 1]
    cu_tot = sum(float(r["cuda_wall_s"]) for r in e2_rows)
    labels = ["GECKO 2\nCUDA (T4)", "GECKO v1\n(2015 cached)", "MUMmer\n(2015 HPC)", "Gepard\n(2015 HPC)", "LASTZ\n(2015)"]
    vals   = [cu_tot, 142954.0, 558360.0, 576889.0, 4752315.0]
    cols   = ["#27ae60", "#7f8c8d", "#95a5a6", "#aeb6bf", "#566573"]
    bars = ax.bar(labels, vals, color=cols, edgecolor="black", linewidth=0.6, width=0.62)
    ax.set_yscale("log")
    ax.set_ylim(30, 3.5e7)
    ax.set_ylabel("Total Wall Time for 45 Pairs (s, log scale)")
    ax.set_title("(B) 10-Species Mammalian Chr1 All-vs-All (45 Pairs)")
    ax.grid(True, axis="y", which="both", alpha=0.22)
    for bar, v in zip(bars, vals):
        if v < 3600:
            txt = f"{v:.1f} s\n({v/60:.2f} min)"
        elif v < 172800:
            txt = f"{v/3600:.1f} h\n({v/cu_tot:,.0f}×)"
        else:
            txt = f"{v/86400:.1f} d\n({v/cu_tot:,.0f}×)"
        ax.text(bar.get_x() + bar.get_width()/2, v * 1.32, txt, ha="center", va="bottom",
                fontsize=7.8, fontweight="bold", color="#1e8449" if v == cu_tot else "#2c3e50")

    # --- Panel C: Exp 3 (Full 24-Chromosome Human x Chimp Scaling) ---
    ax = axes[1, 0]
    xc = np.arange(len(chr24))
    ax.plot(xc, [r["cpu_wall"] for r in chr24], marker="o", color="#2980b9", lw=1.7, ms=4.5,
            label="GECKO 2 Rust CPU (669.7 s total)")
    ax.plot(xc, [r["c_fast_wall"] for r in chr24], marker="^", color="#e67e22", lw=1.5, ms=4.5, ls="--",
            label="C gecko_fast (619.3 s total)")
    ax.plot(xc, [r["cuda_wall"] for r in chr24], marker="*", color="#27ae60", lw=1.9, ms=7.0,
            label="GECKO 2 CUDA T4 (67.2 s total, 10.0×)")
    ax.set_yscale("log")
    ax.set_ylim(0.6, 120)
    ax.set_xticks(xc)
    ax.set_xticklabels([r["chr"].replace("chr", "") for r in chr24], rotation=45)
    ax.set_xlabel("Chromosome (Human GRCh38.p14 × Chimp Clint_PTRv2)")
    ax.set_ylabel("Wall Time (s, log scale)")
    ax.set_title("(C) Full 24-Chromosome Whole-Genome Scaling (~3.1 Gb × 2.9 Gb)")
    ax.grid(True, which="both", alpha=0.22)
    ax.legend(loc="lower left", frameon=True, framealpha=0.92)

    # --- Panel D: Exp 2 Subset (HS chr1 x 9 Mammals: Same-Machine CPU vs CUDA) ---
    ax = axes[1, 1]
    hs_rows = [r for r in e2_rows if r["q"] == "HS"]
    sp_labels = [r["r"] for r in hs_rows]
    xd = np.arange(len(sp_labels))
    w = 0.37
    cpu_v = [float(r["cpu_wall_s"]) for r in hs_rows]
    cu_v  = [float(r["cuda_wall_s"]) for r in hs_rows]
    ax.bar(xd - w/2, cpu_v, width=w, color="#2980b9", edgecolor="black", lw=0.5, label="GECKO 2 Rust CPU (2 vCPUs)")
    ax.bar(xd + w/2, cu_v,  width=w, color="#27ae60", edgecolor="black", lw=0.5, label="GECKO 2 CUDA (Tesla T4)")
    for i, (cv, gv) in enumerate(zip(cpu_v, cu_v)):
        ax.text(i, max(cv, gv) + 2.2, f"{cv/gv:.1f}×", ha="center", va="bottom",
                fontsize=7.8, fontweight="bold", color="#1e8449")
    ax.set_xticks(xd); ax.set_xticklabels(sp_labels)
    ax.set_ylim(0, max(cpu_v) * 1.18)
    ax.set_xlabel("Target Mammalian Chromosome 1 (vs. Homo sapiens chr1)")
    ax.set_ylabel("Wall Time (s)")
    ax.set_title("(D) Same-Machine CPU vs. CUDA Across 9 Mammals (2015 Fig. 2)")
    ax.grid(True, axis="y", alpha=0.22)
    ax.legend(loc="upper right", frameon=True, framealpha=0.92)

    plt.tight_layout(pad=1.4)
    plt.savefig(FIG_DIR / "gecko2_paper_figures.pdf", bbox_inches="tight")
    plt.savefig(FIG_DIR / "gecko2_paper_figures.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # =========================================================================
    # FIGURE 3: 2-Panel (10x10 All-vs-All CUDA Matrix + Disk/RAM Resource Comparison)
    # =========================================================================
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.2, 4.9), gridspec_kw={"width_ratios": [1.15, 1.0]})

    # --- Panel A: 10x10 All-vs-All CUDA Runtime Heatmap ---
    codes = ["HS", "PT", "MM", "PA", "GG", "MMU", "RN", "BT", "CF", "SS"]
    mat = np.full((10, 10), np.nan)
    lookup = {(r["q"], r["r"]): float(r["cuda_wall_s"]) for r in e2_rows}
    for i, a in enumerate(codes):
        for j, b in enumerate(codes):
            if j > i and (a, b) in lookup:
                mat[i, j] = lookup[(a, b)]
    im = ax1.imshow(mat, cmap="YlGnBu", vmin=1.0, vmax=5.2)
    for i in range(10):
        for j in range(10):
            if not np.isnan(mat[i, j]):
                val = mat[i, j]
                ax1.text(j, i, f"{val:.1f}", ha="center", va="center", fontsize=7.5,
                         color="white" if val > 3.6 else "#1c2833", fontweight="bold")
            elif i == j:
                ax1.text(j, i, "—", ha="center", va="center", fontsize=8, color="#bdc3c7")
    ax1.set_xticks(range(10)); ax1.set_xticklabels(codes)
    ax1.set_yticks(range(10)); ax1.set_yticklabels(codes)
    ax1.spines["top"].set_visible(True); ax1.spines["right"].set_visible(True)
    ax1.set_title("(A) 10-Species Chr1 All-vs-All CUDA Wall Time (s) — Total: 108.5 s")
    cbar = plt.colorbar(im, ax=ax1, fraction=0.046, pad=0.04)
    cbar.set_label("Wall Time (s)")

    # --- Panel B: Legacy v1 Temp Disk Explosion vs GECKO 2 Host RSS ---
    ds_lbls = ["TYLCV\n(4 kb)", "Buchnera\n(0.64 Mb)", "E. coli\n(5 Mb)", "Drosophila\n(20 Mb)", "HS×PT chr1\n(249 Mb)"]
    xb = np.arange(len(ds_lbls))
    w = 0.25
    v1_wr  = [e1[i]["legacy_v1"]["total_written_mib"] if "legacy_v1" in e1[i] and e1[i]["legacy_v1"].get("ok") else np.nan for i in ids]
    cpu_r  = [e1[i]["cpu"]["rss_mib"] for i in ids]
    cu_r   = [e1[i]["cuda"]["rss_mib"] for i in ids]

    ax2.bar(xb - w, v1_wr, width=w, color="#c0392b", alpha=0.85, edgecolor="black", lw=0.5, label="Legacy v1 Total Disk Written (MiB)")
    ax2.bar(xb,     cpu_r, width=w, color="#2980b9", alpha=0.85, edgecolor="black", lw=0.5, label="GECKO 2 Rust CPU Peak RSS (MiB)")
    ax2.bar(xb + w, cu_r,  width=w, color="#27ae60", alpha=0.85, edgecolor="black", lw=0.5, label="GECKO 2 CUDA Peak Host RSS (MiB)")
    ax2.set_yscale("log")
    ax2.set_ylim(0.4, 25000)
    ax2.set_xticks(xb); ax2.set_xticklabels(ds_lbls)
    ax2.set_ylabel("Memory / Temp Disk Footprint (MiB, log)")
    ax2.set_title("(B) Legacy v1 Disk I/O Amplification vs. GECKO 2 Host RAM")
    ax2.grid(True, axis="y", which="both", alpha=0.22)
    ax2.legend(loc="upper left", frameon=True, framealpha=0.92, fontsize=7.6)

    ax2.set_xlim(-0.55, 5.05)
    # Annotate key numbers on Panel B
    ax2.text(3 - w, v1_wr[3] * 1.25, "7.3 GiB\ndisk", ha="center", va="bottom", fontsize=7.2, fontweight="bold", color="#922b21")
    ax2.text(4, cpu_r[4] * 1.25, "3.0 GiB\n0 B disk", ha="center", va="bottom", fontsize=7.2, fontweight="bold", color="#1f618d")
    ax2.annotate("494 MiB\n(6.3× less)", xy=(4 + w, cu_r[4]), xytext=(4 + w + 0.30, cu_r[4] * 2.1),
                 ha="center", va="bottom", fontsize=7.2, fontweight="bold", color="#1e8449",
                 arrowprops=dict(arrowstyle="->", color="#1e8449", lw=0.9))

    plt.tight_layout(pad=1.4)
    plt.savefig(FIG_DIR / "gecko2_allvsall_resources.pdf", bbox_inches="tight")
    plt.savefig(FIG_DIR / "gecko2_allvsall_resources.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    print("Generated publication figures in paper/figures/:")
    for p in sorted(FIG_DIR.glob("gecko2_*")):
        print(f"  {p.name:36s} {p.stat().st_size:8,d} B")

if __name__ == "__main__":
    main()
