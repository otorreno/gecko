#!/usr/bin/env python3
"""Populate paper/main.tex result macros & figures from Colab's gecko2_paper_results.zip."""
import csv
import json
import re
import sys
import zipfile
from pathlib import Path

PAPER_DIR = Path(__file__).resolve().parent
MAIN_TEX = PAPER_DIR / "main.tex"
FIG_DIR = PAPER_DIR / "figures"

PUB_2015_CHR1 = 11848.15
PUB_2015_T4 = {
    "cached": 142954.0,
    "nocache": 318591.0,
    "mummer": 558360.0,
    "gepard": 576889.0,
    "lastz": 4752315.0,
}

def fmt_s(v, digits=2):
    return f"{v:.{digits}f}"

def fmt_int(v):
    return f"{int(v):,}"

def fmt_sp(num, den, digits=1):
    if not num or not den or den <= 0:
        return "---"
    val = num / den
    if val >= 100:
        return f"{val:,.0f}$\\times$"
    return f"{val:.{digits}f}$\\times$"

def set_macro(tex: str, name: str, value: str) -> str:
    pattern = re.compile(r"(\\newcommand\{\\" + re.escape(name) + r"\}\{)(.*?)(\})")
    new_tex, n = pattern.subn(lambda m: m.group(1) + value + m.group(3), tex, count=1)
    if n == 0:
        print(f"WARNING: macro \\{name} not found in main.tex", file=sys.stderr)
    return new_tex

def main():
    if len(sys.argv) != 2:
        print("Usage: python3 paper/populate_results.py <path/to/gecko2_paper_results.zip>")
        sys.exit(1)

    zip_path = Path(sys.argv[1])
    if not zip_path.exists():
        sys.exit(f"File not found: {zip_path}")

    FIG_DIR.mkdir(exist_ok=True)
    data = {}
    with zipfile.ZipFile(zip_path) as zf:
        for fn in zf.namelist():
            if fn.endswith(".json"):
                data[fn] = json.loads(zf.read(fn).decode("utf-8"))
            elif fn.endswith(".csv"):
                rows = list(csv.DictReader(zf.read(fn).decode("utf-8").splitlines()))
                data[fn] = rows
            elif fn.endswith((".pdf", ".png")):
                (FIG_DIR / Path(fn).name).write_bytes(zf.read(fn))
                print(f"Extracted figure -> paper/figures/{Path(fn).name}")

    updates = {}

    # Hardware provenance
    if "hardware_info.json" in data:
        hw = data["hardware_info.json"]
        cpu_line = next((l.split(":", 1)[1].strip() for l in hw.get("lscpu", "").splitlines() if "Model name" in l), "Intel Xeon")
        updates["hwCpu"] = f"{cpu_line} ({hw.get('cpu_count', '?')} vCPUs)"
        mem_lines = [l for l in hw.get("mem", "").splitlines() if l.startswith("Mem:")]
        if mem_lines:
            updates["hwRam"] = mem_lines[0].split()[1] + " RAM"
        if hw.get("gpu"):
            updates["hwGpu"] = hw["gpu"].splitlines()[0].replace("_", r"\_")
        if hw.get("gcc"):
            updates["hwGcc"] = hw["gcc"].split()[0] + " " + hw["gcc"].split()[-1]
        if hw.get("nvcc"):
            m = re.search(r"release\s+([0-9.]+)", hw["nvcc"])
            if m:
                updates["hwNvcc"] = f"CUDA {m.group(1)}"

    # Experiment 1
    if "exp1_original_paper_pairwise.json" in data:
        e1 = data["exp1_original_paper_pairwise.json"]
        key_map = {
            "tylcv": "Tylcv", "buchnera": "Buch", "ecoli": "Ecoli",
            "drosophila": "Dros", "hs_pt_chr1": "ChrOne",
        }
        for bid, prefix in key_map.items():
            row = e1.get(bid, {})
            for b, suf in [("legacy_v1", "Vone"), ("c_fast", "Cfast"), ("cpu", "Cpu"), ("wgpu", "Wgpu"), ("cuda", "Cuda")]:
                r = row.get(b, {})
                if r.get("ok"):
                    d = 3 if r["wall_s"] < 1.0 else 2
                    updates[f"eOne{prefix}{suf}"] = fmt_s(r["wall_s"], d)
            v1 = row.get("legacy_v1", {})
            cu = row.get("cuda", {})
            cp = row.get("cpu", {})
            if v1.get("ok") and cu.get("ok"):
                updates[f"eOne{prefix}Sp"] = fmt_sp(v1["wall_s"], cu["wall_s"])
            if v1.get("ok"):
                updates[f"eOne{prefix}Disk"] = f"{v1.get('temp_disk_mib', 0.0):,.1f}"
            if cp.get("ok"):
                updates[f"eOne{prefix}Frags"] = fmt_int(cp["fragments"])
            if bid == "hs_pt_chr1" and v1.get("ok") and cp.get("ok") and cu.get("ok"):
                updates["eOneChrOneTotWr"] = f"{v1.get('total_written_mib', 0.0):,.0f}"
                updates["eOneChrOneVoneFrags"] = fmt_int(v1["fragments"])
                if "cpu_recovers_v1_pct" in v1:
                    updates["eOneChrOneVtwoCovVone"] = f"{v1['cpu_recovers_v1_pct']:.2f}\\%"
                    updates["hlVtwoRecoversVonePct"] = f"{v1['cpu_recovers_v1_pct']:.2f}\\%"
                updates["hlChrOneSpeedupCpuVsVone"] = fmt_sp(v1["wall_s"], cp["wall_s"])
                updates["hlChrOneSpeedupCudaVsVone"] = fmt_sp(v1["wall_s"], cu["wall_s"])
                updates["hlChrOneSpeedupCudaVsPub"] = fmt_sp(PUB_2015_CHR1, cu["wall_s"])
                updates["hlChrOneVoneDiskGiB"] = f"{v1.get('total_written_mib', 0.0)/1024.0:.1f}~GiB"

    # Experiment 2
    if "exp2_all_vs_all_45pairs.json" in data:
        e2 = data["exp2_all_vs_all_45pairs.json"]
        cu_ok = [v["cuda"] for v in e2.values() if v.get("cuda", {}).get("ok")]
        if cu_ok:
            tot_w = sum(r["wall_s"] for r in cu_ok)
            tot_c = sum(r["compute_s"] for r in cu_ok)
            max_r = max(r["rss_mib"] for r in cu_ok)
            tot_f = sum(r["fragments"] for r in cu_ok)
            updates["eTwoCudaWall"] = f"{tot_w:,.1f}"
            updates["eTwoCudaComp"] = f"{tot_c:,.1f}"
            updates["eTwoCudaRss"] = f"{max_r:,.1f}"
            updates["eTwoTotalFrags"] = fmt_int(tot_f)
            updates["hlAllVsAllCudaMin"] = f"{tot_w/60.0:.1f}~min"
            updates["hlAllVsAllCudaSec"] = f"{tot_w:,.1f}~s"
            updates["eTwoSpVoneCached"] = fmt_sp(PUB_2015_T4["cached"], tot_w)
            updates["hlAllVsAllSpeedupVsPubVone"] = fmt_sp(PUB_2015_T4["cached"], tot_w)
            updates["eTwoSpVoneNocache"] = fmt_sp(PUB_2015_T4["nocache"], tot_w)
            updates["eTwoSpMummer"] = fmt_sp(PUB_2015_T4["mummer"], tot_w)
            updates["hlAllVsAllSpeedupVsMummer"] = fmt_sp(PUB_2015_T4["mummer"], tot_w)
            updates["eTwoSpGepard"] = fmt_sp(PUB_2015_T4["gepard"], tot_w)
            updates["eTwoSpLastz"] = fmt_sp(PUB_2015_T4["lastz"], tot_w)
        for code, tag in [("PT","Pt"),("MM","Mm"),("PA","Pa"),("GG","Gg"),("MMU","Mmu"),("RN","Rn"),("BT","Bt"),("CF","Cf"),("SS","Ss")]:
            row = e2.get(f"HS-{code}", {})
            if row.get("cpu", {}).get("ok"):
                updates[f"eTwoHs{tag}Cpu"] = fmt_s(row["cpu"]["wall_s"])
            if row.get("cuda", {}).get("ok"):
                updates[f"eTwoHs{tag}Cuda"] = fmt_s(row["cuda"]["wall_s"])
                updates[f"eTwoHs{tag}Frags"] = fmt_int(row["cuda"]["fragments"])

    # Experiment 3
    if "exp3_human_vs_chimp_24chr.json" in data:
        e3 = data["exp3_human_vs_chimp_24chr.json"]
        totals = {}
        for b, tag in [("c_fast","Cfast"),("cpu","Cpu"),("wgpu","Wgpu"),("cuda","Cuda")]:
            recs = [v[b] for v in e3.values() if v.get(b, {}).get("ok")]
            if recs:
                tw = sum(r["wall_s"] for r in recs)
                tc = sum(r["compute_s"] for r in recs)
                mr = max(r["rss_mib"] for r in recs)
                tf = sum(r["fragments"] for r in recs)
                rc = min(r.get("recall_pct", 100.0) for r in recs)
                totals[b] = (tw, tc, mr, tf, rc)
                updates[f"eThree{tag}Wall"] = f"{tw:,.2f}"
                updates[f"eThree{tag}Comp"] = f"{tc:,.2f}"
                updates[f"eThree{tag}Rss"] = f"{mr:,.1f}"
                if b == "cpu":
                    updates["eThreeCpuFrags"] = fmt_int(tf)
                    updates["hlWholeGenomeCpuSec"] = f"{tw:,.1f}~s"
                else:
                    updates[f"eThree{tag}Recall"] = f"{rc:.3f}\\%"
                if b == "cuda":
                    updates["hlWholeGenomeCudaSec"] = f"{tw:,.1f}~s"
        if "cpu" in totals and "cuda" in totals:
            updates["eThreeCudaSpWall"] = fmt_sp(totals["cpu"][0], totals["cuda"][0], 2)
            updates["eThreeCudaSpComp"] = fmt_sp(totals["cpu"][1], totals["cuda"][1], 2)
        if "cpu" in totals and "wgpu" in totals:
            updates["eThreeWgpuSpWall"] = fmt_sp(totals["cpu"][0], totals["wgpu"][0], 2)
        if "cpu" in totals and "c_fast" in totals:
            pct = 100.0 * (totals["c_fast"][0] - totals["cpu"][0]) / totals["c_fast"][0]
            updates["hlRustVsCfastPct"] = f"{pct:.1f}\\%"
            updates["eThreeCpuVsCfastSp"] = fmt_sp(totals["c_fast"][0], totals["cpu"][0], 2)

    # Experiment 4A
    if "exp4a_slab_memory.csv" in data:
        bmap = {"512M": "FiveTwelve", "1G": "OneG", "2G": "TwoG", "4G": "FourG"}
        for r in data["exp4a_slab_memory.csv"]:
            bt = "Cpu" if r["backend"] == "cpu" else "Cuda"
            mt = bmap.get(r["mem_budget"])
            if mt:
                updates[f"eFourA{bt}{mt}W"] = fmt_s(float(r["wall_s"]))
                updates[f"eFourA{bt}{mt}R"] = f"{float(r['rss_mib']):,.1f}"

    # Experiment 4B
    if "exp4b_wl_fl_sensitivity.csv" in data:
        wmap = {"12": "Twelve", "16": "Sixteen", "24": "TwentyFour", "32": "ThirtyTwo"}
        fmap = {"20": "Twenty", "1000": "Thou"}
        for r in data["exp4b_wl_fl_sensitivity.csv"]:
            wt = wmap.get(r["wl"]); ft = fmap.get(r["fl"])
            if wt and ft:
                base = f"eFourBW{wt}F{ft}"
                updates[base + "Cpu"] = fmt_s(float(r["cpu_wall_s"]))
                updates[base + "Cuda"] = fmt_s(float(r["cuda_wall_s"]))
                updates[base + "Sp"] = f"{float(r['wall_speedup']):.2f}$\\times$"
                updates[base + "Frags"] = fmt_int(int(r["fragments"]))

    tex = MAIN_TEX.read_text(encoding="utf-8")
    for k, v in updates.items():
        tex = set_macro(tex, k, v)
    MAIN_TEX.write_text(tex, encoding="utf-8")
    print(f"Updated {len(updates)} result macros in {MAIN_TEX}")

if __name__ == "__main__":
    main()
