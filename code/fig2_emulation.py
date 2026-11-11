#!/usr/bin/env python3
"""Manuscript Figure 2 from results/emu_mhc_summary.tsv (as_emulate.py output).

A  keep-rate under the primary-hit rule, per MHC locus, for four capture references
   (alt-free; ALT haplotypes only; HLA contigs only; both = hs38DH)
B  CFB keep-rate under each decision rule, alt-free vs hs38DH
C  calibration on the cohort: predicted ejected fraction (1 - keep) under each rule versus the
   observed fraction of reads < 1 kb (median per sample, Additional file 5), CFB and C4A-C4B,
   under each capture reference

Usage: python3 fig2_emulation.py results/emu_mhc_summary.tsv results/fig2_emulation
"""
import csv, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

src, out = sys.argv[1], sys.argv[2]
rows = list(csv.DictReader(open(src), delimiter="\t"))
PL = "400"


def k(locus, ref, rule, bn="5"):
    for r in rows:
        if (r["locus"], r["ref"], r["rule"], r["best_n"], r["prefix_len"]) == (locus, ref, rule, bn, PL):
            return float(r["keep_rate"])
    raise KeyError((locus, ref, rule, bn))


# observed: median per-sample % of reads < 1 kb (Additional file 5, n = 108 / 24)
OBS = {("CFB", "hs38DH"): 75.2, ("C4A_C4B", "hs38DH"): 75.9,
       ("CFB", "noalt"): 2.55, ("C4A_C4B", "noalt"): 4.35}

BLUE, ORANGE, PURPLE, INK = "#1f6fb2", "#c75a34", "#7d5ba6", "#3a3a3a"
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": "#888", "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK})
fig = plt.figure(figsize=(7.2, 2.7))
gs = fig.add_gridspec(1, 3, width_ratios=[1.5, 1.0, 1.3], wspace=0.42)

# ---------------- A
ax = fig.add_subplot(gs[0])
loci = [("CFB", "CFB"), ("C2_CFB", "C2–CFB"), ("C4A_C4B", "C4A–C4B"), ("HLA-A", "HLA-A"),
        ("HLA-DRB1", "HLA-DRB1"), ("ctrl_chr6_43Mb", "chr6 ctrl")]
refs = [("noalt", "alt-free", BLUE), ("ALT_only", "+ ALT only", ORANGE),
        ("HLA_only", "+ HLA only", PURPLE), ("hs38DH", "hs38DH (ALT + HLA)", INK)]
w = 0.19
x = np.arange(len(loci))
for j, (ref, lab, col) in enumerate(refs):
    v = [k(l, ref, "primary") for l, _ in loci]
    ax.bar(x + (j - 1.5) * w, v, w * 0.88, color=col, label=lab)
ax.set_xticks(x)
ax.set_xticklabels([l for _, l in loci], rotation=35, ha="right")
ax.set_ylim(0, 1.05)
ax.set_ylabel("keep-rate (primary-hit rule)")
ax.yaxis.grid(True, color="#e6e6e6", lw=0.6)
ax.set_axisbelow(True)
ax.legend(fontsize=6.5, frameon=False, loc="upper center", bbox_to_anchor=(0.5, 1.22), ncol=2)
ax.set_title("A", loc="left", fontweight="bold", x=-0.18)

# ---------------- B
ax = fig.add_subplot(gs[1])
rules = [("primary", "5", "primary hit"), ("readfish", "5", "readfish, best_n 5"),
         ("readfish", "20", "readfish, best_n 20"), ("mapq1", "5", "primary, MAPQ ≥ 1")]
x = np.arange(len(rules))
for j, (ref, col, lab) in enumerate([("noalt", BLUE, "alt-free"), ("hs38DH", ORANGE, "hs38DH")]):
    v = [k("CFB", ref, r, bn) for r, bn, _ in rules]
    b = ax.bar(x + (j - 0.5) * 0.36, v, 0.32, color=col, label=lab)
    if ref == "hs38DH":
        for xi, vi in zip(x, v):
            ax.text(xi + 0.18, vi + 0.02, f"{vi:.2f}", ha="center", fontsize=6.5, color=INK)
ax.set_xticks(x)
ax.set_xticklabels([l for _, _, l in rules], fontsize=6.5, rotation=35, ha="right")
ax.set_ylim(0, 1.12)
ax.set_ylabel("CFB keep-rate")
ax.yaxis.grid(True, color="#e6e6e6", lw=0.6)
ax.set_axisbelow(True)
ax.legend(fontsize=6.5, frameon=False, loc="upper center", bbox_to_anchor=(0.5, 1.22), ncol=2)
ax.set_title("B", loc="left", fontweight="bold", x=-0.3)

# ---------------- C
ax = fig.add_subplot(gs[2])
cases = [("CFB", "hs38DH", "CFB\nhs38DH"), ("C4A_C4B", "hs38DH", "C4\nhs38DH"),
         ("CFB", "noalt", "CFB\nalt-free"), ("C4A_C4B", "noalt", "C4\nalt-free")]
marks = [("primary", "5", "o", "primary hit"), ("readfish", "5", "s", "readfish"),
         ("mapq1", "5", "^", "MAPQ ≥ 1")]
x = np.arange(len(cases))
for j, (r, bn, m, lab) in enumerate(marks):
    v = [100 * (1 - k(l, ref, r, bn)) for l, ref, _ in cases]
    ax.scatter(x + (j - 1.5) * 0.15, v, marker=m, s=26, facecolor="white", edgecolor=INK, lw=1, label=lab, zorder=3)
ax.scatter(x + 0.3, [OBS[(l, ref)] for l, ref, _ in cases], marker="D", s=34, color=ORANGE,
           edgecolor="white", lw=0.8, label="observed (% reads < 1 kb)", zorder=4)
ax.set_xticks(x)
ax.set_xticklabels([l for _, _, l in cases], fontsize=6.5)
ax.set_ylim(-5, 105)
ax.set_ylabel("% molecules ejected")
ax.yaxis.grid(True, color="#e6e6e6", lw=0.6)
ax.set_axisbelow(True)
ax.legend(fontsize=6.2, frameon=False, loc="upper center", bbox_to_anchor=(0.5, 1.3), ncol=2,
          handletextpad=0.3, columnspacing=0.8)
ax.set_title("C", loc="left", fontweight="bold", x=-0.25)

for ext in ("png", "pdf"):
    fig.savefig(f"{out}.{ext}", dpi=300, bbox_inches="tight")
print(f"wrote {out}.png/.pdf")
