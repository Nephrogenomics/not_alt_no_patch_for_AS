#!/usr/bin/env python3
"""Manuscript Figure 5: is the capture loss caused by patches general, or set by each patch's sequence?

Inputs: patch_divergence.py windows/segments tables and the as_emulate.py summary of run_p14_patches.sh.
Each 10-kb window of every ALT / fix / novel contig alignment was sampled from the primary haplotype and
from the contig haplotype; the keep-rate is read under the reference = primary assembly + that contig's class.

A  keep-rate against local divergence between the contig and its primary region (differences per kb),
   for reads from either haplotype; primary-hit rule (solid) and readfish rule (dashed)
B  per-contig keep-rate by contig class and read haplotype (primary-hit rule)
C  panel genes: keep-rate for reads from either haplotype (primary-hit rule)

Writes <out>.png/.pdf, <out>_windows.tsv, <out>_contigs.tsv and prints the numbers quoted in the text.
"""
import argparse, csv
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OWN = {"alt": "plus_ALT", "fix": "plus_fix", "novel": "plus_novel"}
BINS = [(-0.1, 0.0, "0"), (0.0, 0.5, "≤0.5"), (0.5, 1, "0.5–1"), (1, 2, "1–2"), (2, 5, "2–5"),
        (5, 10, "5–10"), (10, 20, "10–20"), (20, 1e9, ">20")]
PANEL = ["ADAMTS13", "ACTN4", "WDR73", "PTEN", "F13B", "CFH", "CFHR2", "CFHR4", "CFHR5", "C2", "CFB", "C4A"]
BLUE, ORANGE, INK, GRID = "#1f6fb2", "#c75a34", "#3a3a3a", "#e6e6e6"
SHORT = {"alt": "ALT", "fix": "fix", "novel": "novel"}
CLSLAB = {"alt": "GRCh38 ALT", "fix": "fix patch", "novel": "novel patch"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--windows", required=True)
    ap.add_argument("--segments", required=True)
    ap.add_argument("--emu", required=True)
    ap.add_argument("--genes", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    keep = {}
    for r in csv.DictReader(open(a.emu), delimiter="\t"):
        if r["prefix_len"] != "400" or r["best_n"] != "5" or r["rule"] not in ("primary", "readfish"):
            continue
        w, src = r["locus"].split("|src=")
        keep[(w, src, r["ref"], r["rule"])] = (float(r["keep_rate"]), int(r["n"]))
    W = [w for w in csv.DictReader(open(a.windows), delimiter="\t")]
    rows = []
    for w in W:
        rec = dict(w)
        ok = True
        for src in ("primary", "contig"):
            for ref, tag in ((OWN[w["cls"]], "own"), ("primary", "noalt"), ("hg38.p14", "p14")):
                for rule in ("primary", "readfish"):
                    v = keep.get((w["window"], src, ref, rule))
                    if v is None or v[1] < 5:
                        ok = False
                        continue
                    rec[f"keep_{rule}_{tag}_{src}"] = v[0]
        if ok:
            rows.append(rec)
    d = np.array([float(r["diff_per_kb"]) for r in rows])
    cl = np.array([r["cls"] for r in rows])
    K = lambda key: np.array([r[key] for r in rows])
    kp, kc = K("keep_primary_own_primary"), K("keep_primary_own_contig")
    rp, rc = K("keep_readfish_own_primary"), K("keep_readfish_own_contig")

    # ---------------------------------------------------------------- tables
    cols = ["window", "contig", "cls", "chrom", "w_p_start", "w_p_end", "w_c_start", "w_c_end", "strand",
            "aligned", "mm", "ins_ev", "del_ev", "diff_per_kb"] + sorted(k for k in rows[0] if k.startswith("keep_"))
    with open(f"{a.out}_windows.tsv", "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in rows:
            o.write("\t".join(f"{r[c]:.3f}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")
    seg = defaultdict(list)
    for r in rows:
        seg[r["contig"]].append(r)
    segmeta = {}
    for s in csv.DictReader(open(a.segments), delimiter="\t"):
        m = segmeta.setdefault(s["contig"], dict(cls=s["cls"], chrom=s["chrom"], c_len=0, aligned_c=0,
                                                 matches=0, blen=0))
        for k in ("c_len", "aligned_c", "matches", "blen"):
            m[k] += int(s[k])
    crow = []
    for c, m in sorted(segmeta.items()):
        ws = seg.get(c, [])
        mean = lambda k: float(np.mean([w[k] for w in ws])) if ws else float("nan")
        crow.append(dict(contig=c, cls=m["cls"], chrom=m["chrom"], length=m["c_len"],
                         frac_aligned=m["aligned_c"] / m["c_len"] if m["c_len"] else 0,
                         diff_per_kb=1000 * (1 - m["matches"] / m["blen"]) if m["blen"] else float("nan"),
                         n_windows=len(ws),
                         keep_primary_hap=mean("keep_primary_own_primary"),
                         keep_contig_hap=mean("keep_primary_own_contig"),
                         keep_readfish_primary_hap=mean("keep_readfish_own_primary"),
                         keep_readfish_contig_hap=mean("keep_readfish_own_contig"),
                         keep_noalt_contig_hap=mean("keep_primary_noalt_contig")))
    ck = list(crow[0].keys())
    with open(f"{a.out}_contigs.tsv", "w") as o:
        o.write("\t".join(ck) + "\n")
        for r in crow:
            o.write("\t".join(f"{r[k]:.4f}" if isinstance(r[k], float) else str(r[k]) for k in ck) + "\n")

    # ---------------------------------------------------------------- numbers for the text
    print(f"windows analysed: {len(rows)} (" + ", ".join(f"{c} {np.sum(cl == c)}" for c in OWN) + ")")
    print(f"contigs with >=1 window: {sum(1 for r in crow if r['n_windows'])} / {len(crow)}")
    for lo, hi, lab in BINS:
        m = (d > lo) & (d <= hi)
        if m.sum():
            print(f"  diff/kb {lab:>6}: n={m.sum():5d}  keep primary-hap {kp[m].mean():.2f}  contig-hap "
                  f"{kc[m].mean():.2f}  | readfish {rp[m].mean():.2f} / {rc[m].mean():.2f}")
    for c in OWN:
        m = cl == c
        print(f"  {c:5s}: median diff/kb {np.median(d[m]):.2f}; identical windows {np.mean(d[m] == 0):.0%}; "
              f"keep primary-hap {kp[m].mean():.2f}, contig-hap {kc[m].mean():.2f}, "
              f"heterozygote mean {((kp[m] + kc[m]) / 2).mean():.2f}")
    from scipy.stats import spearmanr
    print(f"  Spearman keep vs divergence: primary-hap rho={spearmanr(d, kp)[0]:.2f}, "
          f"contig-hap rho={spearmanr(d, kc)[0]:.2f}")

    # ---------------------------------------------------------------- figure
    plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.edgecolor": "#888", "axes.labelcolor": INK, "xtick.color": INK,
                         "ytick.color": INK})
    fig = plt.figure(figsize=(7.2, 2.8))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.45, 0.95, 1.15], wspace=0.62)

    # A
    ax = fig.add_subplot(gs[0])
    x = np.arange(len(BINS))
    for arr, arr_rf, col, lab in ((kp, rp, BLUE, "reads from primary haplotype"),
                                  (kc, rc, ORANGE, "reads from contig haplotype")):
        mu, lo_, hi_, mu_rf = [], [], [], []
        for lo, hi, _ in BINS:
            m = (d > lo) & (d <= hi)
            v = arr[m]
            mu.append(v.mean() if len(v) else np.nan)
            se = 1.96 * v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else np.nan
            lo_.append(v.mean() - se if len(v) else np.nan)
            hi_.append(v.mean() + se if len(v) else np.nan)
            mu_rf.append(arr_rf[m].mean() if len(v) else np.nan)
        ax.fill_between(x, lo_, hi_, color=col, alpha=0.25, lw=0)
        ax.plot(x, mu, color=col, lw=2, marker="o", ms=4, label=lab)
        ax.plot(x, mu_rf, color=col, lw=1.2, ls="--")
    ax.plot([], [], color="#777", lw=1.2, ls="--", label="readfish rule")
    for i, (lo, hi, _) in enumerate(BINS):
        ax.text(i, 1.05, f"{int(((d > lo) & (d <= hi)).sum())}", ha="center", fontsize=4.8, color="#777", rotation=90, va="bottom")
    ax.text(-0.75, 1.05, "windows", ha="right", fontsize=4.8, color="#777", va="bottom")
    ax.set_xticks(x)
    ax.set_xticklabels([b[2] for b in BINS], rotation=35, ha="right")
    ax.set_xlabel("differences per kb, contig vs primary")
    ax.set_ylim(0, 1.22)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_ylabel("keep-rate (primary-hit rule)")
    ax.axhline(0.5, color="#bbb", lw=0.6, ls=":")
    ax.yaxis.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    fig.legend(*ax.get_legend_handles_labels(), fontsize=6.5, frameon=False, loc="lower center",
               bbox_to_anchor=(0.5, -0.2), ncol=3)
    ax.set_title("a", loc="left", fontweight="bold", x=-0.2)

    # B
    ax = fig.add_subplot(gs[1])
    rng = np.random.default_rng(1)
    pos = 0
    ticks, tl = [], []
    for c in OWN:
        cr = [r for r in crow if r["cls"] == c and r["n_windows"]]
        for key, col, off in (("keep_primary_hap", BLUE, -0.2), ("keep_contig_hap", ORANGE, 0.2)):
            v = np.array([r[key] for r in cr])
            ax.scatter(pos + off + rng.uniform(-0.12, 0.12, len(v)), v, s=4, color=col, alpha=0.45, lw=0)
            ax.plot([pos + off - 0.16, pos + off + 0.16], [np.median(v)] * 2, color=INK, lw=1.6)
        ticks.append(pos)
        tl.append(f"{SHORT[c]}\n{len(cr)}")
        pos += 1
    ax.set_xticks(ticks)
    ax.set_xticklabels(tl, fontsize=6.5)
    ax.set_ylim(-0.02, 1.08)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_ylabel("keep-rate per contig")
    ax.set_xlabel("contig class (n)")
    ax.axhline(0.5, color="#bbb", lw=0.6, ls=":")
    ax.yaxis.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.set_title("b", loc="left", fontweight="bold", x=-0.25)

    # C
    ax = fig.add_subplot(gs[2])
    G = {g["gene"]: g for g in csv.DictReader(open(a.genes), delimiter="\t") if g["counterpart_class"] != "control"}
    names, vp, vc, cls_g = [], [], [], []
    for g in PANEL:
        if g not in G:
            continue
        c, s, e = G[g]["chrom"], int(G[g]["start"]), int(G[g]["end"])
        ws = [r for r in rows if r["chrom"] == c and int(r["w_p_start"]) < e and int(r["w_p_end"]) > s]
        if not ws:
            continue
        names.append(g)
        vp.append(np.mean([r["keep_primary_own_primary"] for r in ws]))
        vc.append(np.mean([r["keep_primary_own_contig"] for r in ws]))
        cls_g.append("/".join(sorted({r["cls"] for r in ws})))
        print(f"  gene {g:9s} {cls_g[-1]:9s} windows={len(ws):3d} primary-hap {vp[-1]:.2f} contig-hap {vc[-1]:.2f}")
    y = np.arange(len(names))[::-1]
    ax.barh(y + 0.19, vp, 0.36, color=BLUE)
    ax.barh(y - 0.19, vc, 0.36, color=ORANGE)
    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=6.5, style="italic")
    for yy, c in zip(y, cls_g):
        ax.text(1.02, yy, c, fontsize=5.5, color="#777", va="center", transform=ax.get_yaxis_transform())
    ax.set_xlim(0, 1.0)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax.axvline(0.5, color="#bbb", lw=0.6, ls=":")
    ax.set_xlabel("keep-rate (primary-hit rule)")
    ax.xaxis.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.set_title("c", loc="left", fontweight="bold", x=-0.38)

    fig.savefig(f"{a.out}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{a.out}.pdf", bbox_inches="tight")


if __name__ == "__main__":
    main()
