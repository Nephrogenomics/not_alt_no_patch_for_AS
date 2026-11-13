#!/usr/bin/env python3
"""Tables of the projected-BED emulation (code/run_projected_bed.sh).

--panel   <panel>_summary.tsv: CKD panel genes with a counterpart contig, refs primary / plus_ALT / hg38.p14,
          each decided against the standard BED and the --extra-bed tags clinical, projALT and proj
--windows <win>_summary.tsv + --div div_windows.tsv: the 10-kb windows of Fig. 5, each contig tested against
          primary + its own class, standard BED versus projected BED
Writes <out>_genes.tsv (one row per gene, source haplotype and rule), <out>_divergence.tsv (by divergence bin)
and <out>_numbers.txt (the values quoted in the README and the manuscript).
"""
import argparse, collections, csv, statistics as st

BINS = [("0", 0, 0), ("<=1", 0, 1), ("1-5", 1, 5), ("5-20", 5, 20), (">20", 20, 1e12)]
OWN = {"alt": "plus_ALT", "fix": "plus_fix", "novel": "plus_novel"}
BEDS = ["", "|clinical", "|projALT", "|proj"]
GENES = ["CFB", "C2", "C4A", "HLA-A", "HLA-DRB1", "ADAMTS13", "ACTN4", "WDR73", "PTEN", "CFH", "CFHR2",
         "CFHR4", "CFHR5", "F13B"]


def dbin(d):
    if d == 0:
        return "0"
    for lab, lo, hi in BINS[1:]:
        if lo < d <= hi:
            return lab


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel", required=True)
    ap.add_argument("--windows", required=True)
    ap.add_argument("--div", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    txt = []

    # ---- panel genes
    R = collections.defaultdict(dict)          # (gene, src, rule) -> ref -> keep
    for r in csv.DictReader(open(a.panel), delimiter="\t"):
        g, src = r["locus"].split("|src=")
        R[(g, src, r["rule"])][r["ref"]] = float(r["keep_rate"])
    refs = ["primary"] + [f"{b}{x}" for b in ("plus_ALT", "hg38.p14") for x in BEDS]
    with open(f"{a.out}_genes.tsv", "w") as o:
        o.write("gene\tsource\trule\t" + "\t".join(r.replace("|", "_") for r in refs) + "\n")
        for (g, src, rule), d in sorted(R.items()):
            o.write(f"{g}\t{src}\t{rule}\t" + "\t".join(f"{d.get(x, float('nan')):.3f}" for x in refs) + "\n")
    txt.append(f"panel genes with a counterpart contig: {len({k[0] for k in R})}")
    txt.append("hg38.p14, primary-hit rule, keep-rate primary / contig haplotype "
               "(standard | clinical | projected ALT only | projected all; primary-only reference):")
    for g in GENES:
        p, c = R.get((g, "primary", "primary")), R.get((g, "contig", "primary"))
        if not p:
            continue
        cells = [f"{p['hg38.p14' + x]:.2f}/{c['hg38.p14' + x]:.2f}" if c else f"{p['hg38.p14' + x]:.2f}/-"
                 for x in BEDS]
        base = f"{p['primary']:.2f}/{c['primary']:.2f}" if c else f"{p['primary']:.2f}/-"
        txt.append(f"  {g:9} " + " | ".join(cells) + f"   ; primary-only {base}")
    txt.append("median keep relative to the primary-only reference, hg38.p14 (primary / contig haplotype):")
    for rule in ("primary", "readfish", "mapq1"):
        cells = []
        for x in BEDS:
            v = {}
            for src in ("primary", "contig"):
                xs = [d["hg38.p14" + x] / d["primary"] for (g, s, ru), d in R.items()
                      if ru == rule and s == src and d.get("primary", 0) > 0]
                v[src] = st.median(xs) if xs else float("nan")
            cells.append(f"{v['primary']:.2f}/{v['contig']:.2f}")
        txt.append(f"  {rule:9} " + " | ".join(cells))

    # ---- genome-wide windows
    win = {r["window"]: (r["cls"], float(r["diff_per_kb"])) for r in csv.DictReader(open(a.div), delimiter="\t")}
    acc = collections.defaultdict(list)
    overall = collections.defaultdict(list)
    for r in csv.DictReader(open(a.windows), delimiter="\t"):
        w, src = r["locus"].split("|src=")
        base, _, tag = r["ref"].partition("|")
        k = float(r["keep_rate"])
        overall[(r["ref"], r["rule"])].append(k)
        if w in win and base == OWN[win[w][0]]:
            acc[(r["rule"], tag, src, dbin(win[w][1]))].append(k)
    with open(f"{a.out}_divergence.tsv", "w") as o:
        o.write("rule\tdiff_per_kb\tn_windows\tstd_primary_hap\tstd_contig_hap\tproj_primary_hap\tproj_contig_hap\n")
        for rule in ("primary", "readfish", "mapq1"):
            for lab, *_ in BINS:
                v = [acc[(rule, t, s, lab)] for t in ("", "proj") for s in ("primary", "contig")]
                o.write(f"{rule}\t{lab}\t{len(v[0])}\t" + "\t".join(f"{st.mean(x):.3f}" for x in v) + "\n")
    txt.append("genome-wide windows, mean keep-rate over molecules (primary-hit / readfish / mapq1), "
               "standard BED -> projected BED:")
    for ref in ("plus_ALT", "plus_fix", "plus_novel", "hg38.p14"):
        s = [f"{st.mean(overall[(ref, ru)]):.3f}" for ru in ("primary", "readfish", "mapq1")]
        p = [f"{st.mean(overall[(ref + '|proj', ru)]):.3f}" for ru in ("primary", "readfish", "mapq1")]
        txt.append(f"  {ref:10} {' / '.join(s)}  ->  {' / '.join(p)}")
    txt.append("by divergence, own class, primary-hit rule (primary / contig haplotype), standard -> projected:")
    for rule in ("primary", "readfish"):
        for lab, *_ in BINS:
            v = [st.mean(acc[(rule, t, s, lab)]) for t in ("", "proj") for s in ("primary", "contig")]
            txt.append(f"  {rule:8} {lab:5} {v[0]:.2f}/{v[1]:.2f} -> {v[2]:.2f}/{v[3]:.2f}")
    open(f"{a.out}_numbers.txt", "w").write("\n".join(txt) + "\n")
    print("\n".join(txt))


if __name__ == "__main__":
    main()
