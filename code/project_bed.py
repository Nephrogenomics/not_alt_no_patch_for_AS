#!/usr/bin/env python3
"""Project a primary-assembly BED onto the ALT / fix / novel contigs that cover it.

Uses the 10-kb alignment windows of code/patch_divergence.py (_windows.tsv): each window links a primary
interval to a contig interval. A BED interval overlapping a window is projected by linear interpolation
inside that window (exact on identical windows, within the window's indel content elsewhere), and the
projected pieces are merged per contig. This is the "projected BED" a user could load on the instrument
so that a read whose primary alignment lands on a contig copy of a target is still on target.

Output BED (to stdout or --out): contig, start, end, name (BED name or chrom:start-end) | contig class.
With --with-primary the input intervals are written first, so the output is the full projected target BED.
With --loci the output keeps one line per (interval, contig), named <name>|src=contig, for sampling.
"""
import argparse, bisect, sys
from collections import defaultdict


def read_bed(p):
    out = []
    for ln in open(p):
        if not ln.strip() or ln.startswith(("#", "track", "browser")):
            continue
        f = ln.rstrip("\n").split("\t")
        out.append((f[0], int(f[1]), int(f[2]), f[3] if len(f) > 3 and f[3] else f"{f[0]}:{f[1]}-{f[2]}"))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bed", required=True)
    ap.add_argument("--windows", required=True, help="_windows.tsv of patch_divergence.py")
    ap.add_argument("--classes", default="alt,fix,novel", help="contig classes to project onto")
    ap.add_argument("--contigs", help="optional file of contig names to restrict to (e.g. those in a FASTA .fai)")
    ap.add_argument("--with-primary", action="store_true")
    ap.add_argument("--loci", action="store_true")
    ap.add_argument("--out")
    a = ap.parse_args()
    keep_cls = set(a.classes.split(","))
    allow = None
    if a.contigs:
        allow = {ln.split("\t")[0].strip() for ln in open(a.contigs) if ln.strip()}

    win = defaultdict(list)            # chrom -> [(p0, p1, contig, c0, c1, strand, cls)]
    with open(a.windows) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ix = {k: i for i, k in enumerate(hdr)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if f[ix["cls"]] not in keep_cls or (allow is not None and f[ix["contig"]] not in allow):
                continue
            win[f[ix["chrom"]]].append((int(f[ix["w_p_start"]]), int(f[ix["w_p_end"]]), f[ix["contig"]],
                                        int(f[ix["w_c_start"]]), int(f[ix["w_c_end"]]), f[ix["strand"]],
                                        f[ix["cls"]]))
    starts = {}
    for c in win:
        win[c].sort()
        starts[c] = [w[0] for w in win[c]]
    maxlen = max((w[1] - w[0] for c in win for w in win[c]), default=0)

    pieces = defaultdict(list)         # (name, contig, cls) -> [(c0, c1)]
    for chrom, s, e, name in read_bed(a.bed):
        ws = win.get(chrom)
        if not ws:
            continue
        i = max(0, bisect.bisect_left(starts[chrom], s - maxlen))
        while i < len(ws) and ws[i][0] < e:
            p0, p1, ctg, c0, c1, strand, cls = ws[i]
            i += 1
            lo, hi = max(s, p0), min(e, p1)
            if lo >= hi:
                continue
            k = (c1 - c0) / (p1 - p0)
            if strand == "+":
                x0, x1 = c0 + (lo - p0) * k, c0 + (hi - p0) * k
            else:
                x0, x1 = c1 - (hi - p0) * k, c1 - (lo - p0) * k
            pieces[(name, ctg, cls)].append((int(x0), int(round(x1))))

    o = open(a.out, "w") if a.out else sys.stdout
    if a.with_primary:
        for chrom, s, e, name in read_bed(a.bed):
            o.write(f"{chrom}\t{s}\t{e}\t{name}\n")
    n = 0
    for (name, ctg, cls), iv in sorted(pieces.items(), key=lambda x: (x[0][1], min(x[1]))):
        iv.sort()
        merged = [list(iv[0])]
        for x0, x1 in iv[1:]:
            if x0 <= merged[-1][1] + 1000:     # window joins: tolerate the indel slack of a 10-kb window
                merged[-1][1] = max(merged[-1][1], x1)
            else:
                merged.append([x0, x1])
        for x0, x1 in merged:
            n += 1
            lab = f"{name}|src=contig" if a.loci else f"{name}|{cls}"
            o.write(f"{ctg}\t{x0}\t{x1}\t{lab}\n")
    print(f"projected {n} intervals onto {len({k[1] for k in pieces})} contigs", file=sys.stderr)


if __name__ == "__main__":
    main()
