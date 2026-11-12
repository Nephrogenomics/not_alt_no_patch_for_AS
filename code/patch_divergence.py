#!/usr/bin/env python3
"""Align every ALT / fix / novel contig of a patched assembly (UCSC hg38.p14) to its primary region and
cut the alignment into windows of local divergence, for the per-patch "no patch" emulation (Fig. 5).

For each row of the UCSC altLocations / fixLocations tables (primary chrom:start-end <-> contig[:a-b]),
the contig segment is aligned to the primary region with minimap2 (mappy, preset asm20, cs tag). The
alignment is walked and cut every --window aligned primary bases. Each window gets its primary interval,
its contig interval, and its differences: mismatches, insertion and deletion events and bases.

Outputs (prefix --out):
  _segments.tsv   one row per location row: contig, class, lengths, aligned fraction, identity, strand
  _windows.tsv    one row per window: coordinates on both sequences, differences, diffs per kb
  _loci.bed       for as_emulate.py --loci-bed: each window twice, sampled from the primary interval
                  (name W|src=primary) and from the contig interval (name W|src=contig)
  _target.bed     the primary regions of all location rows (target BED of the emulation)
"""
import argparse, gzip, re, sys
import json, subprocess
from concurrent.futures import ThreadPoolExecutor

CHUNK_DEFAULT = 50000
CS = re.compile(r"(:\d+|\*[a-z][a-z]|[+\-][A-Za-z]+)")


def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def read_locations(paths):
    rows = []
    for p in paths:
        for ln in opener(p):
            if not ln.strip() or ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if f[0].isdigit():
                f = f[1:]
            if "_" in f[0]:          # UCSC tables also hold the reciprocal contig -> primary rows
                continue
            name = f[3]
            ctg, rng = (name.split(":", 1) + [None])[:2]
            qa = qb = None
            if rng:
                a, b = rng.split("-")
                qa, qb = int(a) - 1, int(b)
            rows.append((f[0], int(f[1]), int(f[2]), ctg, qa, qb))
    return rows


def cls(ctg, orig):
    if ctg.endswith("_fix"):
        return "fix"
    return "alt" if ctg in orig else "novel"


def one(job, chunk=None):
    import mappy, pysam
    fa_path, (chrom, ps, pe, ctg, qa, qb), klass, W = job
    CHUNK = chunk or CHUNK_DEFAULT
    fa = pysam.FastaFile(fa_path)
    clen = fa.get_reference_length(ctg)
    qa = 0 if qa is None else qa
    qb = clen if qb is None else min(qb, clen)
    ref = fa.fetch(chrom, ps, pe).upper()
    qry = fa.fetch(ctg, qa, qb).upper()
    seg = dict(contig=ctg, cls=klass, chrom=chrom, p_start=ps, p_end=pe, c_start=qa, c_end=qb,
               p_len=pe - ps, c_len=qb - qa, aligned_p=0, aligned_c=0, matches=0, blen=0, strand=".")
    wins = []
    al = mappy.Aligner(seq=ref, preset="asm20")
    if not al:
        return seg, wins
    hits = []                     # map in chunks: long cs alignments crash mappy 2.30 on some patches
    for c0 in range(0, len(qry), CHUNK):
        sub = qry[c0:c0 + CHUNK]
        if len(sub) < 2000 or sub.count("N") > 0.5 * len(sub):
            continue
        hs = [h for h in al.map(sub, cs=True) if h.is_primary]
        if hs:
            h = max(hs, key=lambda h: h.mlen)
            hits.append((c0, h))
    if not hits:
        return seg, wins
    fw = sum(h.r_en - h.r_st for _, h in hits if h.strand > 0)
    seg["strand"] = "+" if fw * 2 >= sum(h.r_en - h.r_st for _, h in hits) else "-"
    for c0, h in hits:
        seg["aligned_p"] += h.r_en - h.r_st
        seg["aligned_c"] += h.q_en - h.q_st
        seg["matches"] += h.mlen
        seg["blen"] += h.blen
        # walk: r forward from r_st; q forward from q_st (+) or backward from q_en (-)
        r = h.r_st
        q = c0 + (h.q_st if h.strand > 0 else h.q_en)
        step = 1 if h.strand > 0 else -1
        cur = None

        def flush(c):
            if c and c["n"] > 0:
                qs, qe = sorted((c["q0"], c["q1"]))
                wins.append(dict(contig=ctg, cls=klass, chrom=chrom,
                                 w_p_start=ps + c["r0"], w_p_end=ps + c["r1"],
                                 w_c_start=qa + qs, w_c_end=qa + qe, strand="+" if h.strand > 0 else "-",
                                 aligned=c["n"], mm=c["mm"], ins_ev=c["ie"], ins_bp=c["ib"],
                                 del_ev=c["de"], del_bp=c["db"]))

        def new():
            return dict(r0=r, r1=r, q0=q, q1=q, n=0, mm=0, ie=0, ib=0, de=0, db=0)

        cur = new()
        for op in CS.findall(h.cs):
            if op[0] == ":":
                n = int(op[1:])
                while n > 0:
                    take = min(n, W - cur["n"])
                    r += take; q += step * take; n -= take
                    cur["n"] += take; cur["r1"] = r; cur["q1"] = q
                    if cur["n"] >= W:
                        flush(cur); cur = new()
            elif op[0] == "*":
                r += 1; q += step
                cur["n"] += 1; cur["mm"] += 1; cur["r1"] = r; cur["q1"] = q
                if cur["n"] >= W:
                    flush(cur); cur = new()
            elif op[0] == "+":            # insertion in the contig (query)
                L = len(op) - 1
                q += step * L
                cur["ie"] += 1; cur["ib"] += L; cur["q1"] = q
            else:                          # deletion from the contig (bases only in primary)
                L = len(op) - 1
                r += L
                cur["de"] += 1; cur["db"] += L; cur["r1"] = r
        if cur["n"] >= W // 2:
            flush(cur)
    return seg, wins


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fa", required=True, help="indexed FASTA holding primary chromosomes and the contigs")
    ap.add_argument("--locations", action="append", required=True)
    ap.add_argument("--orig-alt", required=True, help="list of the original GRCh38 ALT contig names")
    ap.add_argument("--window", type=int, default=10000)
    ap.add_argument("--min-window", type=int, default=2000, help="drop windows shorter than this on either side")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    import pysam
    present = set(pysam.FastaFile(a.fa).references)
    orig = set(open(a.orig_alt).read().split())
    rows = [r for r in read_locations(a.locations) if r[3] in present and r[0] in present]
    jobs = [(a.fa, r, cls(r[3], orig), a.window) for r in rows]
    jobs.sort(key=lambda j: -(j[1][2] - j[1][1]))          # big ones first
    print(f"{len(jobs)} location rows", file=sys.stderr)
    segs, wins = [], []
    def isolated(job):
        """run one job in its own interpreter: a mappy segfault then costs one contig, not the run"""
        for chunk in (50000, 10000, 3000):
            p = subprocess.run([sys.executable, __file__, "--one", json.dumps([job, chunk])],
                               capture_output=True, text=True)
            if p.returncode == 0:
                return json.loads(p.stdout)
            print(f"  [retry] {job[1][3]} chunk {chunk}: exit {p.returncode}", file=sys.stderr)
        _, (chrom, ps, pe, ctg, qa, qb), klass, _ = job
        print(f"  [fail] {ctg}: not aligned", file=sys.stderr)
        return dict(contig=ctg, cls=klass, chrom=chrom, p_start=ps, p_end=pe, c_start=qa or 0, c_end=qb or 0,
                    p_len=pe - ps, c_len=(qb or 0) - (qa or 0), aligned_p=0, aligned_c=0, matches=0, blen=0,
                    strand="failed"), []

    with ThreadPoolExecutor(a.threads) as ex:
        for i, (s, w) in enumerate(ex.map(isolated, jobs)):
            segs.append(s); wins += w
            if (i + 1) % 50 == 0:
                print(f"  {i + 1}/{len(jobs)}", file=sys.stderr)

    sk = ["contig", "cls", "chrom", "p_start", "p_end", "c_start", "c_end", "p_len", "c_len", "aligned_p",
          "aligned_c", "matches", "blen", "strand"]
    with open(f"{a.out}_segments.tsv", "w") as o:
        o.write("\t".join(sk + ["frac_aligned_c", "identity"]) + "\n")
        for s in sorted(segs, key=lambda s: (s["cls"], s["contig"], s["p_start"])):
            fa_c = s["aligned_c"] / s["c_len"] if s["c_len"] else 0
            idt = s["matches"] / s["blen"] if s["blen"] else float("nan")
            o.write("\t".join(str(s[k]) for k in sk) + f"\t{fa_c:.4f}\t{idt:.6f}\n")

    wk = ["contig", "cls", "chrom", "w_p_start", "w_p_end", "w_c_start", "w_c_end", "strand", "aligned", "mm",
          "ins_ev", "ins_bp", "del_ev", "del_bp"]
    n = 0
    with open(f"{a.out}_windows.tsv", "w") as o, open(f"{a.out}_loci.bed", "w") as lb:
        o.write("window\t" + "\t".join(wk) + "\tdiff_per_kb\n")
        for w in wins:
            if w["w_p_end"] - w["w_p_start"] < a.min_window or w["w_c_end"] - w["w_c_start"] < a.min_window:
                continue
            n += 1
            wid = f"W{n:06d}"
            d = 1000 * (w["mm"] + w["ins_ev"] + w["del_ev"]) / w["aligned"]
            o.write(wid + "\t" + "\t".join(str(w[k]) for k in wk) + f"\t{d:.3f}\n")
            lb.write(f"{w['chrom']}\t{w['w_p_start']}\t{w['w_p_end']}\t{wid}|src=primary\n")
            lb.write(f"{w['contig']}\t{w['w_c_start']}\t{w['w_c_end']}\t{wid}|src=contig\n")
    with open(f"{a.out}_target.bed", "w") as t:
        for c, s, e, *_ in sorted(set(rows)):
            t.write(f"{c}\t{s}\t{e}\n")
    print(f"{len(segs)} segments, {n} windows", file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) == 3 and sys.argv[1] == "--one":
        job, chunk = json.loads(sys.argv[2])
        job[1] = tuple(job[1])
        print(json.dumps(one(tuple(job), chunk)))
    else:
        main()
