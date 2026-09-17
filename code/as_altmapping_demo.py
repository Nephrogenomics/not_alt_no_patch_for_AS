#!/usr/bin/env python3
"""Figure 2 — emulate the ONT adaptive-sampling (AS) keep/eject decision to show that ALT/HLA
contigs in the capture reference are what deplete MHC (CFB/C2) reads.

Mechanism under test (no new sequencing needed): MinKNOW/readfish aligns each nascent read's
PREFIX with minimap2 (map-ont, via mappy) against the capture reference and KEEPS the molecule
only if the best alignment falls in a target-BED region. We re-run exactly that decision, per
read, against two references:
  --noalt-ref : alt-free primary assembly (the fix)
  --alt-ref   : ALT/HLA-bearing capture reference (hs38DH.fa)
and compare, for reads that truly originate from the MHC target:
  (A) MAPQ of the best prefix alignment       -> collapse to 0 on the ALT reference
  (B) where the best alignment lands           -> on-target chr6 vs ALT/HLA/decoy off-target
  (C) simulated AS keep-rate (best hit in BED) -> predicts the observed depletion ratio (~0.34)
Optionally --bwa-postalt shows that ALT-aware post-processing rescues the reads (mechanism =
ALT-handling, not aligner identity).

Read source (choose one):
  --reads-bam BAM --region chr6:31,900,000-31,960,000   (extract real MHC reads from a
        WELL-COVERED whole-genome ONT control; these are bona-fide MHC molecules = ground truth)
  --reads-fastq FQ                                        (a prepared read set)
  --simulate N                                            (draw N prefixes from the no-alt CFB/C2
        target sequence with a light ONT-like error model; ground truth = target by construction)

Run on mesu (needs mappy + pysam; matplotlib for the figure). Example:
  python3 as_altmapping_demo.py \
      --alt-ref /scratch/$USER/ref/hs38DH.fa \
      --noalt-ref /scratch/$USER/ref/GRCh38-GIABv3_std.fa \
      --target-bed CKD_SU_V6_no_alt.bed \
      --reads-bam <wgs_ont_control>.bam --region chr6:31900000-31960000 \
      --prefix-len 400 --out fig2_as_altmapping
"""
import argparse, sys, os, random

def log(*a): print(*a, file=sys.stderr, flush=True)

def load_bed(path):
    """dict contig -> list of (start,end); accepts standard 3-col BED."""
    iv = {}
    with open(path) as fh:
        for ln in fh:
            ln = ln.strip()
            if not ln or ln.startswith(("#", "track", "browser")): continue
            c = ln.split("\t")
            if len(c) < 3: continue
            iv.setdefault(c[0], []).append((int(c[1]), int(c[2])))
    return iv

def in_bed(bed, ctg, st, en):
    for (a, b) in bed.get(ctg, ()):
        if st < b and en > a:
            return True
    return False

def is_alt_contig(ctg):
    c = ctg.lower()
    return ("_alt" in c) or c.startswith("hla-") or ("_decoy" in c) or ("chrun" in c) or ("_random" in c)

def get_reads(args):
    """Yield (name, seq) read prefixes."""
    plen = args.prefix_len
    if args.reads_bam:
        import pysam
        reg = args.region.replace(",", "")
        chrom, rng = reg.split(":"); a, b = rng.split("-")
        bam = pysam.AlignmentFile(args.reads_bam, "rb")
        seen = set()
        for r in bam.fetch(chrom, int(a), int(b)):
            if r.is_secondary or r.is_supplementary or r.query_sequence is None: continue
            if r.query_name in seen: continue
            seen.add(r.query_name)
            s = r.query_sequence
            # emulate the nascent prefix: leading bases in sequencing order
            yield r.query_name, (s[:plen] if not r.is_reverse else s[-plen:])
    elif args.reads_fastq:
        import pysam
        with pysam.FastxFile(args.reads_fastq) as fq:
            for e in fq:
                yield e.name, e.sequence[:plen]
    elif args.simulate:
        import pysam
        fa = pysam.FastaFile(args.noalt_ref)
        reg = (args.region or "chr6:31900000-31960000").replace(",", "")
        chrom, rng = reg.split(":"); a, b = (int(x) for x in rng.split("-"))
        random.seed(1)
        for i in range(args.simulate):
            st = random.randint(a, max(a, b - plen))
            seq = fa.fetch(chrom, st, st + plen).upper()
            if len(seq) < plen or set(seq) - set("ACGTN"): continue
            yield f"sim_{i}_{chrom}_{st}", ont_errorize(seq, random)
    else:
        log("ERROR: give --reads-bam+--region, --reads-fastq, or --simulate N"); sys.exit(2)

def ont_errorize(seq, rng, sub=0.03, ins=0.02, dele=0.02):
    """Light ONT-like error model for --simulate (ground-truth prefixes)."""
    out = []
    for ch in seq:
        r = rng.random()
        if r < dele: continue
        if r < dele + ins: out.append(rng.choice("ACGT")); out.append(ch); continue
        if r < dele + ins + sub: out.append(rng.choice([x for x in "ACGT" if x != ch])); continue
        out.append(ch)
    return "".join(out)

def decide(aligner, seq, bed, mapq_min=0):
    """Emulate the AS decision on one prefix. Returns (status, best_ctg, mapq, r_st).
    status: 'keep' (best hit in target BED, MAPQ>=mapq_min), 'eject_ambiguous' (best hit IS in
    target but MAPQ<mapq_min = a multi-mapper the decision cannot confidently place — the
    ALT-in-BED case), 'eject_alt', 'eject_offtarget', 'unmapped'.
    mapq_min=0 emulates vanilla readfish (keep on target-region match, no MAPQ gate)."""
    best = None
    for h in aligner.map(seq):              # mappy yields primary + secondaries
        if best is None or h.mapq > best.mapq or (h.mapq == best.mapq and (h.r_en - h.r_st) > (best.r_en - best.r_st)):
            if getattr(h, "is_primary", True) or best is None:
                best = h
    if best is None:
        return "unmapped", None, 0, -1
    intgt = in_bed(bed, best.ctg, best.r_st, best.r_en)
    if intgt:
        return ("keep" if best.mapq >= mapq_min else "eject_ambiguous"), best.ctg, best.mapq, best.r_st
    return ("eject_alt" if is_alt_contig(best.ctg) else "eject_offtarget"), best.ctg, best.mapq, best.r_st

def run_ref(ref, reads, bed, preset, best_n, mapq_min=0):
    import mappy
    log(f"  indexing {os.path.basename(ref)} (map-ont)...")
    al = mappy.Aligner(fn_idx_in=ref, preset=preset, best_n=best_n)
    if not al: log(f"ERROR: cannot open/index {ref}"); sys.exit(3)
    rows = []
    for name, seq in reads:
        if not seq or len(seq) < 50: continue
        st, ctg, mq, rst = decide(al, seq, bed, mapq_min)
        rows.append((name, st, ctg or "*", mq, rst))
    return rows

def summarize(tag, rows):
    from collections import Counter
    c = Counter(r[1] for r in rows)
    n = len(rows); keep = c["keep"]
    mapped = n - c["unmapped"]
    kr = keep / mapped if mapped else 0.0
    log(f"[{tag}] reads={n} mapped={mapped} keep={keep} eject_alt={c['eject_alt']} "
        f"eject_ambiguous={c['eject_ambiguous']} eject_offtarget={c['eject_offtarget']} "
        f"unmapped={c['unmapped']} => keep-rate={kr:.3f}")
    return {"n": n, "mapped": mapped, "keep": keep, "eject_alt": c["eject_alt"],
            "eject_ambiguous": c["eject_ambiguous"], "eject_offtarget": c["eject_offtarget"],
            "unmapped": c["unmapped"], "keep_rate": kr}

def make_figure(noalt_rows, alt_rows, s_noalt, s_alt, outstem, observed_ratio):
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    except Exception as e:
        log(f"[fig] matplotlib unavailable ({e}); skipping figure, TSV/summary still written."); return
    C_NOALT, C_ALT = "#1f6fb2", "#c75a34"
    fig, ax = plt.subplots(1, 3, figsize=(12, 3.6))
    # A: MAPQ of best prefix alignment
    mq_n = [r[3] for r in noalt_rows if r[1] != "unmapped"]
    mq_a = [r[3] for r in alt_rows if r[1] != "unmapped"]
    bins = range(0, 65, 5)
    ax[0].hist(mq_n, bins=bins, alpha=.75, color=C_NOALT, label="alt-free")
    ax[0].hist(mq_a, bins=bins, alpha=.65, color=C_ALT, label="ALT-bearing")
    ax[0].set_title("A · best-prefix MAPQ"); ax[0].set_xlabel("MAPQ"); ax[0].set_ylabel("reads"); ax[0].legend(fontsize=8)
    # B: destination of the best hit
    cats = ["keep", "eject_ambiguous", "eject_alt", "eject_offtarget", "unmapped"]
    labs = ["kept\n(in target)", "ambiguous\n(MAPQ<min)", "ALT/HLA\n(off-target)", "other\noff-target", "unmapped"]
    xn = [s_noalt[k] for k in cats]; xa = [s_alt[k] for k in cats]
    import numpy as np
    x = np.arange(len(cats)); w = .38
    ax[1].bar(x - w/2, xn, w, color=C_NOALT, label="alt-free")
    ax[1].bar(x + w/2, xa, w, color=C_ALT, label="ALT-bearing")
    ax[1].set_xticks(x); ax[1].set_xticklabels(labs, fontsize=7)
    ax[1].set_title("B · where the best hit lands"); ax[1].set_ylabel("reads"); ax[1].legend(fontsize=8)
    # C: simulated keep-rate vs observed depletion
    ax[2].bar([0, 1], [s_noalt["keep_rate"], s_alt["keep_rate"]], color=[C_NOALT, C_ALT], width=.6)
    ax[2].set_xticks([0, 1]); ax[2].set_xticklabels(["alt-free", "ALT-bearing"], fontsize=8)
    ax[2].set_ylim(0, 1.05); ax[2].set_title("C · simulated AS keep-rate")
    ax[2].set_ylabel("kept / mapped")
    ratio = (s_alt["keep_rate"] / s_noalt["keep_rate"]) if s_noalt["keep_rate"] else float("nan")
    ax[2].text(1, s_alt["keep_rate"] + .03, f"ratio {ratio:.2f}", ha="center", fontsize=8, color=C_ALT)
    if observed_ratio:
        ax[2].axhline(observed_ratio * s_noalt["keep_rate"], ls="--", lw=1, color="#666")
        ax[2].text(0.5, observed_ratio * s_noalt["keep_rate"] + .02,
                   f"observed depletion {observed_ratio:.2f}", ha="center", fontsize=7, color="#666")
    fig.tight_layout()
    for ext in ("pdf", "png", "svg"):
        fig.savefig(f"{outstem}.{ext}", dpi=200, bbox_inches="tight")
    log(f"[fig] wrote {outstem}.pdf/.png/.svg")

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--alt-ref", required=True, help="ALT/HLA-bearing capture reference (hs38DH.fa)")
    ap.add_argument("--noalt-ref", required=True, help="alt-free primary assembly (the fix)")
    ap.add_argument("--target-bed", required=True, help="the AS target BED (primary-assembly coords)")
    ap.add_argument("--reads-bam"); ap.add_argument("--region", default=None)
    ap.add_argument("--reads-fastq"); ap.add_argument("--simulate", type=int, default=0)
    ap.add_argument("--prefix-len", type=int, default=400, help="nascent-read prefix length (MinKNOW ~400)")
    ap.add_argument("--preset", default="map-ont"); ap.add_argument("--best-n", type=int, default=5)
    ap.add_argument("--mapq-min", type=int, default=0,
                    help="MAPQ gate for the keep decision. 0 = vanilla readfish (keep on any target-region "
                         "match). >0 emulates an AS/downstream path that rejects ambiguous multi-mappers: an "
                         "in-target hit with MAPQ<min becomes 'eject_ambiguous' — the ALT-in-BED failure mode.")
    ap.add_argument("--observed-ratio", type=float, default=0.34, help="observed MHC depletion ratio for the dashed reference line")
    ap.add_argument("--misplace-kb", type=float, default=10.0,
                    help="[simulate mode] a read is 'correctly placed' only if kept AND its best hit is on the "
                         "true contig within this many kb of its simulated origin; else mis-placed/lost. Reports "
                         "the locus-integrity metric (the one that matters at segmental-duplication loci like RCA).")
    ap.add_argument("--out", default="fig2_as_altmapping")
    args = ap.parse_args()

    # validate inputs early with a clear message (an empty path usually means a shell $(find ...)
    # returned nothing — set the path explicitly) instead of a downstream pysam/mappy traceback.
    for lbl, p in (("--alt-ref", args.alt_ref), ("--noalt-ref", args.noalt_ref), ("--target-bed", args.target_bed)):
        if not p or not os.path.exists(p):
            log(f"ERROR: {lbl} path is empty or missing: '{p}'. Set it explicitly (a $(find ...) that "
                f"returned nothing yields an empty string).")
            sys.exit(2)
    for lbl, p in (("--reads-bam", args.reads_bam), ("--reads-fastq", args.reads_fastq)):
        if p and not os.path.exists(p):
            log(f"ERROR: {lbl} not found: '{p}'"); sys.exit(2)

    bed = load_bed(args.target_bed)
    log(f"target BED: {sum(len(v) for v in bed.values())} intervals over {len(bed)} contigs")
    reads = list(get_reads(args))
    log(f"prefixes to test: {len(reads)} (prefix_len={args.prefix_len})")
    if not reads: log("ERROR: no reads"); sys.exit(2)

    log(f"keep-rule: best hit in target BED AND MAPQ >= {args.mapq_min}"
        + ("  (vanilla readfish)" if args.mapq_min == 0 else "  (MAPQ-gated: multi-mappers ejected)"))
    log("== mapping against the ALT-FREE reference (the fix) ==")
    noalt_rows = run_ref(args.noalt_ref, reads, bed, args.preset, args.best_n, args.mapq_min)
    log("== mapping against the ALT-BEARING reference (hs38DH) ==")
    alt_rows = run_ref(args.alt_ref, reads, bed, args.preset, args.best_n, args.mapq_min)

    s_noalt = summarize("alt-free", noalt_rows)
    s_alt = summarize("ALT-bearing", alt_rows)
    ratio = (s_alt["keep_rate"] / s_noalt["keep_rate"]) if s_noalt["keep_rate"] else float("nan")
    log(f"\n>>> in-silico keep-rate ratio (ALT / alt-free) = {ratio:.3f}   (observed coverage depletion ~{args.observed_ratio})")

    # [locus-integrity metric] For simulated reads the true origin is encoded in the name
    # (sim_<i>_<chrom>_<pos>). A read "correctly placed" = kept AND best hit on the true contig
    # within --misplace-kb of its origin. Everything else (unmapped, off-target, or kept-but-far =
    # mis-mapped to a paralog) means its CFHR3-CFHR1 signal is lost/scrambled. This is the metric
    # that matters when a broad segmental-duplication BED masks the depletion as apparent 'keep'.
    def integrity(rows, thresh_bp):
        tot = ok = 0
        for name, st, ctg, mq, rst in rows:
            if not name.startswith("sim_"): return None
            parts = name.split("_")
            try: tpos = int(parts[-1]); tchr = parts[-2]
            except Exception: continue
            tot += 1
            if st == "keep" and ctg == tchr and rst is not None and rst >= 0 and abs(rst - tpos) <= thresh_bp:
                ok += 1
        return (ok, tot)
    thr = args.misplace_kb * 1000
    for tag, rows in (("alt-free", noalt_rows), ("ALT-bearing", alt_rows)):
        r = integrity(rows, thr)
        if r:
            ok, tot = r
            log(f"[locus-integrity {tag}] correctly placed at true CFHR3-CFHR1 locus "
                f"(<= {args.misplace_kb} kb): {ok}/{tot} = {100*ok/tot:.1f}%  "
                f"-> mis-placed/lost {100*(tot-ok)/tot:.1f}%")

    # per-read TSV (join the two references by read name); positions let you check whether an
    # ALT-inclusive BED would cover where the ALT-mapped reads land.
    d = {r[0]: r for r in alt_rows}
    with open(f"{args.out}_perread.tsv", "w") as fh:
        fh.write("read\tnoalt_status\tnoalt_ctg\tnoalt_pos\tnoalt_mapq\talt_status\talt_ctg\talt_pos\talt_mapq\n")
        for r in noalt_rows:
            a = d.get(r[0], (r[0], "NA", "NA", 0, -1))
            fh.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[4]}\t{r[3]}\t{a[1]}\t{a[2]}\t{a[4]}\t{a[3]}\n")
    log(f"[tsv] wrote {args.out}_perread.tsv")

    make_figure(noalt_rows, alt_rows, s_noalt, s_alt, args.out, args.observed_ratio)

if __name__ == "__main__":
    main()
