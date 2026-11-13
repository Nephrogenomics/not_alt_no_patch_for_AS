#!/usr/bin/env python3
"""Emulate the nanopore adaptive-sampling (AS) keep/eject decision, per read and per locus,
against one or several capture references and under several decision rules.

Why several rules. The study runs used MinKNOW's built-in adaptive sampling, whose decision
algorithm is not public. readfish (Payne et al. 2021) is public: it classifies each read from ALL
alignments mappy returns (primary + secondary), as single_on / single_off / multi_on / multi_off /
no_map, and its documented configuration keeps multi_on ("at least one alignment in a target
region"); it applies no MAPQ filter. A rule that looks only at the primary hit behaves very
differently at loci with near-identical ALT/HLA copies. We therefore never call a single rule
"the" AS rule; we report all of them and let the observed coverage tell which one MinKNOW follows.

  primary    keep iff the PRIMARY alignment overlaps a target interval            (MinKNOW-like hypothesis)
  readfish   readfish make_decision(): coordinate = r_en (+) / r_st (-) of EVERY returned
             alignment; keep iff single_on or multi_on                            (readfish default TOML)
  mapqN      keep iff the primary alignment overlaps a target AND MAPQ >= N       (strict hypothesis)

readfish depends on how many secondary alignments mappy returns (best_n); --best-n accepts a list
and every rule is evaluated for each value.

Capture versus placement. The keep decision (capture) and the final placement of a read (alignment)
are different questions. A kept read is sequenced to full length and re-aligned downstream on the
analysis reference whatever the capture reference said. Two metrics are therefore reported and must
not be confused:
  keep_rate       fraction of molecules kept                          -> what capture loses (irreversible)
  placement_ok    primary hit of the prefix on the true locus         -> alignment-level only; meaningful
                  only when the tested reference contains the origin locus in the same coordinates

Read sources (one):
  --simulate N --loci-bed BED      N molecules per locus drawn from --source-ref (default: first --ref),
                                   both strands, light ONT-like error model; ground truth by construction
  --reads-bam BAM --loci-bed BED   real reads from a whole-genome ONT BAM aligned to the primary assembly;
                                   prefix = first bases in sequencing order (reverse-strand reads are
                                   reverse-complemented back); origin = alignment start in that orientation
  --reads-fastq FQ                 prepared reads, no ground truth (locus = 'fastq')

Expected relative depth. If a fraction k of molecules is kept and full-length (L_full) while the rest
contributes only the ejected stub (L_ej), the depth relative to a fully kept control locus is
  k + (1 - k) * L_ej / L_full
reported as pred_depth_ratio (unmapped prefixes counted as ejected). Set --len-ejected / --len-full
from the data (median read length at the depleted locus / at the control locus).

Example (MHC, three references, two readfish best_n):
  python3 as_emulate.py \
     --ref noalt=work/mhc_noalt.fa --ref hs38DH=work/mhc_ALT.fa --ref hs38DH_noHLA=work/mhc_ALT_noHLA.fa \
     --target-bed data/CKD_SU_V6_no_alt.bed --loci-bed loci_mhc.bed --simulate 1000 \
     --prefix-len 200,400,800 --best-n 5,20 --threads 8 --out results/mhc
"""
import argparse, bisect, gzip, os, random, sys, threading
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def revcomp(s):
    return s.translate(COMP)[::-1]


def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


# ----------------------------------------------------------------------------------------- BED
class Intervals:
    """Merged, sorted intervals per contig with O(log n) overlap / point queries."""

    def __init__(self, path):
        raw = defaultdict(list)
        with opener(path) as fh:
            for ln in fh:
                if not ln.strip() or ln.startswith(("#", "track", "browser")):
                    continue
                c = ln.rstrip("\n").split("\t")
                if len(c) >= 3:
                    raw[c[0]].append((int(c[1]), int(c[2])))
        self.st, self.en = {}, {}
        for ctg, iv in raw.items():
            iv.sort()
            m = [list(iv[0])]
            for a, b in iv[1:]:
                if a <= m[-1][1]:
                    m[-1][1] = max(m[-1][1], b)
                else:
                    m.append([a, b])
            self.st[ctg] = [x[0] for x in m]
            self.en[ctg] = [x[1] for x in m]
        self.n = sum(len(v) for v in self.st.values())

    def overlaps(self, ctg, a, b):
        s = self.st.get(ctg)
        if not s:
            return False
        i = bisect.bisect_right(s, b - 1) - 1
        return i >= 0 and self.en[ctg][i] > a

    def contains(self, ctg, x):
        s = self.st.get(ctg)
        if not s:
            return False
        i = bisect.bisect_right(s, x) - 1
        return i >= 0 and x <= self.en[ctg][i]


def read_loci(path):
    loci = []
    with opener(path) as fh:
        for ln in fh:
            if not ln.strip() or ln.startswith(("#", "track", "browser")):
                continue
            c = ln.rstrip("\n").split("\t")
            a, b = int(c[1]), int(c[2])
            loci.append((c[3] if len(c) > 3 and c[3] else f"{c[0]}:{a}-{b}", c[0], a, b))
    return loci


# ------------------------------------------------------------------------------ contig classes
def contig_class(ctg):
    c = ctg.lower()
    if c.startswith("hla-"):
        return "hla"
    if c.endswith("_alt"):
        return "alt"           # GRCh38 ALT scaffolds (and UCSC novel patches, also named *_alt)
    if c.endswith("_fix"):
        return "fix"           # GRC fix patches (UCSC hg38.p14 and later)
    if "_decoy" in c or c == "chrebv":
        return "decoy"
    if c.startswith("chrun") or c.endswith("_random"):
        return "unplaced"
    return "primary"


CLASSES = ["primary", "alt", "hla", "fix", "decoy", "unplaced"]


# ---------------------------------------------------------------------------------- read source
def errorize(seq, rng, rate):
    """ONT-like errors: total per-base rate split 1/2 substitution, 1/4 insertion, 1/4 deletion."""
    if rate <= 0:
        return seq
    sub, ins, dele = rate / 2, rate / 4, rate / 4
    out = []
    for ch in seq:
        r = rng.random()
        if r < dele:
            continue
        if r < dele + ins:
            out.append(rng.choice("ACGT"))
            out.append(ch)
            continue
        if r < dele + ins + sub:
            out.append(rng.choice([x for x in "ACGT" if x != ch] or "A"))
            continue
        out.append(ch)
    return "".join(out)


def simulate(args, loci, plens):
    """Yield molecules: dict(id, locus, strand, seq, origins{plen: (ctg, st, en)})."""
    import pysam
    fa = pysam.FastaFile(args.source_ref)
    lens = dict(zip(fa.references, fa.lengths))
    rng = random.Random(args.seed)
    pmax, pmin = max(plens), min(plens)
    k = 0
    for name, ctg, a, b in loci:
        if ctg not in lens:
            log(f"  [warn] {name}: {ctg} not in source reference, skipped")
            continue
        clen = lens[ctg]
        for _ in range(args.simulate):
            strand = "+" if rng.random() < 0.5 else "-"
            if strand == "+":
                g = rng.randint(a, max(a, b - pmin))
                gs, ge = g, min(clen, g + pmax)
            else:
                e = rng.randint(min(b, a + pmin), b)
                gs, ge = max(0, e - pmax), e
            seq = fa.fetch(ctg, gs, ge).upper()
            if len(seq) < pmin or seq.count("N") > 0.1 * len(seq):
                continue
            mol = seq if strand == "+" else revcomp(seq)
            org = {p: ((ctg, gs, gs + p) if strand == "+" else (ctg, ge - p, ge)) for p in plens}
            k += 1
            yield {"id": f"sim{k}", "locus": name, "strand": strand,
                   "seq": errorize(mol, rng, args.error_rate), "origins": org}


def from_bam(args, loci, plens):
    import pysam
    bam = pysam.AlignmentFile(args.reads_bam, "rb")
    pmin = min(plens)
    seen = set()
    for name, ctg, a, b in loci:
        n = 0
        for r in bam.fetch(ctg, a, b):
            if r.is_secondary or r.is_supplementary or r.is_unmapped or r.query_sequence is None:
                continue
            if r.query_name in seen:
                continue
            # sequencing start of the molecule on the genome, in its own orientation
            start_g = r.reference_end if r.is_reverse else r.reference_start
            if args.require_prefix_in_locus and not (a <= start_g <= b):
                continue
            seq = revcomp(r.query_sequence) if r.is_reverse else r.query_sequence
            if len(seq) < pmin:
                continue
            seen.add(r.query_name)
            org = {p: ((ctg, start_g - p, start_g) if r.is_reverse else (ctg, start_g, start_g + p))
                   for p in plens}
            yield {"id": r.query_name, "locus": name, "strand": "-" if r.is_reverse else "+",
                   "seq": seq, "origins": org}
            n += 1
            if args.per_locus and n >= args.per_locus:
                break
        log(f"  {name}: {n} reads")


def from_fastq(args, plens):
    import pysam
    with pysam.FastxFile(args.reads_fastq) as fq:
        for e in fq:
            yield {"id": e.name, "locus": "fastq", "strand": ".", "seq": e.sequence,
                   "origins": {p: None for p in plens}}


# ------------------------------------------------------------------------------------ decisions
def decide(hits, bed, mapq_min):
    """Return the per-read record for one prefix under every rule."""
    prim = next((h for h in hits if h.is_primary), hits[0] if hits else None)
    if prim is None:
        return {"n_hits": 0, "prim": None, "rf_class": "no_map",
                "keep_primary": False, "keep_readfish": False, "keep_mapq": False}
    prim_on = bed.overlaps(prim.ctg, prim.r_st, prim.r_en)
    # readfish make_decision(): end coordinate in sequencing direction, any() over all alignments
    on = any(bed.contains(h.ctg, h.r_st if h.strand == -1 else h.r_en) for h in hits)
    rf = ("single_" if len(hits) == 1 else "multi_") + ("on" if on else "off")
    return {"n_hits": len(hits), "prim": prim, "rf_class": rf,
            "keep_primary": prim_on, "keep_readfish": on,
            "keep_mapq": prim_on and prim.mapq >= mapq_min}


def placement_ok(prim, origin, tol):
    if prim is None or origin is None:
        return None
    ctg, a, b = origin
    return prim.ctg == ctg and prim.r_st < b + tol and prim.r_en > a - tol


# ------------------------------------------------------------------------------------- mapping
def load_aligner(path, preset, best_n, threads_idx):
    import mappy
    al = mappy.Aligner(fn_idx_in=path, preset=preset, best_n=best_n, n_threads=threads_idx)
    if not al:
        log(f"ERROR: cannot open/index {path}")
        sys.exit(3)
    return al


def ensure_index(name, path, preset, index_dir):
    """Build a .mmi once so that several best_n values do not re-index the FASTA."""
    if path.endswith(".mmi"):
        return path
    import mappy
    mmi = os.path.join(index_dir, f"{name}.{preset}.mmi")
    if not os.path.exists(mmi):
        log(f"  indexing {path} -> {mmi}")
        mappy.Aligner(fn_idx_in=path, preset=preset, fn_idx_out=mmi)
    return mmi


def map_all(al, seqs, threads):
    import mappy
    local = threading.local()

    def one(seq):
        if not hasattr(local, "buf"):
            local.buf = mappy.ThreadBuffer()
        return list(al.map(seq, buf=local.buf))

    with ThreadPoolExecutor(max_workers=threads) as ex:
        return list(ex.map(one, seqs, chunksize=64))


# ------------------------------------------------------------------------------------- summary
def new_acc():
    d = {"n": 0, "unmapped": 0, "keep_primary": 0, "keep_readfish": 0, "keep_mapq": 0,
         "mapq0": 0, "placement_ok": 0, "placement_n": 0, "hits": 0}
    for c in CLASSES:
        d["prim_" + c] = 0
    d["prim_primary_target"] = 0
    for c in ("single_on", "single_off", "multi_on", "multi_off", "no_map"):
        d["rf_" + c] = 0
    return d


def add(acc, rec, pl_ok):
    acc["n"] += 1
    acc["hits"] += rec["n_hits"]
    acc["rf_" + rec["rf_class"]] += 1
    for k in ("keep_primary", "keep_readfish", "keep_mapq"):
        acc[k] += rec[k]
    p = rec["prim"]
    if p is None:
        acc["unmapped"] += 1
    else:
        acc["prim_" + contig_class(p.ctg)] += 1
        acc["prim_primary_target"] += rec["keep_primary"] and contig_class(p.ctg) == "primary"
        acc["mapq0"] += p.mapq == 0
    if pl_ok is not None:
        acc["placement_n"] += 1
        acc["placement_ok"] += pl_ok


def write_summary(path, accs, args, refnames):
    rules = [("primary", "keep_primary"), ("readfish", "keep_readfish"), (f"mapq{args.mapq_min}", "keep_mapq")]
    r = args.len_ejected / args.len_full
    base = refnames[0] if args.baseline is None else args.baseline
    rows = []
    for (locus, ref, bn, pl), a in sorted(accs.items()):
        n = a["n"] or 1
        for rule, key in rules:
            k = a[key] / n
            b = accs.get((locus, base, bn, pl))
            kb = (b[key] / (b["n"] or 1)) if b else float("nan")
            rows.append([locus, ref, rule, bn, pl, a["n"], f"{k:.4f}",
                         f"{(k / kb) if kb else float('nan'):.4f}",
                         f"{k + (1 - k) * r:.4f}", f"{a['unmapped'] / n:.4f}",
                         f"{a['mapq0'] / n:.4f}", f"{a['hits'] / n:.2f}"]
                        + [f"{a['prim_' + c] / n:.4f}" for c in CLASSES]
                        + [f"{a['rf_multi_on'] / n:.4f}", f"{a['rf_multi_off'] / n:.4f}",
                           (f"{a['placement_ok'] / a['placement_n']:.4f}" if a["placement_n"] else "NA")])
    hdr = (["locus", "ref", "rule", "best_n", "prefix_len", "n", "keep_rate", f"keep_vs_{base}",
            "pred_depth_ratio", "frac_unmapped", "frac_mapq0", "mean_hits"]
           + [f"prim_on_{c}" for c in CLASSES] + ["rf_multi_on", "rf_multi_off", "placement_ok"])
    with open(path, "w") as fh:
        fh.write("\t".join(hdr) + "\n")
        for row in rows:
            fh.write("\t".join(str(x) for x in row) + "\n")
    return hdr, rows


# -------------------------------------------------------------------------------------- figure
def make_figure(hdr, rows, out, refnames, max_loci=40):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except Exception as e:
        log(f"[fig] matplotlib/numpy unavailable ({e}); skipped")
        return
    H = {h: i for i, h in enumerate(hdr)}
    loci = sorted({r[0] for r in rows})
    if len(loci) > max_loci:
        log(f"[fig] {len(loci)} loci > {max_loci}: heatmap skipped (use the summary TSV)")
        return
    bns, pls = sorted({r[3] for r in rows}), sorted({r[4] for r in rows})
    bn, pl = bns[0], (400 if 400 in pls else pls[0])
    rules = []
    for r in rows:
        if r[2] not in rules:
            rules.append(r[2])
    cols = [(ref, rule) for ref in refnames for rule in rules]
    M = np.full((len(loci), len(cols)), np.nan)
    for r in rows:
        if r[3] == bn and r[4] == pl:
            M[loci.index(r[0]), cols.index((r[1], r[2]))] = float(r[H["keep_rate"]])
    fig, ax = plt.subplots(1, 2, figsize=(4 + 0.55 * len(cols) + 4, 1.2 + 0.32 * len(loci)),
                           gridspec_kw={"width_ratios": [len(cols), 4]})
    im = ax[0].imshow(M, vmin=0, vmax=1, cmap="RdYlBu", aspect="auto")
    ax[0].set_xticks(range(len(cols)))
    ax[0].set_xticklabels([f"{a}\n{b}" for a, b in cols], fontsize=7, rotation=90)
    ax[0].set_yticks(range(len(loci)))
    ax[0].set_yticklabels(loci, fontsize=7)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if not np.isnan(M[i, j]):
                ax[0].text(j, i, f"{M[i, j]:.2f}", ha="center", va="center", fontsize=6)
    ax[0].set_title(f"A · keep-rate (best_n={bn}, prefix {pl} b)", fontsize=9)
    fig.colorbar(im, ax=ax[0], fraction=0.03)
    # B: where the primary hit lands, per reference, pooled over loci (rule-independent)
    seen, pooled = set(), {ref: np.zeros(len(CLASSES) + 1) for ref in refnames}
    for r in rows:
        key = (r[0], r[1])
        if r[3] == bn and r[4] == pl and key not in seen:
            seen.add(key)
            n = r[H["n"]]
            pooled[r[1]] += n * np.array([float(r[H[f"prim_on_{c}"]]) for c in CLASSES] + [float(r[H["frac_unmapped"]])])
    labs = CLASSES + ["unmapped"]
    colors = ["#1f6fb2", "#c75a34", "#e0a43a", "#7d5ba6", "#999999", "#cccccc", "#333333"]
    left = np.zeros(len(refnames))
    for k, lab in enumerate(labs):
        v = np.array([pooled[ref][k] / max(1, pooled[ref].sum()) for ref in refnames])
        ax[1].barh(range(len(refnames)), v, left=left, color=colors[k], label=lab)
        left += v
    ax[1].set_yticks(range(len(refnames)))
    ax[1].set_yticklabels(refnames, fontsize=8)
    ax[1].set_xlim(0, 1)
    ax[1].set_xlabel("fraction of prefixes", fontsize=8)
    ax[1].set_title("B · contig class of the primary hit", fontsize=9)
    ax[1].legend(fontsize=6, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=4, frameon=False)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(f"{out}.{ext}", dpi=200, bbox_inches="tight")
    log(f"[fig] wrote {out}.pdf/.png")


# ---------------------------------------------------------------------------------------- main
def parse_list(s, typ=int):
    return [typ(x) for x in str(s).split(",") if x != ""]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref", action="append", required=True, metavar="NAME=PATH",
                    help="capture reference to test (FASTA or .mmi); repeat. The first is the baseline.")
    ap.add_argument("--target-bed", required=True, help="AS target BED as loaded on the instrument")
    ap.add_argument("--ref-bed", action="append", default=[], metavar="NAME=BED",
                    help="target BED for one reference in its own coordinates (e.g. CHM13)")
    ap.add_argument("--extra-bed", action="append", default=[], metavar="TAG=BED",
                    help="also decide every read against this target BED, on the same alignments; "
                         "reported as reference '<ref>|<TAG>' (e.g. a BED projected onto ALT/patch contigs)")
    ap.add_argument("--loci-bed", help="loci to sample and summarise (name in column 4)")
    ap.add_argument("--simulate", type=int, default=0, help="molecules simulated per locus")
    ap.add_argument("--source-ref", help="FASTA for simulated reads (default: first --ref, must be FASTA)")
    ap.add_argument("--reads-bam")
    ap.add_argument("--reads-fastq")
    ap.add_argument("--per-locus", type=int, default=500, help="max real reads per locus (0 = all)")
    ap.add_argument("--require-prefix-in-locus", type=int, default=1)
    ap.add_argument("--prefix-len", default="400", help="comma list, e.g. 200,400,800 (MinKNOW chunk ~0.4-1 s at 400 b/s)")
    ap.add_argument("--best-n", default="5", help="comma list of mappy best_n values (affects readfish)")
    ap.add_argument("--mapq-min", type=int, default=1)
    ap.add_argument("--preset", default="map-ont")
    ap.add_argument("--error-rate", type=float, default=0.05, help="simulated per-base error (live basecaller ~0.03-0.08)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--placement-tol", type=int, default=1000)
    ap.add_argument("--len-ejected", type=float, default=600.0)
    ap.add_argument("--len-full", type=float, default=5700.0)
    ap.add_argument("--baseline", help="reference name used for keep_vs_<baseline> (default: first --ref)")
    ap.add_argument("--threads", type=int, default=os.cpu_count() or 4)
    ap.add_argument("--index-dir", help="where to store .mmi indexes (default: next to --out)")
    ap.add_argument("--perread", action="store_true", help="also write the per-read long table (.tsv.gz)")
    ap.add_argument("--no-figure", action="store_true")
    ap.add_argument("--out", required=True, help="output prefix")
    args = ap.parse_args()

    refs = []
    for s in args.ref:
        name, _, path = s.partition("=")
        if not path or not os.path.exists(path):
            log(f"ERROR: --ref {s!r}: give NAME=PATH to an existing file")
            sys.exit(2)
        refs.append((name, path))
    refnames = [n for n, _ in refs]
    for p, lbl in ((args.target_bed, "--target-bed"), (args.loci_bed, "--loci-bed"),
                   (args.reads_bam, "--reads-bam"), (args.reads_fastq, "--reads-fastq")):
        if p and not os.path.exists(p):
            log(f"ERROR: {lbl} not found: {p}")
            sys.exit(2)
    plens, bns = parse_list(args.prefix_len), parse_list(args.best_n)
    default_bed = Intervals(args.target_bed)
    beds = {n: default_bed for n in refnames}
    for s in args.ref_bed:
        n, _, p = s.partition("=")
        beds[n] = Intervals(p)
    extra = []
    for s_ in args.extra_bed:
        t, _, p = s_.partition("=")
        extra.append((t, Intervals(p)))
    log(f"target BED: {default_bed.n} merged intervals; refs: {', '.join(refnames)}; "
        f"prefix {plens}; best_n {bns}; mapq gate {args.mapq_min}")

    if args.simulate:
        if not args.loci_bed:
            log("ERROR: --simulate needs --loci-bed")
            sys.exit(2)
        args.source_ref = args.source_ref or refs[0][1]
        if args.source_ref.endswith(".mmi"):
            log("ERROR: --source-ref must be a FASTA")
            sys.exit(2)
        mols = list(simulate(args, read_loci(args.loci_bed), plens))
    elif args.reads_bam:
        if not args.loci_bed:
            log("ERROR: --reads-bam needs --loci-bed")
            sys.exit(2)
        mols = list(from_bam(args, read_loci(args.loci_bed), plens))
    elif args.reads_fastq:
        mols = list(from_fastq(args, plens))
    else:
        log("ERROR: give --simulate, --reads-bam or --reads-fastq")
        sys.exit(2)
    log(f"molecules: {len(mols)} over {len({m['locus'] for m in mols})} loci")
    if not mols:
        sys.exit(2)

    index_dir = args.index_dir or (os.path.dirname(os.path.abspath(args.out)) or ".")
    os.makedirs(index_dir, exist_ok=True)
    accs = defaultdict(new_acc)
    pr = gzip.open(f"{args.out}_perread.tsv.gz", "wt") if args.perread else None
    if pr:
        pr.write("read\tlocus\tstrand\tprefix_len\tref\tbest_n\tn_hits\tprim_ctg\tprim_class\tprim_st\t"
                 "prim_mapq\treadfish_class\tkeep_primary\tkeep_readfish\tkeep_mapq\tplacement_ok\n")
    for name, path in refs:
        idx = ensure_index(name, path, args.preset, index_dir) if len(bns) > 1 else path
        for bn in bns:
            log(f"== {name} (best_n={bn}) ==")
            al = load_aligner(idx, args.preset, bn, min(args.threads, 8))
            for pl in plens:
                seqs = [m["seq"][:pl] for m in mols]
                hits = map_all(al, seqs, args.threads)
                for m, h in zip(mols, hits):
                    rec = decide(h, beds[name], args.mapq_min)
                    ok = placement_ok(rec["prim"], m["origins"][pl], args.placement_tol)
                    add(accs[(m["locus"], name, bn, pl)], rec, ok)
                    for t, xb in extra:
                        add(accs[(m["locus"], f"{name}|{t}", bn, pl)], decide(h, xb, args.mapq_min), ok)
                    if pr:
                        p = rec["prim"]
                        pr.write("\t".join(str(x) for x in (
                            m["id"], m["locus"], m["strand"], pl, name, bn, rec["n_hits"],
                            p.ctg if p else "*", contig_class(p.ctg) if p else "unmapped",
                            p.r_st if p else -1, p.mapq if p else 0, rec["rf_class"],
                            int(rec["keep_primary"]), int(rec["keep_readfish"]), int(rec["keep_mapq"]),
                            "NA" if ok is None else int(ok))) + "\n")
                for rn in [name] + [f"{name}|{t}" for t, _ in extra]:
                    tot = [a for (l, r, b, p), a in accs.items() if r == rn and b == bn and p == pl]
                    n = sum(a["n"] for a in tot) or 1
                    log(f"  {rn:<24} prefix {pl:>5}: keep primary {sum(a['keep_primary'] for a in tot) / n:.3f} | "
                        f"readfish {sum(a['keep_readfish'] for a in tot) / n:.3f} | "
                        f"mapq>={args.mapq_min} {sum(a['keep_mapq'] for a in tot) / n:.3f} | "
                        f"unmapped {sum(a['unmapped'] for a in tot) / n:.3f}")
            del al
    if pr:
        pr.close()
        log(f"[tsv] wrote {args.out}_perread.tsv.gz")
    refnames = refnames + [f"{n}|{t}" for n in refnames for t, _ in extra]
    hdr, rows = write_summary(f"{args.out}_summary.tsv", accs, args, refnames)
    log(f"[tsv] wrote {args.out}_summary.tsv ({len(rows)} rows)")
    if not args.no_figure:
        make_figure(hdr, rows, args.out, refnames)


if __name__ == "__main__":
    main()
